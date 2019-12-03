using System;
using System.Runtime.InteropServices;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Graph.FindMaxDeg;
using static ChacoSharp.Graph.Subgraph;
using static ChacoSharp.Utilities.MergeSort;
using static ChacoSharp.Graph.CountWeights;
using static ChacoSharp.Coarsening.KlSpiff;

namespace ChacoSharp.RefinePartition
{
    public static unsafe class RefinePartitionHelper
    {
        /* Construct a graph representing the inter-set communication. */
        public static bool refine_part(vtx_data** graph, /* graph data structure */
            int nvtxs, /* number of vertices in graph */
            bool useEdgeWeights, /* are edge weights being used? */
            int* assign, /* current assignment */
            int architecture, /* 0 => hypercube, d => d-dimensional mesh */
            int ndims_tot, /* if hypercube, number of dimensions */
            int[] mesh_dims, /* if mesh, size in each direction */
            double[] goal /* desired set sizes */
        )
        {
            bilist* set_list = null; /* lists of vtxs in each set */
            bilist* vtx_elems = null; /* space for all vtxs in set_lists */
            bilist* ptr = null; /* loops through set_lists */
            ipairs* pairs = null; /* ordered list of edges in comm graph */
            double* comm_vals = null; /* edge wgts of comm graph for sorting */
            float*[] term_wgts = new float*[2]; /* terminal propagation vector */
            int[][] hops = new int[MAXSETS][]; /* preference weighting */
            for (var i = 0; i < hops.Length; i++)
            {
                hops[i] = new int[MAXSETS];
            }

            IntPtr temp; /* return argument from srealloc_ret() */
            int* indices = null; /* sorted order for communication edges */
            int* space = null; /* space for mergesort */
            int* sizes = null; /* sizes of the different sets */
            int* sub_assign = null; /* new assignment for subgraph */
            int* old_sub_assign = null; /* room for current sub assignment */
            int** edges_list = null; /* lists of comm graph edges */
            int** ewgts_list = null; /* lists of comm graph edge wgts */
            int* ewgts = null; /* loops through ewgts_list */
            int* edges = null; /* edges in communication graph */
            int* adj_sets = null; /* weights connecting sets */
            int* eptr = null; /* loop through edges and edge weights */
            int* ewptr = null; /* loop through edges and edge weights */
            int ewgt; /* weight of an edge */
            vtx_data** subgraph = null; /* subgraph data structure */
            int* nedges = null; /* space for saving graph data */
            int* degrees = null; /* # neighbors of vertices */
            int* glob2loc = null; /* maps full to reduced numbering */
            int* loc2glob = null; /* maps reduced to full numbering */
            int nmax; /* largest subgraph I expect to encounter */
            int set, set1, set2; /* sets vertices belong to */
            int vertex; /* vertex in graph */
            int ncomm; /* # edges in communication graph */
            int dist = -1; /* architectural distance between two sets */
            int nsets_tot = 0; /* total number of processors */
            bool change; /* did change occur in this pass? */
            bool any_change = false; /* has any change occurred? */
            int error; /* out of space? */
            int size; /* array spacing */

            error = 1;
            term_wgts[1] = null;

            if (architecture == 0)
            {
                nsets_tot = 1 << ndims_tot;
            }
            else if (architecture > 0)
            {
                nsets_tot = mesh_dims[0] * mesh_dims[1] * mesh_dims[2];
            }

            hops[0][0] = hops[1][1] = 0;
            if (!TERM_PROP)
            {
                hops[0][1] = hops[1][0] = 1;
            }

            /* Set up convenient data structure for navigating through sets. */
            set_list = (bilist*) Marshal.AllocHGlobal(nsets_tot * sizeof(bilist));
            vtx_elems = (bilist*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(bilist));
            sizes = (int*) Marshal.AllocHGlobal(nsets_tot * sizeof(int));
            if (set_list == null || vtx_elems == null || sizes == null)
            {
                goto skip;
            }

            for (var i = 0; i < nsets_tot; i++)
            {
                set_list[i].next = null;
                sizes[i] = 0;
            }

            for (var i = 1; i <= nvtxs; i++)
            {
                set = assign[i];
                ++sizes[set];
                vtx_elems[i].next = set_list[set].next;
                if (vtx_elems[i].next != null)
                {
                    vtx_elems[i].next->prev = &(vtx_elems[i]);
                }

                vtx_elems[i].prev = &(set_list[set]);
                set_list[set].next = &(vtx_elems[i]);
            }

            /* For each set, find connections to all set neighbors. */
            edges_list = (int**) Marshal.AllocHGlobal(nsets_tot * sizeof(int*));
            if (edges_list == null)
            {
                goto skip;
            }

            for (set = 0; set < nsets_tot - 1; set++)
            {
                edges_list[set] = null;
            }

            ewgts_list = (int**) Marshal.AllocHGlobal(nsets_tot * sizeof(int*));
            if (ewgts_list == null)
            {
                goto skip;
            }

            for (set = 0; set < nsets_tot - 1; set++)
            {
                ewgts_list[set] = null;
            }

            nedges = (int*) Marshal.AllocHGlobal(nsets_tot * sizeof(int));
            adj_sets = (int*) Marshal.AllocHGlobal(nsets_tot * sizeof(int));
            if (nedges == null || adj_sets == null)
            {
                goto skip;
            }

            size = (int) (&(vtx_elems[1]) - &(vtx_elems[0]));
            ncomm = 0;
            ewgt = 1;
            nmax = 0;
            for (set = 0; set < nsets_tot - 1; set++)
            {
                if (sizes[set] > nmax)
                {
                    nmax = sizes[set];
                }

                for (var i = 0; i < nsets_tot; i++)
                {
                    adj_sets[i] = 0;
                }

                for (ptr = set_list[set].next; ptr != null; ptr = ptr->next)
                {
                    vertex = ((int) (ptr - vtx_elems)) / size;
                    for (var j = 1; j < graph[vertex]->nedges; j++)
                    {
                        set2 = assign[graph[vertex]->edges[j]];
                        if (useEdgeWeights)
                        {
                            ewgt = (int) graph[vertex]->ewgts[j];
                        }

                        adj_sets[set2] += ewgt;
                    }
                }

                /* Now save adj_sets data to later construct graph. */
                var t = 0;
                for (var i = set + 1; i < nsets_tot; i++)
                {
                    if (adj_sets[i] != 0)
                    {
                        t++;
                    }
                }

                nedges[set] = t;
                if (t != 0)
                {
                    edges_list[set] = edges = (int*) Marshal.AllocHGlobal(t * sizeof(int));
                    ewgts_list[set] = ewgts = (int*) Marshal.AllocHGlobal(t * sizeof(int));
                    if (edges == null || ewgts == null)
                    {
                        goto skip;
                    }
                }

                t = 0;
                for (var i = set + 1; i < nsets_tot; i++)
                {
                    if (adj_sets[i] != 0)
                    {
                        edges[t] = i;
                        ewgts[t] = adj_sets[i];
                        t++;
                    }
                }

                ncomm += t;
            }

            Marshal.FreeHGlobal((IntPtr) adj_sets);
            adj_sets = null;

            /* Now compact all communication weight information into single */
            /* vector for sorting. */

            pairs = (ipairs*) Marshal.AllocHGlobal((ncomm + 1) * sizeof(ipairs));
            comm_vals = (double*) Marshal.AllocHGlobal((ncomm + 1) * sizeof(double));
            if (pairs == null || comm_vals == null)
            {
                goto skip;
            }

            var u = 0;
            for (set = 0; set < nsets_tot - 1; set++)
            {
                eptr = edges_list[set];
                ewptr = ewgts_list[set];
                for (var k = 0; k < nedges[set]; k++)
                {
                    set2 = eptr[k];
                    pairs[u].val1 = set;
                    pairs[u].val2 = set2;
                    comm_vals[u] = ewptr[k];
                    u++;
                }
            }

            Marshal.FreeHGlobal((IntPtr) nedges);
            nedges = null;

            indices = (int*) Marshal.AllocHGlobal((ncomm + 1) * sizeof(int));
            space = (int*) Marshal.AllocHGlobal((ncomm + 1) * sizeof(int));
            if (indices == null || space == null)
            {
                goto skip;
            }

            ch_mergesort(comm_vals, ncomm, indices, space);
            Marshal.FreeHGlobal((IntPtr) space);
            Marshal.FreeHGlobal((IntPtr) comm_vals);
            space = null;
            comm_vals = null;

            for (set = 0; set < nsets_tot - 1; set++)
            {
                if (edges_list[set] != null)
                {
                    Marshal.FreeHGlobal((IntPtr) edges_list[set]);
                }

                if (ewgts_list[set] != null)
                {
                    Marshal.FreeHGlobal((IntPtr) ewgts_list[set]);
                }
            }

            Marshal.FreeHGlobal((IntPtr) ewgts_list);
            Marshal.FreeHGlobal((IntPtr) edges_list);
            ewgts_list = null;
            edges_list = null;

            /* 2 for 2 subsets, 20 for safety margin. Should check this at run time. */
            nmax = 2 * nmax + 20;

            subgraph = (vtx_data**) Marshal.AllocHGlobal((nmax + 1) * sizeof(vtx_data*));
            degrees = (int*) Marshal.AllocHGlobal((nmax + 1) * sizeof(int));
            glob2loc = (int*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));
            loc2glob = (int*) Marshal.AllocHGlobal((nmax + 1) * sizeof(int));
            sub_assign = (int*) Marshal.AllocHGlobal((nmax + 1) * sizeof(int));
            old_sub_assign = (int*) Marshal.AllocHGlobal((nmax + 1) * sizeof(int));

            if (subgraph == null || degrees == null || glob2loc == null || loc2glob == null ||
                sub_assign == null || old_sub_assign == null)
            {
                goto skip;
            }

            if (TERM_PROP)
            {
                term_wgts[1] = (float*) Marshal.AllocHGlobal((nmax + 1) * sizeof(float));
                if (term_wgts[1] == null)
                {
                    goto skip;
                }
            }
            else
            {
                term_wgts[1] = null;
            }

            /* Do large boundaries first to encourage convergence. */
            any_change = false;
            for (var i = ncomm - 1; i >= 0; i--)
            {
                var t = indices[i];
                set1 = pairs[t].val1;
                set2 = pairs[t].val2;

                /* Make sure subgraphs aren't too big. */
                if (sizes[set1] + sizes[set2] > nmax)
                {
                    nmax = sizes[set1] + sizes[set2];

                    temp = Marshal.ReAllocHGlobal((IntPtr) subgraph, new IntPtr((nmax + 1) * sizeof(vtx_data*)));
                    if (temp == null)
                    {
                        goto skip;
                    }
                    else
                    {
                        subgraph = (vtx_data**) temp;
                    }

                    temp = Marshal.ReAllocHGlobal((IntPtr) degrees, new IntPtr((nmax + 1) * sizeof(int)));
                    if (temp == null)
                    {
                        goto skip;
                    }
                    else
                    {
                        degrees = (int*) temp;
                    }

                    temp = Marshal.ReAllocHGlobal((IntPtr) loc2glob, new IntPtr((nmax + 1) * sizeof(int)));
                    if (temp == null)
                    {
                        goto skip;
                    }
                    else
                    {
                        loc2glob = (int*) temp;
                    }

                    temp = Marshal.ReAllocHGlobal((IntPtr) sub_assign, new IntPtr((nmax + 1) * sizeof(int)));
                    if (temp == null)
                    {
                        goto skip;
                    }
                    else
                    {
                        sub_assign = (int*) temp;
                    }

                    temp = Marshal.ReAllocHGlobal((IntPtr) old_sub_assign, new IntPtr((nmax + 1) * sizeof(int)));
                    if (temp == null)
                    {
                        goto skip;
                    }
                    else
                    {
                        old_sub_assign = (int*) temp;
                    }

                    if (TERM_PROP)
                    {
                        temp = Marshal.ReAllocHGlobal((IntPtr) term_wgts[1], new IntPtr((nmax + 1) * sizeof(float)));
                        if (temp == null)
                        {
                            goto skip;
                        }
                        else
                        {
                            term_wgts[1] = (float*) temp;
                        }
                    }
                }

                if (TERM_PROP)
                {
                    if (architecture == 0)
                    {
                        var j = set1 ^ set2;
                        dist = 0;
                        while (j != 0)
                        {
                            if ((j & 1) != 0)
                            {
                                dist++;
                            }

                            j >>= 1;
                        }
                    }
                    else if (architecture > 0)
                    {
                        dist = Math.Abs((set1 % mesh_dims[0]) - (set2 % mesh_dims[0]));
                        dist +=
                            Math.Abs(((set1 / mesh_dims[0]) % mesh_dims[1]) - ((set2 / mesh_dims[0]) % mesh_dims[1]));
                        dist +=
                            Math.Abs((set1 / (mesh_dims[0] * mesh_dims[1])) - (set2 / (mesh_dims[0] * mesh_dims[1])));
                    }

                    hops[0][1] = hops[1][0] = dist;
                }

                change = kl_refine(graph, subgraph, set_list, vtx_elems, assign, set1, set2, glob2loc, loc2glob,
                    sub_assign, old_sub_assign, degrees, useEdgeWeights, hops, goal, sizes,
                    term_wgts, architecture, mesh_dims);

                any_change |= change;
            }

            error = 0;

            skip:
            if (error != 0)
            {
                Console.WriteLine("\nWARNING: No space to refine partition.");
                Console.WriteLine("         NO PARTITION REFINEMENT PERFORMED.");
            }

            if (edges_list != null)
            {
                for (set = 0; set < nsets_tot - 1; set++)
                {
                    if (edges_list[set] != null)
                    {
                        Marshal.FreeHGlobal((IntPtr) edges_list[set]);
                    }
                }

                Marshal.FreeHGlobal((IntPtr) edges_list);
            }

            if (ewgts_list != null)
            {
                for (set = 0; set < nsets_tot - 1; set++)
                {
                    if (ewgts_list[set] != null)
                    {
                        Marshal.FreeHGlobal((IntPtr) ewgts_list[set]);
                    }
                }

                Marshal.FreeHGlobal((IntPtr) ewgts_list);
            }

            Marshal.FreeHGlobal((IntPtr) space);
            Marshal.FreeHGlobal((IntPtr) comm_vals);
            Marshal.FreeHGlobal((IntPtr) nedges);
            Marshal.FreeHGlobal((IntPtr) adj_sets);
            Marshal.FreeHGlobal((IntPtr) term_wgts[1]);
            Marshal.FreeHGlobal((IntPtr) old_sub_assign);
            Marshal.FreeHGlobal((IntPtr) sub_assign);
            Marshal.FreeHGlobal((IntPtr) loc2glob);
            Marshal.FreeHGlobal((IntPtr) glob2loc);
            Marshal.FreeHGlobal((IntPtr) degrees);
            Marshal.FreeHGlobal((IntPtr) subgraph);

            Marshal.FreeHGlobal((IntPtr) indices);
            Marshal.FreeHGlobal((IntPtr) pairs);

            Marshal.FreeHGlobal((IntPtr) sizes);
            Marshal.FreeHGlobal((IntPtr) vtx_elems);
            Marshal.FreeHGlobal((IntPtr) set_list);

            return (any_change);
        }

/* Perform KL between two sets. */
        private static bool kl_refine(vtx_data** graph, /* graph data structure */
            vtx_data** subgraph, /* space for subgraph to refine */
            bilist* set_list, /* lists of vtxs in each set */
            bilist* vtx_elems, /* start of storage for lists */
            int* new_assign, /* set assignments for all vertices */
            int set1, int set2, /* two sets being refined */
            int* glob2loc, /* maps vertices to subgraph vertices */
            int* loc2glob, /* maps subgraph vertices to vertices */
            int* sub_assign, /* new assignment for subgraphs */
            int* old_sub_assign, /* current assignment for subgraphs */
            int* degrees, /* space for forming subgraphs */
            bool useEdgeWeights, /* are edge weights being used? */
            int[][] hops, /* KL set preferences */
            double[] goal, /* desired set sizes */
            int* sizes, /* number of vertices in different sets */
            float*[] term_wgts, /* space for terminal propagation weights */
            int architecture, /* 0 => hypercube, d => d-dimensional mesh */
            int[] mesh_dims /*[3]*/ /* if mesh, how big is it? */
        )
        {
            bilist* ptr; /* element in set_list */
            double[] subgoal = new double[2]; /* goal within two subgraphs */
            double[] weights = {0.0, 0.0}; /* weights for each set */
            double maxdeg; /* largest degree of a vertex */
            double ratio; /* set sizes / goals */
            int* null_ptr; /* argument to klspiff */
            int vwgt_max; /* largest vertex weight */
            int max_dev; /* largest set deviation allowed in KL */
            int subnvtxs; /* number of vtxs in subgraph */
            int vwgt_sum1; /* sum of vertex wgts in first set */
            int vwgt_sum2; /* sum of vertex wgts in second set */
            int subnedges; /* number of edges in subgraph */
            int setA, setB; /* two sets being refined */
            int nsame; /* number of vertices not moved */
            int vtx; /* vertex in subgraph */
            int i; /* loop counter */

            /* Compute all the quantities I'll need. */
            null_ptr = null;
            make_maps_ref(graph, set_list, vtx_elems, new_assign, sub_assign, set1, set2, glob2loc, loc2glob,
                &subnvtxs, &vwgt_max, &vwgt_sum1, &vwgt_sum2);

            for (i = 1; i <= subnvtxs; i++)
            {
                old_sub_assign[i] = sub_assign[i];
            }

            /* Set up goals for this KL invocation. */
            ratio = (vwgt_sum1 + vwgt_sum2) / (goal[set1] + goal[set2]);
            subgoal[0] = ratio * goal[set1];
            subgoal[1] = ratio * goal[set2];

            if (TERM_PROP)
            {
                make_terms_ref(graph, useEdgeWeights, subnvtxs, loc2glob, set1, set2, new_assign, architecture,
                    mesh_dims, term_wgts);
            }

            /* New_assign has overwritten set2 with set1. */
            make_subgraph(graph, subgraph, subnvtxs, &subnedges, new_assign, set1, glob2loc, loc2glob,
                degrees, useEdgeWeights);

            maxdeg = find_maxdeg(subgraph, subnvtxs, useEdgeWeights, (float*) null);

            count_weights(subgraph, subnvtxs, sub_assign, 2, weights, (vwgt_max != 1));

            max_dev = vwgt_max;
            ratio = (subgoal[0] + subgoal[1]) * KL_IMBALANCE / 2;
            if (ratio > max_dev)
            {
                max_dev = (int) ratio;
            }

            klspiff(subgraph, subnvtxs, sub_assign, 2, hops, subgoal, term_wgts, max_dev, maxdeg, useEdgeWeights,
                &null_ptr, weights);

            /* Figure out which modification leaves most vertices intact. */
            nsame = 0;
            for (i = 1; i <= subnvtxs; i++)
            {
                if (old_sub_assign[i] == sub_assign[i])
                {
                    nsame++;
                }
            }

            if (2 * nsame > subnvtxs)
            {
                setA = set1;
                setB = set2;
            }
            else
            {
                setA = set2;
                setB = set1;
            }

            /* Now update the assignments. */
            sizes[setA] = sizes[setB] = 0;
            for (i = 1; i <= subnvtxs; i++)
            {
                vtx = loc2glob[i];
                /* Update the set_lists. */
                ptr = &(vtx_elems[vtx]);
                if (ptr->next != null)
                {
                    ptr->next->prev = ptr->prev;
                }

                if (ptr->prev != null)
                {
                    ptr->prev->next = ptr->next;
                }

                if (sub_assign[i] == 0)
                {
                    new_assign[vtx] = setA;
                    ++sizes[setA];
                    ptr->next = set_list[setA].next;
                    if (ptr->next != null)
                    {
                        ptr->next->prev = ptr;
                    }

                    ptr->prev = &(set_list[setA]);
                    set_list[setA].next = ptr;
                }
                else
                {
                    new_assign[vtx] = setB;
                    ++sizes[setB];
                    ptr->next = set_list[setB].next;
                    if (ptr->next != null)
                    {
                        ptr->next->prev = ptr;
                    }

                    ptr->prev = &(set_list[setB]);
                    set_list[setB].next = ptr;
                }
            }

            remake_graph(subgraph, subnvtxs, loc2glob, degrees, useEdgeWeights);

            return (nsame != subnvtxs);
        }

/* Set up data structures for refine_part. */
        private static void make_maps_ref(vtx_data** graph, /* graph data structure */
            bilist* set_list, /* lists of vertices in each set */
            bilist* vtx_elems, /* start of storage for vertices */
            int* assignment, /* set assignments for graph */
            int* sub_assign, /* assignment file for subgraph */
            int set1, int set2, /* set value denoting subgraph */
            int* glob2loc, /* graph -> subgraph numbering map */
            int* loc2glob, /* subgraph -> graph numbering map */
            int* psub_nvtxs, /* number of vtxs in subgraph */
            int* pvwgt_max, /* returned largest vwgt */
            int* pvwgt_sum1, int* pvwgt_sum2 /* returned set sizes */
        )
        {
            bilist* ptr; /* loops through set lists */
            int vwgt_max; /* largest vertex weight in subgraph */
            int vwgt_sum1, vwgt_sum2; /* sum of vertex weights in sets */
            int vtx; /* vertex in subgraph */
            int size; /* array spacing */
            int j; /* loop counter */

            size = (int) (&(vtx_elems[1]) - &(vtx_elems[0]));
            j = 1;
            vwgt_max = vwgt_sum1 = vwgt_sum2 = 0;
            for (ptr = set_list[set1].next; ptr != null; ptr = ptr->next)
            {
                vtx = ((int) (ptr - vtx_elems)) / size;
                sub_assign[j] = 0;
                glob2loc[vtx] = j;
                loc2glob[j] = vtx;
                if (graph[vtx]->vwgt > vwgt_max)
                {
                    vwgt_max = graph[vtx]->vwgt;
                }

                vwgt_sum1 += graph[vtx]->vwgt;
                j++;
            }

            for (ptr = set_list[set2].next; ptr != null; ptr = ptr->next)
            {
                vtx = ((int) (ptr - vtx_elems)) / size;
                sub_assign[j] = 1;
                glob2loc[vtx] = j;
                loc2glob[j] = vtx;
                if (graph[vtx]->vwgt > vwgt_max)
                {
                    vwgt_max = graph[vtx]->vwgt;
                }

                vwgt_sum2 += graph[vtx]->vwgt;
                assignment[vtx] = set1;
                j++;
            }

            *pvwgt_sum1 = vwgt_sum1;
            *pvwgt_sum2 = vwgt_sum2;
            *pvwgt_max = vwgt_max;
            *psub_nvtxs = j - 1;
        }

/* Compute the terminal constraints for next partition. */
        private static void make_terms_ref(vtx_data** graph, /* data structure for graph */
            bool useEdgeWeights, /* are edge weights being used? */
            int subnvtxs, /* number of vtxs in subgraph */
            int* loc2glob, /* mapping from subgraph to graph */
            int set0, int set1, /* two processors I'm choosing between */
            int* assignment, /* set for each vertex */
            int architecture, /* 0 => hypercube, 1 => mesh */
            int[] mesh_dims /*[3]*/, /* if mesh, size of mesh */
            float*[] term_wgts /* terminal weights for each vertex */
        )
        {
            double term_wgt; /* terminal weight */
            float edge_wgt; /* weight of an edge */
            int dist0 = 0, dist1 = 0; /* distance from set to set0 and set1 */
            int set; /* set neighbor vtx belongs to */
            int vtx; /* vertex number */
            int neighbor; /* neighboring vertex number */
            int x; /* bitwise difference between sets */
            int i, j; /* loop counters */

            /* NOTE: CURRENTLY ONLY WORKING FOR BISECTION. */

            edge_wgt = 1;
            for (i = 1; i <= subnvtxs; i++)
            {
                term_wgt = 0;
                vtx = loc2glob[i];

                for (j = 1; j < graph[vtx]->nedges; j++)
                {
                    neighbor = graph[vtx]->edges[j];
                    set = assignment[neighbor];
                    if (set != set0 && set != set1)
                    {
                        if (architecture == 0)
                        {
                            dist0 = 0;
                            x = set ^ set0;
                            while (x != 0)
                            {
                                if ((x & 1) != 0)
                                {
                                    ++dist0;
                                }

                                x >>= 1;
                            }

                            dist1 = 0;
                            x = set ^ set1;
                            while (x != 0)
                            {
                                if ((x & 1) != 0)
                                {
                                    ++dist1;
                                }

                                x >>= 1;
                            }
                        }

                        else if (architecture > 0)
                        {
                            dist0 = Math.Abs((set % mesh_dims[0]) - (set0 % mesh_dims[0]));
                            dist0 +=
                                Math.Abs(((set / mesh_dims[0]) % mesh_dims[1]) - ((set0 / mesh_dims[0]) % mesh_dims[1]));
                            dist0 +=
                                Math.Abs((set / (mesh_dims[0] * mesh_dims[1])) - (set0 / (mesh_dims[0] * mesh_dims[1])));

                            dist1 = Math.Abs((set % mesh_dims[0]) - (set1 % mesh_dims[0]));
                            dist1 +=
                                Math.Abs(((set / mesh_dims[0]) % mesh_dims[1]) - ((set1 / mesh_dims[0]) % mesh_dims[1]));
                            dist1 +=
                                Math.Abs((set / (mesh_dims[0] * mesh_dims[1])) - (set1 / (mesh_dims[0] * mesh_dims[1])));
                        }

                        if (useEdgeWeights)
                        {
                            edge_wgt = graph[vtx]->ewgts[j];
                        }

                        term_wgt += edge_wgt * (dist0 - dist1);
                    }
                }

                (term_wgts[1])[i] = (float) term_wgt;
            }
        }
    }
}
