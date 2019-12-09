#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
#pragma warning disable HAA0101 // Array allocation for params parameter

using System;
using System.Diagnostics;
using System.Runtime.InteropServices;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Utilities.Timer;
using static ChacoSharp.Coarsening.NWayKlv;
using static ChacoSharp.Coarsening.MakeFGraph;

namespace ChacoSharp.Coarsening
{
    public static unsafe class KlvSpiff
    {
        /* Improve (weighted) vertex separator.  Two sets are 0/1; separator = 2. */


        public static void klvspiff(vtx_data** graph, /* list of graph info for each vertex */
            int nvtxs, /* number of vertices in graph */
            int* sets, /* local partitioning of vtxs */
            double[] goal, /* desired set sizes */
            int max_dev, /* largest deviation from balance allowed */
            int** bndy_list, /* list of vertices on boundary (0 ends) */
            double[] weights /* vertex weights in each set */
        )
        {
            bilist** lbuckets; /* space for bucket sorts for left moves */
            bilist** rbuckets; /* space for bucket sorts for right moves */
            bilist* llistspace; /* space for all left bidirectional elements */
            bilist* rlistspace; /* space for all right bidirectional elements */
            int* ldvals; /* change in penalty for each possible move */
            int* rdvals; /* change in penalty for each possible move */
            int* edges; /* loops through neighbor lists */
            double time, time1; /* timing variables */
            int dval; /* largest transition cost for a vertex */
            int maxdval; /* largest transition cost for all vertices */
            bool error; /* out of space? */
            int i, j; /* loop counters */

            time = seconds();

            if (DEBUG_TRACE)
            {
                Trace.WriteLine($"<Entering klvspiff, {nameof(nvtxs)} = {nvtxs:d}>");
            }

            /* Find largest possible change. */
            maxdval = 0;
            for (i = 1; i <= nvtxs; i++)
            {
                if (graph[i]->vwgt > maxdval)
                {
                    maxdval = graph[i]->vwgt;
                }

                dval = -graph[i]->vwgt;
                edges = graph[i]->edges;
                for (j = graph[i]->nedges - 1; j != 0; j--)
                {
                    dval += graph[*(++edges)]->vwgt;
                }

                if (dval > maxdval)
                {
                    maxdval = dval;
                }
            }

            /* Allocate a bunch of space for KLV. */
            time1 = seconds();
            error =
                klv_init(&lbuckets, &rbuckets, &llistspace, &rlistspace, &ldvals, &rdvals, nvtxs, maxdval);
            kl_init_time += seconds() - time1;

            if (!error)
            {
                if (DEBUG_KL != DebugFlagKL.NoDebugging)
                {
                    Trace.WriteLine(" Before KLV: ");
                    countup_vtx_sep(graph, nvtxs, sets);
                }

                time1 = seconds();
                error = nway_klv(graph, nvtxs, lbuckets, rbuckets, llistspace, rlistspace, ldvals, rdvals, sets,
                    maxdval, goal, max_dev, bndy_list, weights);
                nway_kl_time += seconds() - time1;

                if (DEBUG_KL == DebugFlagKL.MoreInfo || DEBUG_KL == DebugFlagKL.PrintBucket)
                {
                    Trace.WriteLine(" After KLV: ");
                    countup_vtx_sep(graph, nvtxs, sets);
                }
            }

            if (error)
            {
                Trace.WriteLine($"\nWARNING: No space to perform KLV on graph with {nvtxs:d} vertices.");
                Trace.WriteLine("         NO LOCAL REFINEMENT PERFORMED.\n");
            }

            free_klv(lbuckets, rbuckets, llistspace, rlistspace, ldvals, rdvals);

            kl_total_time += seconds() - time;
        }

        public static void countup_vtx_sep(vtx_data** graph, /* list of graph info for each vertex */
            int nvtxs, /* number of vertices in graph */
            int* sets /* local partitioning of vtxs */
        )
        {
            int vtx, set; /* vertex and set in graph */
            int sep_size; /* size of the separator */
            int i, j, k; /* loop counters */

            sep_size = 0;
            j = k = 0;
            for (i = 1; i <= nvtxs; i++)
            {
                if (sets[i] == 0)
                {
                    j += graph[i]->vwgt;
                }

                if (sets[i] == 1)
                {
                    k += graph[i]->vwgt;
                }

                if (sets[i] == 2)
                {
                    sep_size += graph[i]->vwgt;
                }
            }

            Trace.WriteLine($"Set sizes = {j:d}/{k:d}, Separator size = {sep_size:d}\n");

            /* Now check that it really is a separator. */
            for (i = 1; i <= nvtxs; i++)
            {
                set = sets[i];
                if (set != 2)
                {
                    for (j = 1; j < graph[i]->nedges; j++)
                    {
                        vtx = graph[i]->edges[j];
                        if (sets[vtx] != 2 && sets[vtx] != set)
                        {
                            Trace.WriteLine($"Error: {i:d} (set {set:d}) adjacent to {vtx:d} (set {sets[vtx]:d})");
                        }
                    }
                }
            }
        }

        public static void free_klv(
            /* Free everything malloc'd for KLV. */
            bilist** lbuckets, /* space for bucket sorts */
            bilist** rbuckets, /* space for bucket sorts */
            bilist* llistspace, /* space for all bidirectional elements */
            bilist* rlistspace, /* space for all bidirectional elements */
            int* ldvals, /* change in penalty for each possible move */
            int* rdvals /* change in penalty for each possible move */
        )
        {
            Marshal.FreeHGlobal((IntPtr) rlistspace);
            Marshal.FreeHGlobal((IntPtr) llistspace);
            Marshal.FreeHGlobal((IntPtr) rdvals);
            Marshal.FreeHGlobal((IntPtr) ldvals);
            Marshal.FreeHGlobal((IntPtr) rbuckets);
            Marshal.FreeHGlobal((IntPtr) lbuckets);
        }

        public static bool klv_init(bilist*** lbucket_ptr, /* space for left bucket sorts */
            bilist*** rbucket_ptr, /* space for right bucket sorts */
            bilist** llistspace, /* space for elements of linked lists */
            bilist** rlistspace, /* space for elements of linked lists */
            int** ldvals, /* change in separator for left moves */
            int** rdvals, /* change in separator for right moves */
            int nvtxs, /* number of vertices in the graph */
            int maxchange /* maximum change by moving a vertex */
        )
        {
            int sizeb; /* size of set of buckets */
            int sizel; /* size of set of pointers for all vertices */

            /* Allocate appropriate data structures for buckets, and listspace. */

            sizeb = (2 * maxchange + 1) * sizeof(bilist*);
            *lbucket_ptr = (bilist**) Marshal.AllocHGlobal(sizeb);
            *rbucket_ptr = (bilist**) Marshal.AllocHGlobal(sizeb);

            *ldvals = (int*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));
            *rdvals = (int*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));

            sizel = (nvtxs + 1) * sizeof(bilist);
            *llistspace = (bilist*) Marshal.AllocHGlobal(sizel);
            *rlistspace = (bilist*) Marshal.AllocHGlobal(sizel);

            return *lbucket_ptr == null || *rbucket_ptr == null || *ldvals == null || *rdvals == null ||
                   *llistspace == null || *rlistspace == null;
        }

        /* Find vertices on boundary of partition, and change their assignments. */

        public static int find_bndy(vtx_data** graph, /* array of vtx data for graph */
            int nvtxs, /* number of vertices in graph */
            int* assignment, /* processor each vertex gets assigned to */
            int new_val, /* assignment value for boundary vtxs */
            int** pbndy_list /* returned list, end with zero */
        )
        {
            int* bndy_list; /* returned list, end with zero */
            int* edges; /* loops through edge list */
            int list_length; /* returned number of vtxs on boundary */
            int set, set2; /* set a vertex is in */
            int i, j; /* loop counters */

            bndy_list = (int*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));

            list_length = 0;
            for (i = 1; i <= nvtxs; i++)
            {
                set = assignment[i];
                edges = graph[i]->edges;
                for (j = graph[i]->nedges - 1; j != 0; j--)
                {
                    set2 = assignment[*(++edges)];
                    if (set2 != set)
                    {
                        bndy_list[list_length++] = i;
                        break;
                    }
                }
            }

            bndy_list[list_length] = 0;

            for (i = 0; i < list_length; i++)
            {
                assignment[bndy_list[i]] = new_val;
            }

            /* Shrink out unnecessary space */
            *pbndy_list = (int*) Marshal.ReAllocHGlobal((IntPtr) bndy_list, new IntPtr((list_length + 1) * sizeof(int)));

            return (list_length);
        }

/* Find a vertex separator on one side of an edge separator. */

        public static int find_side_bndy(vtx_data** graph, /* array of vtx data for graph */
            int nvtxs, /* number of vertices in graph */
            int* assignment, /* processor each vertex gets assigned to */
            int side, /* side to take vertices from */
            int new_val, /* assignment value for boundary vtxs */
            int** pbndy_list /* returned list, end with zero */
        )

        {
            int* edges; /* loops through edge list */
            int* bndy_list; /* returned list, end with zero */
            int list_length; /* returned number of vtxs on boundary */
            int set, set2; /* set a vertex is in */
            int i, j; /* loop counters */

            if (*pbndy_list != null)
            {
                /* Contains list of all vertices on boundary. */
                bndy_list = *pbndy_list;
                i = list_length = 0;
                while (bndy_list[i] != 0)
                {
                    if (assignment[bndy_list[i]] == side)
                    {
                        bndy_list[list_length++] = bndy_list[i];
                    }

                    ++i;
                }
            }

            else
            {
                bndy_list = (int*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));

                list_length = 0;
                for (i = 1; i <= nvtxs; i++)
                {
                    set = assignment[i];
                    if (set == side)
                    {
                        edges = graph[i]->edges;
                        for (j = graph[i]->nedges - 1; j != 0; j--)
                        {
                            set2 = assignment[*(++edges)];
                            if (set2 != set)
                            {
                                bndy_list[list_length++] = i;
                                break;
                            }
                        }
                    }
                }
            }

            bndy_list[list_length] = 0;

            for (i = 0; i < list_length; i++)
            {
                assignment[bndy_list[i]] = new_val;
            }

            /* Shrink out unnecessary space */
            *pbndy_list = (int*) Marshal.ReAllocHGlobal((IntPtr) bndy_list, new IntPtr((list_length + 1) * sizeof(int)));

            return (list_length);
        }

        public static bool flatten(vtx_data ** graph,       /* array of vtx data for graph */
            int                nvtxs,       /* number of vertices in graph */
            int                nedges,      /* number of edges in graph */
            vtx_data ***pcgraph,     /* coarsened version of graph */
            int *              pcnvtxs,     /* number of vtxs in coarsened graph */
            int *              pcnedges,    /* number of edges in coarsened graph */
            int **             pv2cv,       /* pointer to v2cv */
            bool                useEdgeWeights, /* are edge weights being used? */
            int                igeom,       /* dimensions of geometric data */
            float **           coords,      /* coordinates for vertices */
            float **           ccoords      /* coordinates for coarsened vertices */
)
{
  double Thresh; /* minimal acceptable size reduction */
  int *  v2cv;   /* map from vtxs to coarse vtxs */
  int    cnvtxs; /* number of vertices in flattened graph */

  Thresh = .9;

  v2cv = (int*)Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));

  find_flat(graph, nvtxs, &cnvtxs, v2cv);

  if (cnvtxs <= Thresh * nvtxs) { /* Sufficient shrinkage? */
    makefgraph(graph, nvtxs, nedges, pcgraph, cnvtxs, pcnedges, v2cv, useEdgeWeights, igeom, coords,
               ccoords);

    *pcnvtxs = cnvtxs;
    *pv2cv   = v2cv;
    return true;
  }

  /* Not worth bothering */
  Marshal.FreeHGlobal((IntPtr)v2cv);
  return false;
}

private static void find_flat(vtx_data **graph,   /* data structure for storing graph */
               int               nvtxs,   /* number of vertices in graph */
               int *             pcnvtxs, /* number of coarse vertices */
               int *             v2cv     /* map from vtxs to coarse vtxs */
)
{
  /* Look for cliques with the same neighbor set.  These are matrix */
  /* rows corresponding to multiple degrees of freedom on a node. */
  /* They can be flattened out, generating a smaller graph. */

  int *scatter;   /* for checking neighbor list identity */
  int *hash;      /* hash value for each vertex */
  int  this_hash; /* particular hash value */
  int  neighbor;  /* neighbor of a vertex */
  int  cnvtxs;    /* number of distinct vertices */
  int  i, j;      /* loop counters */

  hash    = (int*)Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));
  scatter = (int*)Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));

  /* compute hash values */

  for (i = 1; i <= nvtxs; i++) {
    this_hash = i;
    for (j = 1; j < graph[i]->nedges; j++) {
      this_hash += graph[i]->edges[j];
    }
    hash[i] = this_hash;
  }

  for (i = 1; i <= nvtxs; i++) {
    v2cv[i]    = 0;
    scatter[i] = 0;
  }

  /* Find supernodes. */
  cnvtxs = 0;

  for (i = 1; i <= nvtxs; i++) {
    if (v2cv[i] == 0) { /* Not yet flattened. */
      v2cv[i] = ++cnvtxs;
      for (j = 1; j < graph[i]->nedges; j++) {
        neighbor = graph[i]->edges[j];
        if (neighbor > i && hash[neighbor] == hash[i] &&   /* same hash value */
            graph[i]->nedges == graph[neighbor]->nedges && /* same degree */
            v2cv[neighbor] == 0 &&                         /* neighbor not flattened */
            SameStructure(i, neighbor, graph, scatter)) {
          v2cv[neighbor] = cnvtxs;
        }
      }
    }
  }

  *pcnvtxs = cnvtxs;

  Marshal.FreeHGlobal((IntPtr)scatter);
  Marshal.FreeHGlobal((IntPtr)hash);
}

private static bool SameStructure(int node1, int node2,     /* two vertices which might have same nonzeros */
                  vtx_data **graph,  /* data structure for storing graph */
                  int *             scatter /* array for checking vertex labelling */
)
{
  bool same; /* are two vertices indistinguisable? */
  int i;    /* loop counter */

  scatter[node1] = node1;
  for (i = 1; i < graph[node1]->nedges; i++) {
    scatter[graph[node1]->edges[i]] = node1;
  }

  for (i = 1; i < graph[node2]->nedges; i++) {
    if (scatter[graph[node2]->edges[i]] != node1) {
      break;
    }
  }
  same = (i == graph[node2]->nedges && scatter[node2] == node1);

  return (same);
}
    }
}
