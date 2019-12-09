using System;
using System.Diagnostics;
using System.Runtime.InteropServices;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Coarsening.MakeCoords;
using static ChacoSharp.Utilities.Timer;

namespace ChacoSharp.Coarsening
{
    public static unsafe class MakeFGraph
    {
        public static void makefgraph(vtx_data** graph, /* array of vtx data for graph */
            int nvtxs, /* number of vertices in graph */
            int nedges, /* number of edges in graph */
            vtx_data*** pcgraph, /* coarsened version of graph */
            int cnvtxs, /* number of vtxs in coarsened graph */
            int* pcnedges, /* number of edges in coarsened graph */
            int* v2cv, /* mapping from vtxs to coarsened vtxs */
            bool useEdgeWeights, /* are edge weights being used? */
            int igeom, /* dimensions of geometric data */
            float** coords, /* coordinates for vertices */
            float** ccoords /* coordinates for coarsened vertices */
        )
        {
            vtx_data** cgraph = null; /* coarsened version of graph */
            vtx_data* links = null; /* space for all the vertex data */
            vtx_data** gptr = null; /* loops through cgraph */
            vtx_data* cgptr = null; /* loops through cgraph */
            int* iptr = null; /* loops through integer arrays */
            int* seenflag = null; /* flags for vtxs already put in edge list */
            int* sptr = null; /* loops through seenflags */
            int* cv2v_vals = null; /* vtxs corresponding to each cvtx */
            int* cv2v_ptrs = null; /* indices into cv2v_vals */
            float* eweights = null; /* space for edge weights in coarsened graph */
            float* ewptr = null; /* loops through eweights */
            float* fptr = null; /* loops through eweights */
            float ewgt; /* edge weight */
            double ewgt_sum; /* sum of edge weights */
            double time; /* timing parameters */
            int nseen; /* number of edges of coarse graph seen so far */
            int vtx; /* vertex in original graph */
            int cvtx; /* vertex in coarse graph */
            int cnedges; /* twice number of edges in coarsened graph */
            int neighbor; /* neighboring vertex */
            int size; /* space needed for coarsened graph */
            int* edges = null; /* space for edges in coarsened graph */
            int* eptr = null; /* loops through edges data structure */
            int cneighbor; /* neighboring vertex number in coarsened graph */
            int i, j; /* loop counters */

            /* Compute the number of vertices and edges in the coarsened graph, */
            /* and construct start pointers into coarsened edge array. */
            time = seconds();

            /* Construct mapping from original graph vtxs to coarsened graph vtxs. */
            cv2v_vals = (int*) Marshal.AllocHGlobal(nvtxs * sizeof(int));
            cv2v_ptrs = (int*) Marshal.AllocHGlobal((cnvtxs + 2) * sizeof(int));
            makecv2v(nvtxs, cnvtxs, v2cv, cv2v_vals, cv2v_ptrs);

            /* Compute an upper bound on the number of coarse graph edges. */
            cnedges = nedges - (nvtxs - cnvtxs);

            /* Now allocate space for the new graph.  Overallocate and realloc later. */
            *pcgraph = cgraph = (vtx_data**) Marshal.AllocHGlobal((cnvtxs + 1) * sizeof(vtx_data*));
            links = (vtx_data*) Marshal.AllocHGlobal(cnvtxs * sizeof(vtx_data));

            size = 2 * cnedges + cnvtxs;
            edges = (int*) Marshal.AllocHGlobal(size * sizeof(int));
            if (COARSEN_EWGTS)
            {
                ewptr = eweights = (float*) Marshal.AllocHGlobal(size * sizeof(float));
            }

            /* Zero all the seen flags. */
            seenflag = (int*) Marshal.AllocHGlobal((cnvtxs + 1) * sizeof(int));
            sptr = seenflag;
            for (i = cnvtxs; i != 0; i--)
            {
                *(++sptr) = 0;
            }

            /* Use the renumbering to fill in the edge lists for the new graph. */
            cnedges = 0;
            eptr = edges;
            ewgt = 1;

            sptr = cv2v_vals;
            for (cvtx = 1; cvtx <= cnvtxs; cvtx++)
            {
                nseen = 1;

                cgptr = cgraph[cvtx] = links++;

                if (COARSEN_VWGTS)
                {
                    cgptr->vwgt = 0;
                }
                else
                {
                    cgptr->vwgt = 1;
                }

                eptr[0] = cvtx;
                cgptr->edges = eptr;
                if (COARSEN_EWGTS)
                {
                    cgptr->ewgts = ewptr;
                }
                else
                {
                    cgptr->ewgts = null;
                }

                ewgt_sum = 0;
                for (i = cv2v_ptrs[cvtx + 1] - cv2v_ptrs[cvtx]; i != 0; i--)
                {
                    vtx = *sptr++;

                    iptr = graph[vtx]->edges;
                    if (useEdgeWeights)
                    {
                        fptr = graph[vtx]->ewgts;
                    }

                    for (j = graph[vtx]->nedges - 1; j != 0; j--)
                    {
                        neighbor = *(++iptr);
                        cneighbor = v2cv[neighbor];
                        if (cneighbor != cvtx)
                        {
                            if (useEdgeWeights)
                            {
                                ewgt = *(++fptr);
                            }

                            ewgt_sum += ewgt;

                            /* Seenflags being used as map from cvtx to index. */
                            if (seenflag[cneighbor] == 0)
                            {
                                /* New neighbor. */
                                cgptr->edges[nseen] = cneighbor;
                                if (COARSEN_EWGTS)
                                {
                                    cgptr->ewgts[nseen] = ewgt;
                                }

                                seenflag[cneighbor] = nseen++;
                            }
                            else
                            {
                                /* Already seen neighbor. */
                                if (COARSEN_EWGTS)
                                {
                                    cgptr->ewgts[seenflag[cneighbor]] += ewgt;
                                }
                            }
                        }
                        else if (useEdgeWeights)
                        {
                            ++fptr;
                        }
                    }
                }

                /* Now clear the seenflag values. */
                iptr = cgptr->edges;
                for (j = nseen - 1; j != 0; j--)
                {
                    seenflag[*(++iptr)] = 0;
                }

                if (COARSEN_EWGTS)
                {
                    cgptr->ewgts[0] = (float) -ewgt_sum;
                }

                /* Increment pointers into edges list. */
                cgptr->nedges = nseen;
                eptr += nseen;
                if (COARSEN_EWGTS)
                {
                    ewptr += nseen;
                }

                cnedges += nseen - 1;
            }

            Marshal.FreeHGlobal((IntPtr) seenflag);

            /* Form new vertex weights by adding those from contracted edges. */
            if (COARSEN_VWGTS)
            {
                gptr = graph;
                for (i = 1; i <= nvtxs; i++)
                {
                    cgraph[v2cv[i]]->vwgt += (*(++gptr))->vwgt;
                }
            }

            /* Reduce arrays to actual sizes */
            cnedges /= 2;
            size = 2 * cnedges + cnvtxs;
            eptr = edges;
            edges = (int*) Marshal.ReAllocHGlobal((IntPtr) edges, (IntPtr) (size * sizeof(int)));
            if (eptr != edges)
            {
                /* Need to reset pointers in graph. */
                for (i = 1; i <= cnvtxs; i++)
                {
                    cgraph[i]->edges = edges;
                    edges += cgraph[i]->nedges;
                }
            }

            if (COARSEN_EWGTS)
            {
                ewptr = eweights;
                eweights = (float*) Marshal.ReAllocHGlobal((IntPtr) eweights, (IntPtr) (size * sizeof(float)));
                if (ewptr != eweights)
                {
                    /* Need to reset pointers in graph. */
                    for (i = 1; i <= cnvtxs; i++)
                    {
                        cgraph[i]->ewgts = eweights;
                        eweights += cgraph[i]->nedges;
                    }
                }
            }

            /* If desired, make new vtx coordinates = center-of-mass of their parents. */
            if (coords != null && ccoords != null && igeom > 0)
            {
                makeccoords(graph, cnvtxs, cv2v_ptrs, cv2v_vals, igeom, coords, ccoords);
            }

            *pcnedges = cnedges;

            Marshal.FreeHGlobal((IntPtr) cv2v_ptrs);
            Marshal.FreeHGlobal((IntPtr) cv2v_vals);

            if (DEBUG_COARSEN)
            {
                Trace.WriteLine($" Coarse graph has {cnvtxs:D} vertices and {cnedges:D} edges");
            }

            make_cgraph_time += seconds() - time;
        }

        static void makecv2v(int nvtxs, /* number of vertices in graph */
            int cnvtxs, /* number of vtxs in coarsened graph */
            int* v2cv, /* mapping from vtxs to coarsened vtxs */
            int* cv2v_vals, /* vtxs corresponding to each cvtx */
            int* cv2v_ptrs /* indices into cv2c_vals */
        )

        {
            int sum; /* cumulative offests into vals array */
            int i; /* loop counter */

            /* First find number of vtxs associated with each coarse graph vtx. */

            for (i = 1; i <= cnvtxs + 1; i++)
            {
                cv2v_ptrs[i] = 0;
            }

            for (i = 1; i <= nvtxs; i++)
            {
                ++cv2v_ptrs[v2cv[i] + 1]; /* +1 offsets and simplifies next loop. */
            }

            cv2v_ptrs[1] = 0;

            /* Now make this a cumulative total to index into cv2v_vals. */
            sum = 0;
            for (i = 2; i <= cnvtxs + 1; i++)
            {
                cv2v_ptrs[i] += sum;
                sum = cv2v_ptrs[i];
            }

            /* Now go ahead and set the cv2v_vals. */
            for (i = 1; i <= nvtxs; i++)
            {
                cv2v_vals[cv2v_ptrs[v2cv[i]]] = i;
                ++cv2v_ptrs[v2cv[i]];
            }

            /* Finally, reset the cv2v_ptrs values. */
            for (i = cnvtxs; i != 0; i--)
            {
                cv2v_ptrs[i] = cv2v_ptrs[i - 1];
            }

            cv2v_ptrs[1] = 0;
        }
    }
}
