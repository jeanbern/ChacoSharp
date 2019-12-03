using System;
using System.Runtime.InteropServices;
using static ChacoSharp.StaticConstants;

namespace ChacoSharp.Graph
{
    public static unsafe class Reformat
    {
        /* Change from a FORTRAN graph style to our graph data structure. */

        public static bool reformat(int* start, /* start of edge list for each vertex */
            int* adjacency, /* edge list data */
            int nvtxs, /* number of vertices in graph */
            int* pnedges, /* ptr to number of edges in graph */
            int* vwgts, /* weights for all vertices */
            float* ewgts, /* weights for all edges */
            vtx_data*** pgraph /* ptr to array of vtx data for graph */
        )
        {
            if (FullTrace)
            {
                Console.WriteLine($"<Entering {nameof(reformat)}>");
            }

            vtx_data** graph = null; /* array of vtx data for graph */
            vtx_data* links = null; /* space for data for all vtxs */
            int* edges = null; /* space for all adjacency lists */
            float* eweights = null; /* space for all edge weights */
            int* eptr = null; /* steps through adjacency list */
            int* eptr_save = null; /* saved index into adjacency list */
            float* wptr = null; /* steps through edge weights list */
            int self_edge; /* number of self loops detected */
            int size; /* length of all edge lists */
            double sum; /* sum of edge weights for a vtx */
            int i, j; /* loop counters */

            var useEdgeWeights = (ewgts != null); /* are edge weights being used? */
            var using_vwgts = (vwgts != null); /* are vertex weights being used? */

            graph = (vtx_data**) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(vtx_data*));
            *pgraph = graph;
            if (graph == null)
            {
                return true;
            }

            graph[1] = null;

            /* Set up all the basic data structure for the vertices. */
            /* Replace many small mallocs by a few large ones. */
            links = (vtx_data*) Marshal.AllocHGlobal((nvtxs) * sizeof(vtx_data));
            if (links == null)
            {
                return true;
            }

            for (i = 1; i <= nvtxs; i++)
            {
                graph[i] = links++;
            }

            graph[1]->edges = null;
            graph[1]->ewgts = null;

            /* Now fill in all the data fields. */
            if (start != null)
            {
                *pnedges = start[nvtxs] / 2;
            }
            else
            {
                *pnedges = 0;
            }

            size = 2 * (*pnedges) + nvtxs;
            edges = (int*) Marshal.AllocHGlobal(size * sizeof(int));
            if (edges == null)
            {
                return true;
            }

            if (useEdgeWeights)
            {
                eweights = (float*) Marshal.AllocHGlobal(size * sizeof(float));
                if (eweights == null)
                {
                    return true;
                }
            }

            if (start != null)
            {
                eptr = adjacency + start[0];
                wptr = ewgts;
            }

            self_edge = 0;

            for (i = 1; i <= nvtxs; i++)
            {
                if (using_vwgts)
                {
                    graph[i]->vwgt = *(vwgts++);
                }
                else
                {
                    graph[i]->vwgt = 1;
                }

                if (start != null)
                {
                    size = start[i] - start[i - 1];
                }
                else
                {
                    size = 0;
                }

                graph[i]->nedges = size + 1;
                graph[i]->edges = edges;
                *edges++ = i;
                eptr_save = eptr;
                for (j = size; j != 0; j--)
                {
                    if (*eptr != i)
                    {
                        *edges++ = *eptr++;
                    }
                    else
                    {
                        /* Self edge, skip it. */
                        //if (self_edge == 0)
                        {
                            Console.WriteLine("WARNING: Self edge ({0:d},{1:d}) being ignored", i, i);
                        }

                        ++self_edge;
                        eptr++;
                        --(graph[i]->nedges);
                        --(*pnedges);
                    }
                }

                if (useEdgeWeights)
                {
                    graph[i]->ewgts = eweights;
                    eweights++;
                    sum = 0;
                    for (j = size; j != 0; j--)
                    {
                        if (*eptr_save++ != i)
                        {
                            sum += *wptr;
                            *eweights++ = *wptr++;
                        }
                        else
                        {
                            wptr++;
                        }
                    }

                    graph[i]->ewgts[0] = (float) -sum;
                }
                else
                {
                    graph[i]->ewgts = null;
                }
            }

            if (self_edge > 1)
            {
                Console.WriteLine("WARNING: {0:d} self edges were detected and ignored", self_edge);
            }

            return false;
        }
    }
}
