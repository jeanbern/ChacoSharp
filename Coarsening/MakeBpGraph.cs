﻿using System;
using System.Diagnostics;
using System.IO;
using System.Runtime.InteropServices;

namespace ChacoSharp.Coarsening
{
    public static unsafe class MakeBpGraph
    {
        /* Make a bipartite graph from vertex separator and neighbors. */

        public static void make_bpgraph(vtx_data** graph, /* list of graph info for each vertex */
            int* sets, /* local partitioning of vtxs */
            int* bndy_list, /* list of vertices on boundary (0 ends) */
            int sep_size, /* length of bndy_list */
            int set_match, /* side to match against */
            int** ppointers, /* start/stop of adjacency lists */
            int** pindices, /* adjacency list for each vertex */
            int** pvweight, /* weight of each vertex */
            int** ploc2glob, /* maps bp number to full graph */
            int* pnleft, /* number of nodes in left half */
            int* pnright, /* number of nodes in right half */
            bool using_vwgts /* are vertices weighted? */
        )
        {
            int* loc2glob = null; /* maps bp number to full graph */
            int* pointers = null; /* start/stop of adjacency lists */
            int* indices = null; /* adjacency list for each vertex */
            int* vwgts = null; /* saves vwgts so they can be overwritten */
            int nleft, nright; /* # vertices in halves of bipartite graph */
            int nedges; /* number of edges in bipartite graph */
            int vtx; /* vertex in graph */
            int neighbor; /* neighbor of vertex */
            int i, j, k; /* loop counters */

            /* First count everything that needs to be counted. */
            nleft = sep_size;
            nright = 0;
            nedges = 0;
            for (i = 0; i < sep_size; i++)
            {
                vtx = bndy_list[i];
                for (j = 1; j < graph[vtx]->nedges; j++)
                {
                    neighbor = graph[vtx]->edges[j];
                    if (sets[neighbor] == set_match)
                    {
                        ++nedges;
                        if (graph[neighbor]->edges[0] > 0)
                        {
                            /* Not yet seen */
                            ++nright;
                            /* Flag him as seen already. */
                            graph[neighbor]->edges[0] = -1;
                        }
                    }
                }
            }

            pointers = (int*) Marshal.AllocHGlobal((nleft + nright + 1) * sizeof(int));
            indices = (int*) Marshal.AllocHGlobal((2 * nedges + 1) * sizeof(int));

            /* Now set up data structures to make construction easier */
            loc2glob = (int*) Marshal.AllocHGlobal((nleft + nright) * sizeof(int));

            if (!using_vwgts)
            {
                for (i = 0; i < nleft; i++)
                {
                    vtx = bndy_list[i];
                    loc2glob[i] = vtx;
                    graph[vtx]->edges[0] = i;
                }

                k = nleft;
                for (i = 0; i < nleft; i++)
                {
                    vtx = bndy_list[i];
                    for (j = 1; j < graph[vtx]->nedges; j++)
                    {
                        neighbor = graph[vtx]->edges[j];
                        if (sets[neighbor] == set_match)
                        {
                            if (graph[neighbor]->edges[0] == -1)
                            {
                                loc2glob[k] = neighbor;
                                /* Reflag him as seen already with glob2loc value. */
                                graph[neighbor]->edges[0] = k;
                                k++;
                            }
                        }
                    }
                }
            }
            else
            {
                vwgts = (int*) Marshal.AllocHGlobal((nleft + nright) * sizeof(int));

                for (i = 0; i < nleft; i++)
                {
                    vtx = bndy_list[i];
                    loc2glob[i] = vtx;
                    vwgts[i] = graph[vtx]->vwgt;
                    /* Use edges[0] as a seen flag and as a glob2loc value. */
                    graph[vtx]->edges[0] = i;
                }

                k = nleft;
                for (i = 0; i < nleft; i++)
                {
                    vtx = bndy_list[i];
                    for (j = 1; j < graph[vtx]->nedges; j++)
                    {
                        neighbor = graph[vtx]->edges[j];
                        if (sets[neighbor] == set_match)
                        {
                            if (graph[neighbor]->edges[0] == -1)
                            {
                                /* First occurrence. */
                                loc2glob[k] = neighbor;
                                vwgts[k] = graph[neighbor]->vwgt;
                                /* Use edges[0] as a seen flag and as a glob2loc value. */
                                graph[neighbor]->edges[0] = k;
                                k++;
                            }
                        }
                    }
                }
            }

            /* I can now construct graph directly */
            nedges = 0;
            pointers[0] = 0;
            for (i = 0; i < nleft; i++)
            {
                vtx = loc2glob[i];
                for (j = 1; j < graph[vtx]->nedges; j++)
                {
                    neighbor = graph[vtx]->edges[j];
                    if (sets[neighbor] == set_match)
                    {
                        indices[nedges++] = graph[neighbor]->edges[0];
                    }
                }

                pointers[i + 1] = nedges;
            }

            for (i = nleft; i < nleft + nright; i++)
            {
                vtx = loc2glob[i];
                for (j = 1; j < graph[vtx]->nedges; j++)
                {
                    neighbor = graph[vtx]->edges[j];
                    if (sets[neighbor] == 2)
                    {
                        indices[nedges++] = graph[neighbor]->edges[0];
                    }
                }

                pointers[i + 1] = nedges;
            }

            /* Now restore the edges[0] values. */
            for (i = 0; i < nleft + nright; i++)
            {
                graph[loc2glob[i]]->edges[0] = loc2glob[i];
            }

            /*
            check_bpgraph(nleft, nright, pointers, indices);
            */

            if (using_vwgts)
            {
                *pvweight = vwgts;
            }
            else
            {
                *pvweight = null;
            }

            *ploc2glob = loc2glob;
            *ppointers = pointers;
            *pindices = indices;
            *pnleft = nleft;
            *pnright = nright;
        }

        public static void check_bpgraph(int n_left, int n_right, int* pointers, int* indices)
        {
            int i, j, k, neighbor;

            for (i = 0; i < n_left; i++)
            {
                for (j = pointers[i]; j < pointers[i + 1]; j++)
                {
                    neighbor = indices[j];
                    if (neighbor < n_left || neighbor >= n_left + n_right)
                    {
                        Trace.WriteLine($"Bad edge ({i:d}, {neighbor:d})");
                    }

                    /* Check for counter-edge */
                    for (k = pointers[neighbor]; k < pointers[neighbor + 1]; k++)
                    {
                        if (indices[k] == i)
                        {
                            break;
                        }
                    }

                    if (k == pointers[neighbor + 1])
                    {
                        Trace.WriteLine($"Flip edge ({k:d}, {i:d}) not found");
                    }
                }
            }

            for (i = n_left; i < n_left + n_right; i++)
            {
                for (j = pointers[i]; j < pointers[i + 1]; j++)
                {
                    neighbor = indices[j];
                    if (neighbor < 0 || neighbor >= n_left)
                    {
                        Trace.WriteLine($"Bad edge ({i:d}, {neighbor:d})");
                    }

                    /* Check for counter-edge */
                    for (k = pointers[neighbor]; k < pointers[neighbor + 1]; k++)
                    {
                        if (indices[k] == i)
                        {
                            break;
                        }
                    }

                    if (k == pointers[neighbor + 1])
                    {
                        Console.WriteLine($"Flip edge ({k:d}, {i:d}) not found");
                    }
                }
            }
        }

        public static void print_bpgraph(int nleft, int nright, int* pointers, int* indices, int* vwgts)
        {
            int i, j, nedges, nvtxs;

            nvtxs = nleft + nright;
            nedges = (pointers[nvtxs] - pointers[0]) / 2;

            using (var file = File.OpenWrite("BPGRAPH"))
            using (var stream = new StreamWriter(file))
            {
                stream.WriteLine("{0:d} {1:d}", nvtxs, nedges);

                for (i = 0; i < nvtxs; i++)
                {
                    if (vwgts != null)
                    {
                        stream.WriteLine("{0:d}     ", vwgts[i]);
                    }

                    for (j = pointers[i]; j < pointers[i + 1]; j++)
                    {
                        stream.WriteLine("{0:d} ", indices[j]);
                    }

                    stream.WriteLine();
                }
            }
        }
    }
}
