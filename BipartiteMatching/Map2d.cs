#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
using System;
using System.Diagnostics;
using System.Runtime.InteropServices;
using ChacoSharp.Utilities;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.BipartiteMatching.CheckBp;
using static ChacoSharp.BipartiteMatching.FindIndex;
using static ChacoSharp.BipartiteMatching.MoveVtxs;

namespace ChacoSharp.BipartiteMatching
{
    public static unsafe class Map2d
    {
        private const int nsection = 2;
        private const int nlists = 4; /* number of lists to generate */
        private const int nsets = 4;

        public static void map2d(vtx_data** graph, /* data structure with vertex weights */
            double*[] xvecs, /* vectors to partition */
            int vertexCount, /* number of vertices */
            int* sets, /* set each vertex gets assigned to */
            double[] goal, /* desired set sizes */
            int maxVertexWeight /* largest vertex weight */
        )
        {
            double*[][] vals = new double*[4][]; /* values in sorted lists */
            for (var i = 0; i < vals.Length; i++)
            {
                vals[i] = new double*[MAXSETS];
            }

            double[] dist = new double[4]; /* trial separation point */
            double[] size = new double[4]; /* sizes of each set being modified */
            int*[][] indices = new int*[4][]; /* indices sorting lists */
            for (var i = 0; i < indices.Length; i++)
            {
                indices[i] = new int*[MAXSETS];
            }

            int[][] startvtx = new int[4][]; /* indices defining separation */
            for (var i = 0; i < startvtx.Length; i++)
            {
                startvtx[i] = new int[MAXSETS];
            }

            N_VTX_CHECKS = N_VTX_MOVES = 0;

            /* Generate all the lists of values that need to be sorted. */
            genvals2d(xvecs, vals, vertexCount);

            /* Sort the lists of values. */
            sorts2d(vals, indices, vertexCount);

            /* Now initialize dists and assign to sets. */
            inits2d(graph, xvecs, vals, indices, vertexCount, dist, startvtx, size, sets);

            /* Determine the largest and smallest allowed set sizes. */
            /* (For now, assume all sets must be same size, but can easily change.) */

            if (DEBUG_BPMATCH == DebugFlagBP.ErrorChecking)
            {
                Trace.WriteLine(" Calling check before movevtxs");
                checkbp(graph, xvecs, sets, dist, vertexCount, nsection);
            }

            movevtxs(graph, vertexCount, nsets, dist, indices, vals, startvtx, sets, size, goal, maxVertexWeight);

            if (DEBUG_BPMATCH != DebugFlagBP.NoDebugging)
            {
                Trace.WriteLine($" {nameof(N_VTX_CHECKS)} = {N_VTX_CHECKS:d}, {nameof(N_VTX_MOVES)} = {N_VTX_MOVES:d}");
                checkbp(graph, xvecs, sets, dist, vertexCount, nsection);
            }

            free2d(vals, indices);
        }

        private static void genvals2d(
            /* Create lists of sets of values to be sorted. */
            double*[] xvecs, /* vectors to partition */
            double*[][] vals /*[4][MAXSETS]*/, /* ptrs to lists of values */
            int nvtxs /* number of values */
        )
        {
            double*[] temp = new double*[nlists]; /* place holders for vals */
            int i; /* loop counter */

            for (i = 0; i < nlists; i++)
            {
                temp[i] = (double*) Marshal.AllocHGlobal(nvtxs * sizeof(double));
            }

            for (i = 1; i <= nvtxs; i++)
            {
                temp[0][i - 1] = 4 * xvecs[1][i];
                temp[1][i - 1] = 4 * xvecs[2][i];
                temp[2][i - 1] = 4 * (xvecs[1][i] + xvecs[2][i]);
                temp[3][i - 1] = 4 * (xvecs[2][i] - xvecs[1][i]);
            }

            vals[0][1] = vals[1][0] = vals[2][3] = vals[3][2] = temp[0];
            vals[0][2] = vals[2][0] = vals[1][3] = vals[3][1] = temp[1];
            vals[0][3] = vals[3][0] = temp[2];
            vals[1][2] = vals[2][1] = temp[3];
        }

        private static void inits2d(vtx_data** graph, /* graph data structure for vertex weights */
            double*[] xvecs, /* values to partition with */
            double*[][] vals /*[4][MAXSETS]*/, /* values in sorted lists */
            int*[][] indices /*[4][MAXSETS]*/, /* indices sorting lists */
            int nvtxs, /* number of vertices */
            double[] dist, /* trial separation point */
            int[][] startvtx /*[4][MAXSETS]*/, /* indices defining separation */
            double[] size, /* size of each set being modified */
            int* sets /* set each vertex gets assigned to */
        )
        {
            /*
                xmid = .25 * (vals[0][1][indices[0][1][vertexCount / 2]] + vals[0][1][indices[0][1][vertexCount / 2 - 1]]);
                ymid = .25 * (vals[0][2][indices[0][2][vertexCount / 2]] + vals[0][2][indices[0][2][vertexCount / 2 - 1]]);
            */
            var xmid = .5 * vals[0][1][indices[0][1][nvtxs / 2]]; /* median x and y values */
            var ymid = .5 * vals[0][2][indices[0][2][nvtxs / 2]]; /* median x and y values */

            dist[0] = -xmid - ymid;
            dist[1] = xmid - ymid;
            dist[2] = -xmid + ymid;
            dist[3] = xmid + ymid;

            /* Now initialize startvtx. */
            startvtx[0][1] = startvtx[2][3] = nvtxs / 2;
            startvtx[0][2] = startvtx[1][3] = nvtxs / 2;
            startvtx[1][2] = findindex(indices[1][2], vals[1][2], dist[2] - dist[1], nvtxs);
            startvtx[0][3] = findindex(indices[0][3], vals[0][3], dist[3] - dist[0], nvtxs);

            for (var i = 0; i < nsets; i++)
            {
                size[i] = 0;
            }

            var bestset = 0; /* set vertex wants to be in */
            for (var i = 1; i <= nvtxs; i++)
            {
                /* Which set is this vertex in? */
                int signy; /* sign values for different target points */
                var signx = signy = -1; /* sign values for different target points */
                double bestval = 0; /* values for determining set preferences */
                for (var j = 0; j < nsets; j++)
                {
                    var val = -dist[j] + 2 * (signx * xvecs[1][i] + signy * xvecs[2][i]); /* values for determining set preferences */
                    if (j == 0 || val < bestval)
                    {
                        bestval = val;
                        bestset = j;
                    }

                    if (signx == 1)
                    {
                        signy *= -1;
                    }

                    signx *= -1;
                }

                sets[i] = bestset;
                size[bestset] += graph[i]->vwgt;
            }
        }

        private static void sorts2d(
            /* Sort the lists needed to find the splitter. */
            double *[][] vals/*[4][MAXSETS]*/,    /* lists of values to sort */
            int *[][]   indices/*[4][MAXSETS]*/, /* indices of sorted lists */
            int     nvtxs                /* number of vertices */
        )
        {
            int *[]temp = new int*[4];    /* place holders for indices */
            int  i;          /* loop counter */

            var space = (int*)Marshal.AllocHGlobal(nvtxs * sizeof(int)); /* space for mergesort routine */

            for (i = 0; i < nlists; i++) {
                temp[i] = (int*)Marshal.AllocHGlobal(nvtxs * sizeof(int));
            }

            MergeSort.ch_mergesort(vals[0][1], nvtxs, temp[0], space);
            MergeSort.ch_mergesort(vals[0][2], nvtxs, temp[1], space);
            MergeSort.ch_mergesort(vals[0][3], nvtxs, temp[2], space);
            MergeSort.ch_mergesort(vals[1][2], nvtxs, temp[3], space);

            Marshal.FreeHGlobal((IntPtr)space);

            indices[0][1] = indices[1][0] = indices[2][3] = indices[3][2] = temp[0];
            indices[0][2] = indices[2][0] = indices[1][3] = indices[3][1] = temp[1];
            indices[0][3] = indices[3][0] = temp[2];
            indices[1][2] = indices[2][1] = temp[3];
        }

        /* Free the space used in the bpmatch routines. */
        private static void free2d(double*[][] vals, int*[][] indices)
        {
            Marshal.FreeHGlobal((IntPtr) vals[0][1]);
            Marshal.FreeHGlobal((IntPtr) vals[0][2]);
            Marshal.FreeHGlobal((IntPtr) vals[0][3]);
            Marshal.FreeHGlobal((IntPtr) vals[1][2]);

            Marshal.FreeHGlobal((IntPtr) indices[0][1]);
            Marshal.FreeHGlobal((IntPtr) indices[0][2]);
            Marshal.FreeHGlobal((IntPtr) indices[0][3]);
            Marshal.FreeHGlobal((IntPtr) indices[1][2]);
        }
    }
}
