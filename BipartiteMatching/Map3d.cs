#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
#pragma warning disable HAA0101 // Array allocation for params parameter
using System;
using System.Runtime.InteropServices;
using static ChacoSharp.BipartiteMatching.CheckBp;
using static ChacoSharp.BipartiteMatching.MoveVtxs;
using static ChacoSharp.BipartiteMatching.FindIndex;
using static ChacoSharp.Utilities.MergeSort;
using static ChacoSharp.StaticConstants;

namespace ChacoSharp.BipartiteMatching
{
    public static unsafe class Map3d
    {
        public static void map3d(vtx_data** graph, /* graph data structure */
            double*[] xvecs, /* vectors to partition */
            int nvtxs, /* number of vertices */
            int* sets, /* set each vertex gets assigned to */
            double[] goal, /* desired set sizes */
            int vwgt_max /* largest vertex weight */
        )
        {
            double*[][] vals = new double*[8][]; //[8][MAXSETS];     /* values in sorted lists */
            for (var i = 0; i < vals.Length; i++)
            {
                vals[i] = new double*[MAXSETS];
            }

            double[] dist = new double[8]; /* trial separation point */
            double[] size = new double[8]; /* sizes of each set being modified */
            int*[][] indices = new int*[8][]; /* indices sorting lists */
            for (var i = 0; i < indices.Length; i++)
            {
                indices[i] = new int*[MAXSETS];
            }

            int[][] startvtx = new int[8][]; /* indices defining separation */
            for (var i = 0; i < startvtx.Length; i++)
            {
                startvtx[i] = new int[MAXSETS];
            }

            int nsection = 3; /* number of xvectors */
            int nsets = 8; /* number of sets being divided into */

            N_VTX_CHECKS = N_VTX_MOVES = 0;

            /* Generate all the lists of values that need to be sorted. */
            genvals3d(xvecs, vals, nvtxs);

            /* Sort the lists of values. */
            sorts3d(vals, indices, nvtxs);

            /* Now initialize distances using median values, and assign to sets. */
            inits3d(graph, xvecs, vals, indices, nvtxs, dist, startvtx, size, sets);

            /* Determine largest and smallest allowed sizes for each set. */
            /* (For now all are the same, but easily changed.) */

            if (DEBUG_BPMATCH == DebugFlagBP.ErrorChecking)
            {
                Console.WriteLine(" Calling check before movevtxs");
                checkbp(graph, xvecs, sets, dist, nvtxs, nsection);
            }

            movevtxs(graph, nvtxs, nsets, dist, indices, vals, startvtx, sets, size, goal, vwgt_max);

            if (DEBUG_BPMATCH != DebugFlagBP.NoDebugging)
            {
                Console.WriteLine(" N_VTX_CHECKS = {0:d}, N_VTX_MOVES = {1:d}", N_VTX_CHECKS, N_VTX_MOVES);
                checkbp(graph, xvecs, sets, dist, nvtxs, nsection);
            }

            free3d(vals, indices);
        }

        private static void free3d(double*[][] vals /*[8][MAXSETS]*/, int*[][] indices /*[8][MAXSETS]*/)
        {

            Marshal.FreeHGlobal((IntPtr) vals[0][1]);
            Marshal.FreeHGlobal((IntPtr) vals[0][2]);
            Marshal.FreeHGlobal((IntPtr) vals[0][4]);
            Marshal.FreeHGlobal((IntPtr) vals[0][3]);
            Marshal.FreeHGlobal((IntPtr) vals[1][2]);
            Marshal.FreeHGlobal((IntPtr) vals[0][5]);
            Marshal.FreeHGlobal((IntPtr) vals[1][4]);
            Marshal.FreeHGlobal((IntPtr) vals[0][6]);
            Marshal.FreeHGlobal((IntPtr) vals[2][4]);
            Marshal.FreeHGlobal((IntPtr) vals[0][7]);
            Marshal.FreeHGlobal((IntPtr) vals[1][6]);
            Marshal.FreeHGlobal((IntPtr) vals[2][5]);
            Marshal.FreeHGlobal((IntPtr) vals[3][4]);

            Marshal.FreeHGlobal((IntPtr) indices[0][1]);
            Marshal.FreeHGlobal((IntPtr) indices[0][2]);
            Marshal.FreeHGlobal((IntPtr) indices[0][4]);
            Marshal.FreeHGlobal((IntPtr) indices[0][3]);
            Marshal.FreeHGlobal((IntPtr) indices[1][2]);
            Marshal.FreeHGlobal((IntPtr) indices[0][5]);
            Marshal.FreeHGlobal((IntPtr) indices[1][4]);
            Marshal.FreeHGlobal((IntPtr) indices[0][6]);
            Marshal.FreeHGlobal((IntPtr) indices[2][4]);
            Marshal.FreeHGlobal((IntPtr) indices[0][7]);
            Marshal.FreeHGlobal((IntPtr) indices[1][6]);
            Marshal.FreeHGlobal((IntPtr) indices[2][5]);
            Marshal.FreeHGlobal((IntPtr) indices[3][4]);
        }

        private static void inits3d(vtx_data** graph, /* graph data structure for vertex weights */
            double*[] xvecs, /* values to partition with */
            double*[][] vals /*[8][MAXSETS]*/, /* values in sorted lists */
            int*[][] indices /*[8][MAXSETS]*/, /* indices sorting lists */
            int nvtxs, /* number of vertices */
            double[] dist, /* trial separation point */
            int[][] startvtx /*[8][MAXSETS]*/, /* indices defining separation */
            double[] size, /* size of each set being modified */
            int* sets /* set each vertex gets assigned to */
        )
        {
            double xmid, ymid, zmid; /* median x, y and z values */
            double val, bestval; /* values for determining set preferences */
            int bestset = 0; /* set vertex wants to be in */
            int signx, signy, signz; /* sign values for different target points */
            int nsets = 8; /* number of different sets */
            int i, j; /* loop counters */

            /*
                xmid = .25 * (vals[0][1][indices[0][1][nvtxs / 2]] +
                              vals[0][1][indices[0][1][nvtxs / 2 - 1]]);
                ymid = .25 * (vals[0][2][indices[0][2][nvtxs / 2]] +
                              vals[0][2][indices[0][2][nvtxs / 2 - 1]]);
                zmid = .25 * (vals[0][4][indices[0][4][nvtxs / 2]] +
                              vals[0][4][indices[0][4][nvtxs / 2 - 1]]);
            */
            xmid = .5 * vals[0][1][indices[0][1][nvtxs / 2]];
            ymid = .5 * vals[0][2][indices[0][2][nvtxs / 2]];
            zmid = .5 * vals[0][4][indices[0][4][nvtxs / 2]];

            dist[0] = -xmid - ymid - zmid;
            dist[1] = xmid - ymid - zmid;
            dist[2] = -xmid + ymid - zmid;
            dist[3] = xmid + ymid - zmid;
            dist[4] = -xmid - ymid + zmid;
            dist[5] = xmid - ymid + zmid;
            dist[6] = -xmid + ymid + zmid;
            dist[7] = xmid + ymid + zmid;

            /* Now initialize startvtxs. */
            startvtx[0][1] = startvtx[2][3] = startvtx[4][5] = startvtx[6][7] = nvtxs / 2;
            startvtx[0][2] = startvtx[1][3] = startvtx[4][6] = startvtx[5][7] = nvtxs / 2;
            startvtx[0][4] = startvtx[1][5] = startvtx[2][6] = startvtx[3][7] = nvtxs / 2;
            startvtx[0][3] = startvtx[4][7] = findindex(indices[0][3], vals[0][3], dist[3] - dist[0], nvtxs);
            startvtx[1][2] = startvtx[5][6] = findindex(indices[1][2], vals[1][2], dist[2] - dist[1], nvtxs);
            startvtx[0][5] = startvtx[2][7] = findindex(indices[0][5], vals[0][5], dist[5] - dist[0], nvtxs);
            startvtx[1][4] = startvtx[3][6] = findindex(indices[1][4], vals[1][4], dist[4] - dist[1], nvtxs);
            startvtx[0][6] = startvtx[1][7] = findindex(indices[0][6], vals[0][6], dist[6] - dist[0], nvtxs);
            startvtx[2][4] = startvtx[3][5] = findindex(indices[2][4], vals[2][4], dist[4] - dist[2], nvtxs);
            startvtx[0][7] = findindex(indices[0][7], vals[0][7], dist[7] - dist[0], nvtxs);
            startvtx[1][6] = findindex(indices[1][6], vals[1][6], dist[6] - dist[1], nvtxs);
            startvtx[2][5] = findindex(indices[2][5], vals[2][5], dist[5] - dist[2], nvtxs);
            startvtx[3][4] = findindex(indices[3][4], vals[3][4], dist[4] - dist[3], nvtxs);

            /* Finally, determine the set sizes based on this splitter. */

            for (i = 0; i < nsets; i++)
            {
                size[i] = 0;
            }

            for (i = 1; i <= nvtxs; i++)
            {
                /* Which set is this vertex in? */
                signx = signy = signz = -1;
                bestval = 0;
                for (j = 0; j < nsets; j++)
                {
                    val = -dist[j] + 2 * (signx * xvecs[1][i] + signy * xvecs[2][i] + signz * xvecs[3][i]);
                    if (j == 0 || val < bestval)
                    {
                        bestval = val;
                        bestset = j;
                    }

                    if (signx == 1 && signy == 1)
                    {
                        signz *= -1;
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

        private static void genvals3d(
            /* Create lists of sets of values to be sorted. */
            double*[] xvecs, /* vectors to partition */
            double*[][] vals /*[8][MAXSETS]*/, /* ptrs to lists of values */
            int nvtxs /* number of values */
        )
        {
            int nsets = 8; /* number of sets */
            int nlists = 13; /* number of lists to generate */
            double*[] temp = new double*[13]; /* place holders for vals */
            int i, j; /* loop counter */

            for (i = 0; i < nlists; i++)
            {
                temp[i] = (double*) Marshal.AllocHGlobal(nvtxs * sizeof(double));
            }

            for (i = 1; i <= nvtxs; i++)
            {
                temp[0][i - 1] = 4 * xvecs[1][i];
                temp[1][i - 1] = 4 * xvecs[2][i];
                temp[2][i - 1] = 4 * xvecs[3][i];
                temp[3][i - 1] = 4 * (xvecs[1][i] + xvecs[2][i]);
                temp[4][i - 1] = 4 * (-xvecs[1][i] + xvecs[2][i]);
                temp[5][i - 1] = 4 * (xvecs[1][i] + xvecs[3][i]);
                temp[6][i - 1] = 4 * (-xvecs[1][i] + xvecs[3][i]);
                temp[7][i - 1] = 4 * (xvecs[2][i] + xvecs[3][i]);
                temp[8][i - 1] = 4 * (-xvecs[2][i] + xvecs[3][i]);
                temp[9][i - 1] = 4 * (xvecs[1][i] + xvecs[2][i] + xvecs[3][i]);
                temp[10][i - 1] = 4 * (-xvecs[1][i] + xvecs[2][i] + xvecs[3][i]);
                temp[11][i - 1] = 4 * (xvecs[1][i] - xvecs[2][i] + xvecs[3][i]);
                temp[12][i - 1] = 4 * (-xvecs[1][i] - xvecs[2][i] + xvecs[3][i]);
            }

            vals[0][1] = vals[2][3] = vals[4][5] = vals[6][7] = temp[0];
            vals[0][2] = vals[1][3] = vals[4][6] = vals[5][7] = temp[1];
            vals[0][4] = vals[1][5] = vals[2][6] = vals[3][7] = temp[2];
            vals[0][3] = vals[4][7] = temp[3];
            vals[1][2] = vals[5][6] = temp[4];
            vals[0][5] = vals[2][7] = temp[5];
            vals[1][4] = vals[3][6] = temp[6];
            vals[0][6] = vals[1][7] = temp[7];
            vals[2][4] = vals[3][5] = temp[8];
            vals[0][7] = temp[9];
            vals[1][6] = temp[10];
            vals[2][5] = temp[11];
            vals[3][4] = temp[12];

            for (i = 0; i < nsets; i++)
            {
                for (j = i + 1; j < nsets; j++)
                {
                    vals[j][i] = vals[i][j];
                }
            }
        }

        private static void sorts3d(
            /* Sort the lists needed to find the splitter. */
            double*[][] vals /*[8][MAXSETS]*/, /* lists of values to sort */
            int*[][] indices /*[8][MAXSETS]*/, /* indices of sorted lists */
            int nvtxs /* number of vertices */
        )
        {
            int* space; /* space for mergesort routine */
            int nsets = 8; /* number of sets */
            int nlists = 13; /* number of directions to sort */
            int*[] temp = new int*[13]; /* place holders for indices */
            int i, j; /* loop counter */

            space = (int*) Marshal.AllocHGlobal(nvtxs * sizeof(int));

            for (i = 0; i < nlists; i++)
            {
                temp[i] = (int*) Marshal.AllocHGlobal(nvtxs * sizeof(int));
            }

            ch_mergesort(vals[0][1], nvtxs, temp[0], space);
            ch_mergesort(vals[0][2], nvtxs, temp[1], space);
            ch_mergesort(vals[0][4], nvtxs, temp[2], space);
            ch_mergesort(vals[0][3], nvtxs, temp[3], space);
            ch_mergesort(vals[1][2], nvtxs, temp[4], space);
            ch_mergesort(vals[0][5], nvtxs, temp[5], space);
            ch_mergesort(vals[1][4], nvtxs, temp[6], space);
            ch_mergesort(vals[0][6], nvtxs, temp[7], space);
            ch_mergesort(vals[2][4], nvtxs, temp[8], space);
            ch_mergesort(vals[0][7], nvtxs, temp[9], space);
            ch_mergesort(vals[1][6], nvtxs, temp[10], space);
            ch_mergesort(vals[2][5], nvtxs, temp[11], space);
            ch_mergesort(vals[3][4], nvtxs, temp[12], space);

            Marshal.FreeHGlobal((IntPtr) space);

            indices[0][1] = indices[2][3] = indices[4][5] = indices[6][7] = temp[0];
            indices[0][2] = indices[1][3] = indices[4][6] = indices[5][7] = temp[1];
            indices[0][4] = indices[1][5] = indices[2][6] = indices[3][7] = temp[2];
            indices[0][3] = indices[4][7] = temp[3];
            indices[1][2] = indices[5][6] = temp[4];
            indices[0][5] = indices[2][7] = temp[5];
            indices[1][4] = indices[3][6] = temp[6];
            indices[0][6] = indices[1][7] = temp[7];
            indices[2][4] = indices[3][5] = temp[8];
            indices[0][7] = temp[9];
            indices[1][6] = temp[10];
            indices[2][5] = temp[11];
            indices[3][4] = temp[12];

            for (i = 0; i < nsets; i++)
            {
                for (j = i + 1; j < nsets; j++)
                {
                    indices[j][i] = indices[i][j];
                }
            }
        }
    }
}
