#pragma warning disable HAA0101 // Array allocation for params parameter
#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
using System;
using static ChacoSharp.StaticConstants;

namespace ChacoSharp.BipartiteMatching
{
    public static unsafe class CheckBp
    {
        /* Confirm that the bipartite match algorithm did the right thing. */
        public static void checkbp(
            vtx_data** graph, /* graph data structure for vertex weights */
            double*[] xvecs, /* values to partition */
            int* sets, /* set assignments returned */
            double[] dists, /* distances that separate sets */
            int nvtxs, /* number of vertices */
            int ndims /* number of dimensions for division */
        )
        {
            int[] signs = new int[MAXDIMS]; /* signs for coordinates of target points */
            int[] sizes = new int[MAXSETS]; /* size of each set */
            int[] weights = new int[MAXSETS]; /* size of each set */
            double setval = 0.0; /* value from assigned set */
            double val, bestval = 0.0; /* value to decide set assignment */
            const double tolerance = 1.0e-8; /* numerical tolerance */
            bool error = false; /* are errors encountered? */
            int bestset = -1; /* set vtx should be assigned to */

            var nsets = 1 << ndims; /* number of sets */

            for (var i = 0; i < nsets; i++)
            {
                sizes[i] = 0;
                weights[i] = 0;
            }

            for (var i = 1; i <= nvtxs; i++)
            {
                /* Is vertex closest to the set it is assigned to? */
                for (var j = 0; j < MAXDIMS; j++)
                {
                    signs[j] = -1;
                }

                bestval = 0;
                for (var j = 0; j < nsets; j++)
                {
                    val = -dists[j];
                    for (var k = 1; k <= ndims; k++)
                    {
                        val += 2 * signs[k - 1] * xvecs[k][i];
                    }

                    if (j == sets[i])
                    {
                        setval = val;
                    }

                    if (j == 0 || val < bestval)
                    {
                        bestval = val;
                        bestset = j;
                    }

                    if (signs[0] == 1 && signs[1] == 1)
                    {
                        signs[2] *= -1;
                    }

                    if (signs[0] == 1)
                    {
                        signs[1] *= -1;
                    }

                    signs[0] *= -1;
                }

                if (Math.Abs(setval - bestval) >= tolerance * (Math.Abs(setval) + Math.Abs(bestval)))
                {
                    Console.WriteLine(" Vtx {0:d} in set {1:d} ({2:e}), but should be in {3:d} ({4:e})", i, sets[i], setval, bestset, bestval);
                    error = true;
                }

                ++sizes[sets[i]];
                weights[sets[i]] += graph[i]->vwgt;
            }

            Console.Write(" Sizes:");
            for (var i = 0; i < nsets; i++)
            {
                Console.Write(" {0:d}({1:d})", sizes[i], weights[i]);
            }

            Console.WriteLine();

            if (error)
            {
                throw new InvalidOperationException("ERROR in checkbp");
            }
        }
    }
}
