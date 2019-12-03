using System;
using System.Runtime.InteropServices;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Assignment.Median;
using static ChacoSharp.Assignment.RecMedian;
using static ChacoSharp.BipartiteMatching.Map2d;
using static ChacoSharp.BipartiteMatching.Map3d;

namespace ChacoSharp.Assignment
{
    public static unsafe class Mapper
    {
        /// <summary>
        ///
        /// </summary>
        /// <param name="graph">data structure with vertex weights</param>
        /// <param name="xvecs">continuous indicator vectors</param>
        /// <param name="nvtxs">number of vtxs in graph</param>
        /// <param name="active">space for nmyvals ints</param>
        /// <param name="sets">returned processor assignment for my vtxs</param>
        /// <param name="ndims">number of dimensions being divided into</param>
        /// <param name="cube_or_mesh">0 => hypercube, d => d-dimensional mesh</param>
        /// <param name="nsets">number of sets to divide into</param>
        /// <param name="mediantype">type of eigenvector partitioning to use</param>
        /// <param name="goal">desired set sizes</param>
        /// <param name="vertexWeightMax">largest vertex weight</param>
        public static void mapper(vtx_data** graph, double*[] xvecs, int nvtxs, int* active, int* sets, int ndims, bool cube_or_mesh, int nsets,MappingType mediantype, double[] goal, int vertexWeightMax)
        {
            double[] temp_goal = new double[2]; /* combined set goals if using option 1. */
            double wbelow; /* weight of vertices with negative values */
            double wabove; /* weight of vertices with positive values */
            int* temp_sets; /* sets vertices get assigned to */
            int bits; /* bits for assigning set numbers */
            int i, j; /* loop counters */

            /* NOTE: THIS EXPECTS XVECS, NOT YVECS! */

            var useVertexWeights = vertexWeightMax != 1;

            if (ndims == 1 && mediantype == MappingType.MinCost)
            {
                mediantype = MappingType.IndependantMedians; /* simpler call than normal option 1. */
            }

            switch (mediantype)
            {
                case MappingType.CutAtOrigin:
                {
                    /* Divide at zero instead of median. */
                    bits = 1;
                    temp_sets = (int*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));
                    for (j = 1; j <= nvtxs; j++)
                    {
                        sets[j] = 0;
                    }

                    for (i = 1; i <= ndims; i++)
                    {
                        temp_goal[0] = temp_goal[1] = 0;
                        for (j = 0; j < (1 << ndims); j++)
                        {
                            if ((bits & j) != 0)
                            {
                                temp_goal[1] += goal[j];
                            }
                            else
                            {
                                temp_goal[0] += goal[j];
                            }
                        }

                        bits <<= 1;

                        wbelow = wabove = 0;
                        var vertexWeight = 1; /* weight of a vertex */
                        for (j = 1; j <= nvtxs; j++)
                        {
                            if (useVertexWeights)
                            {
                                vertexWeight = graph[j]->vwgt;
                            }

                            if (xvecs[i][j] < 0)
                            {
                                wbelow += vertexWeight;
                            }
                            else if (xvecs[i][j] > 0)
                            {
                                wabove += vertexWeight;
                            }
                        }

                        median_assign(graph, xvecs[i], nvtxs, temp_goal, useVertexWeights, temp_sets, wbelow, wabove, 0.0);

                        for (j = 1; j <= nvtxs; j++)
                        {
                            sets[j] = (sets[j] << 1) + temp_sets[j];
                        }
                    }

                    Marshal.FreeHGlobal((IntPtr) temp_sets);
                    break;
                }
                /* Divide using min-cost assignment. */
                case MappingType.MinCost:
                {
                    if (ndims == 2)
                    {
                        map2d(graph, xvecs, nvtxs, sets, goal, vertexWeightMax);

                    }
                    else if (ndims == 3)
                    {
                        map3d(graph, xvecs, nvtxs, sets, goal, vertexWeightMax);
                    }
                    else
                    {
                        throw new ArgumentOutOfRangeException(nameof(ndims));
                    }

                    break;
                }
                case MappingType.RecursiveMedian:
                {
                    /* Divide recursively using medians. */
                    rec_median_k(graph, xvecs, nvtxs, active, ndims, cube_or_mesh, goal, useVertexWeights, sets);
                    break;
                }
                case MappingType.IndependantMedians:
                {
                    /* Cut with independent medians => unbalanced. */
                    bits = 1;
                    temp_sets = (int*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));
                    for (j = 1; j <= nvtxs; j++)
                    {
                        sets[j] = 0;
                    }

                    for (i = 1; i <= ndims; i++)
                    {
                        temp_goal[0] = temp_goal[1] = 0;
                        for (j = 0; j < (1 << ndims); j++)
                        {
                            if ((bits & j) != 0)
                            {
                                temp_goal[1] += goal[j];
                            }
                            else
                            {
                                temp_goal[0] += goal[j];
                            }
                        }

                        bits <<= 1;

                        median(graph, xvecs[i], nvtxs, active, temp_goal, useVertexWeights, temp_sets);
                        for (j = 1; j <= nvtxs; j++)
                        {
                            sets[j] = (sets[j] << 1) + temp_sets[j];
                        }
                    }

                    Marshal.FreeHGlobal((IntPtr) temp_sets);
                    break;
                }
                case MappingType.Striped:
                {
                    /* Stripe the domain. */
                    rec_median_1(graph, xvecs[1], nvtxs, active, cube_or_mesh, nsets, goal, useVertexWeights, sets, true);
                    break;
                }
            }
        }
    }
}
