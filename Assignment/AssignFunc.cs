#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
#pragma warning disable HAA0101 // Array allocation for params parameter
using System;
using System.Diagnostics;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Assignment.Y2X;
using static ChacoSharp.Assignment.Mapper;
using static ChacoSharp.Assignment.Rotate;
using static ChacoSharp.Optimize.Opt3d;
using static ChacoSharp.Optimize.Opt2d;
using static ChacoSharp.Utilities.TriProd;

namespace ChacoSharp.Assignment
{
    public static unsafe class AssignFunc
    {
        /// <summary>
        ///
        /// </summary>
        /// <param name="graph">data structure with vtx weights</param>
        /// <param name="yvecs">ptr to list of y-vectors (lengths nvtxs+1)</param>
        /// <param name="nvtxs">number of vertices in graph</param>
        /// <param name="ndims">number of vectors for dividing</param>
        /// <param name="cubeOrMesh">0 => hypercube, d => d-dimensional mesh</param>
        /// <param name="nsets">number of sets to divide into</param>
        /// <param name="wsqrt">sqrt of vertex weights</param>
        /// <param name="sets">processor assignment for my vtxs</param>
        /// <param name="active">space for nvtxs integers</param>
        /// <param name="mediantype">which partitioning strategy to use</param>
        /// <param name="goal">desired set sizes</param>
        /// <param name="vertexWeightMax">largest vertex weight</param>
        public static void Assign(vtx_data** graph, double*[] yvecs, int nvtxs, int ndims, bool cubeOrMesh, int nsets, double* wsqrt, int* sets, int* active, MappingType mediantype, double[] goal, int vertexWeightMax)
        {
            if (yvecs == null)
            {
                throw new ArgumentNullException(nameof(yvecs));
            }

            if (DEBUG_TRACE)
            {
                Trace.WriteLine($"<Entering {nameof(Assign)}, {nameof(nvtxs)} = {nvtxs:D}, {nameof(ndims)} = {ndims:D}>");
            }

            var useVertexWeights = vertexWeightMax != 1;

            switch (ndims)
            {
                case 1:
                {
                    break;
                }
                case 2:
                {
                    var theta = opt2d(graph, yvecs, nvtxs, nvtxs);
                    Rotate2d(yvecs, nvtxs, theta);
                    break;
                }
                case 3:
                {
                    if (DEBUG_ASSIGN)
                    {
                        var temp = tri_prod(yvecs[1], yvecs[2], yvecs[3], wsqrt, nvtxs);
                        Trace.WriteLine($"Before rotation, 3-way orthogonality = {temp:e}");
                    }

                    double phi, gamma, theta; /* angles for optimal rotation */
                    opt3d(graph, yvecs, nvtxs, nvtxs, wsqrt, &theta, &phi, &gamma, useVertexWeights);
                    Rotate3d(yvecs, nvtxs, theta, phi, gamma);

                    if (DEBUG_ASSIGN)
                    {
                        var temp = tri_prod(yvecs[1], yvecs[2], yvecs[3], wsqrt, nvtxs);
                        Trace.WriteLine($"After rotation ({theta:f},{phi:f},{gamma:f}), 3-way orthogonality = {temp:e}");
                    }

                    break;
                }
                default:
                {
                    throw new ArgumentOutOfRangeException(nameof(ndims));
                }
            }

            /* Unscale yvecs to get xvecs. */
            y2x(yvecs, ndims, nvtxs, wsqrt);
            mapper(graph, yvecs, nvtxs, active, sets, ndims, cubeOrMesh, nsets, mediantype, goal, vertexWeightMax);
        }
    }
}