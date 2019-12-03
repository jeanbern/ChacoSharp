#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
using System;
using System.Runtime.InteropServices;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Eigen.EigenSolve;
using static ChacoSharp.Utilities.MergeSort;
using static ChacoSharp.Graph.Subgraph;
using static ChacoSharp.Graph.FindMaxDeg;
using static ChacoSharp.Connect.FindComps;
using static ChacoSharp.Connect.Connected;
using static ChacoSharp.Utilities.MakeMaps;
using static ChacoSharp.Utilities.MakeVwSqrt;
using static ChacoSharp.Assignment.Y2X;

namespace ChacoSharp
{
    public static unsafe class ChacoSequenceJPetit
    {
/*
 	Aixo es una versio modifica del fitxer chaco_sequence.c de la llibreria Chaco per treure tota la brossa que diua i que no cal.

    Google translate:
        This is a modified version of the chaco_sequence.c file of the library Chaco to remove all the waste that he says and that is not necessary.
*/

/* JPS START */
        public static int[] jps_eigord;

        public static float[] jps_eigvec;
/* JPS END */

        public static void sequence(
            vtx_data** graph, /* graph data structure */
            int nvtxs, /* number of vertices in graph */
            int nedges, /* number of edges in graph */
            bool useEdgeWeights, /* are edge weights being used? */
            double* sqrtVertexWeights, /* sqrt of vertex weights (length nvtxs+1) */
            LanczosType solver_flag, /* which eigensolver should I use? */
            bool rqi_flag, /* use multilevel eigensolver? */
            int vmax, /* if so, how many vtxs to coarsen down to? */
            double eigtol /* tolerance on eigenvectors */
        )
        {
            if (FullTrace)
            {
                Console.WriteLine($"<Entering {nameof(sequence)}>");
            }
/* JPS
            extern char SEQ_FILENAME[];
*/
            vtx_data** subgraph = null; /* subgraph data structure */
            edgeslist* edgeslist; /* edges added for connectivity */
            double*[] yvecs = new double*[MAXDIMS + 1]; /* space for pointing to eigenvectors */
            double[] evals = new double[MAXDIMS + 1]; /* corresponding eigenvalues */
            double[] goal = new double[2]; /* needed for eigen convergence mode = 1 */
            float*[] terminalWeights = new float*[2]; /* dummy vector for terminal weights */
            int* graphToSubgraphMap = null; /* maps graph vtxs to subgraph vtxs */
            int* subtraphToGraphMap = null; /* maps subgraph vtxs to graph vtxs */
            int* degree = null; /* degrees of vertices in subgraph */
            var useVertexWeights = sqrtVertexWeights != null;/* are vertex weights being used? */

            /* Sort each connected component seperately. */
            // Stores the component number for each vertex.
            var vertexComponentMap = (int*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));
/* JPS
            orderfile = fopen(SEQ_FILENAME, "w");
*/
            var space = (int*) Marshal.AllocHGlobal(nvtxs * sizeof(int));
            var connectedComponentCount = find_edges(graph, nvtxs, vertexComponentMap, space, &edgeslist);
            ++connectedComponentCount;

            Console.WriteLine("Found connected components: " + connectedComponentCount);

            free_edgeslist(edgeslist);
            yvecs[1] = (double*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(double));

            terminalWeights[1] = null;
            double* subvwsqrt = null; /* vwsqrt vector for subgraphs */

            if (RQI_CONVERGENCE_MODE != 0)
            {
                throw new InvalidOperationException(nameof(RQI_CONVERGENCE_MODE) + " should be 0 when using " + nameof(SEQUENCE) + " = true");
            }

            if (LANCZOS_CONVERGENCE_MODE != 0)
            {
                throw new InvalidOperationException(nameof(LANCZOS_CONVERGENCE_MODE) + " should be 0 when using " + nameof(SEQUENCE) + " = true");
            }

            var largestConnectedComponentSize = nvtxs;
            int* setsize = null; /* size of each connected component */
            if (connectedComponentCount > 1)
            {
                /* Find size of largest set. */
                setsize = (int*) Marshal.AllocHGlobal(connectedComponentCount * sizeof(int));
                for (var comp = 0; comp < connectedComponentCount; comp++)
                {
                    setsize[comp] = 0;
                }

                for (var i = 1; i <= nvtxs; i++)
                {
                    ++setsize[vertexComponentMap[i]];
                }

                largestConnectedComponentSize = 0;
                for (var comp = 0; comp < connectedComponentCount; comp++)
                {
                    if (setsize[comp] > largestConnectedComponentSize)
                    {
                        largestConnectedComponentSize = setsize[comp];
                    }
                }

                graphToSubgraphMap = (int*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));
                subtraphToGraphMap = (int*) Marshal.AllocHGlobal((largestConnectedComponentSize + 1) * sizeof(int));
                subgraph = (vtx_data**) Marshal.AllocHGlobal((largestConnectedComponentSize + 1) * sizeof(vtx_data*));
                degree = (int*) Marshal.AllocHGlobal((largestConnectedComponentSize + 1) * sizeof(int));
                if (useVertexWeights)
                {
                    subvwsqrt = (double*) Marshal.AllocHGlobal((largestConnectedComponentSize + 1) * sizeof(double));
                }
            }

            var indices = (int*) Marshal.AllocHGlobal(largestConnectedComponentSize * sizeof(int));

            for (var comp = 0; comp < connectedComponentCount; comp++)
            {
                int subgraphVertexCount; /* number of vertices in subgraph */
                int subgraphEdgeCount; /* number of edges in subgraph */
                if (connectedComponentCount > 1)
                {
                    make_maps2(vertexComponentMap, nvtxs, comp, graphToSubgraphMap, subtraphToGraphMap);
                    subgraphVertexCount = setsize[comp];
                    make_subgraph(graph, subgraph, subgraphVertexCount, &subgraphEdgeCount, vertexComponentMap, comp, graphToSubgraphMap, subtraphToGraphMap, degree, useEdgeWeights);

                    if (useVertexWeights)
                    {
                        make_subvector(sqrtVertexWeights, subvwsqrt, subgraphVertexCount, subtraphToGraphMap);
                    }
                }
                else
                {
                    subgraph = graph;
                    subgraphVertexCount = nvtxs;
                    subgraphEdgeCount = nedges;
                    subvwsqrt = sqrtVertexWeights;
                }

                var maxWeightedVertexDegree = find_maxdeg(subgraph, subgraphVertexCount, useEdgeWeights, (float*) null);
                Console.WriteLine("maxWeightedVertexDegree: " + maxWeightedVertexDegree);

                double totalVertexWeight; /* sum of all vertex weights */
                int componentVertexWeightMax; /* largest vertex weight in component */
                if (useVertexWeights)
                {
                    componentVertexWeightMax = 0;
                    totalVertexWeight = 0;
                    for (var i = 1; i <= subgraphVertexCount; i++)
                    {
                        componentVertexWeightMax = Math.Max(subgraph[i]->vwgt, componentVertexWeightMax);
                        totalVertexWeight += subgraph[i]->vwgt;
                    }
                }
                else
                {
                    componentVertexWeightMax = 1;
                    totalVertexWeight = subgraphVertexCount;
                }

                goal[0] = goal[1] = totalVertexWeight / 2;

                Console.WriteLine("useVertexWeights: " + useVertexWeights);
                Console.WriteLine("useEdgeWeights: " + useEdgeWeights);

                if (subgraphVertexCount == 1)
                {
                    yvecs[1][1] = 0;
                }
                else
                {
                    eigensolve(subgraph, subgraphVertexCount, subgraphEdgeCount, maxWeightedVertexDegree, componentVertexWeightMax, subvwsqrt, useVertexWeights, useEdgeWeights, terminalWeights, 0, null, yvecs, evals, false, space, goal, solver_flag, rqi_flag, vmax, 1, MappingType.IndependantMedians, eigtol);
                }

                if (connectedComponentCount > 1)
                {
                    remake_graph(subgraph, subgraphVertexCount, subtraphToGraphMap, degree, useEdgeWeights);
                }

                /* Sort values in eigenvector */
                if (useVertexWeights)
                {
                    y2x(yvecs, 1, subgraphVertexCount, subvwsqrt);
                }

                ch_mergesort(&(yvecs[1][1]), subgraphVertexCount, indices, space);

                /* Print out the order and the corresponding component of the eigenvector */
                if (connectedComponentCount == 1)
                {
                    for (var i = 0; i < subgraphVertexCount; i++)
                    {
/* JPS START */
                        jps_eigord[i] = indices[i] + 1;
                        jps_eigvec[i] = (float) yvecs[1][indices[i] + 1];
/* JPS END */
/* JPS
			fprintf(orderfile, "%-7d   %9.6f\n", indices[i] + 1, yvecs[1][indices[i] + 1]);
*/
                    }
                }
                else
                {
                    for (var i = 0; i < subgraphVertexCount; i++)
                    {
/* JPS START */
                        jps_eigord[i] = subtraphToGraphMap[indices[i] + 1];
                        jps_eigvec[i] = (float) yvecs[1][indices[i] + 1];
/* JPS END */
/* JPS
		fprintf(orderfile, "%-7d   %9.6f\n", loc2glob[indices[i] + 1], yvecs[1][indices[i] + 1]);
*/
                    }
                }
            }

            if (connectedComponentCount > 1)
            {
                Marshal.FreeHGlobal((IntPtr) degree);
                Marshal.FreeHGlobal((IntPtr) subgraph);
                Marshal.FreeHGlobal((IntPtr) subtraphToGraphMap);
                Marshal.FreeHGlobal((IntPtr) graphToSubgraphMap);
                Marshal.FreeHGlobal((IntPtr) setsize);
                if (useVertexWeights)
                {
                    Marshal.FreeHGlobal((IntPtr) subvwsqrt);
                }
            }

            Marshal.FreeHGlobal((IntPtr) yvecs[1]);
            Marshal.FreeHGlobal((IntPtr) indices);
            Marshal.FreeHGlobal((IntPtr) space);
            Marshal.FreeHGlobal((IntPtr) vertexComponentMap);

            Console.WriteLine("{0:d} connected components found.", connectedComponentCount);
        }
    }
}
