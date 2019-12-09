#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
#pragma warning disable HAA0101 // Array allocation for params parameter
using System;
using System.Diagnostics;

namespace ChacoSharp.Graph
{
    public static unsafe class CheckGraph
    {
        /// <summary>
        /// Check graph for errors
        /// </summary>
        /// <param name="graph">graph data structure</param>
        /// <param name="nvtxs">number of vertices</param>
        /// <param name="nedges">number of edges</param>
        /// <returns>A flag representing if inconsistencies were found.</returns>
        public static bool check_graph(vtx_data** graph, int nvtxs, int nedges)
        {
            var edgeWeightTotal = 0.0d; /* sum of edge weights */

            var flag = false;
            var isolatedVertexCount = 0;
            var badWeightVertexCount = 0;
            var useEdgeWeights = graph[1]->ewgts != null;
            var neighborCount = 0;
            for (var i = 1; i <= nvtxs; i++)
            {
                neighborCount += graph[i]->nedges - 1;

                if (graph[i]->edges[0] != i)
                {
                    Trace.WriteLine($" Self edge wrong for vtx {i:d}");
                    flag = true;
                }

                if (graph[i]->nedges == 1)
                {
                    if (isolatedVertexCount == 0)
                    {
                        Trace.WriteLine($"WARNING: Vertex {i:d} has no neighbors");
                    }

                    ++isolatedVertexCount;
                }

                if (graph[i]->vwgt <= 0)
                {
                    if (badWeightVertexCount == 0)
                    {
                        Trace.WriteLine($"Vertex {i:d} has bad vertex weight {graph[i]->vwgt:d}.");
                    }

                    ++badWeightVertexCount;
                    flag = true;
                }

                if (useEdgeWeights)
                {
                    edgeWeightTotal = graph[i]->ewgts[0];
                }

                for (var j = 1; j < graph[i]->nedges; j++)
                {
                    var neighbor = graph[i]->edges[j]; /* neighbor of a vertex */
                    if (useEdgeWeights)
                    {
                        edgeWeightTotal += graph[i]->ewgts[j];
                    }

                    /* Move it to the end and delete instead? */
                    if (neighbor == i)
                    {
                        Trace.WriteLine($"Self edge ({i:d},{neighbor:d}) not allowed");
                        flag = true;
                    }

                    if (neighbor < 1 || neighbor > nvtxs)
                    {
                        Trace.WriteLine($"Edge ({i:d},{neighbor:d}) included, but nvtxs = {nvtxs:d}");
                        flag = true;
                    }

                    /* Move it to the end and delete instead? */
                    if (useEdgeWeights && graph[i]->ewgts[j] <= 0)
                    {
                        Trace.WriteLine($"Bad edge weight {graph[i]->ewgts[j]:g} for edge ({i:d}, {neighbor:d})");
                        flag = true;
                    }

                    float eweight; /* edge weight */
                    if (!is_an_edge(graph[neighbor], i, &eweight))
                    {
                        Trace.WriteLine($"Edge ({i:d},{neighbor:d}) included but not ({neighbor:d},{i:d})");
                        flag = true;
                    }
                    // ReSharper disable once CompareOfFloatsByEqualityOperator
                    else if (useEdgeWeights && eweight != graph[i]->ewgts[j])
                    {
                        Trace.WriteLine($"Weight of ({i:d},{neighbor:d})={graph[i]->ewgts[j]:g}, but weight of ({neighbor:d},{i:d})={eweight:g}");
                        flag = true;
                    }
                }

                // ReSharper disable once InvertIf
                if (useEdgeWeights && Math.Abs(edgeWeightTotal) > 1.0e-7 * Math.Abs(graph[i]->ewgts[0]))
                {
                    Trace.WriteLine($"Sum of edge weights for vertex {i:d} = {edgeWeightTotal:g}");
                    flag = true;
                }
            }

            if (isolatedVertexCount > 1)
            {
                Trace.WriteLine($"WARNING: {isolatedVertexCount:d} vertices have no neighbors");
            }

            if (badWeightVertexCount > 1)
            {
                Trace.WriteLine($"{badWeightVertexCount:d} vertices have bad vertex weights");
            }

            // ReSharper disable once InvertIf
            if (neighborCount != 2 * nedges)
            {
                Trace.WriteLine($" twice nedges = {2 * nedges:d}, but I count {neighborCount:d}");
                flag = true;
            }

            return flag;
        }

        public static bool is_an_edge(vtx_data* vertex, /* data for a vertex */
            int v2, /* neighbor to look for */
            float* weight2 /* weight of edge if found */
        )
        {
            int i; /* loop counter */

            for (i = 1; i < vertex->nedges; i++)
            {
                // ReSharper disable once InvertIf
                if (vertex->edges[i] == v2)
                {
                    if (vertex->ewgts != null)
                    {
                        *weight2 = vertex->ewgts[i];
                    }
                    else
                    {
                        *weight2 = 1;
                    }

                    return true;
                }
            }

            return false;
        }
    }
}
