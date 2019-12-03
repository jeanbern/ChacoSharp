#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
#pragma warning disable HAA0101 // Array allocation for params parameter
using System;

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
                    Console.WriteLine(" Self edge wrong for vtx {0:d}", i);
                    flag = true;
                }

                if (graph[i]->nedges == 1)
                {
                    if (isolatedVertexCount == 0)
                    {
                        Console.WriteLine("WARNING: Vertex {0:d} has no neighbors", i);
                    }

                    ++isolatedVertexCount;
                }

                if (graph[i]->vwgt <= 0)
                {
                    if (badWeightVertexCount == 0)
                    {
                        Console.WriteLine("Vertex {0:d} has bad vertex weight {1:d}.", i, graph[i]->vwgt);
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
                        Console.WriteLine("Self edge ({0:d},{1:d}) not allowed", i, neighbor);
                        flag = true;
                    }

                    if (neighbor < 1 || neighbor > nvtxs)
                    {
                        Console.WriteLine("Edge ({0:d},{1:d}) included, but nvtxs = {2:d}", i, neighbor, nvtxs);
                        flag = true;
                    }

                    /* Move it to the end and delete instead? */
                    if (useEdgeWeights && graph[i]->ewgts[j] <= 0)
                    {
                        Console.WriteLine("Bad edge weight {0:g} for edge ({1:d}, {2:d})", graph[i]->ewgts[j], i, neighbor);
                        flag = true;
                    }

                    float eweight; /* edge weight */
                    if (!is_an_edge(graph[neighbor], i, &eweight))
                    {
                        Console.WriteLine("Edge ({0:d},{1:d}) included but not ({2:d},{3:d})", i, neighbor, neighbor, i);
                        flag = true;
                    }
                    // ReSharper disable once CompareOfFloatsByEqualityOperator
                    else if (useEdgeWeights && eweight != graph[i]->ewgts[j])
                    {
                        Console.WriteLine("Weight of ({0:d},{1:d})={2:g}, but weight of ({3:d},{4:d})={5:g}", i, neighbor, graph[i]->ewgts[j], neighbor, i, eweight);
                        flag = true;
                    }
                }

                // ReSharper disable once InvertIf
                if (useEdgeWeights && Math.Abs(edgeWeightTotal) > 1.0e-7 * Math.Abs(graph[i]->ewgts[0]))
                {
                    Console.WriteLine("Sum of edge weights for vertex {0:d} = {1:g}", i, edgeWeightTotal);
                    flag = true;
                }
            }

            if (isolatedVertexCount > 1)
            {
                Console.WriteLine("WARNING: {0:d} vertices have no neighbors", isolatedVertexCount);
            }

            if (badWeightVertexCount > 1)
            {
                Console.WriteLine("{0:d} vertices have bad vertex weights", badWeightVertexCount);
            }

            // ReSharper disable once InvertIf
            if (neighborCount != 2 * nedges)
            {
                Console.WriteLine(" twice nedges = {0:d}, but I count {1:d}", 2 * nedges, neighborCount);
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
