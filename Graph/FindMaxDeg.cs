using System;
using static ChacoSharp.StaticConstants;

namespace ChacoSharp.Graph
{
    public static unsafe class FindMaxDeg
    {
        /// <summary>
        /// Find the maximum weighted degree of a vertex.
        /// </summary>
        /// <param name="graph">graph data structure</param>
        /// <param name="nvtxs">number of vertices</param>
        /// <param name="useEdgeWeights">are edge weights being used?</param>
        /// <param name="pmax_ewgt">returns largest edge weight if not null</param>
        /// <returns>The maximum weighted degree of a vertex.</returns>
        public static double find_maxdeg(vtx_data** graph, int nvtxs, bool useEdgeWeights, float* pmax_ewgt)
        {
            if (FullTrace)
            {
                Console.WriteLine($"<Entering {nameof(find_maxdeg)}>");
            }

            /* Find the maximum weighted degree of a vertex. */
            double maxdeg = 0;
            if (useEdgeWeights)
            {
                if (pmax_ewgt != null)
                {
                    float largestEdgeWeight = 0; /* largest edge weight in graph */
                    for (var i = 1; i <= nvtxs; i++)
                    {
                        if (-graph[i]->ewgts[0] > maxdeg)
                        {
                            maxdeg = -graph[i]->ewgts[0];
                        }

                        var eptr = graph[i]->ewgts; /* loops through edge weights */
                        for (var j = graph[i]->nedges - 1; j != 0; j--)
                        {
                            // ReSharper disable once ArrangeRedundantParentheses
                            var ewgt = *(++eptr); /* edge weight value */
                            if (ewgt > largestEdgeWeight)
                            {
                                largestEdgeWeight = ewgt;
                            }
                        }
                    }

                    *pmax_ewgt = largestEdgeWeight;
                }
                else
                {
                    for (var i = 1; i <= nvtxs; i++)
                    {
                        if (-graph[i]->ewgts[0] > maxdeg)
                        {
                            maxdeg = -graph[i]->ewgts[0];
                        }
                    }
                }
            }
            else
            {
                for (var i = 1; i <= nvtxs; i++)
                {
                    if (graph[i]->nedges > maxdeg)
                    {
                        maxdeg = graph[i]->nedges - 1;
                    }
                }

                if (pmax_ewgt != null)
                {
                    *pmax_ewgt = 1;
                }
            }

            return maxdeg;
        }
    }
}
