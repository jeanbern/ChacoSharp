using System;
using System.Runtime.InteropServices;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Utilities.Randomize;

namespace ChacoSharp.Graph
{
    public static unsafe class SimplePartition
    {
        /* Partition vertices into sets in one of several simplistic ways. */
        public static void simple_part(vtx_data** graph, /* data structure for graph */
            int nvtxs, /* total number of vtxs in graph */
            int* sets, /* sets vertices get assigned to */
            int nsets, /* number of sets at each division */
            in MappingType mappingType, /* type of decomposition */
            double[] goal /* desired set sizes */
        )
        {
            var usingVertexWeights = graph != null;

            switch (mappingType)
            {
                /* Scattered initial decomposition. */
                /* Random initial decomposition. */
                case MappingType.MinCost:
                {
                    ScatteredPartition(graph, nvtxs, sets, nsets, goal, usingVertexWeights);
                    return;
                }
                /* Linearly ordered initial decomposition. */
                case MappingType.RecursiveMedian:
                {
                    RandomPartition(graph, nvtxs, sets, nsets, goal, usingVertexWeights);
                    return;
                }
                case MappingType.IndependantMedians:
                {
                    LinearPartition(graph, nvtxs, sets, nsets, goal, usingVertexWeights);
                    return;
                }
            }
        }

        private static void ScatteredPartition(vtx_data** graph, int nvtxs, int* sets, int nsets, double[] goal, bool usingVertexWeights)
        {
            var weights = new int[MAXSETS]; /* weight assigned to given set so far */
            if (DEBUG_TRACE)
            {
                Console.WriteLine("Generating scattered partition, nvtxs = {0:d}", nvtxs);
            }

            for (var j = 0; j < nsets; j++)
            {
                weights[j] = 0;
            }

            var set = 0;
            for (var i = 1; i <= nvtxs; i++)
            {
                double bestRatio = 2; /* lowest ratio of weight/goal */
                for (var j = 0; j < nsets; j++)
                {
                    var ratio = weights[j] / goal[j]; /* weight/goal */
                    // ReSharper disable once InvertIf
                    if (ratio < bestRatio)
                    {
                        bestRatio = ratio;
                        set = j;
                    }
                }

                if (usingVertexWeights)
                {
                    weights[set] += graph[i]->vwgt;
                }
                else
                {
                    weights[set]++;
                }

                sets[i] = set;
            }
        }

        private static void LinearPartition(vtx_data** graph, int nvtxs, int* sets, int nsets, double[] goal, bool usingVertexWeights)
        {
            if (DEBUG_TRACE)
            {
                Console.WriteLine("Generating linear partition, nvtxs = {0:d}", nvtxs);
            }

            var weight = 0.0d;
            var cutoff = goal[0];
            var set = 0;
            var vwgt = 1;
            var sum = 0.0d;
            for (var i = 1; i <= nvtxs; i++)
            {
                if (usingVertexWeights)
                {
                    vwgt = graph[i]->vwgt;
                }

                if (set < nsets - 1 &&
                    (weight >= cutoff ||
                     (weight + vwgt >= cutoff && sum + vwgt - goal[set] > vwgt - goal[set + 1])))
                {
                    cutoff += goal[++set];
                    sum = 0;
                }

                weight += vwgt;
                sum += vwgt;

                sets[i] = set;
            }
        }

        private static void RandomPartition(vtx_data** graph, int nvtxs, int* sets, int nsets, double[] goal, bool usingVertexWeights)
        {
            if (DEBUG_TRACE)
            {
                Console.WriteLine("Generating random partition, nvtxs = {0:d}", nvtxs);
            }

            /* Construct random order in which to step through graph. */
            int* order = (int*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int)); /* random ordering of vertices */
            for (var i = 1; i <= nvtxs; i++)
            {
                order[i] = i;
            }

            randomize(order, nvtxs);

            var weight = 0.0d;
            var cutoff = goal[0];
            var set = 0;
            var vwgt = 1;
            var sum = 0.0d;
            for (var i = 1; i <= nvtxs; i++)
            {
                if (usingVertexWeights)
                {
                    vwgt = graph[order[i]]->vwgt;
                }

                if (set < nsets - 1 &&
                    (weight >= cutoff ||
                     (weight + vwgt >= cutoff && sum + vwgt - goal[set] > vwgt - goal[set + 1])))
                {
                    cutoff += goal[++set];
                    sum = 0;
                }

                weight += vwgt;
                sum += vwgt;

                sets[order[i]] = set;
            }

            Marshal.FreeHGlobal((IntPtr) order);
        }
    }
}
