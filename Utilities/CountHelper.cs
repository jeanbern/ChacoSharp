#pragma warning disable HAA0101 // Array allocation for params parameter
#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
using System.Diagnostics;

namespace ChacoSharp.Utilities
{
    public static unsafe class CountHelper
    {
        public static void count(vtx_data** graph, /* graph data structure */
            int nvtxs, /* number of vtxs in graph */
            int* sets, /* processor each vertex is assigned to */
            int nsets, /* number of sets partitioned into */
            int[][] hops, /* hops metric between sets */
            bool dump, /* flag for extended output */
            bool useEdgeWeights /* are edge weights being used? */
        )
        {
            /* number of vtxs in each set */
            var nguys = new int[nsets];

            for (var i = 0; i < nsets; i++)
            {
                nguys[i] = 0;
            }

            var ncross = 0;
            var nhops = 0;
            for (var i = 1; i <= nvtxs; i++)
            {
                nguys[sets[i]] += graph[i]->vwgt;

                for (var j = 1; j < graph[i]->nedges; j++)
                {
                    var neighbor = graph[i]->edges[j]; /* neighbor of a vertex */
                    if (sets[neighbor] == sets[i])
                    {
                        continue;
                    }

                    if (useEdgeWeights)
                    {
                        ncross += (int) graph[i]->ewgts[j];
                        nhops += (int) graph[i]->ewgts[j] * hops[sets[i]][sets[neighbor]];
                    }
                    else
                    {
                        ncross++;
                        nhops += hops[sets[i]][sets[neighbor]];
                    }
                }
            }

            ncross /= 2;
            nhops /= 2;

            var nmax = nguys[0];
            var nmin = nguys[0];
            for (var i = 1; i < nsets; i++)
            {
                if (nguys[i] > nmax)
                {
                    nmax = nguys[i];
                }

                if (nguys[i] < nmin)
                {
                    nmin = nguys[i];
                }
            }

            Trace.WriteLine($"In subgraph: Cuts={ncross:d}, Hops={nhops:d}; Max={nmax:d}, Min={nmin:d} (nvtxs={nvtxs:d}).");

            if (dump)
            {
                for (var i = 0; i < nsets; i++)
                {
                    Trace.WriteLine($" Size of {i:d} = {nguys[i]:d}");
                }

                for (var i = 0; i < nvtxs; i++)
                {
                    Trace.WriteLine($"{sets[i]:d}");
                }

                Trace.WriteLine("\n");
            }
        }
    }
}
