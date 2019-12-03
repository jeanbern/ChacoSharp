#pragma warning disable HAA0101 // Array allocation for params parameter
#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
using System;

namespace ChacoSharp.Utilities
{
    public static unsafe class CountHelper
    {
        public static void count(vtx_data **graph, /* graph data structure */
        int               nvtxs, /* number of vtxs in graph */
        int *             sets,  /* processor each vertex is assigned to */
        int               nsets, /* number of sets partitioned into */
        int [][] hops,    /* hops metric between sets */
        bool dump,                /* flag for extended output */
        bool useEdgeWeights          /* are edge weights being used? */
        )
        {
            int  ncross;     /* number of outgoing edges */
            int  nhops;      /* number of hops */
            int  neighbor;   /* neighbor of a vertex */
            int  nmax, nmin; /* largest and smallest set sizes */
            int  i, j;       /* loop counters */

            /* number of vtxs in each set */
            var nguys = new int[nsets];

            for (i = 0; i < nsets; i++) {
                nguys[i] = 0;
            }

            ncross = nhops = 0;
            for (i = 1; i <= nvtxs; i++) {
                nguys[sets[i]] += graph[i]->vwgt;

                for (j = 1; j < graph[i]->nedges; j++) {
                    neighbor = graph[i]->edges[j];
                    if (sets[neighbor] != sets[i]) {
                        if (useEdgeWeights) {
                            ncross += (int)graph[i]->ewgts[j];
                            nhops += (int)graph[i]->ewgts[j] * hops[sets[i]][sets[neighbor]];
                        }
                        else {
                            ncross++;
                            nhops += hops[sets[i]][sets[neighbor]];
                        }
                    }
                }
            }

            ncross /= 2;
            nhops /= 2;

            nmax = nguys[0];
            nmin = nguys[0];
            for (i = 1; i < nsets; i++) {
                if (nguys[i] > nmax) {
                    nmax = nguys[i];
                }
                if (nguys[i] < nmin) {
                    nmin = nguys[i];
                }
            }
            Console.WriteLine("In subgraph: Cuts={0:d}, Hops={1:d}; Max={2:d}, Min={3:d} (nvtxs={4:d}).",
                ncross, nhops, nmax, nmin, nvtxs);

            if (dump) {
                for (i = 0; i < nsets; i++) {
                    Console.WriteLine(" Size of {0:d} = {1:d}", i, nguys[i]);
                }

                for (i = 0; i < nvtxs; i++) {
                    Console.WriteLine("{0:d}", sets[i]);
                }
                Console.WriteLine("\n");
            }
        }
    }
}
