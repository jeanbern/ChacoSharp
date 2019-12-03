using System;
using static ChacoSharp.StaticConstants;

namespace ChacoSharp.Utilities
{
    public static unsafe class MakeVwSqrt
    {
        public static void makevwsqrt(double* vwsqrt, vtx_data** graph, int nvtxs)
/* Make vector of square roots of vertex weights. */
/* vector returned */
/* graph data structure */
/* number of vertices in graph */
        {
            int vwgt; /* vertex weight */
            int i; /* loop counter */

            for (i = 1; i <= nvtxs; i++)
            {
                vwgt = graph[i]->vwgt;
                if (vwgt <= NSQRTS)
                {
                    vwsqrt[i] = SQRTS[vwgt];
                }
                else
                {
                    vwsqrt[i] = Math.Sqrt(vwgt);
                }
            }
        }

/* Extract the subgraph vwsqrt values */
        public static void make_subvector(double* vec, double* subvec, int subnvtxs, int* loc2glob)
/* vector for all vertices */
/* vector for vertices in subgraph */
/* number of vtxs in subgraph */
/* subgraph -> graph numbering map */
        {
            int i;

            for (i = 1; i <= subnvtxs; i++)
            {
                ++subvec;
                (*subvec) = vec[loc2glob[i]];
            }
        }
    }
}
