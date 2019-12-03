namespace ChacoSharp.Graph
{
    public static unsafe class CountWeights
    {
        public static void count_weights(vtx_data **graph,      /* data structure for graph */
        int               nvtxs,      /* number of vtxs in graph */
        int *             sets,       /* set each vertex is assigned to */
        int               nsets,      /* number of sets in this division */
        double []          weights,    /* vertex weights in each set */
        bool               using_vwgts /* are vertex weights being used? */
        )

        {
            int i; /* loop counters */

            /* Compute the sum of vertex weights for each set. */
            for (i = 0; i < nsets; i++) {
                weights[i] = 0;
            }

            if (using_vwgts) {
                for (i = 1; i <= nvtxs; i++) {
                    weights[sets[i]] += graph[i]->vwgt;
                }
            }

            else {
                for (i = 1; i <= nvtxs; i++) {
                    weights[sets[i]]++;
                }
            }
        }
    }
}
