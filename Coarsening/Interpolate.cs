namespace ChacoSharp.Coarsening
{
    public static unsafe class Interpolate
    {
        /// <summary>
        /// Interpolate the eigenvector from the coarsened to the original graph.
        /// This may require regenerating mflag and v2cv arrays from merged array.
        /// I also need to undo the merged edges in the reverse order from that in which they were collapsed.
        /// </summary>
        public static void ch_interpolate(double*[] vecs, /* approximate eigenvectors for graph */
            double*[] cvecs, /* exact eigenvectors for coarse graph */
            int ndims, /* number of vectors to interpolate */
            vtx_data** graph, /* array of vtx data for graph */
            int nvtxs, /* number of vertices in graph */
            int* v2cv, /* mapping from vtxs to cvtxs */
            bool useEdgeWeights /* are edge weights being used in fine graph? */
        )
        {
            double* vec; /* pointers into vecs and vecs */
            double* cvec; /* pointers into vecs and vecs */
            int* eptr; /* loops through edge lists */
            float* ewptr; /* loops through edge weights */
            float ewgt; /* value for edge weight */
            double ewsum; /* sum of incident edge weights */
            double sum; /* sum of values of neighbors */
            int neighbor; /* neighboring vertex */
            int degree; /* number of neighbors */
            int npasses = 1; /* number of Gauss-Seidel iterations */
            int pass; /* loops through Gauss-Seidel iterations */
            int i, j, k; /* loop counters */

            /* Uncompress the coarse eigenvector by replicating matched values. */
            for (j = 1; j <= ndims; j++)
            {
                vec = vecs[j];
                cvec = cvecs[j];
                for (i = 1; i <= nvtxs; i++)
                {
                    vec[i] = cvec[v2cv[i]];
                }
            }

            /* Now do a single pass of Gauss-Seidel to smooth eigenvectors. */
            for (pass = 1; pass <= npasses; pass++)
            {
                if (useEdgeWeights)
                {
                    for (j = 1; j <= ndims; j++)
                    {
                        vec = vecs[j];
                        for (i = 1; i <= nvtxs; i++)
                        {
                            eptr = graph[i]->edges;
                            ewptr = graph[i]->ewgts;
                            sum = 0;
                            ewsum = 0;
                            degree = graph[i]->nedges - 1;
                            for (k = degree; k != 0; k--)
                            {
                                neighbor = *(++eptr);
                                ewgt = *(++ewptr);
                                sum += ewgt * vec[neighbor];
                                ewsum += ewgt;
                            }

                            vec[i] = sum / ewsum;
                        }
                    }
                }
                else
                {
                    for (j = 1; j <= ndims; j++)
                    {
                        vec = vecs[j];
                        for (i = 1; i <= nvtxs; i++)
                        {
                            eptr = graph[i]->edges;
                            sum = 0;
                            degree = graph[i]->nedges - 1;
                            for (k = degree; k != 0; k--)
                            {
                                neighbor = *(++eptr);
                                sum += vec[neighbor];
                            }

                            vec[i] = sum / degree;
                        }
                    }
                }
            }
        }
    }
}
