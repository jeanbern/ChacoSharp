using static ChacoSharp.Utilities.Perturb;
using static ChacoSharp.StaticConstants;

namespace ChacoSharp.Eigen
{
    public static unsafe class Splarax
    {
        /* Sparse linked A(matrix) times x(vector), double precision. */
        public static void splarax(double* result, /* result of matrix vector multiplication */
            vtx_data** mat, /* graph data structure */
            int n, /* number of rows/columns in matrix */
            double* vec, /* vector being multiplied by matrix */
            double* vwsqrt, /* square roots of vertex weights */
            double* work /* work vector from 1-n */
        )
        {
            vtx_data* mat_i; /* an entry in "mat" */
            double sum; /* sums inner product of matrix-row & vector */
            int* colpntr; /* loops through indices of nonzeros in a row */
            float* wgtpntr; /* loops through values of nonzeros */
            int i, j; /* loop counters */
            double* wrkpntr; /* loops through indices of work vector */
            double* vwsqpntr; /* loops through indices of vwsqrt */
            double* vecpntr; /* loops through indices of vec */
            double* respntr; /* loops through indices of result */
            int last_edge; /* last edge in edge list */

            if (vwsqrt == null)
            {
                /* No vertex weights */
                if (mat[1]->ewgts == null)
                {
                    /* No edge weights */
                    respntr = result;
                    for (i = 1; i <= n; i++)
                    {
                        mat_i = mat[i];
                        colpntr = mat_i->edges;
                        last_edge = mat_i->nedges - 1;
                        sum = last_edge * vec[*colpntr++];
                        for (j = last_edge; j != 0; j--)
                        {
                            sum -= vec[*colpntr++];
                        }

                        *(++respntr) = sum;
                    }
                }
                else
                {
                    /* Edge weights */
                    respntr = result;
                    for (i = 1; i <= n; i++)
                    {
                        mat_i = mat[i];
                        colpntr = mat_i->edges;
                        wgtpntr = mat_i->ewgts;
                        sum = 0.0;
                        for (j = mat_i->nedges; j != 0; j--)
                        {
                            sum -= *wgtpntr++ * vec[*colpntr++];
                        }

                        *(++respntr) = sum; /* -sum if want -Ax */
                    }
                }

                if (PERTURB && NPERTURB > 0 && PERTURB_MAX > 0.0)
                {
                    perturb(result, vec);
                }
            }
            else
            {
                /* Vertex weights */
                if (vwsqrt != null)
                {
                    wrkpntr = work;
                    vecpntr = vec;
                    vwsqpntr = vwsqrt;
                    for (i = n; i != 0; i--)
                    {
                        *(++wrkpntr) = *(++vecpntr) / *(++vwsqpntr);
                    }
                }

                if (mat[1]->ewgts == null)
                {
                    /* No edge weights. */
                    respntr = result;
                    for (i = 1; i <= n; i++)
                    {
                        mat_i = mat[i];
                        colpntr = mat_i->edges;
                        last_edge = mat_i->nedges - 1;
                        sum = (last_edge) * work[*colpntr++];
                        for (j = last_edge; j!= 0; j--)
                        {
                            sum -= work[*colpntr++];
                        }

                        *(++respntr) = sum;
                    }
                }
                else
                {
                    /* Edge weights. */
                    respntr = result;
                    for (i = 1; i <= n; i++)
                    {
                        mat_i = mat[i];
                        colpntr = mat_i->edges;
                        wgtpntr = mat_i->ewgts;
                        sum = 0.0;
                        for (j = mat_i->nedges; j != 0; j--)
                        {
                            sum -= *wgtpntr++ * work[*colpntr++];
                        }

                        *(++respntr) = sum; /* -sum if want -Ax */
                    }
                }

                if (PERTURB && NPERTURB > 0 && PERTURB_MAX > 0.0)
                {
                    perturb(result, work);
                }

                if (vwsqrt != null)
                {
                    respntr = result;
                    vwsqpntr = vwsqrt;
                    for (i = n; i != 0; i--)
                    {
                        *(++respntr) /= *(++vwsqpntr);
                    }
                }
            }
        }

/* Sparse linked A(matrix) times x(vector)i, float version. */
        public static void splarax_float(float* result, /* result of matrix vector multiplication */
            vtx_data** mat, /* graph data structure */
            int n, /* number of rows/columns in matrix */
            float* vec, /* vector being multiplied by matrix */
            float* vwsqrt, /* square roots of vertex weights */
            float* work /* work vector from 1-n */
        )
        {
            vtx_data* mat_i; /* an entry in "mat" */
            double sum; /* sums inner product of matrix-row & vector */
            int* colpntr; /* loops through indices of nonzeros in a row */
            float* wgtpntr; /* loops through values of nonzeros */
            int i, j; /* loop counters */
            float* wrkpntr; /* loops through indices of work vector */
            float* vwsqpntr; /* loops through indices of vwsqrt */
            float* vecpntr; /* loops through indices of vec */
            float* respntr; /* loops through indices of result */
            int last_edge; /* last edge in edge list */

            if (vwsqrt == null)
            {
                /* No vertex weights */
                if (mat[1]->ewgts == null)
                {
                    /* No edge weights */
                    respntr = result;
                    for (i = 1; i <= n; i++)
                    {
                        mat_i = mat[i];
                        colpntr = mat_i->edges;
                        last_edge = mat_i->nedges - 1;
                        sum = (last_edge) * vec[*colpntr++];
                        for (j = last_edge; j != 0; j--)
                        {
                            sum -= vec[*colpntr++];
                        }

                        *(++respntr) = (float)sum;
                    }
                }
                else
                {
                    /* Edge weights */
                    respntr = result;
                    for (i = 1; i <= n; i++)
                    {
                        mat_i = mat[i];
                        colpntr = mat_i->edges;
                        wgtpntr = mat_i->ewgts;
                        sum = 0.0;
                        for (j = mat_i->nedges; j != 0; j--)
                        {
                            sum -= *wgtpntr++ * vec[*colpntr++];
                        }

                        *(++respntr) = (float)sum; /* -sum if want -Ax */
                    }
                }

                if (PERTURB && NPERTURB > 0 && PERTURB_MAX > 0.0)
                {
                    perturb_float(result, vec);
                }
            }
            else
            {
                /* Vertex weights */
                if (vwsqrt != null)
                {
                    wrkpntr = work;
                    vecpntr = vec;
                    vwsqpntr = vwsqrt;
                    for (i = n; i != 0; i--)
                    {
                        *(++wrkpntr) = *(++vecpntr) / *(++vwsqpntr);
                    }
                }

                if (mat[1]->ewgts == null)
                {
                    /* No edge weights. */
                    respntr = result;
                    for (i = 1; i <= n; i++)
                    {
                        mat_i = mat[i];
                        colpntr = mat_i->edges;
                        last_edge = mat_i->nedges - 1;
                        sum = (last_edge) * work[*colpntr++];
                        for (j = last_edge; j!= 0; j--)
                        {
                            sum -= work[*colpntr++];
                        }

                        *(++respntr) = (float)sum;
                    }
                }
                else
                {
                    /* Edge weights. */
                    respntr = result;
                    for (i = 1; i <= n; i++)
                    {
                        mat_i = mat[i];
                        colpntr = mat_i->edges;
                        wgtpntr = mat_i->ewgts;
                        sum = 0.0;
                        for (j = mat_i->nedges; j!= 0; j--)
                        {
                            sum -= *wgtpntr++ * work[*colpntr++];
                        }

                        *(++respntr) = (float)sum; /* -sum if want -Ax */
                    }
                }

                if (PERTURB && NPERTURB > 0 && PERTURB_MAX > 0.0)
                {
                    perturb_float(result, work);
                }

                if (vwsqrt != null)
                {
                    respntr = result;
                    vwsqpntr = vwsqrt;
                    for (i = n; i != 0; i--)
                    {
                        *(++respntr) /= *(++vwsqpntr);
                    }
                }
            }
        }
    }
}
