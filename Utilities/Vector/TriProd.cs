namespace ChacoSharp.Utilities
{
    public static unsafe class TriProd
    {
        /* Form inner product of three vectors. */
        public static double tri_prod(double* v1, double* v2, double* v3, double* wsqrt, int n)

        {
            double dot = 0;
            int i;

            if (wsqrt == null)
            {
                /* Unweighted case.  Use weights = 1. */
                for (i = 1; i <= n; i++)
                {
                    dot += v1[i] * v2[i] * v3[i];
                }
            }
            else
            {
                /* Unweighted case. */
                for (i = 1; i <= n; i++)
                {
                    dot += v1[i] * v2[i] * v3[i] / wsqrt[i];
                }
            }

            return (dot);
        }
    }
}
