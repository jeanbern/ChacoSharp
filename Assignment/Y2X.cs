namespace ChacoSharp.Assignment
{
    public static unsafe class Y2X
    {
        public static void y2x(double*[] xvecs, /* pointer to list of x-vectors */
            int ndims, /* number of divisions to make (# xvecs) */
            int nmyvtxs, /* number of vertices I own (length xvecs) */
            double* wsqrt /* sqrt of vertex weights */
        )

/* Convert from y to x by dividing by wsqrt. */
        {
            double* wptr; /* loops through wsqrt */
            double* xptr; /* loops through elements of a xvec */
            int i, j; /* loop counters */

            if (wsqrt == null)
            {
                return;
            }

            for (i = 1; i <= ndims; i++)
            {
                xptr = xvecs[i];
                wptr = wsqrt;
                for (j = nmyvtxs; j != 0; j--)
                {
                    *(++xptr) /= *(++wptr);
                }
            }
        }

        public static void x2y(double*[] yvecs, /* pointer to list of y-vectors */
            int ndims, /* number of divisions to make (# yvecs) */
            int nmyvtxs, /* number of vertices I own (length yvecs) */
            double* wsqrt /* sqrt of vertex weights */
        )

/* Convert from x to y by multiplying by wsqrt. */
        {
            double* wptr; /* loops through wsqrt */
            double* yptr; /* loops through elements of a yvec */
            int i, j; /* loop counters */

            if (wsqrt == null)
            {
                return;
            }

            for (i = 1; i <= ndims; i++)
            {
                yptr = yvecs[i];
                wptr = wsqrt;
                for (j = nmyvtxs; j != 0; j--)
                {
                    *(++yptr) *= *(++wptr);
                }
            }
        }
    }
}
