namespace ChacoSharp.Utilities
{
    public static unsafe class VecScale
    {
        /* Scales vector by diagonal matrix (passed as vector) over range. */
        public static void scale_diag(double* vec, /* the vector to scale */
            int beg, int end, /* specify the range to norm over */
            double* diag /* vector to scale by */
        )
        {
            /* otherwise return vec unchanged */
            if (diag == null)
            {
                return;
            }

            vec += beg;
            diag += beg;
            for (var i = end - beg + 1; i != 0; i--)
            {
                *vec++ *= *diag++;
            }
        }

        /* Scales vector by diagonal matrix (passed as vector) over range. */
        public static void scale_diag_float(float* vec, /* the vector to scale */
            int beg, int end, /* specify the range to norm over */
            float* diag /* vector to scale by */
        )
        {
            if (diag == null)
            {
                return;
            }

            vec += beg;
            diag += beg;
            for (var i = end - beg + 1; i != 0; i--)
            {
                *vec++ *= *diag++;
            }
        }

        /// <summary>
        /// Scale - fills vec1 with alpha*vec2 over range, double version.
        /// </summary>
        public static void vecscale(double* vec1, int beg, int end, double alpha, double* vec2)
        {
            int i;

            vec1 += beg;
            vec2 += beg;
            for (i = end - beg + 1; i != 0; i--)
            {
                (*vec1++) = alpha * (*vec2++);
            }
        }

        /// <summary>
        /// Scale - fills vec1 with alpha*vec2 over range, float version.
        /// </summary>
        public static void vecscale_float(float* vec1, int beg, int end, float alpha, float* vec2)
        {
            int i;

            vec1 += beg;
            vec2 += beg;
            for (i = end - beg + 1; i != 0; i--)
            {
                (*vec1++) = alpha * (*vec2++);
            }
        }
    }
}
