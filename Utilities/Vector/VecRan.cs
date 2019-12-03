using static ChacoSharp.Utilities.Randomize;
using static ChacoSharp.Utilities.Norm;

namespace ChacoSharp.Utilities
{
    public static unsafe class VecRan
    {
        /* Fill double vector with random numbers over a range. */
        public static void vecran(double* vec, int beg, int end)
        {
            var pntr = vec + beg;
            for (var i = end - beg + 1; i != 0; i--)
            {
                *pntr++ = drandom();
            }

            ch_normalize(vec, beg, end);
        }

        /* Fill float vector with random numbers over a range. */
        public static void vecran_float(float* vec, int beg, int end)
        {
            var pntr = vec + beg;
            for (var i = end - beg + 1; i != 0; i--)
            {
                *pntr++ = (float) drandom();
            }

            normalize_float(vec, beg, end);
        }
    }
}
