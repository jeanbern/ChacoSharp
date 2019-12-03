namespace ChacoSharp.Utilities
{
    public static unsafe class CpVec
    {
        /* Copy a range of a double vector to a double vector */
        public static void cpvec(double* copy, int beg, int end, double* vec)
        {
            int i;

            copy = copy + beg;
            vec = vec + beg;
            for (i = end - beg + 1; i != 0; i--)
            {
                *copy++ = *vec++;
            }
        }

        /* Copy a range of a float vector to a double vector */
        public static void float_to_double(double* copy, int beg, int end, float* vec)
        {
            int i;

            copy = copy + beg;
            vec = vec + beg;
            for (i = end - beg + 1; i != 0; i--)
            {
                *copy++ = *vec++;
            }
        }

        /* Copy a range of a double vector to a float vector */
        public static void double_to_float(float* copy, int beg, int end, double* vec)
        {
            int i;

            copy = copy + beg;
            vec = vec + beg;
            for (i = end - beg + 1; i != 0; i--)
            {
                *copy++ = (float) *vec++;
            }
        }
    }
}
