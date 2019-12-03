namespace ChacoSharp.Utilities
{
    public static unsafe class SetVec
    {
        /* Set a double precision vector to constant over range. */
        public static void setvec(double* vec, int beg, int end, double setval)
        {
            int i;

            vec = vec + beg;
            for (i = end - beg + 1; i != 0; i--)
            {
                (*vec++) = setval;
            }
        }

        /* Set a float precision vector to constant over range. */
        public static void setvec_float(float* vec, int beg, int end, float setval)
        {
            int i;

            vec = vec + beg;
            for (i = end - beg + 1; i !=  0; i--)
            {
                (*vec++) = setval;
            }
        }
    }
}
