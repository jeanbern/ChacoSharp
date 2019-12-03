namespace ChacoSharp.Utilities
{
    public static unsafe class Update
    {
        /* update - fills double vec1 with vec2 + alpha*vec3 over range*/
        public static void update(double* vec1, int beg, int end, double* vec2, double fac, double* vec3)
        {
            int i;

            vec1 += beg;
            vec2 += beg;
            vec3 += beg;
            for (i = end - beg + 1; i != 0; i--)
            {
                (*vec1++) = (*vec2++) + fac * (*vec3++);
            }
        }

/* update - fills float vec1 with vec2 + alpha*vec3 over range*/
        public static void update_float(float* vec1, int beg, int end, float* vec2, float fac, float* vec3)
        {
            int i;

            vec1 += beg;
            vec2 += beg;
            vec3 += beg;
            for (i = end - beg + 1; i != 0; i--)
            {
                (*vec1++) = (*vec2++) + fac * (*vec3++);
            }
        }
    }
}
