// ReSharper disable InconsistentNaming
namespace ChacoSharp.Utilities
{
    public static unsafe class Scadd
    {
        /* Scaled add - fills double vec1 with vec1 + alpha*vec2 over range*/
        public static void scadd(double *vec1, int beg, int end, double fac, double *vec2)
        {
            int i;

            vec1 = vec1 + beg;
            vec2 = vec2 + beg;
            for (i = end - beg + 1; i != 0; i--) {
                (*vec1++) += fac * (*vec2++);
            }
        }

/* Scaled add - fills float vec1 with vec1 + alpha*vec2 over range*/
        public static void scadd_float(float *vec1, int beg, int end, float fac, float *vec2)
        {
            int i;

            vec1 = vec1 + beg;
            vec2 = vec2 + beg;
            for (i = end - beg + 1; i != 0; i--) {
                (*vec1++) += fac * (*vec2++);
            }
        }

/* Scaled add - fills double vec1 with vec1 + alpha*vec2 where vec2 is float */
        public static void scadd_mixed(double *vec1, int beg, int end, double fac, float *vec2)
        {
            int i;

            vec1 = vec1 + beg;
            vec2 = vec2 + beg;
            for (i = end - beg + 1; i != 0; i--) {
                (*vec1++) += fac * (*vec2++);
            }
        }

    }
}
