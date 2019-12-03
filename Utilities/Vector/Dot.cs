namespace ChacoSharp.Utilities
{
    public static unsafe class Dot
    {
        /* Returns scalar product of two double n-vectors. */
        public static double dot(double *vec1, int beg, int end, double *vec2)
        {
            int    i;
            double sum;

            sum  = 0.0;
            vec1 = vec1 + beg;
            vec2 = vec2 + beg;
            for (i = end - beg + 1; i != 0; i--) {
                sum += (*vec1++) * (*vec2++);
            }
            return (sum);
        }

/* Returns scalar product of two float n-vectors. */
        public static double dot_float(float *vec1, int beg, int end, float *vec2)
        {
            int   i;
            double sum;

            sum  = 0.0f;
            vec1 = vec1 + beg;
            vec2 = vec2 + beg;
            for (i = end - beg + 1; i != 0; i--) {
                sum += (*vec1++) * (*vec2++);
            }
            return sum;
        }
    }
}
