﻿namespace ChacoSharp.Eigen
{
    public static unsafe class CkSturmCount
    {
        /* Check sturmcnt */
        public static void cksturmcnt(double* vec, int beg, int end, double x1, double x2, int* x1ck, int* x2ck, int* numck)
        {
            int i, count;

            count = 0;
            for (i = beg; i <= end; i++)
            {
                if (vec[i] > x1)
                {
                    count += 1;
                }
            }

            *x1ck = end - count;

            count = 0;
            for (i = beg; i <= end; i++)
            {
                if (vec[i] > x2)
                {
                    count += 1;
                }
            }

            *x2ck = end - count;

            count = 0;
            for (i = beg; i <= end; i++)
            {
                if (vec[i] > x1 && vec[i] < x2)
                {
                    count += 1;
                }
            }

            *numck = count;
        }
    }
}
