using System;
using static ChacoSharp.StaticConstants;

namespace ChacoSharp.Eigen
{
    public static unsafe class SoListOut
    {
        /* Print out the orthogonalization set, double version */
        public static void solistout(orthlink** solist, /* vector of pntrs to orthlnks */
            int ngood, /* number of good vecs on list */
            int j /* current number of Lanczos steps */
        )
        {
            int i; /* index */

            for (i = 1; i <= ngood; i++)
            {
                if ((solist[i])->index <= (j / 2))
                {
                    Console.Write(".");
                }
                else
                {
                    Console.Write("+");
                }

                /* Really detailed output: printf("\n"); printf("depth
                   %d\n",(solist[i])->depth); printf("index %d\n",(solist[i])->index);
                   printf("ritzval %g\n",(solist[i])->ritzval); printf("betaji
                   %g\n",(solist[i])->betaji); printf("tau %g\n",(solist[i])->tau);
                   printf("prevtau %g\n",(solist[i])->prevtau);
                   vecout((solist[i])->vec,1,n,"vec", null); */
            }

            Console.WriteLine("{0:d}", ngood);

            if (DEBUG_EVECS > 2)
            {
                Console.Write("  actual indices: ");
                for (i = 1; i <= ngood; i++)
                {
                    Console.Write(" {0:d}", solist[i]->index);
                }

                Console.WriteLine();
            }
        }

/* Print out the orthogonalization set, float version */
        public static void solistout_float(orthlink_float** solist, /* vector of pntrs to orthlnks */
            int ngood, /* number of good vecs on list */
            int j /* current number of Lanczos steps */
        )
        {
            int i; /* index */

            for (i = 1; i <= ngood; i++)
            {
                if ((solist[i])->index <= (j / 2))
                {
                    Console.Write(".");
                }
                else
                {
                    Console.Write("+");
                }
            }

            Console.WriteLine("{0:d}", ngood);
        }
    }
}
