using System;
using System.Runtime.InteropServices;
using static ChacoSharp.Utilities.Scadd;
using static ChacoSharp.Utilities.Dot;

namespace ChacoSharp.Eigen
{
    public static unsafe class Orthogonalization
    {
        /* Check orthogonality of vector set */
        public static void checkorth(double** mat, int n, int dim)
        {
            int i, j; /* loop idices */
            double measure; /* Froebenius norm */
            double prod; /* value of dot product */
            double worst; /* greatest off-diagonal dot product */
            int lim; /* index of last vec to check against */
            int screenlim; /* value of lim that will fit on screen */
            int option; /* which option to use */

            /* The T/F argument in the conditionals is just a convenient option: */

            screenlim = 20;
            option = 3;

            /* Check orthogonality over whole set. */
            if (option == 1)
            {
                Console.WriteLine("Orthogonality check:");
                for (i = 1; i <= dim; i++)
                {
                    Console.Write("{0:d})", i);
                    for (j = 1; j <= i; j++)
                    {
                        prod = dot(mat[i], 1, n, mat[j]);
                        /* printf(" %g ",prod); */
                        /* printf(" %4.2e ",prod); */
                        /* printf(" %4.2e ",fabs(prod)); */
                        Console.Write(" {0:d}", -(int) Math.Log10(prod));
                    }

                    Console.WriteLine();
                }
            }

            if (option == 2)
            {
                Console.Write("Frobenius orthogonality measure:");
                measure = 0;
                for (i = 1; i <= dim; i++)
                {
                    for (j = i; j <= dim; j++)
                    {
                        prod = dot(mat[i], 1, n, mat[j]);
                        if (i == j)
                        {
                            measure += Math.Abs(1.0 - prod);
                        }
                        else
                        {
                            measure += 2.0 * Math.Abs(prod);
                        }
                    }
                }

                Console.WriteLine("{0:g} ", measure);
            }

            /* Check orthogonality against last vector. Allows you to build up orthogonality
               matrix much faster if previous columns stay the same when add a new column,
               but may interact with other debug output to give a confusing presentation. */
            if (option == 3)
            {
                Console.Write("{0:d}) ", dim);
                lim = Math.Min(dim, screenlim);
                worst = 0;
                for (i = 1; i <= dim; i++)
                {
                    prod = dot(mat[i], 1, n, mat[dim]);
                    if (i <= lim)
                    {
                        Console.Write(" {0:d}", -(int) Math.Log10(Math.Abs(prod)));
                    }

                    if ((i != dim) && (Math.Abs(prod) > Math.Abs(worst)))
                    {
                        worst = prod;
                    }
                }

                Console.WriteLine(" worst {0:e}", worst);
            }
        }

/* Check orthogonality of vector set */
        public static void checkorth_float(float** mat, int n, int dim)
        {
            int i, j; /* loop idices */
            double measure; /* Froebenius norm */
            double prod; /* value of dot product */
            double worst; /* greatest off-diagonal dot product */
            int lim; /* index of last vec to check against */
            int screenlim; /* value of lim that will fit on screen */
            int option; /* which option to use */

            /* The T/F argument in the conditionals is just a convenient option: */

            screenlim = 20;
            option = 3;

            /* Check orthogonality over whole set. */
            if (option == 1)
            {
                Console.WriteLine("Orthogonality check:");
                for (i = 1; i <= dim; i++)
                {
                    Console.Write("{0:d})", i);
                    for (j = 1; j <= i; j++)
                    {
                        prod = dot_float(mat[i], 1, n, mat[j]);
                        /* printf(" %g ",prod); */
                        /* printf(" %4.2e ",prod); */
                        /* printf(" %4.2e ",fabs(prod)); */
                        Console.Write(" {0:d}", -(int) Math.Log10(prod));
                    }

                    Console.WriteLine();
                }
            }

            if (option == 2)
            {
                Console.Write("Frobenius orthogonality measure:");
                measure = 0;
                for (i = 1; i <= dim; i++)
                {
                    for (j = i; j <= dim; j++)
                    {
                        prod = dot_float(mat[i], 1, n, mat[j]);
                        if (i == j)
                        {
                            measure += Math.Abs(1.0 - prod);
                        }
                        else
                        {
                            measure += 2.0 * Math.Abs(prod);
                        }
                    }
                }

                Console.WriteLine("{0:g} ", measure);
            }

            /* Check orthogonality against last vector. Allows you to build up orthogonality
               matrix much faster if previous columns stay the same when add a new column,
               but may interact with other debug output to give a confusing presentation. */
            if (option == 3)
            {
                Console.Write("{0:d}) ", dim);
                lim = Math.Min(dim, screenlim);
                worst = 0;
                for (i = 1; i <= dim; i++)
                {
                    prod = dot_float(mat[i], 1, n, mat[dim]);
                    if (i <= lim)
                    {
                        Console.Write(" {0:d}", -(int) Math.Log10(Math.Abs(prod)));
                    }

                    if ((i != dim) && (Math.Abs(prod) > Math.Abs(worst)))
                    {
                        worst = prod;
                    }
                }

                Console.WriteLine(" worst {0:e}", worst);
            }
        }

        public static void sorthog(double* vec, /* vector to be orthogonalized */
            int n, /* length of the columns of orth */
            orthlink** solist, /* set of vecs to orth. against */
            int ngood /* number of vecs in solist */
        )
        {
            double alpha;
            double* dir;
            int i;

            for (i = 1; i <= ngood; i++)
            {
                dir = (solist[i])->vec;
                alpha = -dot(vec, 1, n, dir) / dot(dir, 1, n, dir);
                scadd(vec, 1, n, alpha, dir);
            }
        }

        public static void sorthog_float(float* vec, /* vector to be orthogonalized */
            int n, /* length of the columns of orth */
            orthlink_float** solist, /* set of vecs to orth. against */
            int ngood /* number of vecs in solist */
        )
        {
            float alpha;
            float* dir;
            int i;

            for (i = 1; i <= ngood; i++)
            {
                dir = (solist[i])->vec;
                alpha = (float) (-dot_float(vec, 1, n, dir) / dot_float(dir, 1, n, dir));
                scadd_float(vec, 1, n, alpha, dir);
            }
        }

        /* Allocate space for new orthlink, double version. */
        public static orthlink* makeorthlnk()
        {
            return (orthlink*) Marshal.AllocHGlobal(sizeof(orthlink));
        }

        /* Allocate space for new orthlink, float version. */
        public static orthlink_float* makeorthlnk_float()
        {
            return (orthlink_float*) Marshal.AllocHGlobal(sizeof(orthlink_float));
        }

        public static void orthogonalize(double* vec, /* vector to be orthogonalized */
            int n, /* length of the columns of orth */
            orthlink* orthlist /* set of vectors to orthogonalize against */
        )
        {
            orthlink* curlnk = orthlist;
            while (curlnk != (orthlink*) IntPtr.Zero)
            {
                orthogvec(vec, 1, n, curlnk->vec);
                curlnk = curlnk->pntr;
            }
        }

        public static void orthogvec(
            double* vec1, /* vector to be orthogonalized */
            int beg, int end, /* start and stop range for vector */
            double* vec2 /* vector to be orthogonalized against */
        )
        {
            double alpha = -dot(vec1, beg, end, vec2) / dot(vec2, beg, end, vec2);
            scadd(vec1, beg, end, alpha, vec2);
        }

        public static void orthogvec_float(
            float* vec1, /* vector to be orthogonalized */
            int beg, int end, /* start and stop range for vector */
            float* vec2 /* vector to be orthogonalized against */
        )
        {
            float alpha = (float) (-dot_float(vec1, beg, end, vec2) / dot_float(vec2, beg, end, vec2));
            scadd_float(vec1, beg, end, alpha, vec2);
        }

        /* Orthogonalize a double vector to all one's */
        public static void orthog1(double* x, int beg, int end)
        {
            int i;

            var len = end - beg + 1;
            var sum = 0.0;
            var pntr = x + beg;
            for (i = len; i != 0; i--)
            {
                sum += *pntr++;
            }

            sum /= len;
            pntr = x + beg;
            for (i = len; i != 0; i--)
            {
                *pntr++ -= sum;
            }
        }

/* Orthogonalize a float vector to all one's */
        public static void orthog1_float(float* x, int beg, int end)
        {
            var len = end - beg + 1;
            var sum = 0.0f;
            var pntr = x + beg;
            for (var i = len; i != 0; i--)
            {
                sum += *pntr++;
            }

            sum /= len;
            pntr = x + beg;
            for (var i = len; i != 0; i--)
            {
                *pntr++ -= sum;
            }
        }
    }
}
