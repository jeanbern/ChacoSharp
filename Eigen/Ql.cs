using System;

namespace ChacoSharp.Eigen
{
    public static unsafe class Ql
    {
        private static double SIGN(double a, double b)
        {
            return b < 0.0 ? -Math.Abs(a) : Math.Abs(a);
        }

        /*  */
        /// <summary>
        /// Eigensolution of real symmetric tridiagonal matrix using the algorithm of Numerical Recipes p. 380.
        /// Removed eigenvector calculation and added return codes:
        ///     1 if maximum number of iterations is exceeded, 0 otherwise.
        /// NOTE CAREFULLY: the vector e is used as workspace, the eigenvals are returned in the vector d.
        /// </summary>
        /// <returns>True if there was an error.</returns>
        public static  bool ql(double* d, double* e, int n)

        {
            int    m, l, iter, i;
            double s, r, p, g, f, dd, c, b;

            e[n] = 0.0;

            for (l = 1; l <= n; l++) {
                iter = 0;
                do {
                    for (m = l; m <= n - 1; m++) {
                        dd = Math.Abs(d[m]) + Math.Abs(d[m + 1]);
                        if (Math.Abs(e[m]) + dd == dd) {
                            break;
                        }
                    }
                    if (m != l) {
                        if (iter++ == 50) {
                            return true;
                            /* ... not converging; bail out with error code. */
                        }
                        g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                        r = Math.Abs((g * g) + 1.0);
                        g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
                        s = c = 1.0;
                        p     = 0.0;
                        for (i = m - 1; i >= l; i--) {
                            f = s * e[i];
                            b = c * e[i];
                            if (Math.Abs(f) >= Math.Abs(g)) {
                                c        = g / f;
                                r        = Math.Sqrt((c * c) + 1.0);
                                e[i + 1] = f * r;
                                c *= (s = 1.0 / r);
                            }
                            else {
                                s        = f / g;
                                r        = Math.Sqrt((s * s) + 1.0);
                                e[i + 1] = g * r;
                                s *= (c = 1.0 / r);
                            }
                            g        = d[i + 1] - p;
                            r        = (d[i] - g) * s + 2.0 * c * b;
                            p        = s * r;
                            d[i + 1] = g + p;
                            g        = c * r - b;
                        }
                        d[l] = d[l] - p;
                        e[l] = g;
                        e[m] = 0.0;
                    }
                } while (m != l);
            }
            return false; /* ... things seem ok */
        }
    }
}
