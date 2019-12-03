using System;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Eigen.CheckEig;

namespace ChacoSharp.Eigen
{
    public static unsafe class Warnings
    {
        public static void doubleout(double number, int mode)
/* argument to print */
/* currently just one */
        {
            if (mode == 1)
            {
                if (Math.Abs(number) < 100)
                {
                    Console.Write("  {0:f}", number);
                }
                else
                {
                    Console.Write("  {0:g}", number);
                }
            }
        }

        /* Post various warnings about computation  */
        public static void warnings(double* workn, /* work vector (1..n) */
            vtx_data** A, /* graph */
            double*[] y, /* eigenvectors */
            int n, /* number of vtxs */
            double[] lambda, /* ritz approximation to eigenvals of A */
            double* vwsqrt, /* square roots of vertex weights */
            double* Ares, /* how well Lanczos calc. eigpair lambda,y */
            double[] bound, /* on ritz pair approximations to eig pairs of A */
            int* index, /* the Ritz index of an eigenpair */
            int d, /* problem dimension = number of eigvecs to find */
            int j, /* number of Lanczos iterations used */
            int maxj, /* maximum number of Lanczos iterations */
            double Sres_max, /* Max value of Sres */
            double eigtol, /* tolerance on eigenvectors */
            double* u, /* Lanczos vector; here used as workspace */
            double Anorm /* Gershgorin bound on eigenvalue */
        )
        {
            bool warning1 = false; /* warning1 cond. (eigtol not achieved) true? */
            bool warning2 = false; /* warning2 cond. (premature orth. loss) true? */
            bool warning3 = false; /* warning3 cond. (suspected misconvergence) true? */
            int i; /* loop index */
            bool hosed; /* flag for serious Lanczos problems */
            int pass; /* which time through we are on */

            hosed = false;
            for (pass = 1; pass <= 2; pass++)
            {
                if (DEBUG_EVECS > 0 || WARNING_EVECS > 0)
                {
                    if (LANCZOS_CONVERGENCE_MODE == 1)
                    {
                        Console.WriteLine("Note about warnings: in partition convergence monitoring mode.");
                    }

                    for (i = 1; i <= d; i++)
                    {
                        Ares[i] = checkeig(workn, A, y[i], n, lambda[i], vwsqrt, u);
                    }
                }

                if (DEBUG_EVECS > 0)
                {
                    if (pass == 1)
                    {
                        Console.WriteLine("Lanczos itns. = {0:d}", j);
                    }

                    Console.WriteLine("          lambda                Ares est.              Ares          index");
                    for (i = 1; i <= d; i++)
                    {
                        Console.Write("{0:d}.", i);
                        doubleout(lambda[i], 1);
                        doubleout(bound[i], 1);
                        doubleout(Ares[i], 1);
                        Console.WriteLine("   {0:d}", index[i]);
                    }

                    Console.WriteLine();
                }

                if (WARNING_EVECS > 0)
                {
                    warning1 = false;
                    warning2 = false;
                    warning3 = A != null;
                    for (i = 1; i <= d; i++)
                    {
                        if (Ares[i] > eigtol)
                        {
                            warning1 = true;
                        }

                        if (Ares[i] > WARNING_ORTHTOL * bound[i] && Ares[i] > .01 * eigtol)
                        {
                            warning2 = true;
                        }

                        if (Ares[i] > WARNING_MISTOL * bound[i] && Ares[i] > .01 * eigtol)
                        {
                            warning3 = true;
                        }
                    }

                    if (j == maxj)
                    {
                        Console.WriteLine("WARNING: Maximum number of Lanczos iterations reached.");
                    }

                    if (warning2 && !warning3)
                    {
                        Console.WriteLine("WARNING: Minor loss of orthogonality (Ares/est. > {0:g}).",
                            WARNING_ORTHTOL);
                    }

                    if (warning3)
                    {
                        Console.WriteLine("WARNING: Substantial loss of orthogonality (Ares/est. > {0:g}).",
                            WARNING_MISTOL);
                    }

                    if (warning1)
                    {
                        Console.WriteLine("WARNING: Eigen pair tolerance ({0:g}) not achieved.", eigtol);
                    }
                }

                if (WARNING_EVECS > 1)
                {
                    if (warning1 || warning2 || warning3)
                    {
                        if (DEBUG_EVECS <= 0)
                        {
                            Console.WriteLine("          lambda                Ares est.              Ares          index");
                            for (i = 1; i <= d; i++)
                            {
                                Console.Write("{0:d}.", i);
                                doubleout(lambda[i], 1);
                                doubleout(bound[i], 1);
                                doubleout(Ares[i], 1);
                                Console.WriteLine("   {0:d}", index[i]);
                            }
                        }

                        /* otherwise gets printed above */
                    }
                }

                if (warning1 || warning2 || warning3 || WARNING_EVECS > 2)
                {
                    if (Sres_max > SRESTOL)
                    {
                        Console.WriteLine("WARNING: Maximum eigen residual of T ({0:g}) exceeds SRESTOL.", Sres_max);
                    }
                }

                if (WARNING_EVECS > 2)
                {
                    if (SRES_SWITCHES > 0)
                    {
                        Console.WriteLine("WARNING: Switched routine for computing evec of T {0:d} times.", SRES_SWITCHES);
                        SRES_SWITCHES = 0;
                    }
                }

                /* Put the best face on things ... */
                for (i = 1; i <= d; i++)
                {
                    if (lambda[i] < 0 || lambda[i] > Anorm + eigtol)
                    {
                        hosed = true;
                    }
                }

                if (hosed)
                {
                    Console.WriteLine("ERROR: Sorry, out-of-bounds eigenvalue indicates serious breakdown.");
                    Console.WriteLine("       Try different parameters or another eigensolver.");
                    if (pass == 2)
                    {
                        throw new InvalidOperationException();
                    }
                }

            } /* Pass loop */
        }
    }
}
