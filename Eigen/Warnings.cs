using System;
using System.Diagnostics;
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
                    Trace.Write($"  {number:f}");
                }
                else
                {
                    Trace.Write($"  {number:g}");
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
                        Trace.WriteLine("Note about warnings: in partition convergence monitoring mode.");
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
                        Trace.WriteLine($"Lanczos itns. = {j:d}");
                    }

                    Trace.WriteLine("          lambda                Ares est.              Ares          index");
                    for (i = 1; i <= d; i++)
                    {
                        Trace.Write($"{i:d}.");
                        doubleout(lambda[i], 1);
                        doubleout(bound[i], 1);
                        doubleout(Ares[i], 1);
                        Trace.WriteLine($"   {index[i]:d}");
                    }

                    Trace.WriteLine("");
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
                        Trace.WriteLine("WARNING: Maximum number of Lanczos iterations reached.");
                    }

                    if (warning2 && !warning3)
                    {
                        Trace.WriteLine($"WARNING: Minor loss of orthogonality (Ares/est. > {WARNING_ORTHTOL:g}).");
                    }

                    if (warning3)
                    {
                        Trace.WriteLine($"WARNING: Substantial loss of orthogonality (Ares/est. > {WARNING_MISTOL:g}).");
                    }

                    if (warning1)
                    {
                        Trace.WriteLine($"WARNING: Eigen pair tolerance ({eigtol:g}) not achieved.");
                    }
                }

                if (WARNING_EVECS > 1)
                {
                    if (warning1 || warning2 || warning3)
                    {
                        if (DEBUG_EVECS <= 0)
                        {
                            Trace.WriteLine("          lambda                Ares est.              Ares          index");
                            for (i = 1; i <= d; i++)
                            {
                                Trace.Write($"{i:d}.");
                                doubleout(lambda[i], 1);
                                doubleout(bound[i], 1);
                                doubleout(Ares[i], 1);
                                Trace.WriteLine($"   {index[i]:d}");
                            }
                        }

                        /* otherwise gets printed above */
                    }
                }

                if (warning1 || warning2 || warning3 || WARNING_EVECS > 2)
                {
                    if (Sres_max > SRESTOL)
                    {
                        Trace.WriteLine($"WARNING: Maximum eigen residual of T ({Sres_max:g}) exceeds SRESTOL.");
                    }
                }

                if (WARNING_EVECS > 2)
                {
                    if (SRES_SWITCHES > 0)
                    {
                        Trace.WriteLine($"WARNING: Switched routine for computing evec of T {SRES_SWITCHES:d} times.");
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
                    Trace.WriteLine("ERROR: Sorry, out-of-bounds eigenvalue indicates serious breakdown.");
                    Trace.WriteLine("       Try different parameters or another eigensolver.");
                    if (pass == 2)
                    {
                        throw new InvalidOperationException();
                    }
                }

            } /* Pass loop */
        }
    }
}
