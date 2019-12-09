using System;
using System.Diagnostics;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Utilities.CpVec;
using static ChacoSharp.Utilities.Norm;
using static ChacoSharp.Eigen.TevecHelper;
using static ChacoSharp.Eigen.TriSolve;

namespace ChacoSharp.Eigen
{
    public static unsafe class GetExtVal
    {
        /* Finds first extended eigenpair of system corresponding to
   tridiagonal T using using Rafael's bisection technique. */

        public static void get_extval(double* alpha, /* j-vector of Lanczos scalars (using elements 1 to j) */
            double* beta, /* (j+1)-vector of " " (has 0 element but using 1 to j-1) */
            int j, /* number of Lanczos iterations taken */
            double ritzval, /* Ritz value */
            double* s, /* Ritz vector (length n, re-computed in this routine) */
            double eigtol, /* tolerance on eigenpair */
            double wnorm_g, /* W-norm of n-vector g, the rhs in the extended eig. problem */
            double sigma, /* the norm constraint on the extended eigenvector */
            double* extval, /* the extended eigenvalue this routine computes */
            double* v, /* the j-vector solving the extended eig problem in T */
            double* work1, /* j-vector of workspace */
            double* work2 /* j-vector of workspace */
        )
        {
            double lambda_low; /* lower bound on extended eval */
            double lambda_high; /* upper bound on extended eval */
            double tol; /* bisection tolerance */
            double norm_v; /* norm of the extended T eigenvector v */
            double lambda; /* the parameter that iterates to extval */
            int cnt; /* debug iteration counter */
            double diff; /* distance between lambda limits */

            /* Compute the Ritz vector */
            Tevec(alpha, beta - 1, j, ritzval, s);

            /* Shouldn't happen, but just in case ... */
            if (wnorm_g == 0.0)
            {
                *extval = ritzval;
                cpvec(v, 1, j, s);
                if (DEBUG_EVECS > 0)
                {
                    Trace.WriteLine("Degenerate extended eigenvector problem (g = 0).");
                }

                return;
                /* ... not really an extended eigenproblem; just return Ritz pair */
            }

            /* Set up the bisection parameters */
            lambda_low = ritzval - wnorm_g / sigma;
            lambda_high = ritzval - (wnorm_g / sigma) * s[1];
            lambda = 0.5 * (lambda_low + lambda_high);
            tol = eigtol * eigtol * (1 + Math.Abs(lambda_low) + Math.Abs(lambda_high));

            if (DEBUG_EVECS > 2)
            {
                Trace.WriteLine("Computing extended eigenpairs of T");
                Trace.WriteLine($"  target norm_v (= sigma) {sigma:g}");
                Trace.WriteLine($"  bisection tolerance {tol:g}");
            }

            if (DEBUG_EVECS > 3)
            {
                Trace.WriteLine("  lambda iterates to the extended eigenvalue");
                Trace.WriteLine("         lambda_low           lambda            lambda_high      norm_v");
            }

            /* Bisection loop - iterate until norm constraint is satisfied */
            cnt = 1;
            diff = 2 * tol;
            while (diff > tol)
            {
                lambda = 0.5 * (lambda_low + lambda_high);
                tri_solve(alpha, beta, j, lambda, v, wnorm_g, work1, work2);
                norm_v = ch_norm(v, 1, j);
                if (DEBUG_EVECS > 3)
                {
                    Trace.WriteLine($"{cnt++:i}   {lambda_low:f}  {lambda:f}  {lambda_high:f}  {norm_v:g}");
                }

                if (norm_v <= sigma)
                {
                    lambda_low = lambda;
                }

                if (norm_v >= sigma)
                {
                    lambda_high = lambda;
                }

                diff = lambda_high - lambda_low;
            }

            /* Return the extended eigenvalue (eigvec is automatically returned) */
            *extval = lambda;
        }

    }
}
