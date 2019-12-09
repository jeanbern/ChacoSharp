using System;
using System.Diagnostics;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Utilities.MkVec;
using static ChacoSharp.Utilities.CpVec;
using static ChacoSharp.Utilities.ShellSort;
using static ChacoSharp.Eigen.Ql;
using static ChacoSharp.Eigen.Bisect;

namespace ChacoSharp.Eigen
{
    public static unsafe class GetRitzVals
    {
        /// <summary>
        /// Finds needed eigenvalues of tridiagonal T using either the QL algorithm or Sturm sequence bisection,
        /// whichever is predicted to be faster based on a simple complexity model.
        /// If one fails (which is rare), the other is tried.
        /// The return value is 0 if one of the routines succeeds.
        /// If they both fail, the return value is 1, and Lanczos should compute the best approximation it can based on previous iterations.
        /// </summary>
        /// <param name="alpha">vector of Lanczos scalars</param>
        /// <param name="beta">vector of Lanczos scalars</param>
        /// <param name="j">number of Lanczos iterations taken</param>
        /// <param name="Anorm">Gershgorin estimate</param>
        /// <param name="workj">work vector for Sturm sequence</param>
        /// <param name="ritz">array holding evals</param>
        /// <param name="d">problem dimension = num. eigenpairs needed</param>
        /// <param name="left_goodlim">number of ritz pairs checked on left end</param>
        /// <param name="right_goodlim">number of ritz pairs checked on right end</param>
        /// <param name="eigtol">tolerance on eigenpair</param>
        /// <param name="bis_safety">bisection tolerance function divisor</param>
        /// <returns>True if an error was encountered.</returns>
        public static bool get_ritzvals(double* alpha,
            double* beta,
            int j,
            double Anorm,
            double* workj,
            double* ritz,
            int d,
            int left_goodlim,
            int right_goodlim,
            double eigtol,
            double bis_safety
        )
        {
            int nvals_left; /* numb. evals to find on left end of spectrum */
            int nvals_right; /* numb. evals to find on right end of spectrum */
            double bisection_tol; /* width of interval bisection should converge to */
            int pred_steps; /* predicts # of required bisection steps per eval */
            int tot_pred_steps; /* predicts total # of required bisection steps */
            double* ritz_sav = null; /* copy of ritzvals for debugging */
            int bisect_flag; /* return status of bisect() */
            bool ql_flag; /* return status of ql() */
            bool local_debug; /* whether to check bisection results with ql */

            /* Determine number of ritzvals to find on left and right ends */
            nvals_left = Math.Max(d, left_goodlim);
            nvals_right = Math.Min(j - nvals_left, right_goodlim);

            /* Estimate work for bisection vs. ql assuming bisection takes 5j flops per
               step, ql takes 30j^2 flops per call. (Ignore sorts, copies, addressing.) */

            bisection_tol = eigtol * eigtol / bis_safety;
            pred_steps = (int) (Math.Log10(Anorm / bisection_tol) / Math.Log10(2.0)) + 1;
            tot_pred_steps = (nvals_left + nvals_right) * pred_steps;

            bisect_flag = 0;
            ql_flag = false;

            if (5 * tot_pred_steps < 30 * j)
            {
                if (DEBUG_EVECS > 2)
                {
                    Trace.WriteLine("  tridiagonal solver: bisection");
                }

                /* Set local_debug = TRUE for a table checking bisection against QL. */
                local_debug = false;
                if (local_debug)
                {
                    ritz_sav = mkvec(1, j);
                    cpvec(ritz_sav, 1, j, alpha);
                    cpvec(workj, 0, j, beta);
                    ql_flag = ql(ritz_sav, workj, j);
                    if (ql_flag)
                    {
                        throw new InvalidOperationException("Aborting debugging procedure in get_ritzvals().");
                    }

                    shell_sort(j, &ritz_sav[1]);
                }

                bisect_flag = bisect(alpha, beta, j, Anorm, workj, ritz, nvals_left, nvals_right, bisection_tol,
                    ritz_sav, pred_steps + 10);

                if (local_debug)
                {
                    frvec(ritz_sav, 1);
                }
            }
            else
            {
                if (DEBUG_EVECS > 2)
                {
                    Trace.WriteLine("  tridiagonal solver: ql");
                }

                cpvec(ritz, 1, j, alpha);
                cpvec(workj, 0, j, beta);
                ql_flag = ql(ritz, workj, j);
                shell_sort(j, &ritz[1]);
            }

            if (bisect_flag != 0 && !ql_flag)
            {
                if (DEBUG_EVECS > 0 || WARNING_EVECS > 0)
                {
                    Trace.WriteLine("WARNING: Sturm bisection of T failed; switching to QL.");
                }

                if (DEBUG_EVECS > 1 || WARNING_EVECS > 1)
                {
                    if (bisect_flag == 1)
                    {
                        Console.WriteLine("         - failure detected in sturmcnt().");
                    }

                    if (bisect_flag == 2)
                    {
                        Trace.WriteLine("         - maximum number of bisection steps reached");
                    }
                }

                cpvec(ritz, 1, j, alpha);
                cpvec(workj, 0, j, beta);
                ql_flag = ql(ritz, workj, j);
                shell_sort(j, &ritz[1]);
            }

            if (ql_flag && bisect_flag == 0)
            {
                if (DEBUG_EVECS > 0 || WARNING_EVECS > 0)
                {
                    Trace.WriteLine("WARNING: QL failed for T; switching to Sturm bisection.");
                }

                bisect_flag = bisect(alpha, beta, j, Anorm, workj, ritz, nvals_left, nvals_right, bisection_tol,
                    ritz_sav, pred_steps + 3);
            }

            if (bisect_flag != 0 && ql_flag)
            {
                if (DEBUG_EVECS > 0 || WARNING_EVECS > 0)
                {
                    return true; /* can't recover; bail out with error code */
                }
            }

            return false; /* ... things seem ok. */
        }
    }
}
