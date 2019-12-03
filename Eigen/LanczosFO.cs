using System;
using System.Runtime.InteropServices;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Eigen.Orthogonalization;
using static ChacoSharp.Eigen.Scan;
using static ChacoSharp.Eigen.MakeEigVec;
using static ChacoSharp.Eigen.Splarax;
using static ChacoSharp.Eigen.LanczosTime;
using static ChacoSharp.Eigen.GetRitzVals;
using static ChacoSharp.Eigen.Warnings;
using static ChacoSharp.Utilities.MkVec;
using static ChacoSharp.Utilities.VecRan;
using static ChacoSharp.Utilities.Dot;
using static ChacoSharp.Utilities.Update;
using static ChacoSharp.Utilities.SetVec;
using static ChacoSharp.Utilities.Norm;
using static ChacoSharp.Utilities.VecScale;
using static ChacoSharp.Symmlq.SymmLq;

namespace ChacoSharp.Eigen
{
    public static unsafe class LanczosFO
    {
        #region longComment

        /* Lanczos iteration with FULL orthogonalization.
   Works in standard (version 1) or inverse operator (version 2) mode. */
/* Finds the lowest (zero) eigenvalue, but does not print it
   or compute corresponding eigenvector. */
/* Does NOT find distinct eigenvectors corresponding to multiple
   eigenvectors and hence will not lead to good partitioners
   for symmetric graphs. Will be useful for graphs known
   (believed) not to have multiple eigenvectors, e.g. graphs
   in which a random set of edges have been added or perturbed. */
/* May fail on small graphs with high multiplicities (e.g. k5) if Ritz
   pairs converge before there have been as many iterations as the number
   of eigenvalues sought. This is rare and a different random number
   seed will generally alleviate the problem. */
/* Convergence check uses Paige bji estimate over the whole
   spectrum of T. This is a lot of work, but we are trying to be
   extra safe. Since we are orthogonalizing fully, we assume the
   bji esitmates are very good and don't provide a contingency for
   when they don't match the residuals. */
/* A lot of the time in this routine (say half) is spent in ql finding
   the evals of T on each iteration. This could be reduced by only using
   ql say every 10 steps. This might require that the orthogonalization
   be done with Householder (rather than Gram-Scmidt as currently)
   to avoid a problem with zero beta values causing subsequent breakdown.
   But this routine is really intended to be used for smallish problems
   where the Lanczos runs will be int, e.g. when we are using the inverse
   operator method. In the inverse operator case, it's really important to
   stop quickly to avoid additional back solves. If the ql work is really
   a problem then we should be using a selective othogonalization algorithm.
   This routine provides a convenient reference point for how well those
   routines are functioning since it has the same basic structure but just
   does more orthogonalizing. */
/* The algorithm orthogonalizes the starting vector and the residual vectors
   against the vector of all ones since we know that is the null space of
   the Laplacian. This generally saves net time because Lanczos tends to
   converge faster. */
/* Replaced call to ql() with call to get_ritzvals(). This starts with ql
   bisection, whichever is predicted to be faster based on a simple complexity
   model. If that fails it switches to the other. This routine should be
   safer and faster than straight ql (which does occasionally fail). */
/* NOTE: This routine indexes beta (and workj) from 1 to maxj+1, whereas
   selective orthogonalization indexes them from 0 to maxj. */

/* Comments for Lanczos with inverted operator: */
/* Used Symmlq for the back solve since already maintaining that code for
   the RQI/Symmlq multilevel method. Straight CG would only be marginally
   faster. */
/* The orthogonalization against the vector of all ones in Symmlq is not
   as efficient as possible - could in principle use orthog1 method, but
   don't want to disrupt the RQI/Symmlq code. Also, probably would only
   need to orthogonalize out a given mode periodically. But we want to be
   extra robust numericlly. */

        #endregion longComment

        public static void lanczos_FO(vtx_data** A, /* graph data structure */
            int n, /* number of rows/columns in matrix */
            int d, /* problem dimension = # evecs to find */
            double*[] y, /* columns of y are eigenvectors of A  */
            double[] lambda, /* ritz approximation to eigenvals of A */
            double[] bound, /* on ritz pair approximations to eig pairs of A */
            double eigtol, /* tolerance on eigenvectors */
            double* vwsqrt, /* square root of vertex weights */
            double maxdeg, /* maximum degree of graph */
            int version /* 1 = standard mode, 2 = inverse operator mode */
        )

        {
            int i, j; /* indices */
            int maxj; /* maximum number of Lanczos iterations */
            double* u;
            double* r; /* Lanczos vectors */
            double* Aq; /* sparse matrix-vector product vector */
            double* alpha;
            double* beta; /* the Lanczos scalars from each step */
            double* ritz; /* copy of alpha for tqli */
            double* workj; /* work vector (eg. for tqli) */
            double* workn; /* work vector (eg. for checkeig) */
            double* s; /* eigenvector of T */
            double** q; /* columns of q = Lanczos basis vectors */
            double* bj; /* beta(j)*(last element of evecs of T) */
            double bis_safety; /* real safety factor for bisection algorithm */
            double Sres; /* how well TevecHelper calculated eigvecs */
            double Sres_max; /* Maximum value of Sres */
            bool inc_bis_safety; /* need to increase bisection safety */
            double* Ares; /* how well Lanczos calculated each eigpair */
            double* inv_lambda; /* eigenvalues of inverse operator */
            int* index; /* the Ritz index of an eigenpair */
            orthlink* orthlist = null; /* vectors to orthogonalize against in Lanczos */
            orthlink* orthlist2 = null; /* vectors to orthogonalize against in Symmlq */
            orthlink* temp; /* for expanding orthogonalization list */
            double* ritzvec = null; /* ritz vector for current iteration */
            double* zeros = null; /* vector of all zeros */
            double* ones = null; /* vector of all ones */
            scanlink* scanlist; /* list of fields for min ritz vals */
            scanlink* curlnk; /* for traversing the scanlist */
            double bji_tol; /* tol on bji estimate of A e-residual */
            bool converged; /* has the iteration converged? */
            double time; /* current clock time */
            double shift, rtol; /* symmlq input */
            bool precon = false, goodb = false; /* symmlq input */
            int nout; /* symmlq input */
            bool checka = false;
            int intlim; /* symmlq input */
            double anorm, acond; /* symmlq output */
            double rnorm, ynorm; /* symmlq output */
            int istop, itn; /* symmlq output */
            double macheps; /* machine precision calculated by symmlq */
            double normxlim; /* a stopping criteria for symmlq */
            int itnmin; /* enforce minimum number of iterations */
            int symmlqitns = 0; /* # symmlq itns */
            double* wv1 = null;
            double* wv2 = null;
            double* wv3 = null; /* Symmlq work space */
            double* wv4 = null;
            double* wv5 = null;
            double* wv6 = null; /* Symmlq work space */
            int long_n = 0; /* long int copy of n for symmlq */
            bool ritzval_flag; /* status flag for ql() */
            double Anorm; /* Norm estimate of the Laplacian matrix */
            int left, right; /* ranges on the search for ritzvals */
            bool memory_ok; /* true as long as don't run out of memory */

            if (DEBUG_TRACE)
            {
                Console.WriteLine("<Entering lanczos_FO>");
            }

            if (DEBUG_EVECS > 0)
            {
                if (version == 1)
                {
                    Console.WriteLine("Full orthogonalization Lanczos, matrix size = {0:d}", n);
                }
                else
                {
                    Console.WriteLine("Full orthogonalization Lanczos, inverted operator, matrix size = {0:d}", n);
                }
            }

            /* Initialize time. */
            time = lanc_seconds();

            if (n < d + 1)
            {
                throw new ArgumentOutOfRangeException(nameof(n), "ERROR: System too small for number of eigenvalues requested.");
                /* d+1 since don't use zero eigenvalue pair */
            }

            /* Allocate Lanczos space. */
            maxj = LANCZOS_MAXITNS;
            u = mkvec(1, n);
            r = mkvec(1, n);
            Aq = mkvec(1, n);
            ritzvec = mkvec(1, n);
            zeros = mkvec(1, n);
            setvec(zeros, 1, n, 0.0);
            workn = mkvec(1, n);
            Ares = mkvec(1, d);
            inv_lambda = mkvec(1, d);
            index = (int*) Marshal.AllocHGlobal((d + 1) * sizeof(int));
            alpha = mkvec(1, maxj);
            beta = mkvec(1, maxj + 1);
            ritz = mkvec(1, maxj);
            s = mkvec(1, maxj);
            bj = mkvec(1, maxj);
            workj = mkvec(1, maxj + 1);
            q = (double**) Marshal.AllocHGlobal((maxj + 1) * sizeof(double*));
            scanlist = mkscanlist(d);

            if (version == 2)
            {
                /* Allocate Symmlq space all in one chunk. */
                wv1 = (double*) Marshal.AllocHGlobal(6 * (n + 1) * sizeof(double));
                wv2 = &wv1[(n + 1)];
                wv3 = &wv1[2 * (n + 1)];
                wv4 = &wv1[3 * (n + 1)];
                wv5 = &wv1[4 * (n + 1)];
                wv6 = &wv1[5 * (n + 1)];

                /* Set invariant symmlq parameters */
                precon = false; /* false until we figure out a good way */
                goodb = false; /* should be false for this application */
                checka = false; /* if don't know by now, too bad */
                intlim = n; /* set to enforce a maximum number of Symmlq itns */
                itnmin = 0; /* set to enforce a minimum number of Symmlq itns */
                shift = 0.0; /* since just solving rather than doing RQI */
                symmlqitns = 0; /* total number of Symmlq iterations */
                nout = 0; /* Effectively disabled - see notes in symmlq.f */
                rtol = 1.0e-5; /* requested residual tolerance */
                normxlim = double.MaxValue; /* Effectively disables ||x|| termination criterion */
                long_n = n; /* copy to long for linting */
            }

            /* Initialize. */
            vecran(r, 1, n);
            if (vwsqrt == null)
            {
                /* whack one's direction from initial vector */
                orthog1(r, 1, n);

                /* list the ones direction for later use in Symmlq */
                if (version == 2)
                {
                    orthlist2 = makeorthlnk();
                    ones = mkvec(1, n);
                    setvec(ones, 1, n, 1.0);
                    orthlist2->vec = ones;
                    orthlist2->pntr = null;
                }
            }
            else
            {
                /* whack vwsqrt direction from initial vector */
                orthogvec(r, 1, n, vwsqrt);

                if (version == 2)
                {
                    /* list the vwsqrt direction for later use in Symmlq */
                    orthlist2 = makeorthlnk();
                    orthlist2->vec = vwsqrt;
                    orthlist2->pntr = null;
                }
            }

            beta[1] = ch_norm(r, 1, n);
            q[0] = zeros;
            bji_tol = eigtol;
            orthlist = null;
            Sres_max = 0.0;
            Anorm = 2 * maxdeg; /* Gershgorin estimate for ||A|| */
            bis_safety = BISECTION_SAFETY;
            inc_bis_safety = false;
            init_time += lanc_seconds() - time;

            /* Main Lanczos loop. */
            j = 1;
            converged = false;
            memory_ok = true;
            while ((j <= maxj) && (converged == false) && memory_ok)
            {
                time = lanc_seconds();

                /* Allocate next Lanczos vector. If fail, back up one step and compute approx. eigvec. */
                q[j] = mkvec_ret(1, n);
                if (q[j] == null)
                {
                    memory_ok = false;
                    if (DEBUG_EVECS > 0 || WARNING_EVECS > 0)
                    {
                        Console.WriteLine("WARNING: Lanczos out of memory; computing best approximation available.");
                    }

                    if (j <= 2)
                    {
                        throw new InvalidOperationException("ERROR: Sorry, can't salvage Lanczos.");
                        /* ... save yourselves, men.  */
                    }

                    j--;
                }

                vecscale(q[j], 1, n, 1.0 / beta[j], r);
                blas_time += lanc_seconds() - time;
                time = lanc_seconds();
                if (version == 1)
                {
                    splarax(Aq, A, n, q[j], vwsqrt, workn);
                }
                else if (version == 2)
                {
                    symmlq(long_n, &(q[j][1]), &wv1[1], &wv2[1], &wv3[1], &wv4[1], &Aq[1], &wv5[1], &wv6[1],
                        checka, goodb, precon, &shift, &nout, &intlim, &rtol, &istop, &itn, &anorm, &acond,
                        &rnorm, &ynorm, (double*) A, vwsqrt, (double*) orthlist2, &macheps, &normxlim,
                        &itnmin);
                    symmlqitns += itn;
                    if (DEBUG_EVECS > 2)
                    {
                        Console.WriteLine("Symmlq report:      rtol {0:g}", rtol);
                        Console.WriteLine("  system norm {0:g}, solution norm {1:g}", anorm, ynorm);
                        Console.WriteLine("  system condition {0:g}, residual {1:g}", acond, rnorm);
                        Console.WriteLine("  termination condition {0:d}, iterations {1:d}", istop, itn);
                    }
                }

                splarax_time += lanc_seconds() - time;
                time = lanc_seconds();
                update(u, 1, n, Aq, -beta[j], q[j - 1]);
                alpha[j] = dot(u, 1, n, q[j]);
                update(r, 1, n, u, -alpha[j], q[j]);
                blas_time += lanc_seconds() - time;
                time = lanc_seconds();
                if (vwsqrt == null)
                {
                    orthog1(r, 1, n);
                }
                else
                {
                    orthogvec(r, 1, n, vwsqrt);
                }

                orthogonalize(r, n, orthlist);
                temp = orthlist;
                orthlist = makeorthlnk();
                orthlist->vec = q[j];
                orthlist->pntr = temp;
                beta[j + 1] = ch_norm(r, 1, n);
                orthog_time += lanc_seconds() - time;

                time = lanc_seconds();
                left = j / 2;
                right = j - left + 1;
                if (inc_bis_safety)
                {
                    bis_safety *= 10;
                    inc_bis_safety = false;
                }

                ritzval_flag = get_ritzvals(alpha, beta + 1, j, Anorm, workj + 1, ritz, d, left, right, eigtol,
                    bis_safety);
                /* ... have to off-set beta and workj since full orthogonalization
                       indexes these from 1 to maxj+1 whereas selective orthog.
                       indexes them from 0 to maxj */

                if (ritzval_flag)
                {
                    throw new InvalidOperationException("ERROR: Both Sturm bisection and QL failed.");
                    /* ... give up. */
                }

                ql_time += lanc_seconds() - time;

                /* Convergence check using Paige bji estimates. */
                time = lanc_seconds();
                for (i = 1; i <= j; i++)
                {
                    Sres = TevecHelper.Tevec(alpha, beta, j, ritz[i], s);
                    if (Sres > Sres_max)
                    {
                        Sres_max = Sres;
                    }

                    if (Sres > SRESTOL)
                    {
                        inc_bis_safety = true;
                    }

                    bj[i] = s[j] * beta[j + 1];
                }

                tevec_time += lanc_seconds() - time;

                time = lanc_seconds();
                if (version == 1)
                {
                    scanmin(ritz, 1, j, &scanlist);
                }
                else
                {
                    scanmax(ritz, 1, j, &scanlist);
                }

                converged = true;
                if (j < d)
                {
                    converged = false;
                }
                else
                {
                    curlnk = scanlist;
                    while (curlnk != null)
                    {
                        if (bj[curlnk->indx] > bji_tol)
                        {
                            converged = false;
                        }

                        curlnk = curlnk->pntr;
                    }
                }

                scan_time += lanc_seconds() - time;
                j++;
            }

            j--;

            /* Collect eigenvalue and bound information. */
            time = lanc_seconds();
                mkeigvecs(scanlist, lambda, bound, index, bj, d, &Sres_max, alpha, beta + 1, j, s, y, n, q);

            evec_time += lanc_seconds() - time;

            /* Analyze computation for and report additional problems */
            time = lanc_seconds();
            if (DEBUG_EVECS > 0 && version == 2)
            {
                Console.WriteLine("\nTotal Symmlq iterations {0:d}", symmlqitns);
            }

            if (version == 2)
            {
                for (i = 1; i <= d; i++)
                {
                    lambda[i] = 1.0 / lambda[i];
                }
            }

            warnings(workn, A, y, n, lambda, vwsqrt, Ares, bound, index, d, j, maxj, Sres_max, eigtol, u, Anorm);
            debug_time += lanc_seconds() - time;

            /* Free any memory allocated in this routine. */
            time = lanc_seconds();
            frvec(u, 1);
            frvec(r, 1);
            frvec(Aq, 1);
            frvec(ritzvec, 1);
            frvec(zeros, 1);
            if (vwsqrt == null && version == 2)
            {
                frvec(ones, 1);
            }

            frvec(workn, 1);
            frvec(Ares, 1);
            frvec(inv_lambda, 1);
            Marshal.FreeHGlobal((IntPtr) index);
            frvec(alpha, 1);
            frvec(beta, 1);
            frvec(ritz, 1);
            frvec(s, 1);
            frvec(bj, 1);
            frvec(workj, 1);
            if (version == 2)
            {
                frvec(wv1, 0);
            }

            while (scanlist != null)
            {
                curlnk = scanlist->pntr;
                Marshal.FreeHGlobal((IntPtr) scanlist);
                scanlist = curlnk;
            }

            for (i = 1; i <= j; i++)
            {
                frvec(q[i], 1);
            }

            while (orthlist != null)
            {
                temp = orthlist->pntr;
                Marshal.FreeHGlobal((IntPtr) orthlist);
                orthlist = temp;
            }

            while (version == 2 && orthlist2 != null)
            {
                temp = orthlist2->pntr;
                Marshal.FreeHGlobal((IntPtr) orthlist2);
                orthlist2 = temp;
            }

            Marshal.FreeHGlobal((IntPtr) q);
            init_time += lanc_seconds() - time;
        }
    }
}
