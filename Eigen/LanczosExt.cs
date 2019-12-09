using System;
using System.Diagnostics;
using System.Runtime.InteropServices;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Eigen.Orthogonalization;
using static ChacoSharp.Eigen.CheckEig;
using static ChacoSharp.Eigen.Scan;
using static ChacoSharp.Eigen.GetExtVal;
using static ChacoSharp.Eigen.SoListOut;
using static ChacoSharp.Eigen.Splarax;
using static ChacoSharp.Eigen.LanczosTime;
using static ChacoSharp.Eigen.GetRitzVals;
using static ChacoSharp.Eigen.Warnings;
using static ChacoSharp.Utilities.MkVec;
using static ChacoSharp.Utilities.CpVec;
using static ChacoSharp.Utilities.Dot;
using static ChacoSharp.Utilities.Update;
using static ChacoSharp.Utilities.SetVec;
using static ChacoSharp.Utilities.Scadd;
using static ChacoSharp.Utilities.Norm;
using static ChacoSharp.Utilities.VecScale;
using static ChacoSharp.Eigen.TevecHelper;

namespace ChacoSharp.Eigen
{
    public static unsafe class LanczosExt
    {
        /* This is based on lanczos_SO but is simplified and tailored to solve
   the extended eigen-problem using Rafael's technique. Solves Nx = lambda*Wx + Wg.
   In the terminal propagation context we know that e^T Wg = 0 by virue of the way
   g is constructed. This code, however, is designed to work in the more general
   case of arbitrary g (but this capability has not been tested). Returns 0 if it
   thinks everything's OK, nonzero otherwise. */

        public static bool lanczos_ext(vtx_data** A, /* sparse matrix in row linked list format */
            int n, /* problem size */
            int d, /* problem dimension = number of eigvecs to find */
            double*[] y, /* columns of y are eigenvectors of A  */
            double eigtol, /* tolerance on eigenvectors */
            double* vwsqrt, /* square roots of vertex weights */
            double maxdeg, /* maximum degree of graph */
            int version, /* flags which version of sel. orth. to use */
            double* gvec, /* the rhs n-vector in the extended eigen problem */
            double sigma /* specifies the norm constraint on extended
                                              eigenvector */
        )
        {
            int i, j, k; /* indices */
            int maxj; /* maximum number of Lanczos iterations */
            double* u;
            double* r; /* Lanczos vectors */
            double* alpha;
            double* beta; /* the Lanczos scalars from each step */
            double* ritz; /* copy of alpha for ql */
            double* workj; /* work vector, e.g. copy of beta for ql */
            double* workn; /* work vector, e.g. product Av for checkeig */
            double* s; /* eigenvector of T */
            double** q; /* columns of q are Lanczos basis vectors */
            double* bj; /* beta(j)*(last el. of corr. eigvec s of T) */
            double Sres; /* how well Tevec calculated eigvec s */
            double Sres_max; /* Max value of Sres */
            bool inc_bis_safety; /* has Sres increased? */
            double* Ares; /* how well Lanczos calc. eigpair lambda,y */
            int* index; /* the Ritz index of an eigenpair */
            orthlink** solist; /* vec. of structs with vecs. to orthog. against */
            scanlink* scanlist; /* linked list of fields to do with min ritz vals */
            scanlink* curlnk; /* for traversing the scanlist */
            double bis_safety; /* real safety for T bisection algorithm */
            bool converged; /* has the iteration converged? */
            double goodtol; /* error tolerance for a good Ritz vector */
            int ngood; /* total number of good Ritz pairs at current step */
            int maxngood; /* biggest val of ngood through current step */
            int left_ngood; /* number of good Ritz pairs on left end */
            int lastpause; /* Most recent step with good ritz vecs */
            bool nopauses; /* Have there been any pauses? */
            int interval; /* number of steps between pauses */
            double time; /* Current clock time */
            int left_goodlim; /* number of ritz pairs checked on left end */
            double Anorm; /* Norm estimate of the Laplacian matrix */
            int pausemode; /* which Lanczos pausing criterion to use */
            bool pause; /* whether to pause */
            int temp; /* used to prevent redundant index computations */
            double* extvec; /* n-vector solving the extended A eigenproblem */
            double* v; /* j-vector solving the extended T eigenproblem */
            double extval = 0.0; /* computed extended eigenvalue (of both A and T) */
            double* work1;
            double* work2; /* work vectors */
            double check; /* to check an orthogonality condition */
            double numerical_zero; /* used for zero in presence of round-off  */
            bool ritzval_flag; /* status flag for get_ritzvals() */
            bool memory_ok; /* true until memory runs out */

            if (DEBUG_TRACE)
            {
                Trace.WriteLine("<Entering lanczos_ext>");
            }

            if (DEBUG_EVECS > 0)
            {
                Trace.WriteLine($"Selective orthogonalization Lanczos for extended eigenproblem, matrix size = {n:d}.");
            }

            /* Initialize time. */
            time = lanc_seconds();

            if (d != 1)
            {
                throw new ArgumentOutOfRangeException(nameof(d), "ERROR: Extended Lanczos only available for bisection.");
                /* ... something must be wrong upstream. */
            }

            if (n < d + 1)
            {
                throw new ArgumentOutOfRangeException(nameof(n), "ERROR: System too small for number of eigenvalues requested.");
                /* ... d+1 since don't use zero eigenvalue pair */
            }

            /* Allocate space. */
            maxj = LANCZOS_MAXITNS;
            u = mkvec(1, n);
            r = mkvec(1, n);
            workn = mkvec(1, n);
            Ares = mkvec(0, d);
            index = (int*) Marshal.AllocHGlobal((d + 1) * sizeof(int));
            alpha = mkvec(1, maxj);
            beta = mkvec(0, maxj);
            ritz = mkvec(1, maxj);
            s = mkvec(1, maxj);
            bj = mkvec(1, maxj);
            workj = mkvec(0, maxj);
            q = (double**) Marshal.AllocHGlobal((maxj + 1) * sizeof(double*));
            solist = (orthlink**) Marshal.AllocHGlobal((maxj + 1) * sizeof(orthlink*));
            scanlist = mkscanlist(d);
            extvec = mkvec(1, n);
            v = mkvec(1, maxj);
            work1 = mkvec(1, maxj);
            work2 = mkvec(1, maxj);

            /* Set some constants governing orthogonalization */
            ngood = 0;
            maxngood = 0;
            Anorm = 2 * maxdeg; /* Gershgorin estimate for ||A|| */
            goodtol = Anorm * Math.Sqrt(DOUBLE_EPSILON); /* Parlett & Scott's bound, p.224 */
            interval = 2 + Math.Min(LANCZOS_SO_INTERVAL - 2, n / (2 * LANCZOS_SO_INTERVAL));
            bis_safety = BISECTION_SAFETY;
            numerical_zero = 1.0e-13;

            if (DEBUG_EVECS > 0)
            {
                Trace.WriteLine($"  maxdeg {maxdeg:g}");
                Trace.WriteLine($"  goodtol {goodtol:g}");
                Trace.WriteLine($"  interval {interval:d}");
                Trace.WriteLine($"  maxj {maxj:d}");
            }

            /* Initialize space. */
            cpvec(r, 1, n, gvec);
            if (vwsqrt != null)
            {
                scale_diag(r, 1, n, vwsqrt);
            }

            check = ch_norm(r, 1, n);
            if (vwsqrt == null)
            {
                orthog1(r, 1, n);
            }
            else
            {
                orthogvec(r, 1, n, vwsqrt);
            }

            check = Math.Abs(check - ch_norm(r, 1, n));
            if (check > 10 * numerical_zero && WARNING_EVECS > 0)
            {
                Trace.WriteLine("WARNING: In terminal propagation, rhs should have no component in the");
                Trace.WriteLine($"         nullspace of the Laplacian, so check val {check:g} should be negligible.");
            }

            beta[0] = ch_norm(r, 1, n);
            q[0] = mkvec(1, n);
            setvec(q[0], 1, n, 0.0);
            setvec(bj, 1, maxj, double.MaxValue);

            if (beta[0] < numerical_zero)
            {
                /* The rhs vector, Dg, of the transformed problem is numerically zero or is
                   in the null space of the Laplacian, so this is not a well posed extended
                   eigenproblem. Set maxj to zero to force a quick exit but still clean-up
                   memory and return(1) to indicate to eigensolve that it should call the
                   default eigensolver routine for the standard eigenproblem. */
                maxj = 0;
            }

            /* Main Lanczos loop. */
            j = 1;
            lastpause = 0;
            pausemode = 1;
            left_ngood = 0;
            left_goodlim = 0;
            converged = false;
            Sres_max = 0.0;
            inc_bis_safety = false;
            nopauses = true;
            memory_ok = true;
            init_time += lanc_seconds() - time;
            while ((j <= maxj) && (!converged) && memory_ok)
            {
                time = lanc_seconds();

                /* Allocate next Lanczos vector. If fail, back up to last pause. */
                q[j] = mkvec_ret(1, n);
                if (q[j] == null)
                {
                    memory_ok = false;
                    if (DEBUG_EVECS > 0 || WARNING_EVECS > 0)
                    {
                        Trace.WriteLine("WARNING: Lanczos_ext out of memory; computing best approximation available.");
                    }

                    if (nopauses)
                    {
                        throw new InvalidOperationException("ERROR: Sorry, can't salvage Lanczos_ext.");
                        /* ... save yourselves, men.  */
                    }

                    for (i = lastpause + 1; i <= j - 1; i++)
                    {
                        frvec(q[i], 1);
                    }

                    j = lastpause;
                }

                /* Basic Lanczos iteration */
                vecscale(q[j], 1, n, 1.0 / beta[j - 1], r);
                blas_time += lanc_seconds() - time;
                time = lanc_seconds();
                splarax(u, A, n, q[j], vwsqrt, workn);
                splarax_time += lanc_seconds() - time;
                time = lanc_seconds();
                update(r, 1, n, u, -beta[j - 1], q[j - 1]);
                alpha[j] = dot(r, 1, n, q[j]);
                update(r, 1, n, r, -alpha[j], q[j]);
                blas_time += lanc_seconds() - time;

                /* Selective orthogonalization */
                time = lanc_seconds();
                if (vwsqrt == null)
                {
                    orthog1(r, 1, n);
                }
                else
                {
                    orthogvec(r, 1, n, vwsqrt);
                }

                if ((j == (lastpause + 1)) || (j == (lastpause + 2)))
                {
                    sorthog(r, n, solist, ngood);
                }

                orthog_time += lanc_seconds() - time;
                beta[j] = ch_norm(r, 1, n);
                time = lanc_seconds();
                pause = lanpause(j, lastpause, interval, q, n, &pausemode, version, beta[j]);
                pause_time += lanc_seconds() - time;
                if (pause)
                {
                    nopauses = false;
                    lastpause = j;

                    /* Compute limits for checking Ritz pair convergence. */
                    if (version == 2)
                    {
                        if (left_ngood + 2 > left_goodlim)
                        {
                            left_goodlim = left_ngood + 2;
                        }
                    }

                    /* Special case: need at least d Ritz vals on left. */
                    left_goodlim = Math.Max(left_goodlim, d);

                    /* Special case: can't find more than j total Ritz vals. */
                    if (left_goodlim > j)
                    {
                        left_goodlim = Math.Min(left_goodlim, j);
                    }

                    /* Find Ritz vals using faster of Sturm bisection or ql. */
                    time = lanc_seconds();
                    if (inc_bis_safety)
                    {
                        bis_safety *= 10;
                        inc_bis_safety = false;
                    }

                    ritzval_flag =
                        get_ritzvals(alpha, beta, j, Anorm, workj, ritz, d, left_goodlim, 0, eigtol, bis_safety);
                    ql_time += lanc_seconds() - time;

                    if (ritzval_flag)
                    {
                        throw new InvalidOperationException("ERROR: Lanczos_ext failed in computing eigenvalues of T.");
                        /* ... we recover from this in lanczos_SO, but don't worry here. */
                    }

                    /* Scan for minimum evals of tridiagonal. */
                    time = lanc_seconds();
                    scanmin(ritz, 1, j, &scanlist);
                    scan_time += lanc_seconds() - time;

                    /* Compute Ritz pair bounds at left end. */
                    time = lanc_seconds();
                    setvec(bj, 1, j, 0.0);
                    for (i = 1; i <= left_goodlim; i++)
                    {
                        Sres = Tevec(alpha, beta - 1, j, ritz[i], s);
                        if (Sres > Sres_max)
                        {
                            Sres_max = Sres;
                        }

                        if (Sres > SRESTOL)
                        {
                            inc_bis_safety = true;
                        }

                        bj[i] = s[j] * beta[j];
                    }

                    ritz_time += lanc_seconds() - time;

                    /* Show portion of spectrum checked for Ritz convergence. */
                    if (DEBUG_EVECS > 2)
                    {
                        time = lanc_seconds();
                        Trace.WriteLine("\nindex         Ritz vals            bji bounds");
                        for (i = 1; i <= left_goodlim; i++)
                        {
                            Trace.Write("  {i:d}");
                            doubleout(ritz[i], 1);
                            doubleout(bj[i], 1);
                            Trace.WriteLine("");
                        }

                        Trace.WriteLine("");
                        curlnk = scanlist;
                        while (curlnk != null)
                        {
                            temp = curlnk->indx;
                            if ((temp > left_goodlim) && (temp < j))
                            {
                                Trace.Write($"  {temp:d}");
                                doubleout(ritz[temp], 1);
                                doubleout(bj[temp], 1);
                                Trace.WriteLine("");
                            }

                            curlnk = curlnk->pntr;
                        }

                        Trace.WriteLine("                            -------------------");
                        Trace.WriteLine("                goodtol:    {goodtol:f}\n");
                        debug_time += lanc_seconds() - time;
                    }

                    get_extval(alpha, beta, j, ritz[1], s, eigtol, beta[0], sigma, &extval, v, work1, work2);

                    /* Check convergence of iteration. */
                    if (Math.Abs(beta[j] * v[j]) < eigtol)
                    {
                        converged = true;
                    }
                    else
                    {
                        converged = false;
                    }

                    if (!converged)
                    {
                        ngood = 0;
                        left_ngood = 0; /* for setting left_goodlim on next loop */

                        /* Compute converged Ritz pairs on left end */
                        time = lanc_seconds();
                        for (i = 1; i <= left_goodlim; i++)
                        {
                            if (bj[i] <= goodtol)
                            {
                                ngood += 1;
                                left_ngood += 1;
                                if (ngood > maxngood)
                                {
                                    maxngood = ngood;
                                    solist[ngood] = makeorthlnk();
                                    (solist[ngood])->vec = mkvec(1, n);
                                }

                                (solist[ngood])->index = i;
                                Sres = Tevec(alpha, beta - 1, j, ritz[i], s);
                                if (Sres > Sres_max)
                                {
                                    Sres_max = Sres;
                                }

                                if (Sres > SRESTOL)
                                {
                                    inc_bis_safety = true;
                                }

                                setvec((solist[ngood])->vec, 1, n, 0.0);
                                for (k = 1; k <= j; k++)
                                {
                                    scadd((solist[ngood])->vec, 1, n, s[k], q[k]);
                                }
                            }
                        }

                        ritz_time += lanc_seconds() - time;

                        if (DEBUG_EVECS > 2)
                        {
                            time = lanc_seconds();

                            /* Show some info on the orthogonalization. */
                            Trace.Write($"  j {j:d}; goodlim lft {left_goodlim:d}, rgt {0:d}; list ");
                            solistout(solist, ngood, j);

                            /* Assemble current approx. eigenvector, check residual directly. */
                            setvec(y[1], 1, n, 0.0);
                            for (k = 1; k <= j; k++)
                            {
                                scadd(y[1], 1, n, v[k], q[k]);
                            }

                            Trace.WriteLine($"  extended eigenvalue {extval:g}");
                            Trace.WriteLine($"  est. extended residual {Math.Abs(v[j] * beta[j]):g}");
                            checkeig_ext(workn, u, A, y[1], n, extval, vwsqrt, gvec, eigtol, false);

                            Trace.WriteLine("---------------------end of iteration---------------------\n");
                            debug_time += lanc_seconds() - time;
                        }
                    }
                }

                j++;
            }

            j--;

            if (DEBUG_EVECS > 0)
            {
                time = lanc_seconds();
                if (maxj == 0)
                {
                    Trace.WriteLine("Not extended eigenproblem -- calling ordinary eigensolver.");
                }
                else
                {
                    Trace.WriteLine($"  Lanczos_ext itns: {j:d}");
                    Trace.WriteLine($"  extended eigenvalue: {extval:g}");
                    if (j == maxj)
                    {
                        Trace.WriteLine("WARNING: Maximum number of Lanczos iterations reached.");
                    }
                }

                debug_time += lanc_seconds() - time;
            }

            if (maxj != 0)
            {
                /* Compute (scaled) extended eigenvector. */
                time = lanc_seconds();
                setvec(y[1], 1, n, 0.0);
                for (k = 1; k <= j; k++)
                {
                    scadd(y[1], 1, n, v[k], q[k]);
                }

                evec_time += lanc_seconds() - time;
                /* Note: assign() will scale this y vector back to x (since y = Dx) */

                /* Compute and check residual directly. */
                if (DEBUG_EVECS > 0 || WARNING_EVECS > 0)
                {
                    time = lanc_seconds();
                    checkeig_ext(workn, u, A, y[1], n, extval, vwsqrt, gvec, eigtol, true);
                    debug_time += lanc_seconds() - time;
                }
            }

            /* free up memory */
            time = lanc_seconds();
            frvec(u, 1);
            frvec(r, 1);
            frvec(workn, 1);
            frvec(Ares, 0);
            Marshal.FreeHGlobal((IntPtr) index);
            frvec(alpha, 1);
            frvec(beta, 0);
            frvec(ritz, 1);
            frvec(s, 1);
            frvec(bj, 1);
            frvec(workj, 0);
            for (i = 0; i <= j; i++)
            {
                frvec(q[i], 1);
            }

            Marshal.FreeHGlobal((IntPtr) q);
            while (scanlist != null)
            {
                curlnk = scanlist->pntr;
                Marshal.FreeHGlobal((IntPtr) scanlist);
                scanlist = curlnk;
            }

            for (i = 1; i <= maxngood; i++)
            {
                frvec((solist[i])->vec, 1);
                Marshal.FreeHGlobal((IntPtr) solist[i]);
            }

            Marshal.FreeHGlobal((IntPtr) solist);
            frvec(extvec, 1);
            frvec(v, 1);
            frvec(work1, 1);
            frvec(work2, 1);
            init_time += lanc_seconds() - time;

            /* see note on beta[0] and maxj above */
            return maxj == 0;
        }

        /* See comments in lanczos_ext(). */

        public static bool lanczos_ext_float(vtx_data** A, /* sparse matrix in row linked list format */
            int n, /* problem size */
            int d, /* problem dimension = number of eigvecs to find */
            double*[] y, /* columns of y are eigenvectors of A  */
            double eigtol, /* tolerance on eigenvectors */
            double* vwsqrt, /* square roots of vertex weights */
            double maxdeg, /* maximum degree of graph */
            int version, /* flags which version of sel. orth. to use */
            double* gvec, /* the rhs n-vector in the extended eigen problem */
            double sigma /* specifies the norm constraint on extended
                                                 eigenvector */
        )
        {
            int i, j, k; /* indices */
            int maxj; /* maximum number of Lanczos iterations */
            float* u;
            float* r; /* Lanczos vectors */
            double* u_double; /* double version of u */
            double* alpha;
            double* beta; /* the Lanczos scalars from each step */
            double* ritz; /* copy of alpha for ql */
            double* workj; /* work vector, e.g. copy of beta for ql */
            float* workn; /* work vector, e.g. product Av for checkeig */
            double* workn_double; /* work vector, e.g. product Av for checkeig */
            double* s; /* eigenvector of T */
            float** q; /* columns of q are Lanczos basis vectors */
            double* bj; /* beta(j)*(last el. of corr. eigvec s of T) */
            double bis_safety; /* real safety factor for T bisection */
            double Sres; /* how well Tevec calculated eigvec s */
            double Sres_max; /* Max value of Sres */
            bool inc_bis_safety; /* need to increase bisection safety */
            double* Ares; /* how well Lanczos calc. eigpair lambda,y */
            int* index; /* the Ritz index of an eigenpair */
            orthlink_float** solist; /* vec. of structs with vecs. to orthog. against */
            scanlink* scanlist; /* linked list of fields to do with min ritz vals */
            scanlink* curlnk; /* for traversing the scanlist */
            double bji_tol; /* tol on bji est. of eigen residual of A */
            bool converged; /* has the iteration converged? */
            double goodtol; /* error tolerance for a good Ritz vector */
            int ngood; /* total number of good Ritz pairs at current step */
            int maxngood; /* biggest val of ngood through current step */
            int left_ngood; /* number of good Ritz pairs on left end */
            int lastpause; /* Most recent step with good ritz vecs */
            bool nopauses; /* Have there been any pauses? */
            int interval; /* number of steps between pauses */
            double time; /* Current clock time */
            int left_goodlim; /* number of ritz pairs checked on left end */
            double Anorm; /* Norm estimate of the Laplacian matrix */
            int pausemode; /* which Lanczos pausing criterion to use */
            bool pause; /* whether to pause */
            int temp; /* used to prevent redundant index computations */
            double* extvec; /* n-vector solving the extended A eigenproblem */
            double* v; /* j-vector solving the extended T eigenproblem */
            double extval = 0.0; /* computed extended eigenvalue (of both A and T) */
            double* work1;
            double* work2; /* work vectors */
            double check; /* to check an orthogonality condition */
            double numerical_zero; /* used for zero in presence of round-off  */
            bool ritzval_flag; /* status flag for get_ritzvals() */
            double resid; /* residual */
            bool memory_ok; /* true until memory runs out */
            float* vwsqrt_float = null; /* float version of vwsqrt */

            if (DEBUG_TRACE)
            {
                Trace.WriteLine("<Entering lanczos_ext_float>");
            }

            if (DEBUG_EVECS > 0)
            {
                Trace.WriteLine($"Selective orthogonalization Lanczos for extended eigenproblem, matrix size = {n:d}.");
            }

            /* Initialize time. */
            time = lanc_seconds();

            if (d != 1)
            {
                throw new ArgumentOutOfRangeException(nameof(d), "ERROR: Extended Lanczos only available for bisection.");
                /* ... something must be wrong upstream. */
            }

            if (n < d + 1)
            {
                throw new ArgumentOutOfRangeException(nameof(n), "ERROR: System too small for number of eigenvalues requested.");
                /* ... d+1 since don't use zero eigenvalue pair */
            }

            /* Allocate space. */
            maxj = LANCZOS_MAXITNS;
            u = mkvec_float(1, n);
            u_double = mkvec(1, n);
            r = mkvec_float(1, n);
            workn = mkvec_float(1, n);
            workn_double = mkvec(1, n);
            Ares = mkvec(0, d);
            index = (int*) Marshal.AllocHGlobal((d + 1) * sizeof(int));
            alpha = mkvec(1, maxj);
            beta = mkvec(0, maxj);
            ritz = mkvec(1, maxj);
            s = mkvec(1, maxj);
            bj = mkvec(1, maxj);
            workj = mkvec(0, maxj);
            q = (float**) Marshal.AllocHGlobal((maxj + 1) * sizeof(float*));
            solist = (orthlink_float**) Marshal.AllocHGlobal((maxj + 1) * sizeof(orthlink_float*));
            scanlist = mkscanlist(d);
            extvec = mkvec(1, n);
            v = mkvec(1, maxj);
            work1 = mkvec(1, maxj);
            work2 = mkvec(1, maxj);

            /* Set some constants governing orthogonalization */
            ngood = 0;
            maxngood = 0;
            bji_tol = eigtol;
            Anorm = 2 * maxdeg; /* Gershgorin estimate for ||A|| */
            goodtol = Anorm * Math.Sqrt(DOUBLE_EPSILON); /* Parlett & Scott's bound, p.224 */
            interval = 2 + Math.Min(LANCZOS_SO_INTERVAL - 2, n / (2 * LANCZOS_SO_INTERVAL));
            bis_safety = BISECTION_SAFETY;
            numerical_zero = 1.0e-6;

            if (DEBUG_EVECS > 0)
            {
                Trace.WriteLine($"  maxdeg {maxdeg:g}");
                Trace.WriteLine($"  goodtol {goodtol:g}");
                Trace.WriteLine($"  interval {interval:d}");
                Trace.WriteLine($"  maxj {maxj:d}");
            }

            /* Make a float copy of vwsqrt */
            if (vwsqrt != null)
            {
                vwsqrt_float = mkvec_float(0, n);
                double_to_float(vwsqrt_float, 1, n, vwsqrt);
            }

            /* Initialize space. */
            double_to_float(r, 1, n, gvec);
            if (vwsqrt_float != null)
            {
                scale_diag_float(r, 1, n, vwsqrt_float);
            }

            check = norm_float(r, 1, n);
            if (vwsqrt_float == null)
            {
                orthog1_float(r, 1, n);
            }
            else
            {
                orthogvec_float(r, 1, n, vwsqrt_float);
            }

            check = Math.Abs(check - norm_float(r, 1, n));
            if (check > 10 * numerical_zero && WARNING_EVECS > 0)
            {
                Trace.WriteLine("WARNING: In terminal propagation, rhs should have no component in the");
                Trace.WriteLine($"         nullspace of the Laplacian, so check val {check:g} should be zero.");
            }

            beta[0] = norm_float(r, 1, n);
            q[0] = mkvec_float(1, n);
            setvec_float(q[0], 1, n, 0.0f);
            setvec(bj, 1, maxj, double.MaxValue);

            if (beta[0] < numerical_zero)
            {
                /* The rhs vector, Dg, of the transformed problem is numerically zero or is
                   in the null space of the Laplacian, so this is not a well posed extended
                   eigenproblem. Set maxj to zero to force a quick exit but still clean-up
                   memory and return(1) to indicate to eigensolve that it should call the
                   default eigensolver routine for the standard eigenproblem. */
                maxj = 0;
            }

            /* Main Lanczos loop. */
            j = 1;
            lastpause = 0;
            pausemode = 1;
            left_ngood = 0;
            left_goodlim = 0;
            converged = false;
            Sres_max = 0.0;
            inc_bis_safety = false;
            nopauses = true;
            memory_ok = true;
            init_time += lanc_seconds() - time;
            while ((j <= maxj) && (!converged) && memory_ok)
            {
                time = lanc_seconds();

                /* Allocate next Lanczos vector. If fail, back up to last pause. */
                q[j] = mkvec_ret_float(1, n);
                if (q[j] == null)
                {
                    memory_ok = false;
                    if (DEBUG_EVECS > 0 || WARNING_EVECS > 0)
                    {
                        Trace.WriteLine("WARNING: Lanczos_ext out of memory; computing best approximation available.");
                    }

                    if (nopauses)
                    {
                        throw new InvalidOperationException("ERROR: Sorry, can't salvage Lanczos_ext.");
                        /* ... save yourselves, men.  */
                    }

                    for (i = lastpause + 1; i <= j - 1; i++)
                    {
                        frvec_float(q[i], 1);
                    }

                    j = lastpause;
                }

                /* Basic Lanczos iteration */
                vecscale_float(q[j], 1, n, (float) (1.0 / beta[j - 1]), r);
                blas_time += lanc_seconds() - time;
                time = lanc_seconds();
                splarax_float(u, A, n, q[j], vwsqrt_float, workn);
                splarax_time += lanc_seconds() - time;
                time = lanc_seconds();
                update_float(r, 1, n, u, (float) (-beta[j - 1]), q[j - 1]);
                alpha[j] = dot_float(r, 1, n, q[j]);
                update_float(r, 1, n, r, (float) (-alpha[j]), q[j]);
                blas_time += lanc_seconds() - time;

                /* Selective orthogonalization */
                time = lanc_seconds();
                if (vwsqrt_float == null)
                {
                    orthog1_float(r, 1, n);
                }
                else
                {
                    orthogvec_float(r, 1, n, vwsqrt_float);
                }

                if ((j == (lastpause + 1)) || (j == (lastpause + 2)))
                {
                    sorthog_float(r, n, solist, ngood);
                }

                orthog_time += lanc_seconds() - time;
                beta[j] = norm_float(r, 1, n);
                time = lanc_seconds();
                pause = lanpause_float(j, lastpause, interval, q, n, &pausemode, version, beta[j]);
                pause_time += lanc_seconds() - time;
                if (pause)
                {
                    nopauses = false;
                    lastpause = j;

                    /* Compute limits for checking Ritz pair convergence. */
                    if (version == 2)
                    {
                        if (left_ngood + 2 > left_goodlim)
                        {
                            left_goodlim = left_ngood + 2;
                        }
                    }

                    /* Special case: need at least d Ritz vals on left. */
                    left_goodlim = Math.Max(left_goodlim, d);

                    /* Special case: can't find more than j total Ritz vals. */
                    if (left_goodlim > j)
                    {
                        left_goodlim = Math.Min(left_goodlim, j);
                    }

                    /* Find Ritz vals using faster of Sturm bisection or ql. */
                    time = lanc_seconds();
                    if (inc_bis_safety)
                    {
                        bis_safety *= 10;
                        inc_bis_safety = false;
                    }

                    ritzval_flag =
                        get_ritzvals(alpha, beta, j, Anorm, workj, ritz, d, left_goodlim, 0, eigtol, bis_safety);
                    ql_time += lanc_seconds() - time;

                    if (ritzval_flag)
                    {
                        throw new InvalidOperationException("ERROR: Lanczos_ext failed in computing eigenvalues of T.");
                        /* ... we recover from this in lanczos_SO, but don't worry here. */
                    }

                    /* Scan for minimum evals of tridiagonal. */
                    time = lanc_seconds();
                    scanmin(ritz, 1, j, &scanlist);
                    scan_time += lanc_seconds() - time;

                    /* Compute Ritz pair bounds at left end. */
                    time = lanc_seconds();
                    setvec(bj, 1, j, 0.0);
                    for (i = 1; i <= left_goodlim; i++)
                    {
                        Sres = Tevec(alpha, beta - 1, j, ritz[i], s);
                        if (Sres > Sres_max)
                        {
                            Sres_max = Sres;
                        }

                        if (Sres > SRESTOL)
                        {
                            inc_bis_safety = true;
                        }

                        bj[i] = s[j] * beta[j];
                    }

                    ritz_time += lanc_seconds() - time;

                    /* Show the portion of the spectrum checked for convergence. */
                    if (DEBUG_EVECS > 2)
                    {
                        time = lanc_seconds();
                        Trace.WriteLine("\nindex         Ritz vals            bji bounds");
                        for (i = 1; i <= left_goodlim; i++)
                        {
                            Trace.Write($"  {i:d}");
                            doubleout(ritz[i], 1);
                            doubleout(bj[i], 1);
                            Trace.WriteLine("");
                        }

                        Trace.WriteLine("");
                        curlnk = scanlist;
                        while (curlnk != null)
                        {
                            temp = curlnk->indx;
                            if ((temp > left_goodlim) && (temp < j))
                            {
                                Trace.Write($"  {temp:d}");
                                doubleout(ritz[temp], 1);
                                doubleout(bj[temp], 1);
                                Trace.WriteLine("");
                            }

                            curlnk = curlnk->pntr;
                        }

                        Trace.WriteLine("                            -------------------");
                        Trace.WriteLine($"                goodtol:    {goodtol:f}\n");
                        debug_time += lanc_seconds() - time;
                    }

                    get_extval(alpha, beta, j, ritz[1], s, eigtol, beta[0], sigma, &extval, v, work1, work2);

                    /* check convergence of Ritz pairs */
                    time = lanc_seconds();
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

                    if (!converged)
                    {
                        ngood = 0;
                        left_ngood = 0; /* for setting left_goodlim on next loop */

                        /* Compute converged Ritz pairs on left end */
                        time = lanc_seconds();
                        for (i = 1; i <= left_goodlim; i++)
                        {
                            if (bj[i] <= goodtol)
                            {
                                ngood += 1;
                                left_ngood += 1;
                                if (ngood > maxngood)
                                {
                                    maxngood = ngood;
                                    solist[ngood] = makeorthlnk_float();
                                    (solist[ngood])->vec = mkvec_float(1, n);
                                }

                                (solist[ngood])->index = i;
                                Sres = Tevec(alpha, beta - 1, j, ritz[i], s);
                                if (Sres > Sres_max)
                                {
                                    Sres_max = Sres;
                                }

                                if (Sres > SRESTOL)
                                {
                                    inc_bis_safety = true;
                                }

                                setvec_float((solist[ngood])->vec, 1, n, 0.0f);
                                for (k = 1; k <= j; k++)
                                {
                                    scadd_float((solist[ngood])->vec, 1, n, (float) s[k], q[k]);
                                }
                            }
                        }

                        ritz_time += lanc_seconds() - time;

                        if (DEBUG_EVECS > 2)
                        {
                            time = lanc_seconds();
                            Trace.Write($"  j {j:d}; goodlim lft {left_goodlim:d}, rgt {0:d}; list ");
                            solistout_float(solist, ngood, j);
                            Trace.WriteLine("---------------------end of iteration---------------------\n");
                            debug_time += lanc_seconds() - time;
                        }
                    }
                }

                j++;
            }

            j--;

            if (DEBUG_EVECS > 0)
            {
                time = lanc_seconds();
                if (maxj == 0)
                {
                    Trace.WriteLine("Not extended eigenproblem -- calling ordinary eigensolver.");
                }
                else
                {
                    Trace.WriteLine($"  Lanczos_ext itns: {j:d}");
                    Trace.WriteLine($"  eigenvalue: {ritz[1]:g}");
                    Trace.WriteLine($"  extended eigenvalue: {extval:g}");
                }

                debug_time += lanc_seconds() - time;
            }

            if (maxj != 0)
            {
                /* Compute (scaled) extended eigenvector. */
                time = lanc_seconds();
                setvec(y[1], 1, n, 0.0);
                for (k = 1; k <= j; k++)
                {
                    scadd_mixed(y[1], 1, n, v[k], q[k]);
                }

                evec_time += lanc_seconds() - time;
                /* Note: assign() will scale this y vector back to x (since y = Dx) */

                /* Compute and check residual directly. Use the Ay = extval*y + Dg version of
                   the problem for convenience. Note that u and v are used here as workspace */
                time = lanc_seconds();
                splarax(workn_double, A, n, y[1], vwsqrt, u_double);
                scadd(workn_double, 1, n, -extval, y[1]);
                scale_diag(gvec, 1, n, vwsqrt);
                scadd(workn_double, 1, n, -1.0, gvec);
                resid = ch_norm(workn_double, 1, n);
                if (DEBUG_EVECS > 0)
                {
                    Trace.WriteLine($"  extended residual: {resid:g}");
                }

                if (WARNING_EVECS > 0 && resid > eigtol)
                {
                    Trace.WriteLine($"WARNING: Extended residual ({resid:g}) greater than tolerance ({eigtol:g}).");
                }

                debug_time += lanc_seconds() - time;
            }

            /* free up memory */
            time = lanc_seconds();
            frvec_float(u, 1);
            frvec(u_double, 1);
            frvec_float(r, 1);
            frvec_float(workn, 1);
            frvec(workn_double, 1);
            frvec(Ares, 0);
            Marshal.FreeHGlobal((IntPtr) index);
            frvec(alpha, 1);
            frvec(beta, 0);
            frvec(ritz, 1);
            frvec(s, 1);
            frvec(bj, 1);
            frvec(workj, 0);
            for (i = 0; i <= j; i++)
            {
                frvec_float(q[i], 1);
            }

            Marshal.FreeHGlobal((IntPtr) q);
            while (scanlist != null)
            {
                curlnk = scanlist->pntr;
                Marshal.FreeHGlobal((IntPtr) scanlist);
                scanlist = curlnk;
            }

            for (i = 1; i <= maxngood; i++)
            {
                frvec_float((solist[i])->vec, 1);
                Marshal.FreeHGlobal((IntPtr) solist[i]);
            }

            Marshal.FreeHGlobal((IntPtr) solist);
            frvec(extvec, 1);
            frvec(v, 1);
            frvec(work1, 1);
            frvec(work2, 1);
            if (vwsqrt != null)
            {
                frvec_float(vwsqrt_float, 1);
            }

            init_time += lanc_seconds() - time;

            /* see note on beta[0] and maxj above */
            return maxj == 0;
        }
    }
}
