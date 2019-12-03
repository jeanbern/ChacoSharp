using System;
using System.Runtime.InteropServices;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Assignment.AssignFunc;
using static ChacoSharp.Eigen.Orthogonalization;
using static ChacoSharp.Eigen.Scan;
using static ChacoSharp.Utilities.CpVec;
using static ChacoSharp.Eigen.CheckEig;
using static ChacoSharp.Eigen.SoListOut;
using static ChacoSharp.Eigen.Splarax;
using static ChacoSharp.Eigen.LanczosTime;
using static ChacoSharp.Eigen.GetRitzVals;
using static ChacoSharp.Eigen.Warnings;
using static ChacoSharp.Utilities.MkVec;
using static ChacoSharp.Utilities.VecRan;
using static ChacoSharp.Utilities.Dot;
using static ChacoSharp.Utilities.Update;
using static ChacoSharp.Utilities.SetVec;
using static ChacoSharp.Utilities.Scadd;
using static ChacoSharp.Utilities.Norm;
using static ChacoSharp.Utilities.VecScale;
using static ChacoSharp.Eigen.TevecHelper;

namespace ChacoSharp.Eigen
{
    public static unsafe class LanczosSoFloat
    {
        /* See comments in lanczos_SO() */

        public static void lanczos_SO_float(vtx_data** A, /* sparse matrix in row linked list format */
            int n, /* problem size */
            int d, /* problem dimension = number of eigvecs to find */
            double*[] y, /* columns of y are eigenvectors of A  */
            double[] lambda, /* ritz approximation to eigenvals of A */
            double[] bound, /* on ritz pair approximations to eig pairs of A */
            double eigtol, /* tolerance on eigenvectors */
            double* vwsqrt, /* square roots of vertex weights */
            double maxdeg, /* maximum degree of graph */
            int version, /* flags which version of sel. orth. to use */
            bool cube_or_mesh, /* 0 => hypercube, d => d-dimensional mesh */
            int nsets, /* number of sets to divide into */
            int* assignment, /* set number of each vtx (length n+1) */
            int* active, /* space for nvtxs integers */
            MappingType mediantype, /* which partitioning strategy to use */
            double[] goal, /* desired set sizes */
            int vwgt_max /* largest vertex weight */
        )
        {
            double bis_safety; /* real safety factor for T bisection */
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
            int right_ngood; /* number of good Ritz pairs on right end */
            int lastpause; /* Most recent step with good ritz vecs */
            bool firstpause; /* Is this the first pause? */
            bool nopauses; /* Have there been any pauses? */
            int interval; /* number of steps between pauses */
            double time; /* Current clock time */
            int left_goodlim; /* number of ritz pairs checked on left end */
            int right_goodlim; /* number of ritz pairs checked on right end */
            double Anorm; /* Norm estimate of the Laplacian matrix */
            int pausemode; /* which Lanczos pausing criterion to use */
            bool pause; /* whether to pause */
            int temp; /* used to prevent redundant index computations */
            int* old_assignment = null; /* set # of each vtx on previous pause, length n+1 */
            int* assgn_pntr; /* pntr to assignment vector */
            int* old_assgn_pntr; /* pntr to previous assignment vector */
            int assigndiff; /* # of differences between old and new assignment */
            int assigntol; /* tolerance on convergence of assignment vector */
            bool ritzval_flag; /* status flag for get_ritzvals() */
            bool memory_ok; /* True until lanczos runs out of memory */
            float* vwsqrt_float; /* float version of vwsqrt */

            if (DEBUG_TRACE)
            {
                Console.WriteLine("<Entering lanczos_so_float>");
            }

            if (DEBUG_EVECS > 0)
            {
                Console.WriteLine("Selective orthogonalization Lanczos (v. {0:d}), matrix size = {1:d}.", version, n);
            }

            /* Initialize time. */
            time = lanc_seconds();

            if (n < d + 1)
            {
                throw new ArgumentOutOfRangeException(nameof(n), "ERROR: System too small for number of eigenvalues requested.");
                /* d+1 since don't use zero eigenvalue pair */
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
            if (LANCZOS_CONVERGENCE_MODE == 1)
            {
                old_assignment = (int*) Marshal.AllocHGlobal((n + 1) * sizeof(int));
            }

            /* Set some constants governing the orthogonalization heuristic. */
            ngood = 0;
            maxngood = 0;
            bji_tol = eigtol;
            assigntol = (int) (eigtol * n);
            Anorm = 2 * maxdeg; /* Gershgorin estimate for ||A|| */
            goodtol = Anorm * Math.Sqrt(DOUBLE_EPSILON); /* Parlett & Scott's bound, p.224 */
            interval = 2 + Math.Min(LANCZOS_SO_INTERVAL - 2, n / (2 * LANCZOS_SO_INTERVAL));
            bis_safety = BISECTION_SAFETY;

            if (DEBUG_EVECS > 0)
            {
                Console.WriteLine("  maxdeg {0:g}", maxdeg);
                Console.WriteLine("  goodtol {0:g}", goodtol);
                Console.WriteLine("  interval {0:d}", interval);
                Console.WriteLine("  maxj {0:d}", maxj);
                if (LANCZOS_CONVERGENCE_MODE == 1)
                {
                    Console.WriteLine("  assigntol {0:d}", assigntol);
                }
            }

            /* Make a float copy of vwsqrt */
            if (vwsqrt == null)
            {
                vwsqrt_float = null;
            }
            else
            {
                vwsqrt_float = mkvec_float(0, n);
                double_to_float(vwsqrt_float, 1, n, vwsqrt);
            }

            /* Initialize space. */
            vecran_float(r, 1, n);
            if (vwsqrt_float == null)
            {
                orthog1_float(r, 1, n);
            }
            else
            {
                orthogvec_float(r, 1, n, vwsqrt_float);
            }

            beta[0] = norm_float(r, 1, n);
            q[0] = mkvec_float(1, n);
            setvec_float(q[0], 1, n, 0.0f);
            setvec(bj, 1, maxj, double.MaxValue);

            /* Main Lanczos loop. */
            j = 1;
            lastpause = 0;
            pausemode = 1;
            left_ngood = 0;
            right_ngood = 0;
            left_goodlim = 0;
            right_goodlim = 0;
            converged = false;
            Sres_max = 0.0;
            inc_bis_safety = false;
            ritzval_flag = false;
            memory_ok = true;
            firstpause = false;
            nopauses = true;
            init_time += lanc_seconds() - time;
            while ((j <= maxj) && (!converged) && (!ritzval_flag) && memory_ok)
            {
                time = lanc_seconds();

                /* Allocate next Lanczos vector. If fail, back up to last pause. */
                q[j] = mkvec_ret_float(1, n);
                if (q[j] == null)
                {
                    memory_ok = false;
                    if (DEBUG_EVECS > 0 || WARNING_EVECS > 0)
                    {
                        Console.WriteLine("WARNING: Lanczos out of memory; computing best approximation available.");
                    }

                    if (nopauses)
                    {
                        throw new InvalidOperationException("ERROR: Sorry, can't salvage Lanczos.");
                        /* ... save yourselves, men.  */
                    }

                    for (i = lastpause + 1; i <= j - 1; i++)
                    {
                        frvec_float(q[i], 1);
                    }

                    j = lastpause;
                }

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
                    if (lastpause == 0)
                    {
                        firstpause = true;
                    }
                    else
                    {
                        firstpause = false;
                    }

                    lastpause = j;

                    /* Compute limits for checking Ritz pair convergence. */
                    if (version == 1)
                    {
                        if (left_ngood + 2 > left_goodlim)
                        {
                            left_goodlim = left_ngood + 2;
                        }

                        if (right_ngood + 3 > right_goodlim)
                        {
                            right_goodlim = right_ngood + 3;
                        }
                    }

                    if (version == 2)
                    {
                        if (left_ngood + 2 > left_goodlim)
                        {
                            left_goodlim = left_ngood + 2;
                        }

                        right_goodlim = 0;
                    }

                    /* Special case: need at least d Ritz vals on left. */
                    left_goodlim = Math.Max(left_goodlim, d);

                    /* Special case: can't find more than j total Ritz vals. */
                    if (left_goodlim + right_goodlim > j)
                    {
                        left_goodlim = Math.Min(left_goodlim, j);
                        right_goodlim = j - left_goodlim;
                    }

                    /* Find Ritz vals using faster of Sturm bisection or QL. */
                    time = lanc_seconds();
                    blas_time += lanc_seconds() - time;
                    time = lanc_seconds();
                    if (inc_bis_safety)
                    {
                        bis_safety *= 10;
                        inc_bis_safety = false;
                    }

                    ritzval_flag = get_ritzvals(alpha, beta, j, Anorm, workj, ritz, d, left_goodlim,
                        right_goodlim, eigtol, bis_safety);
                    ql_time += lanc_seconds() - time;

                    /* If get_ritzvals() fails, back up to last pause point and exit main loop. */
                    if (ritzval_flag)
                    {
                        if (DEBUG_EVECS > 0 || WARNING_EVECS > 0)
                        {
                            Console.WriteLine("ERROR: Lanczos failed in computing eigenvalues of T; computing");
                            Console.WriteLine("       best readily available approximation to eigenvector.\n");
                        }

                        if (firstpause)
                        {
                            throw new InvalidOperationException("ERROR: Sorry, can't salvage Lanczos.");
                            /* ... save yourselves, men.  */
                        }

                        for (i = lastpause + 1; i <= j; i++)
                        {
                            frvec_float(q[i], 1);
                        }

                        j = lastpause;
                        get_ritzvals(alpha, beta, j, Anorm, workj, ritz, d, left_goodlim, right_goodlim, eigtol,
                            bis_safety);
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

                    /* Compute Ritz pair bounds at right end. */
                    for (i = j; i > j - right_goodlim; i--)
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
                        Console.WriteLine("index         Ritz vals            bji bounds   (j = %d)\n", j);
                        for (i = 1; i <= left_goodlim; i++)
                        {
                            Console.Write("  {0:d}", i);
                            doubleout(ritz[i], 1);
                            doubleout(bj[i], 1);
                            Console.WriteLine();
                        }

                        Console.WriteLine();
                        curlnk = scanlist;
                        while (curlnk != null)
                        {
                            temp = curlnk->indx;
                            if ((temp > left_goodlim) && (temp < j - right_goodlim))
                            {
                                Console.Write("  {0:d}", temp);
                                doubleout(ritz[temp], 1);
                                doubleout(bj[temp], 1);
                                Console.WriteLine();
                            }

                            curlnk = curlnk->pntr;
                        }

                        Console.WriteLine();
                        for (i = j - right_goodlim + 1; i <= j; i++)
                        {
                            Console.Write("  {0:d}", i);
                            doubleout(ritz[i], 1);
                            doubleout(bj[i], 1);
                            Console.WriteLine();
                        }

                        Console.WriteLine("                            -------------------");
                        Console.WriteLine("                goodtol:    %19.16f\n", goodtol);
                        debug_time += lanc_seconds() - time;
                    }

                    /* Check for convergence. */
                    time = lanc_seconds();
                    if (LANCZOS_CONVERGENCE_MODE != 1 || d > 1)
                    {
                        /* check convergence of residual bound */
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
                    }

                    if (LANCZOS_CONVERGENCE_MODE == 1 && d == 1)
                    {
                        /* check change in partition */
                        if (firstpause)
                        {
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

                            if (!converged)
                            {
                                /* compute current approx. to eigenvectors */
                                i = d;
                                curlnk = scanlist;
                                while (curlnk != null)
                                {
                                    lambda[i] = curlnk->val;
                                    bound[i] = bj[curlnk->indx];
                                    index[i] = curlnk->indx;
                                    curlnk = curlnk->pntr;
                                    i--;
                                }

                                for (i = 1; i <= d; i++)
                                {
                                    Sres = Tevec(alpha, beta - 1, j, lambda[i], s);
                                    if (Sres > Sres_max)
                                    {
                                        Sres_max = Sres;
                                    }

                                    if (Sres > SRESTOL)
                                    {
                                        inc_bis_safety = true;
                                    }

                                    setvec(y[i], 1, n, 0.0);
                                    for (k = 1; k <= j; k++)
                                    {
                                        scadd_mixed(y[i], 1, n, s[k], q[k]);
                                    }
                                }

                                Assign(A, y, n, d, cube_or_mesh, nsets, vwsqrt, assignment, active, mediantype, goal, vwgt_max);
                            }
                        }
                        else
                        {
                            /* copy assignment to old_assignment */
                            assgn_pntr = assignment;
                            old_assgn_pntr = old_assignment;
                            for (i = n + 1; i != 0; i--)
                            {
                                *old_assgn_pntr++ = *assgn_pntr++;
                            }

                            /* compute current approx. to eigenvectors */
                            i = d;
                            curlnk = scanlist;
                            while (curlnk != null)
                            {
                                lambda[i] = curlnk->val;
                                bound[i] = bj[curlnk->indx];
                                index[i] = curlnk->indx;
                                curlnk = curlnk->pntr;
                                i--;
                            }

                            for (i = 1; i <= d; i++)
                            {
                                Sres = Tevec(alpha, beta - 1, j, lambda[i], s);
                                if (Sres > Sres_max)
                                {
                                    Sres_max = Sres;
                                }

                                if (Sres > SRESTOL)
                                {
                                    inc_bis_safety = true;
                                }

                                setvec(y[i], 1, n, 0.0);
                                for (k = 1; k <= j; k++)
                                {
                                    scadd_mixed(y[i], 1, n, s[k], q[k]);
                                }
                            }

                            /* write new assignment */
                            Assign(A, y, n, d, cube_or_mesh, nsets, vwsqrt, assignment, active, mediantype, goal, vwgt_max);

                            assigndiff = 0;
                            assgn_pntr = assignment;
                            old_assgn_pntr = old_assignment;
                            for (i = n + 1; i != 0; i--)
                            {
                                if (*old_assgn_pntr++ != *assgn_pntr++)
                                {
                                    assigndiff++;
                                }
                            }

                            assigndiff = Math.Min(assigndiff, n - assigndiff);
                            if (DEBUG_EVECS > 1)
                            {
                                Console.WriteLine("  j {0:d},  change from last assignment {1:d}\n", j, assigndiff);
                            }

                            if (assigndiff <= assigntol)
                            {
                                converged = true;
                            }
                            else
                            {
                                converged = false;
                            }
                        }
                    }

                    scan_time += lanc_seconds() - time;

                    /* Show current estimates of evals and bounds (for help in tuning) */
                    if (DEBUG_EVECS > 2 && !converged)
                    {
                        time = lanc_seconds();
                        /* Collect eigenvalue and bound information for display, return. */
                        i = d;
                        curlnk = scanlist;
                        while (curlnk != null)
                        {
                            lambda[i] = curlnk->val;
                            bound[i] = bj[curlnk->indx];
                            index[i] = curlnk->indx;
                            curlnk = curlnk->pntr;
                            i--;
                        }

                        /* Compute eigenvectors and display associated info. */
                        Console.WriteLine("j {0:d;}    lambda                Ares est.             Ares          index", j);
                        for (i = 1; i <= d; i++)
                        {
                            Sres = Tevec(alpha, beta - 1, j, lambda[i], s);
                            if (Sres > Sres_max)
                            {
                                Sres_max = Sres;
                            }

                            if (Sres > SRESTOL)
                            {
                                inc_bis_safety = true;
                            }

                            setvec(y[i], 1, n, 0.0);
                            for (k = 1; k <= j; k++)
                            {
                                scadd_mixed(y[i], 1, n, s[k], q[k]);
                            }

                            float_to_double(u_double, 1, n, u);
                            Ares[i] = checkeig(workn_double, A, y[i], n, lambda[i], vwsqrt, u_double);
                            Console.Write("{0:d}.", i);
                            doubleout(lambda[i], 1);
                            doubleout(bound[i], 1);
                            doubleout(Ares[i], 1);
                            Console.WriteLine("   {0:d}", index[i]);
                        }

                        Console.WriteLine();
                        debug_time += lanc_seconds() - time;
                    }

                    if (!converged)
                    {
                        ngood = 0;
                        left_ngood = 0; /* for setting left_goodlim on next loop */
                        right_ngood = 0; /* for setting right_goodlim on next loop */

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
                                    scadd_float((solist[ngood])->vec, 1, n, (float) (s[k]), q[k]);
                                }
                            }
                        }

                        /* Compute converged Ritz pairs on right end */
                        for (i = j; i > j - right_goodlim; i--)
                        {
                            if (bj[i] <= goodtol)
                            {
                                ngood += 1;
                                right_ngood += 1;
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
                                    scadd_float((solist[ngood])->vec, 1, n, (float) (s[k]), q[k]);
                                }
                            }
                        }

                        ritz_time += lanc_seconds() - time;

                        if (DEBUG_EVECS > 2)
                        {
                            time = lanc_seconds();
                            Console.Write("  j {0:d}; goodlim lft {1:d}, rgt {2:d}; list ", j, left_goodlim, right_goodlim);
                            solistout_float(solist, ngood, j);
                            Console.WriteLine("---------------------end of iteration---------------------\n");
                            debug_time += lanc_seconds() - time;
                        }
                    }
                }

                j++;
            }

            j--;

            /* Collect eigenvalue and bound information. Only compute and display info for
               the eigpairs actually used in the partitioning since don't want to spend the
               time or space to compute the null-space of the Laplacian. */
            time = lanc_seconds();
            i = d;
            curlnk = scanlist;
            while (curlnk != null)
            {
                lambda[i] = curlnk->val;
                bound[i] = bj[curlnk->indx];
                index[i] = curlnk->indx;
                curlnk = curlnk->pntr;
                i--;
            }

            scan_time += lanc_seconds() - time;

            /* Compute eigenvectors. */
            time = lanc_seconds();
            for (i = 1; i <= d; i++)
            {
                Sres = Tevec(alpha, beta - 1, j, lambda[i], s);
                if (Sres > Sres_max)
                {
                    Sres_max = Sres;
                }

                setvec(y[i], 1, n, 0.0);
                for (k = 1; k <= j; k++)
                {
                    scadd_mixed(y[i], 1, n, s[k], q[k]);
                }
            }

            evec_time += lanc_seconds() - time;

            time = lanc_seconds();
            float_to_double(u_double, 1, n, u);
            warnings(workn_double, A, y, n, lambda, vwsqrt, Ares, bound, index, d, j, maxj, Sres_max, eigtol,
                u_double, Anorm);
            debug_time += lanc_seconds() - time;

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
            if (vwsqrt_float != null)
            {
                frvec_float(vwsqrt_float, 0);
            }

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
            if (LANCZOS_CONVERGENCE_MODE == 1)
            {
                Marshal.FreeHGlobal((IntPtr) old_assignment);
            }

            init_time += lanc_seconds() - time;
        }
    }
}
