#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
using System;
using System.Diagnostics;
using System.Runtime.InteropServices;
using static ChacoSharp.Assignment.AssignFunc;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Assignment.Y2X;
using static ChacoSharp.Eigen.Splarax;
using static ChacoSharp.Eigen.Warnings;
using static ChacoSharp.Utilities.Dot;
using static ChacoSharp.Utilities.Norm;
using static ChacoSharp.Utilities.Scadd;
using static ChacoSharp.Utilities.VecScale;
using static ChacoSharp.Symmlq.SymmLq;

namespace ChacoSharp.Eigen
{
    public static unsafe class Rqi
    {
        public static void rqi(vtx_data** A, /* matrix/graph being analyzed */
            double*[] yvecs, /* eigenvectors to be refined */
            int index, /* index of vector in yvecs to be refined */
            int n, /* number of rows/columns in matrix */
            double* r1, double* r2, double* v, double* w, double* x, double* y,
            double* work, /* work space for symmlq */
            double tol, /* error tolerance in eigenpair */
            double initshift, /* initial shift */
            double* evalest, /* returned eigenvalue */
            double* vwsqrt, /* square roots of vertex weights */
            orthlink* orthlist, /* lower evecs to orthogonalize against */
            bool cube_or_mesh, /* 0 => hypercube, d => d-dimensional mesh */
            int nsets, /* number of sets to divide into */
            int* assignment, /* set number of each vtx (length n+1) */
            int* active, /* space for nvtxs integers */
            MappingType mediantype, /* which partitioning strategy to use */
            double[] goal, /* desired set sizes */
            int vwgt_max, /* largest vertex weight */
            int ndims /* dimensionality of partition */
        )
        {
            int rqisteps; /* # rqi rqisteps */
            double res; /* convergence quant for rqi */
            double last_res; /* res on previous rqi step */
            double macheps; /* machine precision calculated by symmlq */
            double normxlim; /* a stopping criteria for symmlq */
            double normx; /* norm of the solution vector */
            int symmlqitns; /* # symmlq itns */
            int inv_it_steps; /* initial steps of inverse iteration */
            int itnmin; /* symmlq input */
            double shift, rtol; /* symmlq input */
            bool precon, goodb; /* symmlq input */
            int nout; /* symmlq input */
            bool checka; /* symmlq input */
            int intlim; /* symmlq input */
            double anorm, acond; /* symmlq output */
            double rnorm, ynorm; /* symmlq output */
            int istop, itn; /* symmlq output */
            int long_n; /* copy of n for passing to symmlq */
            bool warning; /* warning on possible misconvergence */
            double factor; /* ratio between previous res and new tol */
            double minfactor; /* minimum acceptable value of factor */
            bool converged; /* has process converged yet? */
            double* u; /* name of vector being refined */
            int* old_assignment = null; /* previous assignment vector */
            int* assgn_pntr; /* pntr to assignment vector */
            int* old_assgn_pntr; /* pntr to previous assignment vector */
            int assigndiff = 0; /* discrepancies between old and new assignment */
            int assigntol = 0; /* tolerance on convergence of assignment vector */
            bool first; /* is this the first RQI step? */
            int i; /* loop index */

            if (DEBUG_TRACE)
            {
                Trace.WriteLine("<Entering rqi>");
            }

            /* Initialize RQI loop */
            u = yvecs[index];
            splarax(y, A, n, u, vwsqrt, r1);
            shift = dot(u, 1, n, y);
            scadd(y, 1, n, -shift, u);
            res = ch_norm(y, 1, n); /* eigen-residual */
            rqisteps = 0; /* a counter */
            symmlqitns = 0; /* a counter */

            /* Set invariant symmlq parameters */
            precon = false; /* false until we figure out a good way */
            goodb = true; /* should be true for this application */
            nout = 0; /* set to 0 for no Symmlq output; 6 for lots */
            checka = false; /* if don't know by now, too bad */
            intlim = n; /* set to enforce a maximum number of Symmlq itns */
            itnmin = 0; /* set to enforce a minimum number of Symmlq itns */
            long_n = n; /* type change for alint */

            if (DEBUG_EVECS > 0)
            {
                Trace.WriteLine($"Using RQI/Symmlq refinement on graph with {n:D} vertices.");
            }

            if (DEBUG_EVECS > 1)
            {
                Trace.WriteLine("  step      lambda est.            Ares          Symmlq its.   istop  factor  delta");
                Trace.Write("    0");
                doubleout(shift, 1);
                doubleout(res, 1);
                Trace.WriteLine("");
            }

            if (RQI_CONVERGENCE_MODE == 1)
            {
                assigntol = (int) (tol * n);
                old_assignment = (int*) Marshal.AllocHGlobal((n + 1) * sizeof(int));
            }

            /* Perform RQI */
            inv_it_steps = 2;
            warning = false;
            factor = 10;
            minfactor = factor / 2;
            first = true;
            if (res < tol)
            {
                converged = true;
            }
            else
            {
                converged = false;
            }

            while (!converged)
            {
                if (res / tol < 1.2)
                {
                    factor = Math.Max(factor / 2, minfactor);
                }

                rtol = res / factor;

                /* exit Symmlq if iterate is this large */
                normxlim = 1.0 / rtol;

                if (rqisteps < inv_it_steps)
                {
                    shift = initshift;
                }

                symmlq(long_n, &u[1], &r1[1], &r2[1], &v[1], &w[1], &x[1], &y[1], work, checka, goodb,
                    precon, &shift, &nout, &intlim, &rtol, &istop, &itn, &anorm, &acond, &rnorm, &ynorm,
                    (double*) A, vwsqrt, (double*) orthlist, &macheps, &normxlim, &itnmin);
                symmlqitns += itn;
                normx = ch_norm(x, 1, n);
                vecscale(u, 1, n, 1.0 / normx, x);
                splarax(y, A, n, u, vwsqrt, r1);
                shift = dot(u, 1, n, y);
                scadd(y, 1, n, -shift, u);
                last_res = res;
                res = ch_norm(y, 1, n);
                if (res > last_res)
                {
                    warning = true;
                }

                rqisteps++;

                if (res < tol)
                {
                    converged = true;
                }

                if (RQI_CONVERGENCE_MODE == 1 && !converged && ndims == 1)
                {
                    if (first)
                    {
                        Assign(A, yvecs, n, 1, cube_or_mesh, nsets, vwsqrt, assignment, active, mediantype, goal,
                            vwgt_max);
                        x2y(yvecs, ndims, n, vwsqrt);
                        first = false;
                        assigndiff = n; /* dummy value for debug chart */
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

                        Assign(A, yvecs, n, ndims, cube_or_mesh, nsets, vwsqrt, assignment, active, mediantype,
                            goal, vwgt_max);
                        x2y(yvecs, ndims, n, vwsqrt);

                        /* count differences in assignment */
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
                        if (assigndiff <= assigntol)
                        {
                            converged = true;
                        }
                    }
                }

                if (DEBUG_EVECS > 1)
                {
                    Trace.Write($"   {rqisteps:d}");
                    doubleout(shift, 1);
                    doubleout(res, 1);
                    Trace.Write($"     {itn:d}");
                    Trace.Write($"          {istop:d}");
                    Trace.Write($"      {factor:g}");
                    if (RQI_CONVERGENCE_MODE == 1)
                    {
                        Trace.WriteLine($"     {assigndiff:d}");
                    }
                    else
                    {
                        Trace.WriteLine("");
                    }
                }
            }

            *evalest = shift;

            if (WARNING_EVECS > 0 && warning)
            {
                Trace.WriteLine("WARNING: Residual convergence not monotonic; RQI may have misconverged.");
            }

            if (DEBUG_EVECS > 0)
            {
                Trace.Write("Eval ");
                doubleout(*evalest, 1);
                Trace.WriteLine($"   RQI steps {rqisteps:d},  Symmlq iterations {symmlqitns:d}.\n");
            }

            if (RQI_CONVERGENCE_MODE == 1)
            {
                Marshal.FreeHGlobal((IntPtr) old_assignment);
            }
        }

        /* Perform Rayleigh Quotient Iteration for extended eigenproblem. */
        public static void rqi_ext()
        {
            throw new NotImplementedException("ERROR: rqi_ext called, but not yet implemented. Sorry.");
        }
    }
}
