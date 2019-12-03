using System;
using System.Runtime.InteropServices;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Coarsening.Coarsen;
using static ChacoSharp.Eigen.LanczosFO;
using static ChacoSharp.Eigen.LanczosSelectiveOrthogonalization;
using static ChacoSharp.Eigen.LanczosExt;
using static ChacoSharp.Eigen.LanczosSoFloat;
using static ChacoSharp.Utilities.Norm;
using static ChacoSharp.Utilities.Perturb;
using static ChacoSharp.Utilities.Timer;
using static ChacoSharp.Assignment.Y2X;

namespace ChacoSharp.Eigen
{
    public static unsafe class EigenSolve
    {
        /* Invoke the eigenvector calculation */
        public static void eigensolve(vtx_data** graph, /* graph data structure */
            int nvtxs, /* number of vertices in graph */
            int nedges, /* number of edges in graph */
            double maxdeg, /* largest (weighted) degree of a vertex */
            int vwgt_max, /* largest vertex weight */
            double* vwsqrt, /* sqrt of vertex weights (length nvtxs+1) */
            bool using_vwgts, /* are vertex weights being used? */
            bool useEdgeWeights, /* are edge weights being used? */
            float*[] term_wgts, /* terminal propagation weight vector */
            int igeom, /* geometric dimensionality if given coords */
            float** coords, /* coordinates of vertices */
            double*[] yvecs, /* space for pointing to eigenvectors */
            double[] evals, /* eigenvalues associated with eigenvectors */
            bool architecture, /* 0 => hypercube, d => d-dimensional mesh */
            int* assignment, /* set number of each vtx (length n+1) */
            double[] goal, /* desired set sizes */
            LanczosType solver_flag, /* flag indicating which solver to use */
            bool rqi_flag, /* use multi-level techniques? */
            int vmax, /* if so, how many vtxs to coarsen down to? */
            int ndims, /* number of eigenvectors (2^d sets) */
            MappingType mediantype, /* which partitioning strategy to use */
            double eigtol /* tolerance on eigenvectors */
        )
        {
            double[] bound = new double[MAXDIMS + 1]; /* ritz approx bounds to eigenpairs */
            double time; /* time marker */
            float*[] dummy_twgt = new float*[2]; /* turns off terminal propagation */
            float* twptr; /* terminal propagation weight vector */
            int* active; /* space for nvtxs values */
            int step; /* current step in RQI counting */
            int nstep; /* number of uncoarsening levels between RQIs */
            int version; /* which version of sel. orth. to use */
            int nsets = 0; /* number of sets to divide into */
            double* g; /* rhs n-vector in the extended eigenproblem */
            double* ptr; /* loops through yvec */
            double w1, w2; /* desired weights of two sets */
            double term_tot; /* sum of terminal weights */
            double sigma; /* norm constraint on extended eigenvector */
            int i, j; /* loop counter */
            bool normal; /* use normal or extended eigensolver? */
            bool autoset_maxitns; /* set LANCZOS_MAXITNS automatically? */
            int prev_maxitns = 0; /* LANCZOS_MAXITNS value above this routine */
            bool autoset_srestol; /* set SRESTOL automatically? */
            double prev_srestol = 0; /* SRESTOL value above this routine */

            if (DEBUG_TRACE)
            {
                Console.WriteLine("<Entering eigensolve, nvtxs = {0:d}, nedges ={1:d}>\n", nvtxs, nedges);
            }

            if (nvtxs <= ndims)
            {
                /* Pathological special case. */
                for (i = 1; i <= ndims; i++)
                {
                    for (j = 1; j <= nvtxs; j++)
                    {
                        yvecs[i][j] = 0;
                    }
                }

                return;
            }

            active = null;

            /* Autoset (if necessary) some parameters for the eigen calculation */
            autoset_maxitns = false;
            autoset_srestol = false;
            if (LANCZOS_MAXITNS < 0)
            {
                autoset_maxitns = true;
                prev_maxitns = LANCZOS_MAXITNS;
                LANCZOS_MAXITNS = 2 * nvtxs;
            }

            if (SRESTOL < 0)
            {
                autoset_srestol = true;
                prev_srestol = SRESTOL;
                SRESTOL = eigtol * eigtol;
            }

            /* Note: When (if ever) rqi_ext is done, change precedence of eigensolvers. */

            if (term_wgts[1] != null && ndims == 1)
            {
                /* then use lanczos_ext */
                if (PERTURB)
                {
                    if (NPERTURB > 0 && PERTURB_MAX > 0.0)
                    {
                        perturb_init(nvtxs);
                        if (DEBUG_PERTURB)
                        {
                            Console.WriteLine("Matrix being perturbed with scale {0:E}\n", PERTURB_MAX);
                        }
                    }
                    else if (DEBUG_PERTURB)
                    {
                        Console.WriteLine("Matrix not being perturbed\n");
                    }
                }

                version = 2;
                if (LANCZOS_CONVERGENCE_MODE == 1)
                {
                    active = (int*) Marshal.AllocHGlobal(nvtxs * sizeof(int));
                }

                w1 = goal[0];
                w2 = goal[1];
                sigma = Math.Sqrt(4 * w1 * w2 / (w1 + w2));
                g = (double*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(double));

                twptr = term_wgts[1];
                term_tot = 0;
                for (i = 1; i <= nvtxs; i++)
                {
                    term_tot += twptr[i];
                }

                term_tot /= (w1 + w2);
                if (using_vwgts)
                {
                    for (i = 1; i <= nvtxs; i++)
                    {
                        g[i] = twptr[i] / graph[i]->vwgt - term_tot;
                    }
                }
                else
                {
                    for (i = 1; i <= nvtxs; i++)
                    {
                        g[i] = twptr[i] - term_tot;
                    }
                }

                time = seconds();

                if (LANCZOS_SO_PRECISION == 2)
                {
                    /* double precision */
                    normal = lanczos_ext(graph, nvtxs, ndims, yvecs, eigtol, vwsqrt, maxdeg, version, g, sigma);
                }
                else
                {
                    /* single precision */
                    normal = lanczos_ext_float(graph, nvtxs, ndims, yvecs, eigtol, vwsqrt, maxdeg, version, g, sigma);
                }

                Marshal.FreeHGlobal((IntPtr) g);
                if (active != null)
                {
                    Marshal.FreeHGlobal((IntPtr) (active));
                }

                active = null;

                if (normal)
                {
                    if (WARNING_EVECS > 2)
                    {
                        Console.WriteLine("WARNING: Not an extended eigenproblem; switching to standard eigensolver.\n");
                    }
                }
                else
                {
                    if (w2 != w1)
                    {
                        if (using_vwgts)
                        {
                            y2x(yvecs, ndims, nvtxs, vwsqrt);
                        }

                        sigma = (w2 - w1) / (w2 + w1);
                        ptr = yvecs[1];
                        for (i = nvtxs; i != 0; i--)
                        {
                            *(++ptr) += sigma;
                        }

                        /* Note: if assign() could skip scaling, next call unnecessary. */
                        if (using_vwgts)
                        {
                            x2y(yvecs, ndims, nvtxs, vwsqrt);
                        }
                    }
                }

                if (PERTURB && NPERTURB > 0 && PERTURB_MAX > 0.0)
                {
                    perturb_clear();
                }

                lanczos_time += seconds() - time;
            }
            else
            {
                normal = true;
            }

            if (normal)
            {
                if (rqi_flag)
                {
                    /* Solve using multi-level scheme RQI/Symmlq. */
                    time = seconds();
                    nstep = COARSE_NLEVEL_RQI;
                    step = 0;
                    dummy_twgt[1] = null;
                    coarsen(graph, nvtxs, nedges, using_vwgts, useEdgeWeights, dummy_twgt, igeom, coords, yvecs,
                        ndims, solver_flag, vmax, eigtol, nstep, step, false);

                    rqi_symmlq_time += seconds() - time;
                }

                else
                {
                    /* Use standard Lanczos. */
                    if (PERTURB)
                    {
                        if (NPERTURB > 0 && PERTURB_MAX > 0.0)
                        {
                            perturb_init(nvtxs);
                            if (DEBUG_PERTURB)
                            {
                                Console.WriteLine("Matrix being perturbed with scale {0:E}\n", PERTURB_MAX);
                            }
                        }
                        else if (DEBUG_PERTURB)
                        {
                            Console.WriteLine("Matrix not being perturbed\n");
                        }
                    }

                    if (solver_flag == LanczosType.FullOrthogonalization)
                    {
                        time = seconds();
                        version = 1;
                        lanczos_FO(graph, nvtxs, ndims, yvecs, evals, bound, eigtol, vwsqrt, maxdeg, version);
                        lanczos_time += seconds() - time;
                    }

                    if (solver_flag == LanczosType.FullOrthogonalizationInverseOperator)
                    {
                        time = seconds();
                        version = 2;
                        lanczos_FO(graph, nvtxs, ndims, yvecs, evals, bound, eigtol, vwsqrt, maxdeg, version);
                        lanczos_time += seconds() - time;
                    }
                    else if (solver_flag == LanczosType.SelectiveOrthogonalization)
                    {
                        version = 2; /* orthog. against left end only */
                        if (LANCZOS_CONVERGENCE_MODE == 1)
                        {
                            active = (int*) Marshal.AllocHGlobal(nvtxs * sizeof(int));
                        }

                        nsets = 1 << ndims;
                        time = seconds();
                        if (LANCZOS_SO_PRECISION == 2)
                        {
                            /* double precision */
                            lanczos_SO(graph, nvtxs, ndims, yvecs, evals, bound, eigtol, vwsqrt, maxdeg, version,
                                architecture, nsets, assignment, active, mediantype, goal, vwgt_max);
                        }
                        else
                        {
                            /* single precision */
                            lanczos_SO_float(graph, nvtxs, ndims, yvecs, evals, bound, eigtol, vwsqrt, maxdeg,
                                version, architecture, nsets, assignment, active, mediantype, goal,
                                vwgt_max);
                        }

                        lanczos_time += seconds() - time;
                    }
                    else if (solver_flag == LanczosType.SelectiveOrthogonalizationDoubleEnded)
                    {
                        if (EXPERT)
                        {
                            version = 1; /* orthog. against both ends */
                        }
                        else
                        {
                            /* this should have been caught earlier ... */
                            version = 2;
                        }

                        if (LANCZOS_CONVERGENCE_MODE == 1)
                        {
                            active = (int*) Marshal.AllocHGlobal(nvtxs * sizeof(int));
                        }

                        nsets = 1 << ndims;
                        time = seconds();
                        if (LANCZOS_SO_PRECISION == 1)
                        {
                            /* Single precision */
                            lanczos_SO_float(graph, nvtxs, ndims, yvecs, evals, bound, eigtol, vwsqrt, maxdeg,
                                version, architecture, nsets, assignment, active, mediantype, goal,
                                vwgt_max);
                        }
                        else
                        {
                            /* Double precision */
                            lanczos_SO(graph, nvtxs, ndims, yvecs, evals, bound, eigtol, vwsqrt, maxdeg, version,
                                architecture, nsets, assignment, active, mediantype, goal, vwgt_max);
                        }

                        lanczos_time += seconds() - time;
                    }
                }

                /*
                file = fopen("CHACO.EVECS", "w");
                for (i = 1; i <= nvtxs; i++) {
                 for (j = 1; j <= ndims; j++) {
                  fprintf(file, "%g ", (yvecs[j])[i]);
                 }
                  fprintf(file, "\n");
                }
                fclose(file);
                */

                if (PERTURB && NPERTURB > 0 && PERTURB_MAX > 0.0)
                {
                    perturb_clear();
                }
            }

            /* This is an attempt to reduce some machine-to-machine
             * variance. If the first value in the eigenvector is negative,
             * reflect the eigenvector...  This may not be needed following
             * the addition of the standard random number generator in util/random.c
             */
            for (nstep = 1; nstep <= ndims; nstep++)
            {
                vecnorm(yvecs[nstep], 1, nvtxs);
            }

            if (DEBUG_EVECS > 4)
            {
                for (nstep = 1; nstep <= ndims; nstep++)
                {
                    vecout(yvecs[nstep], 1, nvtxs, "Eigenvector", null);
                }
            }

            /* Auto-reset (if necessary) some parameters for the eigen calculation */
            if (autoset_maxitns)
            {
                LANCZOS_MAXITNS = prev_maxitns;
            }

            if (autoset_srestol)
            {
                SRESTOL = prev_srestol;
            }

            if (active != null)
            {
                Marshal.FreeHGlobal((IntPtr) active);
            }

            if (DEBUG_TRACE)
            {
                Console.WriteLine("<Leaving eigensolve>\n");
            }
        }
    }
}
