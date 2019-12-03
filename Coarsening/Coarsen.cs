#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
#pragma warning disable HAA0101 // Array allocation for params parameter
using System;
using System.Runtime.InteropServices;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Coarsening.Coarsen1;
using static ChacoSharp.Coarsening.Interpolate;
using static ChacoSharp.Eigen.EigenSolve;
using static ChacoSharp.Eigen.Orthogonalization;
using static ChacoSharp.Eigen.Rqi;
using static ChacoSharp.Graph.FindMaxDeg;
using static ChacoSharp.Graph.FreeGraph;
using static ChacoSharp.Utilities.Dot;
using static ChacoSharp.Utilities.MakeVwSqrt;
using static ChacoSharp.Utilities.Norm;
using static ChacoSharp.Utilities.Scadd;
using static ChacoSharp.Utilities.Timer;

namespace ChacoSharp.Coarsening
{
    public static unsafe class Coarsen
    {
        public static void coarsen(
            /* Coarsen until nvtxs <= vmax, compute and uncoarsen. */
            vtx_data** graph, /* array of vtx data for graph */
            int nvtxs, /* number of vertices in graph */
            int nedges, /* number of edges in graph */
            bool using_vwgts, /* are vertices weights being used? */
            bool useEdgeWeights, /* are edge weights being used? */
            float*[] term_wgts, /* terminal weights */
            int igeom, /* dimension for geometric information */
            float** coords, /* coordinates for vertices */
            double*[] yvecs, /* eigenvectors returned */
            int ndims, /* number of eigenvectors to calculate */
            LanczosType solver_flag, /* which eigensolver to use */
            int vmax, /* largest subgraph to stop coarsening */
            double eigtol, /* tolerance in eigen calculation */
            int nstep, /* number of coarsenings between RQI steps */
            int step, /* current step number */
            bool give_up /* has coarsening bogged down? */
        )
        {
            vtx_data** cgraph; /* array of vtx data for coarsened graph */
            orthlink* orthlist; /* list of lower evecs to suppress */
            orthlink* newlink; /* lower evec to suppress */
            double*[] cyvecs = new double*[MAXDIMS + 1]; /* eigenvectors for subgraph */
            double[] evals = new double[MAXDIMS + 1]; /* eigenvalues returned */
            double[] goal = new double[MAXSETS]; /* needed for convergence mode = 1 */
            double* r1; /* space needed by symmlq/RQI */
            double* r2; /* space needed by symmlq/RQI */
            double* work; /* space needed by symmlq/RQI */
            double* v; /* space needed by symmlq/RQI */
            double* w; /* space needed by symmlq/RQI */
            double* x; /* space needed by symmlq/RQI */
            double* y; /* space needed by symmlq/RQI */
            double* gvec; /* rhs vector in extended eigenproblem */
            double evalest; /* eigenvalue estimate returned by RQI */
            double maxdeg; /* maximum weighted degree of a vertex */
            float** ccoords; /* coordinates for coarsened graph */
            float*[] cterm_wgts = new float*[MAXSETS]; /* coarse graph terminal weights */
            float*[] new_term_wgts = new float*[MAXSETS]; /* terminal weights for Bui's method*/
            float*[] real_term_wgts; /* one of the above */
            float* twptr = null; /* loops through term_wgts */
            float* twptr_save = null; /* copy of twptr */
            float* ctwptr; /* loops through cterm_wgts */
            double* vwsqrt = null; /* square root of vertex weights */
            double norm, alpha; /* values used for orthogonalization */
            double initshift; /* initial shift for RQI */
            double total_vwgt; /* sum of all the vertex weights */
            double w1, w2; /* weights of two sets */
            double term_tot; /* sum of all terminal weights */
            int* space; /* room for assignment in Lanczos */
            int* morespace; /* room for assignment in Lanczos */
            int* v2cv; /* mapping from vertices to coarse vtxs */
            int vwgt_max; /* largest vertex weight */
            bool oldperturb; /* saves PERTURB value */
            int cnvtxs; /* number of vertices in coarsened graph */
            int cnedges; /* number of edges in coarsened graph */
            int nextstep; /* next step in RQI test */
            int nsets; /* number of sets being created */
            int i, j; /* loop counters */
            double time; /* time marker */

            if (DEBUG_COARSEN)
            {
                Console.WriteLine("<Entering coarsen, step={0:D}, nvtxs={1:D}, nedges={2:D}, vmax={3:D}>\n", step, nvtxs, nedges, vmax);
            }

            nsets = 1 << ndims;

            /* Is problem small enough to solve? */
            if (nvtxs <= vmax || give_up)
            {
                if (using_vwgts)
                {
                    vwsqrt = (double*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(double));
                    makevwsqrt(vwsqrt, graph, nvtxs);
                }
                else
                {
                    vwsqrt = null;
                }

                maxdeg = find_maxdeg(graph, nvtxs, useEdgeWeights, null);

                if (using_vwgts)
                {
                    vwgt_max = 0;
                    total_vwgt = 0;
                    for (i = 1; i <= nvtxs; i++)
                    {
                        if (graph[i]->vwgt > vwgt_max)
                        {
                            vwgt_max = graph[i]->vwgt;
                        }

                        total_vwgt += graph[i]->vwgt;
                    }
                }
                else
                {
                    vwgt_max = 1;
                    total_vwgt = nvtxs;
                }

                for (i = 0; i < nsets; i++)
                {
                    goal[i] = total_vwgt / nsets;
                }

                space = (int*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));

                /* If not coarsening ewgts, then need care with term_wgts. */
                if (!useEdgeWeights && term_wgts[1] != null && step != 0)
                {
                    twptr = (float*) Marshal.AllocHGlobal((nvtxs + 1) * (nsets - 1) * sizeof(float));
                    twptr_save = twptr;
                    for (j = 1; j < nsets; j++)
                    {
                        new_term_wgts[j] = twptr;
                        twptr += nvtxs + 1;
                    }

                    for (j = 1; j < nsets; j++)
                    {
                        twptr = term_wgts[j];
                        ctwptr = new_term_wgts[j];
                        for (i = 1; i <= nvtxs; i++)
                        {
                            if (twptr[i] > .5)
                            {
                                ctwptr[i] = 1;
                            }
                            else if (twptr[i] < -.5)
                            {
                                ctwptr[i] = -1;
                            }
                            else
                            {
                                ctwptr[i] = 0;
                            }
                        }
                    }

                    real_term_wgts = new_term_wgts;
                }
                else
                {
                    real_term_wgts = term_wgts;
                    new_term_wgts[1] = null;
                }

                eigensolve(graph, nvtxs, nedges, maxdeg, vwgt_max, vwsqrt, using_vwgts, useEdgeWeights,
                    real_term_wgts, igeom, coords, yvecs, evals, false, space, goal, solver_flag, false, 0,
                    ndims, MappingType.IndependantMedians, eigtol);

                if (real_term_wgts != term_wgts && new_term_wgts[1] != null)
                {
                    Marshal.FreeHGlobal((IntPtr)real_term_wgts[1]);
                }

                Marshal.FreeHGlobal((IntPtr)space);
                space = null;
                Marshal.FreeHGlobal((IntPtr)vwsqrt);
                vwsqrt = null;
                Marshal.FreeHGlobal((IntPtr)twptr_save);
                twptr_save = null;
                return;
            }

            /* Otherwise I have to coarsen. */
            if (coords != null)
            {
                ccoords = (float**) Marshal.AllocHGlobal(igeom * sizeof(float*));
            }
            else
            {
                ccoords = null;
            }

            coarsen1(graph, nvtxs, nedges, &cgraph, &cnvtxs, &cnedges, &v2cv, igeom, coords, ccoords, useEdgeWeights);

            /* If coarsening isn't working very well, give up and partition. */
            give_up = false;
            if (nvtxs * COARSEN_RATIO_MIN < cnvtxs && cnvtxs > vmax)
            {
                Console.WriteLine("WARNING: Coarsening not making enough progress, nvtxs = {0:d}, cnvtxs = {1:d}.", nvtxs, cnvtxs);
                Console.WriteLine("         Recursive coarsening being stopped prematurely.");
                give_up = true;
            }

            /* Create space for subgraph yvecs. */
            for (i = 1; i <= ndims; i++)
            {
                cyvecs[i] = (double*) Marshal.AllocHGlobal((cnvtxs + 1) * sizeof(double));
            }

            /* Make coarse version of terminal weights. */
            if (term_wgts[1] != null)
            {
                twptr = (float*) Marshal.AllocHGlobal((cnvtxs + 1) * (nsets - 1) * sizeof(float));
                twptr_save = twptr;
                for (i = (cnvtxs + 1) * (nsets - 1); i != 0; i--)
                {
                    *twptr++ = 0;
                }

                twptr = twptr_save;
                for (j = 1; j < nsets; j++)
                {
                    cterm_wgts[j] = twptr;
                    twptr += cnvtxs + 1;
                }

                for (j = 1; j < nsets; j++)
                {
                    ctwptr = cterm_wgts[j];
                    twptr = term_wgts[j];
                    for (i = 1; i < nvtxs; i++)
                    {
                        ctwptr[v2cv[i]] += twptr[i];
                    }
                }
            }
            else
            {
                cterm_wgts[1] = null;
            }

            /* Now recurse on coarse subgraph. */
            nextstep = step + 1;
            coarsen(cgraph, cnvtxs, cnedges, COARSEN_VWGTS, COARSEN_EWGTS, cterm_wgts, igeom, ccoords, cyvecs,
                ndims, solver_flag, vmax, eigtol, nstep, nextstep, give_up);

            ch_interpolate(yvecs, cyvecs, ndims, graph, nvtxs, v2cv, useEdgeWeights);

            Marshal.FreeHGlobal((IntPtr)twptr_save);
            twptr_save = null;
            Marshal.FreeHGlobal((IntPtr)v2cv);
            v2cv = null;

            /* I need to do Rayleigh Quotient Iteration each nstep stages. */
            time = seconds();
            if ((step % nstep) == 0)
            {
                oldperturb = PERTURB;
                PERTURB = false;
                /* Should I do some orthogonalization here against vwsqrt? */
                if (using_vwgts)
                {
                    vwsqrt = (double*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(double));
                    makevwsqrt(vwsqrt, graph, nvtxs);

                    for (i = 1; i <= ndims; i++)
                    {
                        orthogvec(yvecs[i], 1, nvtxs, vwsqrt);
                    }
                }
                else
                {
                    for (i = 1; i <= ndims; i++)
                    {
                        orthog1(yvecs[i], 1, nvtxs);
                    }
                }

                /* Allocate space that will be needed in RQI. */
                r1 = (double*) Marshal.AllocHGlobal(7 * (nvtxs + 1) * sizeof(double));
                r2 = &r1[nvtxs + 1];
                v = &r1[2 * (nvtxs + 1)];
                w = &r1[3 * (nvtxs + 1)];
                x = &r1[4 * (nvtxs + 1)];
                y = &r1[5 * (nvtxs + 1)];
                work = &r1[6 * (nvtxs + 1)];

                if (using_vwgts)
                {
                    vwgt_max = 0;
                    total_vwgt = 0;
                    for (i = 1; i <= nvtxs; i++)
                    {
                        if (graph[i]->vwgt > vwgt_max)
                        {
                            vwgt_max = graph[i]->vwgt;
                        }

                        total_vwgt += graph[i]->vwgt;
                    }
                }
                else
                {
                    vwgt_max = 1;
                    total_vwgt = nvtxs;
                }

                for (i = 0; i < nsets; i++)
                {
                    goal[i] = total_vwgt / nsets;
                }

                space = (int*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));
                morespace = (int*) Marshal.AllocHGlobal((nvtxs) * sizeof(int));

                initshift = 0;
                orthlist = null;
                for (i = 1; i < ndims; i++)
                {
                    ch_normalize(yvecs[i], 1, nvtxs);
                    rqi(graph, yvecs, i, nvtxs, r1, r2, v, w, x, y, work, eigtol, initshift, &evalest, vwsqrt,
                        orthlist, false, nsets, space, morespace, MappingType.IndependantMedians, goal, vwgt_max, ndims);

                    /* Now orthogonalize higher yvecs against this one. */
                    norm = dot(yvecs[i], 1, nvtxs, yvecs[i]);
                    for (j = i + 1; j <= ndims; j++)
                    {
                        alpha = -dot(yvecs[j], 1, nvtxs, yvecs[i]) / norm;
                        scadd(yvecs[j], 1, nvtxs, alpha, yvecs[i]);
                    }

                    /* Now prepare for next pass through loop. */
                    initshift = evalest;
                    newlink = makeorthlnk();
                    newlink->vec = yvecs[i];
                    newlink->pntr = orthlist;
                    orthlist = newlink;
                }

                ch_normalize(yvecs[ndims], 1, nvtxs);

                if (term_wgts[1] != null && ndims == 1)
                {
                    /* Solve extended eigen problem */

                    /* If not coarsening ewgts, then need care with term_wgts. */
                    if (!useEdgeWeights && term_wgts[1] != null && step != 0)
                    {
                        twptr = (float*) Marshal.AllocHGlobal((nvtxs + 1) * (nsets - 1) * sizeof(float));
                        twptr_save = twptr;
                        for (j = 1; j < nsets; j++)
                        {
                            new_term_wgts[j] = twptr;
                            twptr += nvtxs + 1;
                        }

                        for (j = 1; j < nsets; j++)
                        {
                            twptr = term_wgts[j];
                            ctwptr = new_term_wgts[j];
                            for (i = 1; i <= nvtxs; i++)
                            {
                                if (twptr[i] > .5)
                                {
                                    ctwptr[i] = 1;
                                }
                                else if (twptr[i] < -.5)
                                {
                                    ctwptr[i] = -1;
                                }
                                else
                                {
                                    ctwptr[i] = 0;
                                }
                            }
                        }

                        real_term_wgts = new_term_wgts;
                    }
                    else
                    {
                        real_term_wgts = term_wgts;
                        new_term_wgts[1] = null;
                    }

                    /* Following only works for bisection. */
                    w1 = goal[0];
                    w2 = goal[1];
                    gvec = (double*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(double));
                    term_tot = 0;
                    for (j = 1; j <= nvtxs; j++)
                    {
                        term_tot += (real_term_wgts[1])[j];
                    }

                    term_tot /= (w1 + w2);
                    if (using_vwgts)
                    {
                        for (j = 1; j <= nvtxs; j++)
                        {
                            gvec[j] = (real_term_wgts[1])[j] / graph[j]->vwgt - term_tot;
                        }
                    }
                    else
                    {
                        for (j = 1; j <= nvtxs; j++)
                        {
                            gvec[j] = (real_term_wgts[1])[j] - term_tot;
                        }
                    }

                    rqi_ext();

                    Marshal.FreeHGlobal((IntPtr)gvec);
                    gvec = null;
                    if (real_term_wgts != term_wgts && new_term_wgts[1] != null)
                    {
                        Marshal.FreeHGlobal((IntPtr)new_term_wgts[1]);
                        new_term_wgts[1] = null;
                    }
                }
                else
                {
                    rqi(graph, yvecs, ndims, nvtxs, r1, r2, v, w, x, y, work, eigtol, initshift, &evalest, vwsqrt,
                        orthlist, false, nsets, space, morespace, MappingType.IndependantMedians, goal, vwgt_max, ndims);
                }

                refine_time += seconds() - time;

                /* Free the space allocated for RQI. */
                Marshal.FreeHGlobal((IntPtr)morespace);
                Marshal.FreeHGlobal((IntPtr)space);
                while (orthlist != null)
                {
                    newlink = orthlist->pntr;
                    Marshal.FreeHGlobal((IntPtr)orthlist);
                    orthlist = newlink;
                }

                Marshal.FreeHGlobal((IntPtr)r1);
                Marshal.FreeHGlobal((IntPtr)vwsqrt);
                vwsqrt = null;
                PERTURB = oldperturb;
            }

            if (DEBUG_COARSEN)
            {
                Console.WriteLine(" Leaving coarsen, step={0:d}", step);
            }

            Marshal.FreeHGlobal((IntPtr)twptr_save);
            twptr_save = null;

            /* Free the space that was allocated. */
            if (ccoords != null)
            {
                for (i = 0; i < igeom; i++)
                {
                    Marshal.FreeHGlobal((IntPtr)ccoords[i]);
                }

                Marshal.FreeHGlobal((IntPtr)ccoords);
            }

            for (i = ndims; i > 0; i--)
            {
                Marshal.FreeHGlobal((IntPtr)cyvecs[i]);
            }

            free_graph(cgraph);
        }
    }
}
