using System;
using System.Diagnostics;
using System.Runtime.InteropServices;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Assignment.AssignFunc;
using static ChacoSharp.Utilities.MakeVwSqrt;
using static ChacoSharp.Coarsening.Coarsen1;
using static ChacoSharp.Graph.FindMaxDeg;
using static ChacoSharp.Connect.Connected;
using static ChacoSharp.Graph.CountWeights;
using static ChacoSharp.Graph.FreeGraph;
using static ChacoSharp.Eigen.EigenSolve;
using static ChacoSharp.Coarsening.KlSpiff;

namespace ChacoSharp.Coarsening
{
    public static unsafe class CoarsenKl
    {

        public static void coarsen_kl(
            /* Coarsen until nvtxs < vmax, compute and uncoarsen. */
            vtx_data** graph, /* array of vtx data for graph */
            int nvtxs, /* number of vertices in graph */
            int nedges, /* number of edges in graph */
            bool using_vwgts, /* are vertices weights being used? */
            bool useEdgeWeights, /* are edge weights being used? */
            float*[] term_wgts, /* weights for terminal propagation */
            int igeom, /* dimension for geometric information */
            float** coords, /* coordinates for vertices */
            int vwgt_max, /* largest vertex weight */
            int* assignment, /* processor each vertex gets assigned to */
            double[] goal, /* desired set sizes */
            bool architecture, /* 0 => hypercube, d => d-dimensional mesh */
            int[][] hops, /* cost of edge between sets */
            LanczosType solver_flag, /* which eigensolver to use */
            int ndims, /* number of eigenvectors to calculate */
            int nsets, /* number of sets being divided into */
            int vmax, /* largest subgraph to stop coarsening */
            MappingType mediantype, /* flag for different assignment strategies */
            bool mkconnected, /* make graph connected before eigensolver? */
            double eigtol, /* tolerance in eigen calculation */
            int nstep, /* number of coarsenings between RQI steps */
            int step, /* current step number */
            int** pbndy_list, /* pointer to returned boundary list */
            double[] weights, /* weights of vertices in each set */
            bool give_up /* has coarsening bogged down? */
        )
        {
            connect_data* cdata; /* data structure for enforcing connectivity */
            vtx_data** cgraph; /* array of vtx data for coarsened graph */
            double*[] yvecs = new double*[MAXDIMS + 1]; /* eigenvectors for subgraph */
            double[] evals = new double[MAXDIMS + 1]; /* eigenvalues returned */
            double[] new_goal = new double[MAXSETS]; /* new goal if not using vertex weights */
            double[] real_goal; /* chooses between goal and new_goal */
            double goal_weight; /* total weight of vertices in goal */
            double* vwsqrt = null; /* square root of vertex weights */
            double maxdeg; /* maximum weighted degree of a vertex */
            double[] temp_goal = new double[2]; /* goal to simulate bisection while striping */
            double[] fake_goal; /* either goal or temp_goal */
            float*[] cterm_wgts = new float*[MAXSETS]; /* terminal weights for coarse graph */
            float*[] new_term_wgts = new float*[MAXSETS]; /* modified for Bui's method */
            float*[] real_term_wgts; /* which of previous two to use */
            float ewgt_max; /* largest edge weight in graph */
            float* twptr; /* loops through term_wgts */
            float* twptr_save; /* copy of twptr */
            float* ctwptr; /* loops through cterm_wgts */
            float** ccoords; /* coarse graph coordinates */
            int* active; /* space for assign routine */
            int* v2cv; /* mapping from fine to coarse vertices */
            int* cv2v; /* mapping from coarse to fine vertices */
            int* bndy_list; /* list of vertices on boundary */
            int* cbndy_list; /* list of vertices of coarse graph on boundary */
            int* mflag; /* flags indicating matching */
            int* cassignment; /* set assignments for coarsened vertices */
            int clist_length; /* length of coarse graph boundary vtx list */
            int list_length; /* length of boundary vtx list */
            int vtx; /* vertex in graph */
            int cvtx; /* vertex in coarse graph */
            int cnvtxs; /* number of vertices in coarsened graph */
            int cnedges; /* number of edges in coarsened graph */
            int cvwgt_max; /* largest vertex weight in coarsened graph */
            int nextstep; /* next step in RQI test */
            int max_dev; /* largest allowed deviation from balance */
            int set; /* set a vertex is in */
            int i, j; /* loop counters */

            if (DEBUG_COARSEN || DEBUG_TRACE)
            {
                Trace.WriteLine($"<Entering coarsen_kl, {nameof(step)}={step:d}, {nameof(nvtxs)}={nvtxs:d}, {nameof(nedges)}={nedges:d}, {nameof(vmax)}={vmax:d}>");
            }

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

                /* Create space for subgraph yvecs. */
                for (i = 1; i <= ndims; i++)
                {
                    yvecs[i] = (double*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(double));
                }

                active = (int*) Marshal.AllocHGlobal(nvtxs * sizeof(int));
                if (mkconnected)
                {
                    make_connected(graph, nvtxs, &nedges, (int*) &(yvecs[1][0]), active, &cdata, useEdgeWeights);
                    if (DEBUG_CONNECTED)
                    {
                        Trace.WriteLine("Enforcing connectivity on coarse graph");
                        print_connected(cdata);
                    }
                }

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

                if (!COARSEN_VWGTS && step != 0)
                {
                    /* Construct new goal */
                    goal_weight = 0;
                    for (i = 0; i < nsets; i++)
                    {
                        goal_weight += goal[i];
                    }

                    for (i = 0; i < nsets; i++)
                    {
                        new_goal[i] = goal[i] * (nvtxs / goal_weight);
                    }

                    real_goal = new_goal;
                }
                else
                {
                    real_goal = goal;
                }

                if (ndims == 1 && nsets > 2)
                {
                    /* Striping. */
                    mediantype = MappingType.Striped;
                    temp_goal[0] = temp_goal[1] = 0;
                    for (i = 0; 2 * i + 1 < nsets; i++)
                    {
                        temp_goal[0] += real_goal[i];
                        temp_goal[1] += real_goal[nsets - 1 - i];
                    }

                    i = nsets / 2;
                    if (2 * i != nsets)
                    {
                        temp_goal[0] += .5 * real_goal[i];
                        temp_goal[1] += .5 * real_goal[i];
                    }

                    fake_goal = temp_goal;
                }
                else
                {
                    mediantype = MAPPING_TYPE;
                    fake_goal = real_goal;
                }

                /* Partition coarse graph with spectral method */
                maxdeg = find_maxdeg(graph, nvtxs, useEdgeWeights, &ewgt_max);
                eigensolve(graph, nvtxs, nedges, maxdeg, vwgt_max, vwsqrt, using_vwgts, useEdgeWeights,
                    real_term_wgts, 0, (float**) null, yvecs, evals, architecture, assignment, fake_goal,
                    solver_flag, false, 0, ndims, mediantype, eigtol);

                if (mkconnected)
                {
                    make_unconnected(graph, &nedges, &cdata, useEdgeWeights);
                }

                Assign(graph, yvecs, nvtxs, ndims, architecture, nsets, vwsqrt, assignment, active, mediantype,
                    real_goal, vwgt_max);
                /*
                simple_part(graph, nvtxs, assignment, nsets, 1, real_goal);
                */
                Marshal.FreeHGlobal((IntPtr) active);

                count_weights(graph, nvtxs, assignment, nsets, weights, (vwgt_max != 1));

                bndy_list = null;
                if (COARSE_KL_BOTTOM || (step % nstep) == 0)
                {
                    if (LIMIT_KL_EWGTS)
                    {
                        compress_ewgts(graph, nvtxs, nedges, ewgt_max, useEdgeWeights);
                    }

                    max_dev = (step == 0) ? vwgt_max : 5 * vwgt_max;
                    goal_weight = 0;
                    for (i = 0; i < nsets; i++)
                    {
                        goal_weight += real_goal[i];
                    }

                    goal_weight *= KL_IMBALANCE / nsets;
                    if (goal_weight > max_dev)
                    {
                        max_dev = (int) goal_weight;
                    }

                    if (KL_ONLY_BNDY)
                    {
                        /* On small graph, quickest to include everybody in list. */
                        bndy_list = (int*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));
                        for (i = 0; i < nvtxs; i++)
                        {
                            bndy_list[i] = i + 1;
                        }

                        bndy_list[nvtxs] = 0;
                    }

                    klspiff(graph, nvtxs, assignment, nsets, hops, real_goal, real_term_wgts, max_dev, maxdeg,
                        useEdgeWeights, &bndy_list, weights);
                    if (LIMIT_KL_EWGTS)
                    {
                        restore_ewgts(graph, nvtxs);
                    }
                }

                else if (KL_ONLY_BNDY)
                {
                    /* Find boundary vertices directly. */
                    bndy_list = (int*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));
                    list_length = 0;
                    for (i = 1; i <= nvtxs; i++)
                    {
                        set = assignment[i];
                        for (j = 1; j < graph[i]->nedges; j++)
                        {
                            if (assignment[graph[i]->edges[j]] != set)
                            {
                                bndy_list[list_length++] = i;
                                break;
                            }
                        }
                    }

                    bndy_list[list_length] = 0;
                    bndy_list = (int*) Marshal.ReAllocHGlobal((IntPtr) bndy_list, new IntPtr((list_length + 1) * sizeof(int)));
                }

                *pbndy_list = bndy_list;

                if (real_term_wgts != term_wgts && new_term_wgts[1] != null)
                {
                    Marshal.FreeHGlobal((IntPtr) real_term_wgts[1]);
                }

                if (vwsqrt != null)
                {
                    Marshal.FreeHGlobal((IntPtr) vwsqrt);
                }

                for (i = ndims; i > 0; i--)
                {
                    Marshal.FreeHGlobal((IntPtr) yvecs[i]);
                }

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

            coarsen1(graph, nvtxs, nedges, &cgraph, &cnvtxs, &cnedges, &v2cv, igeom, coords, ccoords,
                useEdgeWeights);

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

            /* If coarsening isn't working very well, give up and partition. */
            give_up = false;
            if (nvtxs * COARSEN_RATIO_MIN < cnvtxs && cnvtxs > vmax)
            {
                Trace.WriteLine($"WARNING: Coarsening not making enough progress, {nameof(nvtxs)} = {nvtxs:d}, {nameof(cnvtxs)} = {cnvtxs:d}.");
                Trace.WriteLine("         Recursive coarsening being stopped prematurely.");

                give_up = true;
            }

            /* Now recurse on coarse subgraph. */
            if (COARSEN_VWGTS)
            {
                cvwgt_max = 0;
                for (i = 1; i <= cnvtxs; i++)
                {
                    if (cgraph[i]->vwgt > cvwgt_max)
                    {
                        cvwgt_max = cgraph[i]->vwgt;
                    }
                }
            }
            else
            {
                cvwgt_max = 1;
            }

            cassignment = (int*) Marshal.AllocHGlobal((cnvtxs + 1) * sizeof(int));
            nextstep = step + 1;
            coarsen_kl(cgraph, cnvtxs, cnedges, COARSEN_VWGTS, COARSEN_EWGTS, cterm_wgts, igeom, ccoords,
                cvwgt_max, cassignment, goal, architecture, hops, solver_flag, ndims, nsets, vmax,
                mediantype, mkconnected, eigtol, nstep, nextstep, &cbndy_list, weights, give_up);

            /* Interpolate assignment back to fine graph. */
            for (i = 1; i <= nvtxs; i++)
            {
                assignment[i] = cassignment[v2cv[i]];
            }

            if (KL_ONLY_BNDY)
            {
                /* Construct boundary list from coarse boundary list. */
                cv2v = (int*) Marshal.AllocHGlobal((cnvtxs + 1) * sizeof(int));
                mflag = (int*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));
                for (i = 1; i <= cnvtxs; i++)
                {
                    cv2v[i] = 0;
                }

                for (i = 1; i <= nvtxs; i++)
                {
                    mflag[i] = 0;
                    cvtx = v2cv[i];
                    if (cv2v[cvtx] == 0)
                    {
                        cv2v[cvtx] = i;
                    }
                    else
                    {
                        mflag[cv2v[cvtx]] = i;
                    }
                }

                clist_length = 0;
                while (cbndy_list[clist_length] != 0)
                {
                    clist_length++;
                }

                bndy_list = (int*) Marshal.AllocHGlobal((2 * clist_length + 1) * sizeof(int));

                list_length = 0;
                for (i = 0; i < clist_length; i++)
                {
                    vtx = cv2v[cbndy_list[i]];
                    bndy_list[list_length++] = vtx;
                    if (mflag[vtx] != 0)
                    {
                        bndy_list[list_length++] = mflag[vtx];
                    }
                }

                bndy_list[list_length] = 0;
                bndy_list = (int*) Marshal.ReAllocHGlobal((IntPtr) bndy_list, new IntPtr((list_length + 1) * sizeof(int)));

                Marshal.FreeHGlobal((IntPtr) mflag);
                Marshal.FreeHGlobal((IntPtr) cv2v);
                Marshal.FreeHGlobal((IntPtr) cbndy_list);
            }
            else
            {
                bndy_list = null;
            }

            /* Free the space that was allocated. */
            Marshal.FreeHGlobal((IntPtr) cassignment);
            if (cterm_wgts[1] != null)
            {
                Marshal.FreeHGlobal((IntPtr) cterm_wgts[1]);
            }

            free_graph(cgraph);
            Marshal.FreeHGlobal((IntPtr) v2cv);

            /* Smooth using KL every nstep steps. */
            if ((step % nstep) == 0)
            {
                if (!COARSEN_VWGTS && step != 0)
                {
                    /* Construct new goal */
                    goal_weight = 0;
                    for (i = 0; i < nsets; i++)
                    {
                        goal_weight += goal[i];
                    }

                    for (i = 0; i < nsets; i++)
                    {
                        new_goal[i] = goal[i] * (nvtxs / goal_weight);
                    }

                    real_goal = new_goal;
                }
                else
                {
                    real_goal = goal;
                }

                maxdeg = find_maxdeg(graph, nvtxs, useEdgeWeights, &ewgt_max);
                if (LIMIT_KL_EWGTS)
                {
                    compress_ewgts(graph, nvtxs, nedges, ewgt_max, useEdgeWeights);
                }

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

                max_dev = (step == 0) ? vwgt_max : 5 * vwgt_max;
                goal_weight = 0;
                for (i = 0; i < nsets; i++)
                {
                    goal_weight += real_goal[i];
                }

                goal_weight *= KL_IMBALANCE / nsets;
                if (goal_weight > max_dev)
                {
                    max_dev = (int) goal_weight;
                }

                if (!COARSEN_VWGTS)
                {
                    count_weights(graph, nvtxs, assignment, nsets, weights, (vwgt_max != 1));
                }

                klspiff(graph, nvtxs, assignment, nsets, hops, real_goal, real_term_wgts, max_dev, maxdeg,
                    useEdgeWeights, &bndy_list, weights);

                if (real_term_wgts != term_wgts && new_term_wgts[1] != null)
                {
                    Marshal.FreeHGlobal((IntPtr) real_term_wgts[1]);
                }

                if (LIMIT_KL_EWGTS)
                {
                    restore_ewgts(graph, nvtxs);
                }
            }

            *pbndy_list = bndy_list;

            if (ccoords != null)
            {
                for (i = 0; i < igeom; i++)
                {
                    Marshal.FreeHGlobal((IntPtr) ccoords[i]);
                }

                Marshal.FreeHGlobal((IntPtr) ccoords);
            }

            if (DEBUG_COARSEN)
            {
                Trace.WriteLine($" Leaving coarsen_kl, {nameof(step)}={step:d}");
            }
        }

    }
}
