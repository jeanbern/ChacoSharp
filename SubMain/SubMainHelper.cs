#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
using System;
using System.Diagnostics;
using System.Runtime.InteropServices;
using static ChacoSharp.Assignment.AssignFunc;
using static ChacoSharp.Assignment.MergeAssignments;
using static ChacoSharp.ChacoSequenceJPetit;
using static ChacoSharp.Coarsening.BpmImprove;
using static ChacoSharp.Coarsening.CoarsenKl;
using static ChacoSharp.Coarsening.CoarsenKlv;
using static ChacoSharp.Coarsening.KlSpiff;
using static ChacoSharp.Coarsening.KlvSpiff;
using static ChacoSharp.Connect.Connected;
using static ChacoSharp.Connect.ConnectEnforce;
using static ChacoSharp.Eigen.EigenSolve;
using static ChacoSharp.Graph.CountWeights;
using static ChacoSharp.Graph.FindMaxDeg;
using static ChacoSharp.Graph.SimplePartition;
using static ChacoSharp.Graph.Subgraph;
using static ChacoSharp.Inertial.InertialHelper;
using static ChacoSharp.Input.CheckInput;
using static ChacoSharp.Input.ReflectInput;
using static ChacoSharp.Internal.ForceInternal;
using static ChacoSharp.RefineMap.RefineMapHelper;
using static ChacoSharp.RefinePartition.RefinePartitionHelper;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Utilities.DivideProcs;
using static ChacoSharp.Utilities.MakeMaps;
using static ChacoSharp.Utilities.MakeSetLists;
using static ChacoSharp.Utilities.MakeTermProps;
using static ChacoSharp.Utilities.MakeVwSqrt;
using static ChacoSharp.Utilities.MergeGoals;
using static ChacoSharp.Utilities.Randomize;
using static ChacoSharp.Utilities.Timer;

namespace ChacoSharp.SubMain
{
    public static unsafe class SubMainHelper
    {
        private static double* SQRTS; /* precomputed square roots for efficiency */

        /// <summary>
        ///
        /// </summary>
        /// <param name="graph">data structure for graph</param>
        /// <param name="nvtxs">number of vertices in full graph</param>
        /// <param name="nedges">number of edges in graph</param>
        /// <param name="using_vwgts">are vertex weights being used?</param>
        /// <param name="useEdgeWeights">are edge weights being used?</param>
        /// <param name="igeom">geometry dimension if using inertial method</param>
        /// <param name="coords">coordinates of vertices if used</param>
        /// <param name="assignment">set number of each vtx (length n)</param>
        /// <param name="goal">desired sizes for each set</param>
        /// <param name="architecture">0=> hypercube, d=> d-dimensional mesh</param>
        /// <param name="ndims_tot">total number hypercube dimensions</param>
        /// <param name="mesh_dims">extent of mesh in 3 directions</param>
        /// <param name="global_method">global partitioning algorithm </param>
        /// <param name="localPartitioningStrategy">local partitioning algorithm</param>
        /// <param name="rqi_flag">use RQI/Symmlq eigensolver?</param>
        /// <param name="vmax">if so, how many vtxs to coarsen down to</param>
        /// <param name="ndims">number of eigenvectors (2^d sets)</param>
        /// <param name="eigtol">tolerance on eigenvectors</param>
        /// <param name="seed">for random graph mutations</param>
        /// <returns></returns>
        public static bool submain(vtx_data** graph,
            int nvtxs,
            int nedges,
            bool using_vwgts,
            bool useEdgeWeights,
            int igeom,
            float** coords,
            int* assignment,
            double[] goal,
            int architecture,
            int ndims_tot,
            int[] mesh_dims /*[3]*/,
            PartitioningStrategy global_method,
            LocalPartitioningStrategy localPartitioningStrategy,
            bool rqi_flag,
            int vmax,
            int ndims,
            double eigtol,
            int seed
        )
        {
            vtx_data** graph2; /* data structure for graph */
            int[][] hop_mtx = new int[MAXSETS][] /*[MAXSETS][MAXSETS]*/; /* between-set hop cost for KL */
            int i;
            for (i = 0; i < hop_mtx.Length; i++)
            {
                hop_mtx[i] = new int[MAXSETS];
            }

            double* vwsqrt; /* sqrt of vertex weights (length nvtxs+1) */
            double time, time1; /* timing variables */
            string geomname; /* names of input files */
            string inassignname; /* name of assignment input file */
            int old_nsqrts; /* old value of NSQRTS */
            int nsets; /* number of sets created by each divide */
            int nsets_tot; /* total number of sets */
            int bits; /* used in computing hops */
            bool flag; /* return code from check_input */
            bool old_perturb = false; /* saves original perturbation flag */

            if (DEBUG_TRACE)
            {
                Trace.WriteLine("<Entering submain>");
            }

            /* First check all the input for consistency. */

            if (architecture == 1)
            {
                mesh_dims[1] = mesh_dims[2] = 1;
            }
            else if (architecture == 2)
            {
                mesh_dims[2] = 1;
            }

            /* Check for simple special case of 1 processor. */
            {
                var k = 0;
                if (architecture == 0)
                {
                    k = 1 << ndims_tot;
                }
                else if (architecture > 0)
                {
                    k = mesh_dims[0] * mesh_dims[1] * mesh_dims[2];
                }

                if (k == 1)
                {
                    for (i = 1; i <= nvtxs; i++)
                    {
                        assignment[i] = 0;
                    }

                    return false;
                }
            }
            geomname = Geometry_File_Name;
            inassignname = Assign_In_File_Name;

            /* Turn of perturbation if using bisection */
            if (ndims == 1)
            {
                old_perturb = PERTURB;
                PERTURB = false;
            }

            if (PRINT_HEADERS)
            {
                Trace.WriteLine("\n                    Chaco 2.0");
                Trace.WriteLine("          Sandia National Laboratories\n");
            }

            if (CHECK_INPUT)
            {
                /* Check the input for inconsistencies. */
                time1 = seconds();

                flag = check_input(graph, nvtxs, nedges, igeom, coords, assignment, goal,
                    architecture, ndims_tot, mesh_dims, global_method, localPartitioningStrategy, rqi_flag,
                    &vmax, ndims, eigtol);

                check_input_time += seconds() - time1;

                if (flag)
                {
                    Trace.WriteLine("ERROR IN INPUT.");
                    return (true);
                }
            }

            if (ECHO_INPUT_PARAMETERS)
            {
                reflect_input(nvtxs, nedges, igeom, geomname, inassignname,
                    architecture, ndims_tot, mesh_dims, global_method, localPartitioningStrategy,
                    rqi_flag, vmax, ndims, eigtol, seed);
            }

            if (PRINT_HEADERS)
            {
                Trace.WriteLine("\n\nStarting to partition ...\n");
            }

            time = seconds();

            /* Perform some one-time initializations. */
            setrandom(seed);

            if (DEBUG_MACH_PARAMS)
            {
                Trace.WriteLine("Machine parameters:");
                Trace.WriteLine($"  DOUBLE_EPSILON = {DOUBLE_EPSILON:e}");
                Trace.WriteLine($"  DOUBLE_MAX = {double.MaxValue:e}");
            }

            nsets = (1 << ndims);

            old_nsqrts = NSQRTS;
            if (nvtxs < NSQRTS && !using_vwgts)
            {
                NSQRTS = nvtxs;
            }

            SQRTS = (double*) Marshal.AllocHGlobal((NSQRTS + 1) * sizeof(double));
            if (SQRTS == null)
            {
                Trace.WriteLine("ERROR: No space to allocate sqrts");
                return (true);
            }

            for (i = 1; i <= NSQRTS; i++)
            {
                SQRTS[i] = Math.Sqrt((double) i);
            }

            if (using_vwgts && (global_method == PartitioningStrategy.Multilevel_KL || global_method == PartitioningStrategy.Spectral))
            {
                vwsqrt = (double*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(double));
                if (vwsqrt == null)
                {
                    Trace.WriteLine("ERROR: No space to allocate vwsqrt");
                    Marshal.FreeHGlobal((IntPtr) SQRTS);
                    NSQRTS = old_nsqrts;
                    return (true);
                }

                makevwsqrt(vwsqrt, graph, nvtxs);
            }
            else
            {
                vwsqrt = null;
            }

            if (TIME_KERNELS)
            {
                time1 = seconds();
                time_kernels(graph, nvtxs, vwsqrt);
                kernel_time += seconds() - time1;
            }

            if (SEQUENCE)
            {
                sequence(graph, nvtxs, nedges, useEdgeWeights, vwsqrt, LANCZOS_TYPE, rqi_flag, vmax, eigtol);
                goto End_Label;
            }

            /* Initialize cost function for KL-spiff */
            if (global_method == PartitioningStrategy.Multilevel_KL || global_method == PartitioningStrategy.Spectral)
            {
                for (i = 0; i < nsets; i++)
                {
                    hop_mtx[i][i] = 0;
                    for (var j = 0; j < i; j++)
                    {
                        if (KL_METRIC == KernighanLinMetric.Hops)
                        {
                            /* Count hypercube hops */
                            hop_mtx[i][j] = 0;
                            bits = i ^ j;
                            while (bits != 0)
                            {
                                if ((bits & 1) != 0)
                                {
                                    ++hop_mtx[i][j];
                                }

                                bits >>= 1;
                            }
                        }
                        else if (KL_METRIC == KernighanLinMetric.Cuts)
                        {
                            /* Count cut edges */
                            hop_mtx[i][j] = 1;
                        }

                        hop_mtx[j][i] = hop_mtx[i][j];
                    }
                }
            }

            graph2 = graph;
            if (global_method == PartitioningStrategy.Inertial && localPartitioningStrategy != LocalPartitioningStrategy.KernighanLin && !VERTEX_COVER && !using_vwgts)
            {
                graph2 = null;
            }

            if (!(global_method == PartitioningStrategy.ReadFromFile && localPartitioningStrategy == LocalPartitioningStrategy.None))
            {
                balance(graph2, nvtxs, nedges, using_vwgts, useEdgeWeights, vwsqrt, igeom, coords, assignment,
                    goal, architecture, ndims_tot, mesh_dims, global_method, localPartitioningStrategy, rqi_flag, vmax,
                    ndims, eigtol, hop_mtx);
            }

            partition_time += seconds() - time - kernel_time;

            nsets_tot = 0;
            if (architecture == 0)
            {
                nsets_tot = 1 << ndims_tot;
            }
            else if (architecture > 0)
            {
                nsets_tot = mesh_dims[0] * mesh_dims[1] * mesh_dims[2];
            }

            if (graph != null)
            {
                var j = true;
                for (i = 1; i <= REFINE_PARTITION && j; i++)
                {
                    /* Reduce cuts w/ KL? */
                    if (DEBUG_REFINE_PART)
                    {
                        Trace.WriteLine($"\n\nBefore pass {i:d} to refine partition:");
                        countup(graph, nvtxs, assignment, ndims, architecture, ndims_tot, mesh_dims, useEdgeWeights);
                    }

                    j = refine_part(graph, nvtxs, useEdgeWeights, assignment, architecture, ndims_tot, mesh_dims,
                        goal);
                }
            }

            if (graph != null)
            {
                if (INTERNAL_VERTICES)
                {
                    if (DEBUG_INTERNAL)
                    {
                        Trace.WriteLine("\n\nBefore increasing internal vertices:");
                        countup(graph, nvtxs, assignment, ndims, architecture, ndims_tot, mesh_dims, useEdgeWeights);
                    }

                    force_internal(graph, nvtxs, useEdgeWeights, assignment, goal, nsets_tot, INTERNAL_VERTICES ? 1 : 0);
                }
            }

            if (graph != null)
            {
                if (CONNECTED_DOMAINS)
                {
                    /* Force subdomains to be connected */
                    int res = 0;
                    connect_enforce(graph, nvtxs, useEdgeWeights, assignment, goal, nsets_tot, &i, &res);
                    Trace.WriteLine($"\nConnectivity enforcement moved {i:d} vertices total");
                    Trace.WriteLine($"                         largest moved subset = {res:d}\n");
                }
            }

            if (graph != null)
            {
                if (REFINE_MAP)
                {
                    /* Improve the mapping to processors? */
                    if (DEBUG_REFINE_MAP)
                    {
                        Trace.WriteLine("\n\nBefore refining mapping to processors:");
                        countup(graph, nvtxs, assignment, ndims, architecture, ndims_tot, mesh_dims, useEdgeWeights);
                    }

                    refine_map(graph, nvtxs, useEdgeWeights, assignment, architecture, ndims_tot, mesh_dims);
                }
            }

            if (PRINT_GRAPH_PARTITION_METRICS)
            {
                /* Compute graph metrics of partition. */
                time1 = seconds();
                if (graph != null)
                {
                    if (PRINT_HEADERS)
                    {
                        Trace.WriteLine("\n\n                     Partitioning Results");
                    }

                    countup(graph, nvtxs, assignment, ndims, architecture, ndims_tot, mesh_dims, useEdgeWeights);
                }

                count_time += seconds() - time1;
            }

            /* Invoke communication simulator? */
            /*
            if (graph != null) {
                if (SIMULATOR > 0) {
                    simulate(graph, nvtxs, ndims, architecture, ndims_tot,
                         mesh_dims, assignment, useEdgeWeights, outfile);
                }
            }
            */

            End_Label:

            if (vwsqrt != null)
            {
                Marshal.FreeHGlobal((IntPtr) vwsqrt);
            }

            if (SQRTS != null)
            {
                Marshal.FreeHGlobal((IntPtr) SQRTS);
            }

            /* Turn perturbation back on for next invocation. */
            if (ndims == 1)
            {
                PERTURB = old_perturb;
            }

            NSQRTS = old_nsqrts;

            total_time += seconds() - start_time;

            return (false);
        }

        /// <summary>
        /// Print metrics of partition quality.
        /// </summary>
        /// <param name="graph">graph data structure</param>
        /// <param name="nvtxs">number of vtxs in graph</param>
        /// <param name="assignment">set number of each vtx (length nvtxs+1)</param>
        /// <param name="ndims">number of cuts at each level</param>
        /// <param name="architecture">what's the target parallel machine?</param>
        /// <param name="ndims_tot">total number of hypercube dimensions</param>
        /// <param name="mesh_dims">extent of mesh in each dimension</param>
        /// <param name="useEdgeWeights">are edge weights being used?</param>
        private static void countup(vtx_data** graph, int nvtxs, int* assignment, int ndims, int architecture, int ndims_tot, int[] mesh_dims /*[3]*/, bool useEdgeWeights)
        {
            if (VERTEX_SEPARATOR || VERTEX_COVER)
            {
                countup_vtx_sep(graph, nvtxs, assignment);
            }
            else
            {
                if (architecture == 0)
                {
                    countup_cube(graph, nvtxs, assignment, ndims, ndims_tot, useEdgeWeights);
                }

                else if (architecture > 0)
                {
                    countup_mesh(graph, nvtxs, assignment, mesh_dims, useEdgeWeights);
                }
            }
        }

        /// <summary>
        /// Print metrics of partition quality.
        /// </summary>
        /// <param name="graph">graph data structure</param>
        /// <param name="nvtxs">number of vtxs in graph</param>
        /// <param name="assignment">set number of each vtx (length nvtxs+1)</param>
        /// <param name="ndims">number of cuts at each level</param>
        /// <param name="ndims_tot">total number of divisions of graph</param>
        /// <param name="useEdgeWeights">are edge weights being used?</param>
        private static void countup_cube(vtx_data** graph, int nvtxs, int* assignment, int ndims, int ndims_tot, bool useEdgeWeights)
        {
            double* hopsize; /* number of hops for each set */
            double* cutsize; /* number of cuts for each set */
            int* setsize; /* number of vtxs in each set */
            int* setseen; /* flags for sets adjacent to a particular set */
            int* inorder; /* list of vtxs in each set */
            int* startptr; /* indices into inorder array */
            double ncuts = 0.0d; /* total number of edges connecting sets */
            double nhops = 0.0d; /* total cuts weighted by hypercube hops */
            double ewgt; /* edge weight */
            int nsets; /* number of sets after a level */
            int vtx; /* vertex in graph */
            int set, set2, set3; /* sets neighboring vtxs are assigned to */
            int onbdy; /* counts number of neighboring set for a vtx */
            int internalNodes; /* number of internal nodes in a set */
            int min_internal = 0; /* smallest number of internal vertices */
            int max_internal = 0; /* largest number of internal vertices */
            int total_internal = 0; /* total number of internal vertices */
            int min_size = 0, max_size = 0; /* smallest and largest set sizes */
            int tot_size = 0; /* total of all set sizes */
            double bdyvtx_hops; /* bdyvtxs weighted by wire lengths */
            double bdyvtx_hops_tot = 0.0d; /* total bdyvtx_hops */
            double bdyvtx_hops_max = 0.0d; /* largest value of bdyvtx_hops among all sets */
            double bdyvtx_hops_min = 0.0d; /* smallest value of bdyvtx_hops among all sets */
            int neighbor_sets = 0; /* number of neighboring sets for a set */
            double total_bdyvtxs = 0.0d; /* sum of all onbdy values in whole graph  */
            int total_neighbors = 0; /* number of neighboring sets in graph */
            double maxcuts = 0.0d; /* largest cuts among all processors */
            double mincuts = 0.0d; /* smallest cuts among all processors */
            double maxhops = 0.0d; /* largest hops among all processors */
            double minhops = 0.0d; /* smallest hops among all processors */
            double maxbdy = 0.0d; /* largest bdy_vtxs among all processors */
            double minbdy = 0.0d; /* smallest bdy_vtxs among all processors */
            int maxneighbors = 0; /* largest neighbor_sets among all processors */
            int minneighbors = 0; /* smallest neighbor_sets among all processors */
            int neighbor; /* neighbor of a vertex */
            int mask; /* mask for active bits */
            int bits; /* bit pattern for counting hops */
            int start_dims; /* starting dimension for output loop */
            int level; /* recursion level of partition */
            int i, j, k, l, ll; /* loop counters */

            ewgt = 1;

            nsets = (1 << ndims_tot);
            cutsize = (double*) Marshal.AllocHGlobal(nsets * sizeof(double));
            hopsize = (double*) Marshal.AllocHGlobal(nsets * sizeof(double));
            setsize = (int*) Marshal.AllocHGlobal(nsets * sizeof(int));

            setseen = (int*) Marshal.AllocHGlobal(nsets * sizeof(int));
            startptr = (int*) Marshal.AllocHGlobal((nsets + 1) * sizeof(int));
            inorder = (int*) Marshal.AllocHGlobal(nvtxs * sizeof(int));
            for (j = 0; j < nsets; j++)
            {
                setsize[j] = 0;
            }

            for (i = 1; i <= nvtxs; i++)
            {
                ++setsize[assignment[i]];
            }

            /* Modify setsize to become index into vertex list. */
            for (j = 1; j < nsets; j++)
            {
                setsize[j] += setsize[j - 1];
            }

            for (j = nsets - 1; j > 0; j--)
            {
                startptr[j] = setsize[j] = setsize[j - 1];
            }

            startptr[0] = setsize[0] = 0;
            startptr[nsets] = nvtxs;
            for (i = 1; i <= nvtxs; i++)
            {
                set = assignment[i];
                inorder[setsize[set]] = i;
                setsize[set]++;
            }

            if (PRINT_GRAPH_PARTITION_METRICS_ALL_RECURSION_LEVELS)
            {
                /* Print data from all levels of recursion. */
                start_dims = ndims;
                level = 0;
            }
            else
            {
                /* Only print data from final level. */
                start_dims = ndims_tot;
                level = (ndims_tot + ndims - 1) / ndims - 1;
            }

            k = start_dims;
            while (k <= ndims_tot)
            {
                level++;
                nsets = (1 << k);
                for (j = 0; j < nsets; j++)
                {
                    cutsize[j] = 0;
                    hopsize[j] = 0;
                    setsize[j] = 0;
                }

                mask = 0;
                for (j = 0; j < k; j++)
                {
                    mask = (mask << 1) + 1;
                }

                for (i = 1; i <= nvtxs; i++)
                {
                    set = assignment[i] & mask;
                    setsize[set] += graph[i]->vwgt;
                    for (j = 1; j < graph[i]->nedges; j++)
                    {
                        neighbor = graph[i]->edges[j];
                        set2 = assignment[neighbor] & mask;
                        if (set != set2)
                        {
                            if (useEdgeWeights)
                            {
                                ewgt = graph[i]->ewgts[j];
                            }

                            cutsize[set] += ewgt;
                            bits = set ^ set2;
                            for (l = bits; l != 0; l >>= 1)
                            {
                                if ((l & 1) != 0)
                                {
                                    hopsize[set] += ewgt;
                                }
                            }
                        }
                    }
                }

                tot_size = 0;
                max_size = 0;
                for (set = 0; set < nsets; set++)
                {
                    tot_size += setsize[set];
                    if (setsize[set] > max_size)
                    {
                        max_size = setsize[set];
                    }
                }

                min_size = max_size;
                for (set = 0; set < nsets; set++)
                {
                    if (setsize[set] < min_size)
                    {
                        min_size = setsize[set];
                    }
                }

                ncuts = nhops = 0;
                total_bdyvtxs = total_neighbors = 0;
                bdyvtx_hops_tot = bdyvtx_hops_max = bdyvtx_hops_min = 0;
                maxcuts = mincuts = 0;
                maxhops = minhops = 0;
                total_internal = 0;
                min_internal = max_size;
                max_internal = 0;
                maxbdy = minbdy = 0;
                maxneighbors = minneighbors = 0;

                Trace.WriteLine($"\nAfter level {level:d}  (nsets = {nsets:d}):");
                if (PRINT_GRAPH_PARTITION_METRICS_DETAILED)
                {
                    Trace.WriteLine("    set    size      cuts       hops   bndy_vtxs    adj_sets");
                }

                int bdyvtxs = 0; /* sum of onbdy values for a set */
                for (set = 0; set < nsets; set++)
                {
                    internalNodes = setsize[set];
                    for (i = 0; i < nsets; i++)
                    {
                        setseen[i] = 0;
                    }

                    /* Compute number of set neighbors, and number of vtxs on boundary. */
                    /* Loop through multiple assignments defining current set. */
                    bdyvtxs = 0;
                    bdyvtx_hops = 0;
                    for (l = 0; l < (1 << (ndims_tot - k)); l++)
                    {
                        set2 = (l << k) + set;
                        for (i = startptr[set2]; i < startptr[set2 + 1]; i++)
                        {
                            onbdy = 0;
                            vtx = inorder[i];
                            for (j = 1; j < graph[vtx]->nedges; j++)
                            {
                                neighbor = graph[vtx]->edges[j];
                                set3 = assignment[neighbor] & mask;
                                if (set3 != set)
                                {
                                    /* Is vtx on boundary? */
                                    /* Has this neighboring set been seen already? */
                                    if (setseen[set3] >= 0)
                                    {
                                        bits = set ^ set3;
                                        for (ll = bits; ll != 0; ll >>= 1)
                                        {
                                            if ((ll & 1) != 0)
                                            {
                                                ++bdyvtx_hops;
                                            }
                                        }

                                        ++onbdy;
                                        setseen[set3] = -setseen[set3] - 1;
                                    }
                                }
                            }

                            /* Now reset all the setseen values to be positive. */
                            if (onbdy != 0)
                            {
                                for (j = 1; j < graph[vtx]->nedges; j++)
                                {
                                    neighbor = graph[vtx]->edges[j];
                                    set3 = assignment[neighbor] & mask;
                                    if (setseen[set3] < 0)
                                    {
                                        setseen[set3] = -setseen[set3];
                                    }
                                }

                                internalNodes -= graph[vtx]->vwgt;
                            }

                            bdyvtxs += onbdy;
                        }
                    }

                    total_internal += internalNodes;
                    bdyvtx_hops_tot += bdyvtx_hops;
                    if (bdyvtx_hops > bdyvtx_hops_max)
                    {
                        bdyvtx_hops_max = bdyvtx_hops;
                    }

                    if (set == 0 || bdyvtx_hops < bdyvtx_hops_min)
                    {
                        bdyvtx_hops_min = bdyvtx_hops;
                    }

                    if (internalNodes > max_internal)
                    {
                        max_internal = internalNodes;
                    }

                    if (set == 0 || internalNodes < min_internal)
                    {
                        min_internal = internalNodes;
                    }

                    /* Now count up the number of neighboring sets. */
                    neighbor_sets = 0;
                    for (i = 0; i < nsets; i++)
                    {
                        if (setseen[i] != 0)
                        {
                            ++neighbor_sets;
                        }
                    }

                    if (PRINT_GRAPH_PARTITION_METRICS_DETAILED)
                    {
                        Trace.WriteLine($" {set:d}    {setsize[set]:d}    {cutsize[set]:g}     {hopsize[set]:g}   {bdyvtxs:d}      {neighbor_sets:d}");
                    }

                    if (cutsize[set] > maxcuts)
                    {
                        maxcuts = cutsize[set];
                    }

                    if (set == 0 || cutsize[set] < mincuts)
                    {
                        mincuts = cutsize[set];
                    }

                    if (hopsize[set] > maxhops)
                    {
                        maxhops = hopsize[set];
                    }

                    if (set == 0 || hopsize[set] < minhops)
                    {
                        minhops = hopsize[set];
                    }

                    if (bdyvtxs > maxbdy)
                    {
                        maxbdy = bdyvtxs;
                    }

                    if (set == 0 || bdyvtxs < minbdy)
                    {
                        minbdy = bdyvtxs;
                    }

                    if (neighbor_sets > maxneighbors)
                    {
                        maxneighbors = neighbor_sets;
                    }

                    if (set == 0 || neighbor_sets < minneighbors)
                    {
                        minneighbors = neighbor_sets;
                    }

                    ncuts += cutsize[set];
                    nhops += hopsize[set];
                    total_bdyvtxs += bdyvtxs;
                    total_neighbors += neighbor_sets;
                }

                ncuts /= 2;
                nhops /= 2;

                Trace.WriteLine("\n");
                Trace.WriteLine("                            Total      Max/Set      Min/Set\n");
                Trace.WriteLine("                            -----      -------      -------\n");
                Trace.WriteLine($"Set Size:             {tot_size:d}  {max_size:d}  {min_size:d}");
                Trace.WriteLine($"Edge Cuts:            {ncuts:g}  {maxcuts:g}  {mincuts:g}");
                Trace.WriteLine($"Hypercube Hops:       {nhops:g}  {maxhops:g}  {minhops:g}");
                Trace.WriteLine($"Boundary Vertices:    {total_bdyvtxs:g}  {maxbdy:g}  {minbdy:g}");
                Trace.WriteLine($"Boundary Vertex Hops: {bdyvtx_hops_tot:g}  {bdyvtx_hops_max:g}  {bdyvtx_hops_min:g}");
                Trace.WriteLine($"Adjacent Sets:        {total_neighbors:d}  {maxneighbors:d}  {minneighbors:d}");
                Trace.WriteLine($"Internal Vertices:    {total_internal:d}  {max_internal:d}  {min_internal:d}\n");

                if (k == ndims_tot)
                {
                    k++;
                }
                else
                {
                    k += ndims;
                    if (k > ndims_tot)
                    {
                        k = ndims_tot;
                    }
                }
            }

            Marshal.FreeHGlobal((IntPtr) cutsize);
            Marshal.FreeHGlobal((IntPtr) hopsize);
            Marshal.FreeHGlobal((IntPtr) setsize);
            Marshal.FreeHGlobal((IntPtr) setseen);
            Marshal.FreeHGlobal((IntPtr) startptr);
            Marshal.FreeHGlobal((IntPtr) inorder);
        }

        /// <summary>
        /// Print metrics of partition quality.
        /// </summary>
        /// <param name="graph">graph data structure</param>
        /// <param name="nvtxs">number of vtxs in graph</param>
        /// <param name="assignment">set number of each vtx (length nvtxs+1)</param>
        /// <param name="mesh_dims">extent of mesh in each dimension</param>
        /// <param name="useEdgeWeights">are edge weights being used?</param>
        static void countup_mesh(vtx_data** graph, int nvtxs, int* assignment, int[] mesh_dims /*[3]*/, bool useEdgeWeights)
        {
            double* hopsize; /* number of hops for each set */
            double* cutsize; /* number of cuts for each set */
            int* setsize; /* number or weight of vtxs in each set */
            int* setseen; /* flags for sets adjacent to a particular set */
            int* inorder; /* list of vtxs in each set */
            int* startptr; /* indices into inorder array */
            double ncuts; /* total number of edges connecting sets */
            double nhops; /* total cuts weighted by mesh hops */
            double ewgt; /* edge weight */
            int nsets; /* number of sets after a level */
            int vtx; /* vertex in graph */
            int set, set2, set3; /* sets neighboring vtxs are assigned to */
            int onbdy; /* counts number of neighboring set for a vtx */
            int min_size, max_size; /* min and max set sizes */
            int tot_size; /* total of all set sizes */
            double bdyvtxs; /* sum of onbdy values for a set */
            double bdyvtx_hops; /* weights boundary vertices by wires required */
            double bdyvtx_hops_tot; /* weights boundary vertices by wires required */
            double bdyvtx_hops_max; /* weights boundary vertices by wires required */
            double bdyvtx_hops_min; /* weights boundary vertices by wires required */
            int neighbor_sets; /* number of neighboring sets for a set */
            int internalVertexCount; /* number of internal vertices in a set */
            int min_internal; /* smallest number of internal vertices */
            int max_internal; /* largest number of internal vertices */
            int total_internal; /* total number of internal nodes */
            double maxcuts; /* largest cuts among all processors */
            double mincuts; /* smallest cuts among all processors */
            double maxhops; /* largest hops among all processors */
            double minhops; /* smallest hops among all processors */
            double maxbdy; /* largest bdy_vtxs among all processors */
            double minbdy; /* smallest bdy_vtxs among all processors */
            int maxneighbors; /* largest neighbor_sets among all processors */
            int minneighbors; /* smallest neighbor_sets among all processors */
            double total_bdyvtxs; /* sum of all onbdy values in whole graph  */
            int total_neighbors; /* number of neighboring sets in graph */
            int neighbor; /* neighbor of a vertex */
            int x1, y1, z1; /* mesh location of vertex */
            int x2, y2, z2; /* mesh location of neighboring vertex */
            int i, j; /* loop counters */

            ewgt = 1;

            nsets = mesh_dims[0] * mesh_dims[1] * mesh_dims[2];
            cutsize = (double*) Marshal.AllocHGlobal(nsets * sizeof(double));
            hopsize = (double*) Marshal.AllocHGlobal(nsets * sizeof(double));
            setsize = (int*) Marshal.AllocHGlobal(nsets * sizeof(int));

            setseen = (int*) Marshal.AllocHGlobal(nsets * sizeof(int));
            startptr = (int*) Marshal.AllocHGlobal((nsets + 1) * sizeof(int));
            inorder = (int*) Marshal.AllocHGlobal(nvtxs * sizeof(int));
            for (j = 0; j < nsets; j++)
            {
                setsize[j] = 0;
            }

            for (i = 1; i <= nvtxs; i++)
            {
                ++setsize[assignment[i]];
            }

            /* Modify setsize to become index into vertex list. */
            for (j = 1; j < nsets; j++)
            {
                setsize[j] += setsize[j - 1];
            }

            for (j = nsets - 1; j > 0; j--)
            {
                startptr[j] = setsize[j] = setsize[j - 1];
            }

            startptr[0] = setsize[0] = 0;
            startptr[nsets] = nvtxs;
            for (i = 1; i <= nvtxs; i++)
            {
                set = assignment[i];
                inorder[setsize[set]] = i;
                setsize[set]++;
            }

            for (j = 0; j < nsets; j++)
            {
                cutsize[j] = 0;
                hopsize[j] = 0;
                setsize[j] = 0;
            }

            for (i = 1; i <= nvtxs; i++)
            {
                set = assignment[i];
                setsize[set] += graph[i]->vwgt;
                x1 = set % mesh_dims[0];
                y1 = (set / mesh_dims[0]) % mesh_dims[1];
                z1 = set / (mesh_dims[0] * mesh_dims[1]);
                for (j = 1; j < graph[i]->nedges; j++)
                {
                    neighbor = graph[i]->edges[j];
                    set2 = assignment[neighbor];
                    x2 = set2 % mesh_dims[0];
                    y2 = (set2 / mesh_dims[0]) % mesh_dims[1];
                    z2 = set2 / (mesh_dims[0] * mesh_dims[1]);
                    if (set != set2)
                    {
                        if (useEdgeWeights)
                        {
                            ewgt = graph[i]->ewgts[j];
                        }

                        cutsize[set] += ewgt;
                        hopsize[set] += ewgt * (Math.Abs(x1 - x2) + Math.Abs(y1 - y2) + Math.Abs(z1 - z2));
                    }
                }
            }

            max_size = 0;
            tot_size = 0;
            for (set = 0; set < nsets; set++)
            {
                tot_size += setsize[set];
                if (setsize[set] > max_size)
                {
                    max_size = setsize[set];
                }
            }

            min_size = max_size;
            for (set = 0; set < nsets; set++)
            {
                if (setsize[set] < min_size)
                {
                    min_size = setsize[set];
                }
            }

            ncuts = nhops = 0;
            total_bdyvtxs = total_neighbors = 0;
            bdyvtx_hops_tot = bdyvtx_hops_max = bdyvtx_hops_min = 0;
            maxcuts = mincuts = 0;
            maxhops = minhops = 0;
            maxbdy = minbdy = 0;
            maxneighbors = minneighbors = 0;
            total_internal = 0;
            min_internal = max_size;
            max_internal = 0;
            Trace.WriteLine($"\nAfter full partitioning  (nsets = {nsets:d})");

            if (PRINT_GRAPH_PARTITION_METRICS_DETAILED)
            {
                Trace.WriteLine("    set    size      cuts       hops   bndy_vtxs    adj_sets");
            }

            for (set = 0; set < nsets; set++)
            {
                /* Compute number of set neighbors, and number of vtxs on boundary. */
                internalVertexCount = setsize[set];
                for (i = 0; i < nsets; i++)
                {
                    setseen[i] = 0;
                }

                x1 = set % mesh_dims[0];
                y1 = (set / mesh_dims[0]) % mesh_dims[1];
                z1 = set / (mesh_dims[0] * mesh_dims[1]);

                set2 = set;
                bdyvtxs = 0;
                bdyvtx_hops = 0;
                for (i = startptr[set2]; i < startptr[set2 + 1]; i++)
                {
                    onbdy = 0;
                    vtx = inorder[i];
                    for (j = 1; j < graph[vtx]->nedges; j++)
                    {
                        neighbor = graph[vtx]->edges[j];
                        set3 = assignment[neighbor];
                        if (set3 != set)
                        {
                            /* Is vtx on boundary? */
                            /* Has this neighboring set been seen already? */
                            if (setseen[set3] >= 0)
                            {
                                x2 = set3 % mesh_dims[0];
                                y2 = (set3 / mesh_dims[0]) % mesh_dims[1];
                                z2 = set3 / (mesh_dims[0] * mesh_dims[1]);
                                bdyvtx_hops += Math.Abs(x1 - x2) + Math.Abs(y1 - y2) + Math.Abs(z1 - z2);
                                ++onbdy;
                                setseen[set3] = -setseen[set3] - 1;
                            }
                        }
                    }

                    /* Now reset all the setseen values to be positive. */
                    if (onbdy != 0)
                    {
                        for (j = 1; j < graph[vtx]->nedges; j++)
                        {
                            neighbor = graph[vtx]->edges[j];
                            set3 = assignment[neighbor];
                            if (setseen[set3] < 0)
                            {
                                setseen[set3] = -setseen[set3];
                            }
                        }

                        internalVertexCount -= graph[vtx]->vwgt;
                    }

                    bdyvtxs += onbdy;
                }

                total_internal += internalVertexCount;
                bdyvtx_hops_tot += bdyvtx_hops;
                if (bdyvtx_hops > bdyvtx_hops_max)
                {
                    bdyvtx_hops_max = bdyvtx_hops;
                }

                if (set == 0 || bdyvtx_hops < bdyvtx_hops_min)
                {
                    bdyvtx_hops_min = bdyvtx_hops;
                }

                if (internalVertexCount > max_internal)
                {
                    max_internal = internalVertexCount;
                }

                if (set == 0 || internalVertexCount < min_internal)
                {
                    min_internal = internalVertexCount;
                }

                /* Now count up the number of neighboring sets. */
                neighbor_sets = 0;
                for (i = 0; i < nsets; i++)
                {
                    if (setseen[i] != 0)
                    {
                        ++neighbor_sets;
                    }
                }

                if (PRINT_GRAPH_PARTITION_METRICS_DETAILED)
                {
                    Trace.WriteLine($" {set:d}    {setsize[set]:d}    {cutsize[set]:g}     {hopsize[set]:g}   {bdyvtxs:g}      {neighbor_sets:d}");
                }

                if (cutsize[set] > maxcuts)
                {
                    maxcuts = cutsize[set];
                }

                if (set == 0 || cutsize[set] < mincuts)
                {
                    mincuts = cutsize[set];
                }

                if (hopsize[set] > maxhops)
                {
                    maxhops = hopsize[set];
                }

                if (set == 0 || hopsize[set] < minhops)
                {
                    minhops = hopsize[set];
                }

                if (bdyvtxs > maxbdy)
                {
                    maxbdy = bdyvtxs;
                }

                if (set == 0 || bdyvtxs < minbdy)
                {
                    minbdy = bdyvtxs;
                }

                if (neighbor_sets > maxneighbors)
                {
                    maxneighbors = neighbor_sets;
                }

                if (set == 0 || neighbor_sets < minneighbors)
                {
                    minneighbors = neighbor_sets;
                }

                ncuts += cutsize[set];
                nhops += hopsize[set];
                total_bdyvtxs += bdyvtxs;
                total_neighbors += neighbor_sets;
            }

            ncuts /= 2;
            nhops /= 2;
            Trace.WriteLine("");
            Trace.WriteLine("                            Total      Max/Set      Min/Set");
            Trace.WriteLine("                            -----      -------      -------");
            Trace.WriteLine($"Set Size:             {tot_size:d}  {max_size:d}  {min_size:d}");
            Trace.WriteLine($"Edge Cuts:            {ncuts:g}  {maxcuts:g}  {mincuts:g}");
            Trace.WriteLine($"Mesh Hops:            {nhops:g}  {maxhops:g}  {minhops:g}");
            Trace.WriteLine($"Boundary Vertices:    {total_bdyvtxs:g}  {maxbdy:g}  {minbdy:g}");
            Trace.WriteLine($"Boundary Vertex Hops: {bdyvtx_hops_tot:g}  {bdyvtx_hops_max:g}  {bdyvtx_hops_min:g}");
            Trace.WriteLine($"Adjacent Sets:        {total_neighbors:d}  {maxneighbors:d}  {minneighbors:d}");
            Trace.WriteLine($"Internal Vertices:    {total_internal:d}  {max_internal:d}  {min_internal:d}\n");

            Marshal.FreeHGlobal((IntPtr) cutsize);
            Marshal.FreeHGlobal((IntPtr) hopsize);
            Marshal.FreeHGlobal((IntPtr) setsize);
            Marshal.FreeHGlobal((IntPtr) setseen);
            Marshal.FreeHGlobal((IntPtr) startptr);
            Marshal.FreeHGlobal((IntPtr) inorder);
        }

        /// <summary>
        ///
        /// </summary>
        /// <param name="graph">data structure for graph</param>
        /// <param name="nvtxs">number of vertices in full graph</param>
        /// <param name="nedges">number of edges in graph</param>
        /// <param name="using_vwgts">are vertex weights being used?</param>
        /// <param name="useEdgeWeights">are edge weights being used?</param>
        /// <param name="vwsqrt">sqrt of vertex weights (length nvtxs+1)</param>
        /// <param name="igeom">geometric dimension for inertial method</param>
        /// <param name="coords">coordinates for inertial method</param>
        /// <param name="assignment">set number of each vtx (length n)</param>
        /// <param name="goal">desired set sizes</param>
        /// <param name="architecture">0=> hypercube, d=> d-dimensional mesh</param>
        /// <param name="ndims_tot">number of cuts to make in total</param>
        /// <param name="mesh_dims">shape of mesh</param>
        /// <param name="global_method">global partitioning algorithm</param>
        /// <param name="local_method">local partitioning algorithm</param>
        /// <param name="rqi_flag">should I use multilevel eigensolver?</param>
        /// <param name="vmax">if so, how many vertices to coarsen down to?</param>
        /// <param name="ndims">number of eigenvectors (2^d sets)</param>
        /// <param name="eigtol">tolerance on eigenvectors</param>
        /// <param name="hops">between-set hop cost for KL</param>
        static void balance(vtx_data** graph,
            int nvtxs,
            int nedges,
            bool using_vwgts,
            bool useEdgeWeights,
            double* vwsqrt,
            int igeom,
            float** coords,
            int* assignment,
            double[] goal,
            int architecture,
            int ndims_tot,
            int[] mesh_dims,
            PartitioningStrategy global_method,
            LocalPartitioningStrategy local_method,
            bool rqi_flag,
            int vmax,
            int ndims,
            double eigtol,
            int[][] hops
        )
        {
            vtx_data** subgraph = null; /* data structure for subgraph */
            set_info* set_info = null; /* information about each processor subset */
            set_info** set_buckets = null; /* buckets for sorting processor sets */
            set_info* set = null; /* current processor set information */
            int[][] hops_special; //[MAXSETS][MAXSETS]; /* hop mtx for nonstandard cases */
            hops_special = new int[MAXSETS][];
            for (var i = 0; i < hops_special.Length; i++)
            {
                hops_special[i] = new int[MAXSETS];
            }

            float* term_wgts; /* net pull of terminal propagation */
            float* save_term_wgts; /* saved location of term_wgts */
            float*[] all_term_wgts = new float*[MAXSETS]; /* net pull on all sets */
            int* loc2glob; /* mapping from subgraph to graph numbering */
            int* glob2loc; /* mapping from graph to subgraph numbering */
            int* setlists; /* space for linked lists of vertices in sets */
            int* list_ptrs; /* headers of each linked list */
            int* degree; /* degrees of graph vertices from a subgraph */
            int[] subsets = new int[MAXSETS]; /* subsets being created at current step */
            int[] setsize = new int[MAXSETS]; /* sizes of sets created by division */
            double[] merged_goal = new double[MAXSETS]; /* sizes of sets at this partition level */
            double sub_vwgt_sum; /* sum of subgraph vertex weights */
            int[] cut_dirs = new int[MAXDIMS]; /* directions of processor cuts if mesh */
            bool hops_flag; /* use normal or special hops? */
            int ndims_real; /* actual dimension of partitioning */
            int nsets_real; /* actual # subsets to create */
            int maxsize; /* size of largest subgraph */
            int nsets_tot; /* total sets to divide subgraph into */
            int subnvtxs; /* number of vertices in subgraph */
            int subnedges; /* number of edgess in subgraph */
            double* subvwsqrt = null; /* vwsqrt array for subgraph */
            int* subassign = null; /* set assignments for subgraph */
            float** subcoords = null; /* coordinates for subgraph */
            int striping; /* partition with parallel cuts? */
            int pass; /* counts passes through loop */
            int max_proc_size; /* size of largest processor set */
            int old_max_proc_size; /* previous size of largest processor set */
            bool new_dir; /* new cut direction? => no terminal prop */
            int nsets; /* typical number of sets at each cut */
            bool[] done_dir = new bool[3]; /* mesh directions already cut */

            if (DEBUG_TRACE)
            {
                Trace.WriteLine("<Entering balance>");
            }

            if (global_method != PartitioningStrategy.ReadFromFile)
            {
                /* Not read from file. */
                for (var i = 1; i <= nvtxs; i++)
                {
                    assignment[i] = 0;
                }
            }

            if (nvtxs <= 1)
            {
                return;
            }

            /* Compute some simple parameters. */
            save_term_wgts = null;
            maxsize = 0;
            nsets = 1 << ndims;
            subsets[2] = 2; /* Needed for vertex separators */

            nsets_tot = 0;
            if (architecture > 0)
            {
                nsets_tot = mesh_dims[0] * mesh_dims[1] * mesh_dims[2];
            }
            else if (architecture == 0)
            {
                nsets_tot = 1 << ndims_tot;
            }

            /* Construct data structures for keeping track of processor sets. */
            set_buckets = (set_info**) Marshal.AllocHGlobal((nsets_tot + 1) * sizeof(set_info*));
            set_info = (set_info*) Marshal.AllocHGlobal(nsets_tot * sizeof(set_info));
            for (var i = 0; i < nsets_tot; i++)
            {
                set_info[i].setnum = i;
                set_info[i].ndims = -1;
                set_info[i].span[0] = -1;
                set_info[i].next = null;
                set_buckets[i + 1] = null;
            }

            set_buckets[nsets_tot] = &(set_info[0]);

            if (architecture > 0)
            {
                set_info[0].low[0] = set_info[0].low[1] = set_info[0].low[2] = 0;
                set_info[0].span[0] = mesh_dims[0];
                set_info[0].span[1] = mesh_dims[1];
                set_info[0].span[2] = mesh_dims[2];
            }
            else if (architecture == 0)
            {
                set_info[0].ndims = ndims_tot;
            }

            done_dir[0] = done_dir[1] = done_dir[2] = false;
            pass = 0;
            loc2glob = null;
            degree = null;
            glob2loc = null;
            setlists = null;
            list_ptrs = null;

            old_max_proc_size = 2 * nsets_tot;
            max_proc_size = nsets_tot;

            while (max_proc_size > 1)
            {
                /* Some set still needs to be divided. */

                pass++;
                set = set_buckets[max_proc_size];
                set_buckets[max_proc_size] = set->next;

                /* Divide the processors. */
                hops_flag =
                    divide_procs(architecture, ndims, ndims_tot, set_info, set, subsets, (global_method == PartitioningStrategy.Inertial),
                        &ndims_real, &nsets_real, &striping, cut_dirs, mesh_dims, hops_special);

                /* If new cut direction, turn off terminal propagation. */
                new_dir = false;
                if (architecture == 0)
                {
                    if (old_max_proc_size != max_proc_size)
                    {
                        new_dir = true;
                    }
                }
                else if (architecture > 0)
                {
                    if (!done_dir[cut_dirs[0]])
                    {
                        new_dir = true;
                    }

                    done_dir[cut_dirs[0]] = true;
                }

                old_max_proc_size = max_proc_size;

                /* Now place the new sets into the set_info data structure. */
                /* This is indexed by number of processors in set. */
                for (var i = 0; i < nsets_real; i++)
                {
                    int j = 0;
                    if (architecture > 0)
                    {
                        j = set_info[subsets[i]].span[0] * set_info[subsets[i]].span[1] *
                            set_info[subsets[i]].span[2];
                    }
                    else if (architecture == 0)
                    {
                        j = 1 << set_info[subsets[i]].ndims;
                    }

                    set_info[subsets[i]].next = set_buckets[j];
                    set_buckets[j] = &(set_info[subsets[i]]);
                }

                /* Construct desired set sizes for this division step. */
                if (pass == 1)
                {
                    /* First partition. */
                    subgraph = graph;
                    subnvtxs = nvtxs;
                    subnedges = nedges;
                    subvwsqrt = vwsqrt;
                    subcoords = coords;
                    subassign = assignment;
                    all_term_wgts[1] = null;

                    if (!using_vwgts)
                    {
                        sub_vwgt_sum = subnvtxs;
                    }
                    else
                    {
                        sub_vwgt_sum = 0;
                        for (var i = 1; i <= subnvtxs; i++)
                        {
                            sub_vwgt_sum += subgraph[i]->vwgt;
                        }
                    }

                    merge_goals(goal, merged_goal, set_info, subsets, nsets_real, ndims_tot, architecture != 0,
                        mesh_dims, sub_vwgt_sum);
                }

                else
                {
                    /* Not the first cut. */

                    /* After first cut, allocate all space we'll need. */
                    if (pass == 2)
                    {
                        glob2loc = (int*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));
                        loc2glob = (int*) Marshal.AllocHGlobal((maxsize + 1) * sizeof(int));
                        if (!using_vwgts)
                        {
                            subvwsqrt = null;
                        }
                        else
                        {
                            subvwsqrt = (double*) Marshal.AllocHGlobal((maxsize + 1) * sizeof(double));
                        }

                        if (graph != null)
                        {
                            subgraph = (vtx_data**) Marshal.AllocHGlobal((maxsize + 1) * sizeof(vtx_data*));
                            degree = (int*) Marshal.AllocHGlobal((maxsize + 1) * sizeof(int));
                        }
                        else
                        {
                            subgraph = null;
                        }

                        subassign = (int*) Marshal.AllocHGlobal((maxsize + 1) * sizeof(int));
                        if (global_method == PartitioningStrategy.Inertial ||
                            (MATCH_TYPE == MatchingRoutine.maxmatch5_geometric && (global_method == PartitioningStrategy.Multilevel_KL || (global_method == PartitioningStrategy.Spectral && rqi_flag))))
                        {
                            subcoords = (float**) Marshal.AllocHGlobal(3 * sizeof(float*));
                            subcoords[1] = subcoords[2] = null;
                            subcoords[0] = (float*) Marshal.AllocHGlobal((maxsize + 1) * sizeof(float));
                            if (igeom > 1)
                            {
                                subcoords[1] = (float*) Marshal.AllocHGlobal((maxsize + 1) * sizeof(float));
                            }

                            if (igeom > 2)
                            {
                                subcoords[2] = (float*) Marshal.AllocHGlobal((maxsize + 1) * sizeof(float));
                            }
                        }
                        else
                        {
                            subcoords = null;
                        }

                        if (TERM_PROP && graph != null)
                        {
                            term_wgts = (float*) Marshal.AllocHGlobal((nsets - 1) * (maxsize + 1) * sizeof(float));
                            save_term_wgts = term_wgts;
                            for (var i = 1; i < nsets; i++)
                            {
                                all_term_wgts[i] = term_wgts;
                                term_wgts += maxsize + 1;
                            }
                        }
                    }

                    /* Construct mappings between local and global vertex numbering */
                    subnvtxs = make_maps(setlists, list_ptrs, set->setnum, glob2loc, loc2glob);

                    if (TERM_PROP && !new_dir && graph != null)
                    {
                        all_term_wgts[1] = save_term_wgts;
                        make_term_props(graph, subnvtxs, loc2glob, assignment, architecture, ndims_tot, ndims_real,
                            set_info, set->setnum, nsets_real, nsets_tot, subsets, all_term_wgts,
                            useEdgeWeights);
                    }
                    else
                    {
                        all_term_wgts[1] = null;
                    }

                    /* Form the subgraph in our graph format. */
                    if (graph != null)
                    {
                        make_subgraph(graph, subgraph, subnvtxs, &subnedges, assignment, set->setnum, glob2loc,
                            loc2glob, degree, useEdgeWeights);
                    }
                    else
                    {
                        subnedges = 0; /* Otherwise some output is garbage */
                    }

                    if (!using_vwgts)
                    {
                        sub_vwgt_sum = subnvtxs;
                    }
                    else
                    {
                        sub_vwgt_sum = 0;
                        for (var i = 1; i <= subnvtxs; i++)
                        {
                            sub_vwgt_sum += subgraph[i]->vwgt;
                        }
                    }

                    merge_goals(goal, merged_goal, set_info, subsets, nsets_real, ndims_tot, architecture != 0,
                        mesh_dims, sub_vwgt_sum);

                    /* Condense the relevant vertex weight array. */
                    if (using_vwgts && vwsqrt != null)
                    {
                        make_subvector(vwsqrt, subvwsqrt, subnvtxs, loc2glob);
                    }

                    if (global_method == PartitioningStrategy.Inertial ||
                        (MATCH_TYPE == MatchingRoutine.maxmatch5_geometric && (global_method == PartitioningStrategy.Multilevel_KL || (global_method == PartitioningStrategy.Spectral && rqi_flag))))
                    {
                        make_subgeom(igeom, coords, subcoords, subnvtxs, loc2glob);
                    }
                }

                if (DEBUG_TRACE)
                {
                    Trace.WriteLine(string.Format("About to call divide with nvtxs = %d, nedges = %d, ", subnvtxs, subnedges));
                    if (architecture != 0)
                    {
                        Trace.WriteLine($"ndims = {set->ndims:d}");
                    }
                    else if (architecture == 1)
                    {
                        Trace.WriteLine($"mesh = {set->span[0]:d}");
                    }
                    else if (architecture == 2)
                    {
                        Trace.WriteLine($"mesh = {set->span[0]:d}x{set->span[1]:d}");
                    }
                    else if (architecture == 3)
                    {
                        Trace.WriteLine($"mesh = {set->span[0]:d}x{set->span[1]:d}x{set->span[2]:d}");
                    }
                }

                /* Perform a single division step. */
                divide(subgraph, subnvtxs, subnedges, using_vwgts, useEdgeWeights, subvwsqrt, igeom, subcoords,
                    subassign, merged_goal, architecture != 0, all_term_wgts, global_method, local_method,
                    rqi_flag, vmax, ndims_real, eigtol, (hops_flag ? hops_special : hops), nsets_real,
                    striping);

                /* Undo the subgraph construction. */
                if (pass != 1 && graph != null)
                {
                    remake_graph(subgraph, subnvtxs, loc2glob, degree, useEdgeWeights);
                }

                /* Prepare for next division */
                while (max_proc_size > 1 && set_buckets[max_proc_size] == null)
                {
                    --max_proc_size;
                }

                /* Merge the subgraph partitioning with the graph partitioning. */
                if (pass == 1)
                {
                    subgraph = null;
                    subvwsqrt = null;
                    subcoords = null;
                    subassign = null;

                    /* Find size of largest subgraph for recursing. */
                    if (max_proc_size > 1)
                    {
                        for (var i = 0; i < nsets; i++)
                        {
                            setsize[i] = 0;
                        }

                        for (var i = 1; i <= nvtxs; i++)
                        {
                            ++setsize[assignment[i]];
                        }

                        maxsize = 0;
                        for (var i = 0; i < nsets; i++)
                        {
                            if (setsize[i] > maxsize)
                            {
                                maxsize = setsize[i];
                            }
                        }
                    }

                    for (var i = 1; i <= nvtxs; i++)
                    {
                        assignment[i] = subsets[assignment[i]];
                    }

                    /* Construct list of vertices in sets for make_maps. */
                    if (max_proc_size > 1)
                    {
                        setlists = (int*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));
                        list_ptrs = (int*) Marshal.AllocHGlobal(nsets_tot * sizeof(int));

                        make_setlists(setlists, list_ptrs, nsets_real, subsets, assignment, loc2glob, nvtxs, true);
                    }
                }
                else
                {
                    /* Construct list of vertices in sets for make_maps. */
                    if (max_proc_size > 1)
                    {
                        make_setlists(setlists, list_ptrs, nsets_real, subsets, subassign, loc2glob, subnvtxs,
                            false);
                    }

                    merge_assignments(assignment, subassign, subsets, subnvtxs, loc2glob);
                }
            }

            /* Free everything allocated for subgraphs. */
            Marshal.FreeHGlobal((IntPtr) list_ptrs);
            Marshal.FreeHGlobal((IntPtr) setlists);
            Marshal.FreeHGlobal((IntPtr) save_term_wgts);
            Marshal.FreeHGlobal((IntPtr) subassign);
            Marshal.FreeHGlobal((IntPtr) loc2glob);
            Marshal.FreeHGlobal((IntPtr) glob2loc);
            if (graph != null)
            {
                Marshal.FreeHGlobal((IntPtr) degree);
            }

            if (subgraph != null)
            {
                Marshal.FreeHGlobal((IntPtr) subgraph);
            }

            if (subvwsqrt != null)
            {
                Marshal.FreeHGlobal((IntPtr) subvwsqrt);
            }

            if (subcoords != null)
            {
                if (subcoords[0] != null)
                {
                    Marshal.FreeHGlobal((IntPtr) subcoords[0]);
                }

                if (subcoords[1] != null)
                {
                    Marshal.FreeHGlobal((IntPtr) subcoords[1]);
                }

                if (subcoords[2] != null)
                {
                    Marshal.FreeHGlobal((IntPtr) subcoords[2]);
                }

                Marshal.FreeHGlobal((IntPtr) subcoords);
            }

            Marshal.FreeHGlobal((IntPtr) set_info);
            Marshal.FreeHGlobal((IntPtr) set_buckets);
        }

        /// <summary>
        ///
        /// </summary>
        /// <param name="graph">graph data structure</param>
        /// <param name="nvtxs">number of vertices in graph</param>
        /// <param name="nedges">number of edges in graph</param>
        /// <param name="useVertexWeights">are vertex weights being used?</param>
        /// <param name="useEdgeWeights">are edge weights being used?</param>
        /// <param name="vwsqrt">sqrt of vertex weights (length nvtxs+1)</param>
        /// <param name="igeom">geometry dimension for inertial method</param>
        /// <param name="coords">coordinates for inertial method</param>
        /// <param name="assignment">set number of each vtx (length n)</param>
        /// <param name="goal">desired set sizes</param>
        /// <param name="architecture">0 => hypercube, d => d-dimensional mesh</param>
        /// <param name="term_wgts">weights for terminal propagation</param>
        /// <param name="global_method">global partitioning algorithm</param>
        /// <param name="local_method">local partitioning algorithm</param>
        /// <param name="rqi_flag">should I use multilevel eigensolver?</param>
        /// <param name="vmax">if so, # vertices to coarsen down to</param>
        /// <param name="ndims">number of eigenvectors</param>
        /// <param name="eigtol">tolerance on eigenvectors</param>
        /// <param name="hop_mtx">between-set hop costs for KL</param>
        /// <param name="nsets">number of sets to partition into</param>
        /// <param name="striping">partition by striping into pieces?</param>
        private static void divide(vtx_data** graph,
            int nvtxs,
            int nedges,
            bool useVertexWeights,
            bool useEdgeWeights,
            double* vwsqrt,
            int igeom,
            float** coords,
            int* assignment,
            double[] goal,
            bool architecture,
            float*[] term_wgts,
            PartitioningStrategy global_method,
            LocalPartitioningStrategy local_method,
            bool rqi_flag,
            int vmax,
            int ndims,
            double eigtol,
            int[][] hop_mtx,
            int nsets,
            int striping
        )
        {
            connect_data* cdata; /* data for enforcing connectivity */
            double*[] yvecs = new double*[MAXDIMS + 1]; /* space for pointing to eigenvectors */
            double[] evals = new double[MAXDIMS + 1]; /* corresponding eigenvalues */
            double[] weights = new double[MAXSETS]; /* vertex weight in each set */
            double maxdeg; /* maximum weighted degree of a vertex */
            double total_weight; /* weight of all vertices */
            double[] temp_goal = new double[2]; /* goal to simulate bisection while striping */
            double[] fake_goal; /* either goal or temp_goal */
            int* active; /* keeping track of vtxs in BFS (length nvtxs) */
            int* bndy_list; /* list of boundary vtxs */
            int* null_ptr; /* argument to klspiff */
            int* mark; /* for marking vtxs in BFS (length nvtxs+1) */
            int vwgt_max; /* largest vertex weight */
            int max_dev; /* largest allowed deviation from balance */
            MappingType mediantype; /* how to map from eigenvectors to partition */
            LanczosType solver_flag; /* which Lanczos variant to use */
            bool mkconnected; /* enforce connectivity in spectral method? */
            int i; /* loop counters */
            MappingType simple_type; /* which type of simple partitioning to use */

            maxdeg = 0;

            if (DEBUG_TRACE)
            {
                Trace.WriteLine($"<Entering divide, nvtxs = {nvtxs:d}, nedges = {nedges:d}>");
            }

            if (nvtxs <= 0)
            {
                return;
            }

            if (nedges == 0 && global_method != PartitioningStrategy.Inertial)
            {
                global_method = PartitioningStrategy.Linear;
                local_method = LocalPartitioningStrategy.None;
            }

            for (i = 0; i < MAXSETS; i++)
            {
                weights[i] = 0.0;
            }

            bndy_list = null;

            if (striping != 0)
            {
                ndims = 1;
                mediantype = MappingType.Striped;

                temp_goal[0] = temp_goal[1] = 0;
                for (i = 0; 2 * i + 1 < nsets; i++)
                {
                    temp_goal[0] += goal[i];
                    temp_goal[1] += goal[nsets - 1 - i];
                }

                i = nsets / 2;
                if (2 * i != nsets)
                {
                    temp_goal[0] += .5 * goal[i];
                    temp_goal[1] += .5 * goal[i];
                }

                fake_goal = temp_goal;
            }
            else
            {
                mediantype = MAPPING_TYPE;
                fake_goal = goal;
            }

            if (useVertexWeights)
            {
                vwgt_max = 0;
                for (i = 1; i <= nvtxs; i++)
                {
                    if (graph[i]->vwgt > vwgt_max)
                    {
                        vwgt_max = graph[i]->vwgt;
                    }
                }
            }
            else
            {
                vwgt_max = 1;
            }

            /* Perform one of the global partitionings on this sub-graph. */
            if (global_method == PartitioningStrategy.Multilevel_KL)
            {
                /* Multilevel method. */
                mkconnected = MAKE_CONNECTED;
                solver_flag = LANCZOS_TYPE;

                if (VERTEX_SEPARATOR)
                {
                    coarsen_klv(graph, nvtxs, nedges, useVertexWeights, useEdgeWeights, term_wgts, igeom, coords,
                        vwgt_max, assignment, goal, architecture ? 1 : 0, hop_mtx, solver_flag, ndims, nsets,
                        vmax, mediantype, mkconnected, eigtol, COARSE_NLEVEL_KL, 0, &bndy_list, weights,
                        false);
                    if (VERTEX_COVER && !COARSE_BPM)
                    {
                        max_dev = vwgt_max;
                        total_weight = 0;
                        for (i = 0; i < nsets; i++)
                        {
                            total_weight += goal[i];
                        }

                        total_weight *= KL_IMBALANCE / nsets;
                        if (total_weight > max_dev)
                        {
                            max_dev = (int) total_weight;
                        }

                        bpm_improve(graph, assignment, goal, max_dev, &bndy_list, weights, useVertexWeights);
                    }
                }
                else
                {
                    coarsen_kl(graph, nvtxs, nedges, useVertexWeights, useEdgeWeights, term_wgts, igeom, coords, vwgt_max,
                        assignment, goal, architecture, hop_mtx, solver_flag, ndims, nsets, vmax,
                        mediantype, mkconnected, eigtol, COARSE_NLEVEL_KL, 0, &bndy_list, weights, false);

                    if (VERTEX_COVER)
                    {
                        max_dev = vwgt_max;
                        total_weight = 0;
                        for (i = 0; i < nsets; i++)
                        {
                            total_weight += goal[i];
                        }

                        total_weight *= KL_IMBALANCE / nsets;
                        if (total_weight > max_dev)
                        {
                            max_dev = (int) total_weight;
                        }

                        if (goal[0] - weights[0] <= goal[1] - weights[1])
                        {
                            i = 0;
                        }
                        else
                        {
                            i = 1;
                        }

                        find_side_bndy(graph, nvtxs, assignment, i, 2, &bndy_list);
                        count_weights(graph, nvtxs, assignment, nsets + 1, weights, useVertexWeights);
                        if (DEBUG_KL != DebugFlagKL.NoDebugging)
                        {
                            Trace.WriteLine("After KL, before bpm_improve");
                            countup_vtx_sep(graph, nvtxs, assignment);
                        }

                        bpm_improve(graph, assignment, goal, max_dev, &bndy_list, weights, useVertexWeights);
                    }
                }
            }

            else if (global_method == PartitioningStrategy.Spectral)
            {
                /* Spectral method. */
                mkconnected = MAKE_CONNECTED;
                solver_flag = LANCZOS_TYPE;

                active = (int*) Marshal.AllocHGlobal(nvtxs * sizeof(int));
                for (i = 1; i <= ndims; i++)
                {
                    yvecs[i] = (double*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(double));
                }

                if (mkconnected)
                {
                    /* If doing multi-level KL, this happens on coarse graph. */
                    mark = (int*) &(yvecs[1][0]);
                    make_connected(graph, nvtxs, &nedges, mark, active, &cdata, useEdgeWeights);
                    if (DEBUG_CONNECTED)
                    {
                        Trace.WriteLine("Enforcing connectivity");
                        print_connected(cdata);
                    }
                }
                else if (DEBUG_CONNECTED)
                {
                    Trace.WriteLine("Not enforcing connectivity");
                }

                maxdeg = find_maxdeg(graph, nvtxs, useEdgeWeights, (float*) null);
                eigensolve(graph, nvtxs, nedges, maxdeg, vwgt_max, vwsqrt, useVertexWeights, useEdgeWeights, term_wgts,
                    igeom, coords, yvecs, evals, architecture, assignment, fake_goal, solver_flag,
                    rqi_flag, vmax, ndims, mediantype, eigtol);

                Assign(graph, yvecs, nvtxs, ndims, architecture, nsets, vwsqrt, assignment, active, mediantype,
                    goal, vwgt_max);

                for (i = 1; i <= ndims; i++)
                {
                    Marshal.FreeHGlobal((IntPtr) yvecs[i]);
                }

                if (mkconnected)
                {
                    make_unconnected(graph, &nedges, &cdata, useEdgeWeights);
                }

                Marshal.FreeHGlobal((IntPtr) active);
            }

            else if (global_method == PartitioningStrategy.Inertial)
            {
                /* Inertial method. */
                inertial(graph, nvtxs, architecture, nsets, igeom, coords, assignment, goal, useVertexWeights);
            }

            else if (global_method == PartitioningStrategy.Linear)
            {
                /* Linear ordering. */
                simple_type = MappingType.IndependantMedians;
                simple_part(graph, nvtxs, assignment, nsets, simple_type, goal);
            }

            else if (global_method == PartitioningStrategy.Random)
            {
                /* Random partitioning. */
                simple_type = MappingType.RecursiveMedian;
                simple_part(graph, nvtxs, assignment, nsets, simple_type, goal);
            }

            else if (global_method == PartitioningStrategy.Scattered)
            {
                /* Scattered partitioning. */
                simple_type = MappingType.MinCost;
                simple_part(graph, nvtxs, assignment, nsets, simple_type, goal);
            }

            /* Perform a local refinement, if specified, on the global partitioning. */
            if (local_method == LocalPartitioningStrategy.KernighanLin && global_method != PartitioningStrategy.Multilevel_KL)
            {
                /* If global_method == 1, already did KL as part of multilevel-KL. */
                null_ptr = null;
                max_dev = vwgt_max;
                total_weight = 0;
                for (i = 0; i < nsets; i++)
                {
                    total_weight += goal[i];
                }

                total_weight *= KL_IMBALANCE / nsets;
                if (total_weight > max_dev)
                {
                    max_dev = (int) total_weight;
                }

                if (VERTEX_SEPARATOR)
                {
                    find_bndy(graph, nvtxs, assignment, 2, &bndy_list);
                    count_weights(graph, nvtxs, assignment, nsets + 1, weights, (vwgt_max != 1));
                    klvspiff(graph, nvtxs, assignment, goal, max_dev, &bndy_list, weights);
                    if (VERTEX_COVER)
                    {
                        bpm_improve(graph, assignment, goal, max_dev, &bndy_list, weights, useVertexWeights);
                    }
                }
                else
                {
                    if (global_method != PartitioningStrategy.Spectral)
                    {
                        maxdeg = find_maxdeg(graph, nvtxs, useEdgeWeights, (float*) null);
                    }

                    count_weights(graph, nvtxs, assignment, nsets, weights, (vwgt_max != 1));
                    klspiff(graph, nvtxs, assignment, nsets, hop_mtx, goal, term_wgts, max_dev, maxdeg,
                        useEdgeWeights, &null_ptr, weights);
                    if (VERTEX_COVER)
                    {
                        if (goal[0] - weights[0] <= goal[1] - weights[1])
                        {
                            i = 0;
                        }
                        else
                        {
                            i = 1;
                        }

                        find_side_bndy(graph, nvtxs, assignment, i, 2, &bndy_list);
                        count_weights(graph, nvtxs, assignment, nsets + 1, weights, (vwgt_max != 1));
                        if (DEBUG_KL != DebugFlagKL.NoDebugging)
                        {
                            Trace.WriteLine("After KL, before bpm_improve");
                            countup_vtx_sep(graph, nvtxs, assignment);
                        }

                        bpm_improve(graph, assignment, goal, max_dev, &bndy_list, weights, useVertexWeights);
                    }
                }
            }
            else if (global_method != PartitioningStrategy.Multilevel_KL && VERTEX_COVER)
            {
                find_bndy(graph, nvtxs, assignment, 2, &bndy_list);
                count_weights(graph, nvtxs, assignment, nsets + 1, weights, (vwgt_max != 1));
                if (DEBUG_KL != DebugFlagKL.NoDebugging)
                {
                    Trace.WriteLine("Before bpm_improve");
                    countup_vtx_sep(graph, nvtxs, assignment);
                }

                max_dev = vwgt_max;
                total_weight = 0;
                for (i = 0; i < nsets; i++)
                {
                    total_weight += goal[i];
                }

                total_weight *= KL_IMBALANCE / nsets;
                if (total_weight > max_dev)
                {
                    max_dev = (int) total_weight;
                }

                bpm_improve(graph, assignment, goal, max_dev, &bndy_list, weights, useVertexWeights);
            }

            if (bndy_list != null)
            {
                Marshal.FreeHGlobal((IntPtr) bndy_list);
            }
        }

        /// <summary>
        /// Populates entries in the subcoordinate array with values from the main coordinates accoding to the graph mapping.
        /// </summary>
        /// <param name="igeom">1, 2 or 3 dimensional geometry?</param>
        /// <param name="coords">x, y and z coordinates of vertices</param>
        /// <param name="subcoords">x, y and z coordinates in subgraph</param>
        /// <param name="subnvtxs">number of vertices in subgraph</param>
        /// <param name="loc2glob">maps from subgraph to graph numbering</param>
        private static void make_subgeom(int igeom, float** coords, float** subcoords, int subnvtxs, int* loc2glob)
        {
            if (subcoords == null)
            {
                throw new ArgumentNullException(nameof(subcoords));
            }

            switch (igeom)
            {
                case 1:
                {
                    if (subcoords[0] == null)
                    {
                        throw new ArgumentOutOfRangeException(nameof(subcoords));
                    }
                    if (coords[0] == null)
                    {
                        throw new ArgumentOutOfRangeException(nameof(coords));
                    }

                    for (var i = 1; i <= subnvtxs; i++)
                    {
                        subcoords[0][i] = coords[0][loc2glob[i]];
                    }

                    return;
                }
                case 2:
                {
                    if (subcoords[0] == null || subcoords[1] == null)
                    {
                        throw new ArgumentOutOfRangeException(nameof(subcoords));
                    }
                    if (coords[0] == null || coords[1] == null)
                    {
                        throw new ArgumentOutOfRangeException(nameof(coords));
                    }

                    for (var i = 1; i <= subnvtxs; i++)
                    {
                        subcoords[0][i] = coords[0][loc2glob[i]];
                        subcoords[1][i] = coords[1][loc2glob[i]];
                    }

                    return;
                }
                case 3:
                {
                    if (subcoords[0] == null || subcoords[1] == null || subcoords[2] == null)
                    {
                        throw new ArgumentOutOfRangeException(nameof(subcoords));
                    }
                    if (coords[0] == null || coords[1] == null || coords[2] == null)
                    {
                        throw new ArgumentOutOfRangeException(nameof(coords));
                    }
                    for (var i = 1; i <= subnvtxs; i++)
                    {
                        subcoords[0][i] = coords[0][loc2glob[i]];
                        subcoords[1][i] = coords[1][loc2glob[i]];
                        subcoords[2][i] = coords[2][loc2glob[i]];
                    }

                    return;
                }
                default:
                {
                    throw new ArgumentOutOfRangeException(nameof(igeom), $"Geometry dimension should be in [1,3], but was {igeom}");
                }
            }
        }
    }
}
