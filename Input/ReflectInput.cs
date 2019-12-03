#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
using System;
using static ChacoSharp.StaticConstants;

namespace ChacoSharp.Input
{
    public static class ReflectInput
    {
        /* Print out the input options. */
        public static void reflect_input(int nvtxs, /* number of vertices in graph */
            int nedges, /* number of edges in graph */
            int igeom, /* geometric dimension for inertial method */
            string graphname, /* name of graph input file */
            string geomname, /* name of geometry input file */
            string inassignname, /* name of assignment input file */
            int architecture, /* 0=> hypercube, d=> d-dimensional mesh */
            int ndims_tot, /* total number of cuts to make */
            int[] mesh_dims /*[3]*/, /* size of mesh */
            PartitioningStrategy global_method, /* global partitioning algorithm */
            LocalPartitioningStrategy localPartitioningStrategy, /* local partitioning algorithm */
            bool rqi_flag, /* use RQI/Symmlq eigensolver?  */
            int vmax, /* smallest acceptable coarsened nvtxs */
            int ndims, /* partitioning level */
            double eigtol, /* tolerance on eigenvectors */
            long seed) /* random number seed */
        {

            if (DEBUG_TRACE)
            {
                Console.WriteLine("<Entering reflect_input>");
            }

            Console.WriteLine();

            if (PRINT_HEADERS)
            {
                Console.WriteLine("\n           Input and Parameter Values\n");
            }

            if (graphname != null)
            {
                Console.WriteLine("Graph file: `{0}', ", graphname);
            }

            Console.WriteLine("# vertices = {0:d}, # edges = {1:d}", nvtxs, nedges);

            /* Print global partitioning strategy. */

            Console.WriteLine("Global method: ");
            if (global_method == PartitioningStrategy.Multilevel_KL)
            {
                Console.WriteLine("Multilevel-KL");
            }
            else if (global_method == PartitioningStrategy.Spectral)
            {
                Console.WriteLine("Spectral");
            }
            else if (global_method == PartitioningStrategy.Inertial)
            {
                Console.WriteLine("Inertial");
            }
            else if (global_method == PartitioningStrategy.Linear)
            {
                Console.WriteLine("Linear");
            }
            else if (global_method == PartitioningStrategy.Random)
            {
                Console.WriteLine("Random");
            }
            else if (global_method == PartitioningStrategy.Scattered)
            {
                Console.WriteLine("Scattered");
            }
            else if (global_method == PartitioningStrategy.ReadFromFile)
            {
                Console.WriteLine("Read-From-File");
            }

            if (global_method == PartitioningStrategy.Multilevel_KL)
            {
                Console.WriteLine("Number of vertices to coarsen down to: {0:d}", vmax);
                Console.WriteLine("Eigen tolerance: {0:g}", eigtol);
            }

            else if (global_method == PartitioningStrategy.Spectral)
            {
                if (rqi_flag)
                {
                    Console.WriteLine("Multilevel RQI/Symmlq eigensolver");
                    Console.WriteLine("Number of vertices to coarsen down to: {0:d}", vmax);
                    Console.WriteLine("Eigen tolerance: {0:g}", eigtol);
                }
            }

            else if (global_method == PartitioningStrategy.Inertial)
            {
                if (geomname != null)
                {
                    Console.WriteLine("Geometry input file: `{0}', Dimensionality = {1:d}", geomname, igeom);
                }
            }

            else if (global_method == PartitioningStrategy.ReadFromFile)
            {
                Console.WriteLine("Assignment input file: `{0}'", inassignname);
            }

            /* Now describe local method. */
            if (localPartitioningStrategy == LocalPartitioningStrategy.KernighanLin)
            {
                Console.WriteLine("Local method: Kernighan-Lin");
            }
            else if (localPartitioningStrategy == LocalPartitioningStrategy.None)
            {
                Console.WriteLine("Local method: None");
            }

            /* Now describe target architecture. */
            if (architecture == 0)
            {
                Console.WriteLine("Partitioning target: {0:d}-dimensional hypercube", ndims_tot);
            }
            else if (architecture > 0)
            {
                Console.Write("Partitioning target: {0:d}-dimensional mesh of size ", architecture);
                if (architecture == 1)
                {
                    Console.WriteLine("{0:d}", mesh_dims[0]);
                }
                else if (architecture == 2)
                {
                    Console.WriteLine("{0:d}x{1:d}", mesh_dims[0], mesh_dims[1]);
                }
                else if (architecture == 3)
                {
                    Console.WriteLine("{0:d}x{1:d}x{2:d}", mesh_dims[0], mesh_dims[1], mesh_dims[2]);
                }
            }

            if (ndims == 1)
            {
                Console.WriteLine("Partitioning mode: Bisection");
            }
            else if (ndims == 2)
            {
                Console.WriteLine("Partitioning mode: Quadrisection");
            }
            else if (ndims == 3)
            {
                Console.WriteLine("Partitioning mode: Octasection");
            }

            /* Add stuff about communications simulator. */

            Console.WriteLine("Random seed: {0:d}", seed);

            if (ECHO_USER_PARAMETERS)
            {
                reflect_params(global_method, localPartitioningStrategy, rqi_flag, ndims);
            }

            Console.WriteLine("\n");
        }

        public static void reflect_params(PartitioningStrategy global_method, /* global partitioning algorithm */
            LocalPartitioningStrategy localPartitioningStrategy, /* local partitioning algorithm */
            bool rqi_flag, /* use RQI/SYMMLQ eigensolver? */
            int ndims) /* number of eigenvectors to generate */
        {
            Console.WriteLine("Active Parameters:");

            Console.WriteLine("  CHECK_INPUT = {0}", CHECK_INPUT ? "True" : "False");

            if (global_method == PartitioningStrategy.Multilevel_KL || global_method == PartitioningStrategy.Linear)
            {
                Console.WriteLine("  LANCZOS_TYPE:  ");
                if (LANCZOS_TYPE == LanczosType.FullOrthogonalization)
                {
                    Console.Write(" Full orthogonalization");
                }
                else if (LANCZOS_TYPE == LanczosType.FullOrthogonalizationInverseOperator)
                {
                    Console.Write("Full orthogonalization, inverse operator");
                }
                else if (LANCZOS_TYPE == LanczosType.SelectiveOrthogonalization)
                {
                    Console.Write("Selective orthogonalization");
                }
                else if (LANCZOS_TYPE == LanczosType.SelectiveOrthogonalizationDoubleEnded)
                {
                    if (EXPERT)
                    {
                        Console.Write("Selective orthogonalization against both ends");
                    }
                    else
                    {
                        /* Check input should catch this, but just in case ... */
                        //LANCZOS_TYPE = LanczosType.SelectiveOrthogonalization;
                        throw new InvalidOperationException("Only use " + nameof(LanczosType.SelectiveOrthogonalizationDoubleEnded) + " in " + nameof(EXPERT) + " mode");
                        Console.Write("Selective orthogonalization");
                    }
                }

                Console.WriteLine();

                Console.WriteLine("  EIGEN_TOLERANCE = {0:g}", EIGEN_TOLERANCE);

                if (SRESTOL > 0)
                {
                    Console.WriteLine("  SRESTOL = {0:g}", SRESTOL);
                }
                else
                {
                    Console.WriteLine("  SRESTOL = {0:g} ... autoset to square of eigen tolerance", SRESTOL);
                }

                if (LANCZOS_MAXITNS > 0)
                {
                    Console.WriteLine("  LANCZOS_MAXITNS = {0:d}", LANCZOS_MAXITNS);
                }
                else
                {
                    Console.WriteLine("  LANCZOS_MAXITNS = {0:d} ... autoset to twice # vertices", LANCZOS_MAXITNS);
                }

                if (LANCZOS_SO_PRECISION == 1)
                {
                    Console.WriteLine("  LANCZOS_SO_PRECISION = 1 ... single precision");
                }
                else
                {
                    Console.WriteLine("  LANCZOS_SO_PRECISION = 2 ... double precision");
                }

                Console.WriteLine("  LANCZOS_SO_INTERVAL = {0:d}", LANCZOS_SO_INTERVAL);

                if (LANCZOS_CONVERGENCE_MODE == 1)
                {
                    Console.WriteLine("  LANCZOS_CONVERGENCE_MODE = 1 ... partition tolerance");
                }
                else
                {
                    Console.WriteLine("  LANCZOS_CONVERGENCE_MODE = 0 ... residual tolerance");
                }

                Console.WriteLine("  BISECTION_SAFETY = {0:g}", BISECTION_SAFETY);
                if (LANCZOS_TYPE == LanczosType.SelectiveOrthogonalization || LANCZOS_TYPE == LanczosType.SelectiveOrthogonalizationDoubleEnded)
                {
                    if (!LANCZOS_TIME)
                    {
                        Console.WriteLine("  LANCZOS_TIME = 0 ... no detailed timing");
                    }
                    else
                    {
                        Console.WriteLine("  LANCZOS_TIME = 1 ... detailed timing");
                    }
                }

                if (WARNING_EVECS > 0)
                {
                    Console.WriteLine("  WARNING_EVECS = {0:d}", WARNING_EVECS);
                }

                if (MAPPING_TYPE == MappingType.CutAtOrigin)
                {
                    Console.WriteLine("  MAPPING_TYPE = 0 ... cut at origin");
                }
                else if (MAPPING_TYPE == MappingType.MinCost)
                {
                    Console.WriteLine("  MAPPING_TYPE = 1 ... min-cost assignment");
                }
                else if (MAPPING_TYPE == MappingType.RecursiveMedian)
                {
                    Console.WriteLine("  MAPPING_TYPE = 2 ... recursive median");
                }
                else if (MAPPING_TYPE == MappingType.IndependantMedians)
                {
                    Console.WriteLine("  MAPPING_TYPE = 3 ... independent medians");
                }

                Console.WriteLine("  MAKE_CONNECTED = {0}", MAKE_CONNECTED ? "True" : "False");
                Console.WriteLine("  PERTURB = {0}", PERTURB ? "True" : "False");
                if (PERTURB)
                {
                    Console.WriteLine("    NPERTURB = {0:d}", NPERTURB);
                    Console.WriteLine("    PERTURB_MAX = {0:g}", PERTURB_MAX);
                }

                if (ndims == 3)
                {
                    Console.WriteLine("  OPT3D_NTRIES = {0:d}", OPT3D_NTRIES);
                }
            }

            if (global_method == PartitioningStrategy.Multilevel_KL)
            {
                Console.WriteLine("  COARSEN_RATIO_MIN = {0:g}", COARSEN_RATIO_MIN);
                Console.WriteLine("  COARSE_NLEVEL_KL = {0:d}", COARSE_NLEVEL_KL);
                Console.WriteLine("  MATCH_TYPE = {0:d}", MATCH_TYPE);
                Console.WriteLine("  HEAVY_MATCH = {0}", HEAVY_MATCH ? "True" : "False");
                Console.WriteLine("  COARSE_KL_BOTTOM = {0}", COARSE_KL_BOTTOM ? "True" : "False");
                Console.WriteLine("  COARSEN_VWGTS = {0}", COARSEN_VWGTS ? "True" : "False");
                Console.WriteLine("  COARSEN_EWGTS = {0}", COARSEN_EWGTS ? "True" : "False");
                Console.WriteLine("  KL_ONLY_BNDY = {0}", KL_ONLY_BNDY ? "True" : "False");
            }

            if (global_method == PartitioningStrategy.Spectral && rqi_flag)
            {
                Console.WriteLine("  COARSE_NLEVEL_RQI = {0:d}", COARSE_NLEVEL_RQI);
                if (RQI_CONVERGENCE_MODE == 1)
                {
                    Console.WriteLine("  RQI_CONVERGENCE_MODE = 1 ... partition tolerance");
                }
                else
                {
                    Console.WriteLine("  RQI_CONVERGENCE_MODE = 0 ... residual tolerance");
                }

                Console.WriteLine("  COARSEN_RATIO_MIN = {0:g}", COARSEN_RATIO_MIN);
                Console.WriteLine("  COARSEN_VWGTS = {0}", COARSEN_VWGTS ? "True" : "False");
                Console.WriteLine("  COARSEN_EWGTS = {0}", COARSEN_EWGTS ? "True" : "False");
            }

            if (global_method == PartitioningStrategy.Multilevel_KL || localPartitioningStrategy == LocalPartitioningStrategy.KernighanLin)
            {
                Console.WriteLine("  KL_RANDOM = {0}", KL_RANDOM ? "True" : "False");
                if (KL_METRIC == KernighanLinMetric.Cuts)
                {
                    Console.WriteLine("  KL_METRIC = Cuts");
                }
                else if (KL_METRIC == KernighanLinMetric.Hops)
                {
                    Console.WriteLine("  KL_METRIC = Hops");
                }

                Console.WriteLine("  KL_NTRIES_BAD = {0:d}", KL_NTRIES_BAD);
                Console.WriteLine("  KL_BAD_MOVES = {0:d}", KL_BAD_MOVES);
                Console.WriteLine("  KL_UNDO_LIST = {0}", KL_UNDO_LIST ? "True" : "False");
                Console.WriteLine("  KL_IMBALANCE = {0:g}", KL_IMBALANCE);
            }

            if (global_method == PartitioningStrategy.Multilevel_KL || global_method == PartitioningStrategy.Spectral || localPartitioningStrategy == LocalPartitioningStrategy.KernighanLin)
            {
                Console.WriteLine("  TERM_PROP = {0}", TERM_PROP ? "True" : "False");
                if (TERM_PROP)
                {
                    Console.WriteLine("    CUT_TO_HOP_COST = {0:g}", CUT_TO_HOP_COST);
                }
            }

            if (SEQUENCE)
            {
                Console.WriteLine("  SEQUENCE = {0}", SEQUENCE ? "True" : "False");
            }

            Console.WriteLine("  PRINT_GRAPH_PARTITION_METRICS = {0}", PRINT_GRAPH_PARTITION_METRICS ? "True" : "False");
            Console.WriteLine("  PRINT_GRAPH_PARTITION_METRICS_ALL_RECURSION_LEVELS = {0}", PRINT_GRAPH_PARTITION_METRICS_ALL_RECURSION_LEVELS ? "True" : "False");
            Console.WriteLine("  PRINT_GRAPH_PARTITION_METRICS_DETAILED = {0}", PRINT_GRAPH_PARTITION_METRICS_DETAILED ? "True" : "False");
            Console.WriteLine("  MAKE_VWGTS = {0}", MAKE_VWGTS ? "True" : "False");
            Console.WriteLine("  REFINE_MAP = {0}", REFINE_MAP ? "True" : "False");
            Console.WriteLine("  REFINE_PARTITION = {0:d}", REFINE_PARTITION);
            Console.WriteLine("  INTERNAL_VERTICES = {0}", INTERNAL_VERTICES ? "True" : "False");

            if (SIMULATOR != 0)
            {
                Console.WriteLine("  SIMULATOR = {0:d}", SIMULATOR);
                Console.WriteLine("  SIMULATION_ITNS = {0:d}", SIMULATION_ITNS);
                Console.WriteLine("  CUT_COST = {0:g}", CUT_COST);
                Console.WriteLine("  HOP_COST = {0:g}", HOP_COST);
                Console.WriteLine("  BDY_COST = {0:g}", BDY_COST);
                Console.WriteLine("  BDY_HOP_COST = {0:g}", BDY_HOP_COST);
                Console.WriteLine("  STARTUP_COST = {0:g}", STARTUP_COST);
            }

            /* Now print out all the nonzero debug parameters. */
            if (DEBUG_CONNECTED)
            {
                Console.WriteLine("  DEBUG_CONNECTED = True");
            }

            if (DEBUG_PERTURB)
            {
                Console.WriteLine("  DEBUG_PERTURB = True");
            }

            if (DEBUG_ASSIGN)
            {
                Console.WriteLine("  DEBUG_ASSIGN = True");
            }

            if (DEBUG_INERTIAL)
            {
                Console.WriteLine("  DEBUG_INERTIAL = True");
            }

            if (DEBUG_OPTIMIZE)
            {
                Console.WriteLine("  DEBUG_OPTIMIZE = True");
            }

            if (DEBUG_BPMATCH != DebugFlagBP.NoDebugging)
            {
                Console.WriteLine("  DEBUG_BPMATCH = " + DEBUG_BPMATCH);
            }

            if (DEBUG_COARSEN)
            {
                Console.WriteLine("  DEBUG_COARSEN = True");
            }

            if (DEBUG_EVECS != 0)
            {
                Console.WriteLine("  DEBUG_EVECS = {0:d}", DEBUG_EVECS);
            }

            if (DEBUG_KL != DebugFlagKL.NoDebugging)
            {
                Console.WriteLine("  DEBUG_KL = " + DEBUG_KL);
            }

            if (DEBUG_INTERNAL)
            {
                Console.WriteLine("  DEBUG_INTERNAL = True");
            }

            if (DEBUG_REFINE_PART)
            {
                Console.WriteLine("  DEBUG_REFINE_PART = True");
            }

            if (DEBUG_REFINE_MAP)
            {
                Console.WriteLine("  DEBUG_REFINE_MAP = True");
            }

            if (DEBUG_TRACE)
            {
                Console.WriteLine("  DEBUG_TRACE = {0:d}", DEBUG_TRACE);
            }

            if (DEBUG_MACH_PARAMS)
            {
                Console.WriteLine("  DEBUG_MACH_PARAMS = True");
            }
        }
    }
}
