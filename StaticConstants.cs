// ReSharper disable InconsistentNaming
using System;
using ChacoSharp.Coarsening;

#pragma warning disable S2223 // Non-constant static fields should not be visible
#pragma warning disable S1104 // Fields should not have public accessibility
namespace ChacoSharp
{
    public static unsafe class StaticConstants
    {
        /// <summary>
        /// Just generate spectral ordering?
        /// </summary>
        public const bool SEQUENCE = true;

        public const bool MAKE_VWGTS = false; /* Force vertex weights to be degrees+1 ? */

        public const int LANCZOS_CONVERGENCE_MODE = 0; /* residual or partition mode? */
        public const int RQI_CONVERGENCE_MODE = 0; /* residual or partition mode? */

        public const double EIGEN_TOLERANCE = 0.0001; /* eigen-tolerance convergence criteria */
        public const LanczosType LANCZOS_TYPE = LanczosType.SelectiveOrthogonalization; /* type of Lanczos to use */
        public const int LANCZOS_SO_INTERVAL = 10; /* interval between orthogonalizations */
        public const int LANCZOS_SO_PRECISION = 2; /* controls precision in eigen calc. */
        public static double SRESTOL = -1.0d; /* resid tol for T evec computation */
        public static int LANCZOS_MAXITNS = -1; /* maximum Lanczos iterations allowed */

        //TODO: JB: Is this used in my path?
        public static KernighanLinMetric KL_METRIC = KernighanLinMetric.Cuts; /* KL interset cost: 1=>cuts, 2=>hops */

        //TODO: JB: Is this used in my path?
        public const double BISECTION_SAFETY = 10; /* safety factor for bisection algorithm */

        //TODO: JB: Is this used in my path?
        public const MatchingRoutine MATCH_TYPE = MatchingRoutine.maxmatch1; /* which matching routine to use */

        public enum MatchingRoutine
        {
            maxmatch1 = 1,
            maxmatch2 = 2,
            maxmatch3 = 3,
            maxmatch4_Luby = 4,
            maxmatch5_geometric = 5,
            maxmatch9_minimumVertexDegree = 9
        }

        public static int N_VTX_MOVES; /* total number of vertex moves */
        public static int N_VTX_CHECKS; /* total number of moves contemplated */
        public static int KL_MAX_PASS = -1; /* max KL passes; infinite if <= 0 */

        public static int SRES_SWITCHES = 0; /* # switches to backup routine for computing s */
        public static bool PERTURB = false; /* perturb matrix? */
        public static int NSQRTS = 0; /* number of sqrts already computed */
        public static double* SQRTS = null; /* values computed */

        public const double WARNING_ORTHTOL = 2; /* Warning: modest loss of orthogonality */
        public const double WARNING_MISTOL = 100; /* Warning: serious loss of orthogonality */

        public const int OUTPUT_ASSIGN = 0; /* print assignment to a file? */
        public const string Graph_File_Name = null; /* Input graph file name */
        public const string Geometry_File_Name = null; /* Input coordinate file name */
        public const string Assign_In_File_Name = null; /* Input assignment file name */

        public const bool FREE_GRAPH = false; /* free graph data structure after reformat? */

        public const int PROJECTION_AXIS = 0; /* axis to flatten geometry */

        public const bool VERTEX_SEPARATOR = false;


        public const bool FLATTEN = false; /* Merge indistinguishable vtxs first? */
        public const bool COARSE_KLV = true; /* Use KLV as multilevel refinement? */
        public const bool COARSE_BPM = false; /* Use vertex cover as ML refinement? */
        public const bool COARSE_KL_BOTTOM = false; /* force KL invocation at bottom level? */
        public const bool KL_ONLY_BNDY = false; /* start KL w/ vertices on boundary? */
        public const bool LIMIT_KL_EWGTS = false; /* limit range of edge weights in multilevel-KL? */
        public const bool KL_RANDOM = false; /* use randomness in Kernighan-Lin? */
        public const bool KL_UNDO_LIST = false; /* only resort vtxs affected by last pass? */
        public const double KL_IMBALANCE = 0.0d; /* allowed fractional imbalance in KL */
        public const double EWGT_RATIO_MAX = 0.0d; /* if so, max allowed ewgt/nvtxs */
        public const int KL_BAD_MOVES = 0; /* number of unhelpful moves in a row allowed */
        public const int KL_NTRIES_BAD = 0; /* # unhelpful passes before quitting KL */
        public const int COARSE_NLEVEL_KL = 0; /* # levels between KL calls in uncoarsening */

        public const bool REFINE_MAP = false; /* use greedy strategy to improve mapping? */
        public const bool VERTEX_COVER = false; /* use matching to reduce vertex separator? */
        public const bool CONNECTED_DOMAINS = false; /* force subdomain connectivity at end? */
        public const bool INTERNAL_VERTICES = false; /* greedily increase internal vtxs? */
        public const int REFINE_PARTITION = 0; /* number of calls to pairwise_refine to make */

        public const int OPT3D_NTRIES = 0; /* number of local opts to find global min */

        public const MappingType MAPPING_TYPE = MappingType.MinCost; /* how to map from eigenvectors to partition */
        public const bool MAKE_CONNECTED = false; /* connect graph for spectral method? */

        public const bool COARSEN_VWGTS = false; /* use vertex weights while coarsening? */
        public const bool COARSEN_EWGTS = false; /* use edge weights while coarsening? */
        public const bool HEAVY_MATCH = false; /* choose heavy edges in matching? */
        public const double COARSEN_RATIO_MIN = 0.0d; /* min vtx reduction for coarsening */
        public const double PERTURB_MAX = 0.0d; /* maximum perturbation */
        public const int MAXDIMS = 3; /**< Most cuts allowed at one time */
        public const int MAXSETS = 8; /**< 2^MAXDIMS */
        public const int NPERTURB = 0; /* number of edges to perturb */
        public const int COARSE_NLEVEL_RQI = 0; /* do RQI this often in uncoarsening */

        public const bool TERM_PROP = false; /* invoke terminal propagation? */


        public const int SIMULATOR = 0; /* simulate the communication? */
        public const int SIMULATION_ITNS = 0; /* simulator iterations */
        public const double CUT_TO_HOP_COST = 0.0d; /* ..if so, importance of cuts/hops */
        public const double CUT_COST = 0.0d; /* cost of each cut */
        public const double HOP_COST = 0.0d; /* cost of each hop */
        public const double BDY_COST = 0.0d; /* cost of each boundary vertex */
        public const double BDY_HOP_COST = 0.0d; /* cost of each boundary vertex hop */
        public const double STARTUP_COST = 0.0d; /* initiation cost of each message */

        public const bool EXPERT = false; /* user type */
        public const bool DEBUG_CONNECTED = false; /* Debug flag for connected components */
        public const bool DEBUG_INTERNAL = false; /* debug code about force_internal? */
        public const bool DEBUG_REFINE_PART = false; /* debug code about refine_part? */
        public const bool DEBUG_REFINE_MAP = false; /* debug code about refine_map? */
        public const bool DEBUG_MACH_PARAMS = false; /* print out computed machine params? */
        public const bool DEBUG_ASSIGN = false; /* turn on debugging in assignment */
        public const bool DEBUG_INERTIAL = false; /* Debug flag for inertial method */
        public const bool DEBUG_TRACE = true; /* trace the execution of the code */
        public const int DEBUG_COVER = 0; /* Debug output in min vtx cover routines? */
        public const int DEBUG_EVECS = 2; /* debug flag for eigenvector generation */
        public const int WARNING_EVECS = 5; /* print warning messages? 0 to 5, increasing amount of info */
        public const bool DEBUG_OPTIMIZE = false; /* debug flag for optimization */
        public const DebugFlagBP DEBUG_BPMATCH = DebugFlagBP.ErrorChecking; /* turn on debugging for bipartite matching */
        public const bool DEBUG_PERTURB = false; /* debug flag for matrix perturbation */
        public const bool DEBUG_COARSEN = false; /* debug flag for coarsening */
        public const bool FullTrace = true;
        public const bool PRINT_GRAPH_PARTITION_METRICS = false; // OUTPUT METRICS
        public const bool PRINT_GRAPH_PARTITION_METRICS_ALL_RECURSION_LEVELS = false; // OUTPUT METRICS negative
        public const bool PRINT_GRAPH_PARTITION_METRICS_DETAILED = false; // OUTPUT METRICS false = 1, true = 2
        public const bool PRINT_HEADERS = true; /* print section headings for output? */
        public const bool ECHO_USER_PARAMETERS = false; // ECHO: case != 0
        public const bool ECHO_INPUT_PARAMETERS = false; // ECHO: case 2 or -2
        public const bool CHECK_INPUT = true; /* should I check input for correctness? */

        public const DebugFlagKL DEBUG_KL = DebugFlagKL.ImprovementsPerStep; /* Debug flag for Kernighan-Lin */




        public const bool LANCZOS_TIME = false; /* perform detailed timing on Lanczos_SO? */
        public const bool TIME_KERNELS = false; /* benchmark some numerical kernels? */
        public const int OUTPUT_TIME = 0;

        public static double splarax_time; /* time matvecs */
        public static double orthog_time; /* time orthogonalization work */
        public static double tevec_time; /* time tridiagonal eigvec work */
        public static double evec_time; /* time to generate eigenvectors */
        public static double ql_time; /* time tridiagonal eigval work */
        public static double blas_time; /* time for blas (not assembly coded) */
        public static double init_time; /* time for allocating memory, etc. */
        public static double scan_time; /* time for scanning bounds list */
        public static double debug_time; /* time for debug computations and output */
        public static double ritz_time; /* time to generate ritz vectors */
        public static double pause_time; /* time to compute whether to pause */
        public static double refine_time; /* time for RQI/Symmlq iterative refinement */
        public static double coarsen_time;
        public static double match_time;
        public static double make_cgraph_time;
        public static double lanczos_time; /* time spent in Lanczos algorithm */
        public static double rqi_symmlq_time; /* time spent in RQI/Symmlq method */
        public static double start_time; /* time code was entered */
        public static double total_time; /* (almost) total time spent in code */
        public static double check_input_time; /* time spent checking input */
        public static double partition_time; /* time spent partitioning graph */
        public static double kernel_time; /* time spent benchmarking kernels */
        public static double count_time; /* time spent evaluating the answer */
        public static double kl_bucket_time; /* time spent in KL bucketsort */
        public static double kl_total_time;
        public static double kl_init_time;
        public static double nway_kl_time;
        public static double reformat_time; /* time spent reformatting graph */
        public static double inertial_axis_time; /* time spent finding inertial axis */
        public static double median_time; /* time to find medians */
        public static double inertial_time; /* time spend in inertial calculations */

        /// <summary>
        /// Machine Precision.
        /// </summary>
        /// <remarks>
        /// The value of the Epsilon property is not equivalent to machine epsilon,
        /// which represents the upper bound of the relative error due to rounding in floating-point arithmetic.
        /// </remarks>
        public static readonly double DOUBLE_EPSILON = CalculateMachineEpsilon();

        private static double CalculateMachineEpsilon()
        {
            var eps = 1.0d / 16.0d;
            while (1.0d + eps > 1.0d)
            {
                eps /= 2.0d;
            }

            return eps;
        }

        public enum MappingType
        {
            CutAtOrigin = 0,
            MinCost = 1,
            RecursiveMedian = 2,
            IndependantMedians = 3,
            Striped
        }

        public enum LanczosType
        {
            FullOrthogonalization = 1,
            FullOrthogonalizationInverseOperator = 2,
            SelectiveOrthogonalization = 3,
            SelectiveOrthogonalizationDoubleEnded = 4
        }

        public enum PartitioningStrategy
        {
            Multilevel_KL = 1,
            Spectral = 2,
            Inertial = 3,
            Linear = 4,
            Random = 5,
            Scattered = 6,
            ReadFromFile = 7
        }

        public enum LocalPartitioningStrategy
        {
            KernighanLin = 1,
            None = 2
        }

        public enum KernighanLinMetric
        {
            Cuts = 1,
            Hops = 2
        }

        public enum DebugFlagKL
        {
            NoDebugging = 0,
            ImprovementsPerStep = 1,
            MoreInfo = 2,
            PrintBucket = 3
        }

        public enum DebugFlagBP
        {
            NoDebugging = 0,
            Logging = 1,
            ErrorChecking = 2
        }
    }
}
#pragma warning restore S1104 // Fields should not have public accessibility
#pragma warning restore S2223 // Non-constant static fields should not be visible
