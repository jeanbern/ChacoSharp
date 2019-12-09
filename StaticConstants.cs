// ReSharper disable InconsistentNaming

#pragma warning disable S2223 // Non-constant static fields should not be visible
#pragma warning disable S1104 // Fields should not have public accessibility
#pragma warning disable RCS1057 // Add empty line between declarations.
namespace ChacoSharp
{
    public static unsafe class StaticConstants
    {
        /// <summary>Just generate spectral ordering?</summary>
        public const bool SEQUENCE = true;

        /// <summary>Force vertex weights to be degrees+1 ?</summary>
        public const bool MAKE_VWGTS = false;

        /// <summary>eigen-tolerance convergence criteria</summary>
        public const double EIGEN_TOLERANCE = 0.0001;
        /// <summary>resid tol for T evec computation</summary>
        public static double SRESTOL = -1.0d;

        /// <summary>residual or partition mode?</summary>
        public const int LANCZOS_CONVERGENCE_MODE = 0;
        /// <summary>residual or partition mode?</summary>
        public const int RQI_CONVERGENCE_MODE = 0;
        /// <summary>type of Lanczos to use</summary>
        public const LanczosType LANCZOS_TYPE = LanczosType.SelectiveOrthogonalization;
        /// <summary>interval between orthogonalizations</summary>
        public const int LANCZOS_SO_INTERVAL = 10;
        /// <summary>controls precision in eigen calc.</summary>
        public const int LANCZOS_SO_PRECISION = 2;
        /// <summary>maximum Lanczos iterations allowed</summary>
        public static int LANCZOS_MAXITNS = -1;

        //TODO: JB: Is this used in my path?
        /// <summary>KL interset cost: 1=>cuts, 2=>hops</summary>
        public static KernighanLinMetric KL_METRIC = KernighanLinMetric.Cuts;

        //TODO: JB: Is this used in my path?
        /// <summary>safety factor for bisection algorithm</summary>
        public const double BISECTION_SAFETY = 10;

        //TODO: JB: Is this used in my path?
        /// <summary>which matching routine to use</summary>
        public const MatchingRoutine MATCH_TYPE = MatchingRoutine.maxmatch1;

        public enum MatchingRoutine
        {
            maxmatch1 = 1,
            maxmatch2 = 2,
            maxmatch3 = 3,
            maxmatch4_Luby = 4,
            maxmatch5_geometric = 5,
            maxmatch9_minimumVertexDegree = 9
        }

        /// <summary>total number of vertex moves</summary>
        public static int N_VTX_MOVES;
        /// <summary>total number of moves contemplated</summary>
        public static int N_VTX_CHECKS;
        /// <summary>max KL passes; infinite if &lt;= 0</summary>
        public static int KL_MAX_PASS = -1;
        /// <summary># switches to backup routine for computing s</summary>
        public static int SRES_SWITCHES = 0;
        /// <summary>perturb matrix?</summary>
        public static bool PERTURB = false;
        /// <summary>number of sqrts already computed</summary>
        public static int NSQRTS = 0;
        /// <summary>values computed</summary>
        public static double* SQRTS = null;
        /// <summary>Warning: modest loss of orthogonality</summary>
        public const double WARNING_ORTHTOL = 2;
        /// <summary>Warning: serious loss of orthogonality</summary>
        public const double WARNING_MISTOL = 100;
        /// <summary>Input coordinate file name</summary>
        public const string Geometry_File_Name = null;
        /// <summary>Input assignment file name</summary>
        public const string Assign_In_File_Name = null;
        /// <summary>free graph data structure after reformat?</summary>
        public const bool FREE_GRAPH = false;
        /// <summary>axis to flatten geometry</summary>
        public const int PROJECTION_AXIS = 0;
        // TODO:
        /// <summary></summary>
        public const bool VERTEX_SEPARATOR = false;
        /// <summary>Merge indistinguishable vtxs first?</summary>
        public const bool FLATTEN = false;
        /// <summary>Use KLV as multilevel refinement?</summary>
        public const bool COARSE_KLV = true;
        /// <summary>Use vertex cover as ML refinement?</summary>
        public const bool COARSE_BPM = false;
        /// <summary>force KL invocation at bottom level?</summary>
        public const bool COARSE_KL_BOTTOM = false;
        /// <summary>start KL w/ vertices on boundary?</summary>
        public const bool KL_ONLY_BNDY = false;
        /// <summary>limit range of edge weights in multilevel-KL?</summary>
        public const bool LIMIT_KL_EWGTS = false;
        /// <summary>use randomness in Kernighan-Lin?</summary>
        public const bool KL_RANDOM = false;
        /// <summary>only resort vtxs affected by last pass?</summary>
        public const bool KL_UNDO_LIST = false;
        /// <summary>allowed fractional imbalance in KL</summary>
        public const double KL_IMBALANCE = 0.0d;
        /// <summary>if so, max allowed ewgt/nvtxs</summary>
        public const double EWGT_RATIO_MAX = 0.0d;
        /// <summary>number of unhelpful moves in a row allowed</summary>
        public const int KL_BAD_MOVES = 0;
        /// <summary>unhelpful passes before quitting KL</summary>
        public const int KL_NTRIES_BAD = 0; /* #  */
        /// <summary># levels between KL calls in uncoarsening</summary>
        public const int COARSE_NLEVEL_KL = 0;

        /// <summary>use greedy strategy to improve mapping?</summary>
        public const bool REFINE_MAP = false;
        /// <summary>use matching to reduce vertex separator?</summary>
        public const bool VERTEX_COVER = false;
        /// <summary>force subdomain connectivity at end?</summary>
        public const bool CONNECTED_DOMAINS = false;
        /// <summary>greedily increase internal vtxs?</summary>
        public const bool INTERNAL_VERTICES = false;
        /// <summary>number of calls to pairwise_refine to make</summary>
        public const int REFINE_PARTITION = 0;

        /// <summary>number of local opts to find global min</summary>
        public const int OPT3D_NTRIES = 0;

        /// <summary>how to map from eigenvectors to partition</summary>
        public const MappingType MAPPING_TYPE = MappingType.MinCost;
        /// <summary>connect graph for spectral method?</summary>
        public const bool MAKE_CONNECTED = false;

        /// <summary>use vertex weights while coarsening?</summary>
        public const bool COARSEN_VWGTS = false;
        /// <summary>use edge weights while coarsening?</summary>
        public const bool COARSEN_EWGTS = false;
        /// <summary>choose heavy edges in matching?</summary>
        public const bool HEAVY_MATCH = false;
        /// <summary>min vtx reduction for coarsening</summary>
        public const double COARSEN_RATIO_MIN = 0.0d;
        /// <summary>maximum perturbation</summary>
        public const double PERTURB_MAX = 0.0d;
        /// <summary>&lt; Most cuts allowed at one time</summary>
        public const int MAXDIMS = 3;
        /// <summary>&lt; 2^MAXDIMS</summary>
        public const int MAXSETS = 8;
        /// <summary>number of edges to perturb</summary>
        public const int NPERTURB = 0;
        /// <summary>do RQI this often in uncoarsenin</summary>
        public const int COARSE_NLEVEL_RQI = 0; /* g */

        /// <summary>invoke terminal propagation?</summary>
        public const bool TERM_PROP = false;

        /// <summary>simulate the communication?</summary>
        public const int SIMULATOR = 0;
        /// <summary>simulator iterations</summary>
        public const int SIMULATION_ITNS = 0;
        /// <summary>..if so, importance of cuts/hops</summary>
        public const double CUT_TO_HOP_COST = 0.0d;
        /// <summary>cost of each cut</summary>
        public const double CUT_COST = 0.0d;
        /// <summary>cost of each hop</summary>
        public const double HOP_COST = 0.0d;
        /// <summary>cost of each boundary vertex</summary>
        public const double BDY_COST = 0.0d;
        /// <summary>cost of each boundary vertex hop</summary>
        public const double BDY_HOP_COST = 0.0d;
        /// <summary>initiation cost of each message</summary>
        public const double STARTUP_COST = 0.0d;

        /// <summary>user type</summary>
        public const bool EXPERT = false;
        /// <summary>Debug flag for connected components</summary>
        public const bool DEBUG_CONNECTED = false;
        /// <summary>debug code about force_internal?</summary>
        public const bool DEBUG_INTERNAL = false;
        /// <summary>debug code about refine_part?</summary>
        public const bool DEBUG_REFINE_PART = false;
        /// <summary>debug code about refine_map?</summary>
        public const bool DEBUG_REFINE_MAP = false;
        /// <summary>print out computed machine params?</summary>
        public const bool DEBUG_MACH_PARAMS = false;
        /// <summary>turn on debugging in assignment</summary>
        public const bool DEBUG_ASSIGN = false;
        /// <summary>Debug flag for inertial method</summary>
        public const bool DEBUG_INERTIAL = false;
        /// <summary>trace the execution of the code</summary>
        public const bool DEBUG_TRACE = true;
        /// <summary>Debug output in min vtx cover routines?</summary>
        public const int DEBUG_COVER = 0;
        /// <summary>debug flag for eigenvector generation</summary>
        public const int DEBUG_EVECS = 2;
        /// <summary>print warning messages? 0 to 5, increasing amount of info</summary>
        public const int WARNING_EVECS = 5;
        /// <summary>debug flag for optimization</summary>
        public const bool DEBUG_OPTIMIZE = false;
        /// <summary>turn on debugging for bipartite matching</summary>
        public const DebugFlagBP DEBUG_BPMATCH = DebugFlagBP.ErrorChecking;
        /// <summary>debug flag for matrix perturbation</summary>
        public const bool DEBUG_PERTURB = false;
        /// <summary>debug flag for coarsening</summary>
        public const bool DEBUG_COARSEN = false;
        /// <summary>OUTPUT METRICS</summary>
        public const bool PRINT_GRAPH_PARTITION_METRICS = false;
        /// <summary>OUTPUT METRICS negative</summary>
        public const bool PRINT_GRAPH_PARTITION_METRICS_ALL_RECURSION_LEVELS = false;
        /// <summary>OUTPUT METRICS false = 1, true = 2</summary>
        public const bool PRINT_GRAPH_PARTITION_METRICS_DETAILED = false;
        /// <summary>print section headings for output?</summary>
        public const bool PRINT_HEADERS = true;
        /// <summary>ECHO: case != 0</summary>
        public const bool ECHO_USER_PARAMETERS = false;
        /// <summary>ECHO: case 2 or -2</summary>
        public const bool ECHO_INPUT_PARAMETERS = false;
        /// <summary>should I check input for correctness?</summary>
        public const bool CHECK_INPUT = true;

        /// <summary>Debug flag for Kernighan-Lin</summary>
        public const DebugFlagKL DEBUG_KL = DebugFlagKL.ImprovementsPerStep;

        /// <summary>perform detailed timing on Lanczos_SO?</summary>
        public const bool LANCZOS_TIME = false;
        /// <summary>benchmark some numerical kernels?</summary>
        public const bool TIME_KERNELS = false;
        // TODO:
        /// <summary></summary>
        public const int OUTPUT_TIME = 0;

        // TODO:
        /// <summary></summary>
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

        /// <summary>Machine Precision.</summary>
        /// <remarks>
        /// The value of the Epsilon property is not equivalent to machine epsilon,
        /// which represents the upper bound of the relative error due to rounding in floating-point arithmetic.
        /// </remarks>
        /// <see ref="https://docs.microsoft.com/en-us/dotnet/api/system.double.epsilon?view=netcore-3.0"/>
        public static readonly double DOUBLE_EPSILON = CalculateMachineEpsilon();

        /// <summary>
        /// Calculates the value of Machine Epsilon at the time of execution.
        /// </summary>
        /// <returns>The value of Machine Epsilon value; The upper bound of the relative error due to rounding in floating-point arithmetic.</returns>
        /// <see ref="https://stackoverflow.com/a/9393079/103959"/>
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
#pragma warning restore RCS1057 // Add empty line between declarations.
