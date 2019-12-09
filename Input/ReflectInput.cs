#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
using System;
using System.Diagnostics;
using static ChacoSharp.StaticConstants;

namespace ChacoSharp.Input
{
    public static class ReflectInput
    {
        /* Print out the input options. */
        public static void reflect_input(int nvtxs, /* number of vertices in graph */
            int nedges, /* number of edges in graph */
            int igeom, /* geometric dimension for inertial method */
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
                Trace.WriteLine("<Entering reflect_input>");
            }

            Trace.WriteLine("");

            if (PRINT_HEADERS)
            {
                Trace.WriteLine("\n           Input and Parameter Values\n");
            }

            Trace.WriteLine($"# vertices = {nvtxs:d}, # edges = {nedges:d}");

            /* Print global partitioning strategy. */

            Trace.WriteLine("Global method: ");
            if (global_method == PartitioningStrategy.Multilevel_KL)
            {
                Trace.WriteLine("Multilevel-KL");
            }
            else if (global_method == PartitioningStrategy.Spectral)
            {
                Trace.WriteLine("Spectral");
            }
            else if (global_method == PartitioningStrategy.Inertial)
            {
                Trace.WriteLine("Inertial");
            }
            else if (global_method == PartitioningStrategy.Linear)
            {
                Trace.WriteLine("Linear");
            }
            else if (global_method == PartitioningStrategy.Random)
            {
                Trace.WriteLine("Random");
            }
            else if (global_method == PartitioningStrategy.Scattered)
            {
                Trace.WriteLine("Scattered");
            }
            else if (global_method == PartitioningStrategy.ReadFromFile)
            {
                Trace.WriteLine("Read-From-File");
            }

            if (global_method == PartitioningStrategy.Multilevel_KL)
            {
                Trace.WriteLine($"Number of vertices to coarsen down to: {vmax:d}");
                Trace.WriteLine($"Eigen tolerance: {eigtol:g}");
            }

            else if (global_method == PartitioningStrategy.Spectral)
            {
                if (rqi_flag)
                {
                    Trace.WriteLine("Multilevel RQI/Symmlq eigensolver");
                    Trace.WriteLine($"Number of vertices to coarsen down to: {vmax:d}");
                    Trace.WriteLine($"Eigen tolerance: {eigtol:g}");
                }
            }

            else if (global_method == PartitioningStrategy.Inertial)
            {
                if (geomname != null)
                {
                    Trace.WriteLine($"Geometry input file: `{geomname}', Dimensionality = {igeom:d}");
                }
            }

            else if (global_method == PartitioningStrategy.ReadFromFile)
            {
                Trace.WriteLine(string.Format("Assignment input file: `{0}'", inassignname));
            }

            /* Now describe local method. */
            if (localPartitioningStrategy == LocalPartitioningStrategy.KernighanLin)
            {
                Trace.WriteLine("Local method: Kernighan-Lin");
            }
            else if (localPartitioningStrategy == LocalPartitioningStrategy.None)
            {
                Trace.WriteLine("Local method: None");
            }

            /* Now describe target architecture. */
            if (architecture == 0)
            {
                Trace.WriteLine($"Partitioning target: {ndims_tot:d}-dimensional hypercube");
            }
            else if (architecture > 0)
            {
                Trace.Write($"Partitioning target: {architecture:d}-dimensional mesh of size ");
                if (architecture == 1)
                {
                    Trace.WriteLine($"{mesh_dims[0]:d}");
                }
                else if (architecture == 2)
                {
                    Trace.WriteLine($"{mesh_dims[0]:d}x{mesh_dims[1]:d}");
                }
                else if (architecture == 3)
                {
                    Trace.WriteLine($"{mesh_dims[0]:d}x{mesh_dims[1]:d}x{mesh_dims[2]:d}");
                }
            }

            if (ndims == 1)
            {
                Trace.WriteLine("Partitioning mode: Bisection");
            }
            else if (ndims == 2)
            {
                Trace.WriteLine("Partitioning mode: Quadrisection");
            }
            else if (ndims == 3)
            {
                Trace.WriteLine("Partitioning mode: Octasection");
            }

            /* Add stuff about communications simulator. */

            Trace.WriteLine($"Random seed: {seed:d}");

            if (ECHO_USER_PARAMETERS)
            {
                reflect_params(global_method, localPartitioningStrategy, rqi_flag, ndims);
            }

            Trace.WriteLine("\n");
        }

        public static void reflect_params(PartitioningStrategy global_method, /* global partitioning algorithm */
            LocalPartitioningStrategy localPartitioningStrategy, /* local partitioning algorithm */
            bool rqi_flag, /* use RQI/SYMMLQ eigensolver? */
            int ndims) /* number of eigenvectors to generate */
        {
            Trace.WriteLine("Active Parameters:");

            Trace.WriteLine(string.Format("  CHECK_INPUT = {0}", CHECK_INPUT ? "True" : "False"));

            if (global_method == PartitioningStrategy.Multilevel_KL || global_method == PartitioningStrategy.Linear)
            {
                Trace.WriteLine("  LANCZOS_TYPE:  ");
                if (LANCZOS_TYPE == LanczosType.FullOrthogonalization)
                {
                    Trace.Write(" Full orthogonalization");
                }
                else if (LANCZOS_TYPE == LanczosType.FullOrthogonalizationInverseOperator)
                {
                    Trace.Write("Full orthogonalization, inverse operator");
                }
                else if (LANCZOS_TYPE == LanczosType.SelectiveOrthogonalization)
                {
                    Trace.Write("Selective orthogonalization");
                }
                else if (LANCZOS_TYPE == LanczosType.SelectiveOrthogonalizationDoubleEnded)
                {
                    if (EXPERT)
                    {
                        Trace.Write("Selective orthogonalization against both ends");
                    }
                    else
                    {
                        /* Check input should catch this, but just in case ... */
                        //LANCZOS_TYPE = LanczosType.SelectiveOrthogonalization;
                        throw new InvalidOperationException("Only use " + nameof(LanczosType.SelectiveOrthogonalizationDoubleEnded) + " in " + nameof(EXPERT) + " mode");
                        Trace.Write("Selective orthogonalization");
                    }
                }

                Trace.WriteLine("");

                Trace.WriteLine($"  EIGEN_TOLERANCE = {EIGEN_TOLERANCE:g}");

                if (SRESTOL > 0)
                {
                    Trace.WriteLine($"  SRESTOL = {SRESTOL:g}");
                }
                else
                {
                    Trace.WriteLine($"  SRESTOL = {SRESTOL:g} ... autoset to square of eigen tolerance");
                }

                if (LANCZOS_MAXITNS > 0)
                {
                    Trace.WriteLine($"  LANCZOS_MAXITNS = {LANCZOS_MAXITNS:d}");
                }
                else
                {
                    Trace.WriteLine($"  LANCZOS_MAXITNS = {LANCZOS_MAXITNS:d} ... autoset to twice # vertices");
                }

                if (LANCZOS_SO_PRECISION == 1)
                {
                    Trace.WriteLine("  LANCZOS_SO_PRECISION = 1 ... single precision");
                }
                else
                {
                    Trace.WriteLine("  LANCZOS_SO_PRECISION = 2 ... double precision");
                }

                Trace.WriteLine($"  LANCZOS_SO_INTERVAL = {LANCZOS_SO_INTERVAL:d}");

                if (LANCZOS_CONVERGENCE_MODE == 1)
                {
                    Trace.WriteLine("  LANCZOS_CONVERGENCE_MODE = 1 ... partition tolerance");
                }
                else
                {
                    Trace.WriteLine("  LANCZOS_CONVERGENCE_MODE = 0 ... residual tolerance");
                }

                Trace.WriteLine($"  BISECTION_SAFETY = {BISECTION_SAFETY:g}");
                if (LANCZOS_TYPE == LanczosType.SelectiveOrthogonalization || LANCZOS_TYPE == LanczosType.SelectiveOrthogonalizationDoubleEnded)
                {
                    if (!LANCZOS_TIME)
                    {
                        Trace.WriteLine("  LANCZOS_TIME = 0 ... no detailed timing");
                    }
                    else
                    {
                        Trace.WriteLine("  LANCZOS_TIME = 1 ... detailed timing");
                    }
                }

                if (WARNING_EVECS > 0)
                {
                    Trace.WriteLine($"  WARNING_EVECS = {WARNING_EVECS:d}");
                }

                if (MAPPING_TYPE == MappingType.CutAtOrigin)
                {
                    Trace.WriteLine("  MAPPING_TYPE = 0 ... cut at origin");
                }
                else if (MAPPING_TYPE == MappingType.MinCost)
                {
                    Trace.WriteLine("  MAPPING_TYPE = 1 ... min-cost assignment");
                }
                else if (MAPPING_TYPE == MappingType.RecursiveMedian)
                {
                    Trace.WriteLine("  MAPPING_TYPE = 2 ... recursive median");
                }
                else if (MAPPING_TYPE == MappingType.IndependantMedians)
                {
                    Trace.WriteLine("  MAPPING_TYPE = 3 ... independent medians");
                }

                Trace.WriteLine(string.Format("  MAKE_CONNECTED = {0}", MAKE_CONNECTED ? "True" : "False"));
                Trace.WriteLine(string.Format("  PERTURB = {0}", PERTURB ? "True" : "False"));
                if (PERTURB)
                {
                    Trace.WriteLine($"    NPERTURB = {NPERTURB:d}");
                    Trace.WriteLine($"    PERTURB_MAX = {PERTURB_MAX:g}");
                }

                if (ndims == 3)
                {
                    Trace.WriteLine($"  OPT3D_NTRIES = {OPT3D_NTRIES:d}");
                }
            }

            if (global_method == PartitioningStrategy.Multilevel_KL)
            {
                Trace.WriteLine($"  COARSEN_RATIO_MIN = {COARSEN_RATIO_MIN:g}");
                Trace.WriteLine($"  COARSE_NLEVEL_KL = {COARSE_NLEVEL_KL:d}");
                Trace.WriteLine($"  MATCH_TYPE = {MATCH_TYPE:d}");
                Trace.WriteLine(string.Format("  HEAVY_MATCH = {0}", HEAVY_MATCH ? "True" : "False"));
                Trace.WriteLine(string.Format("  COARSE_KL_BOTTOM = {0}", COARSE_KL_BOTTOM ? "True" : "False"));
                Trace.WriteLine(string.Format("  COARSEN_VWGTS = {0}", COARSEN_VWGTS ? "True" : "False"));
                Trace.WriteLine(string.Format("  COARSEN_EWGTS = {0}", COARSEN_EWGTS ? "True" : "False"));
                Trace.WriteLine(string.Format("  KL_ONLY_BNDY = {0}", KL_ONLY_BNDY ? "True" : "False"));
            }

            if (global_method == PartitioningStrategy.Spectral && rqi_flag)
            {
                Trace.WriteLine($"  COARSE_NLEVEL_RQI = {COARSE_NLEVEL_RQI:d)}");
                if (RQI_CONVERGENCE_MODE == 1)
                {
                    Trace.WriteLine("  RQI_CONVERGENCE_MODE = 1 ... partition tolerance");
                }
                else
                {
                    Trace.WriteLine("  RQI_CONVERGENCE_MODE = 0 ... residual tolerance");
                }

                Trace.WriteLine($"  COARSEN_RATIO_MIN = {COARSEN_RATIO_MIN:g}");
                Trace.WriteLine(string.Format("  COARSEN_VWGTS = {0}", COARSEN_VWGTS ? "True" : "False"));
                Trace.WriteLine(string.Format("  COARSEN_EWGTS = {0}", COARSEN_EWGTS ? "True" : "False"));
            }

            if (global_method == PartitioningStrategy.Multilevel_KL || localPartitioningStrategy == LocalPartitioningStrategy.KernighanLin)
            {
                Trace.WriteLine(string.Format("  KL_RANDOM = {0}", KL_RANDOM ? "True" : "False"));
                if (KL_METRIC == KernighanLinMetric.Cuts)
                {
                    Trace.WriteLine("  KL_METRIC = Cuts");
                }
                else if (KL_METRIC == KernighanLinMetric.Hops)
                {
                    Trace.WriteLine("  KL_METRIC = Hops");
                }

                Trace.WriteLine(string.Format("  KL_NTRIES_BAD = {0:d}", KL_NTRIES_BAD));
                Trace.WriteLine(string.Format("  KL_BAD_MOVES = {0:d}", KL_BAD_MOVES));
                Trace.WriteLine(string.Format("  KL_UNDO_LIST = {0}", KL_UNDO_LIST ? "True" : "False"));
                Trace.WriteLine(string.Format("  KL_IMBALANCE = {0:g}", KL_IMBALANCE));
            }

            if (global_method == PartitioningStrategy.Multilevel_KL || global_method == PartitioningStrategy.Spectral || localPartitioningStrategy == LocalPartitioningStrategy.KernighanLin)
            {
                Trace.WriteLine(string.Format("  TERM_PROP = {0}", TERM_PROP ? "True" : "False"));
                if (TERM_PROP)
                {
                    Trace.WriteLine(string.Format("    CUT_TO_HOP_COST = {0:g}", CUT_TO_HOP_COST));
                }
            }

            if (SEQUENCE)
            {
                Trace.WriteLine(string.Format("  SEQUENCE = {0}", SEQUENCE ? "True" : "False"));
            }

            Trace.WriteLine(string.Format("  PRINT_GRAPH_PARTITION_METRICS = {0}", PRINT_GRAPH_PARTITION_METRICS ? "True" : "False"));
            Trace.WriteLine(string.Format("  PRINT_GRAPH_PARTITION_METRICS_ALL_RECURSION_LEVELS = {0}", PRINT_GRAPH_PARTITION_METRICS_ALL_RECURSION_LEVELS ? "True" : "False"));
            Trace.WriteLine(string.Format("  PRINT_GRAPH_PARTITION_METRICS_DETAILED = {0}", PRINT_GRAPH_PARTITION_METRICS_DETAILED ? "True" : "False"));
            Trace.WriteLine(string.Format("  MAKE_VWGTS = {0}", MAKE_VWGTS ? "True" : "False"));
            Trace.WriteLine(string.Format("  REFINE_MAP = {0}", REFINE_MAP ? "True" : "False"));
            Trace.WriteLine(string.Format("  REFINE_PARTITION = {0:d}", REFINE_PARTITION));
            Trace.WriteLine(string.Format("  INTERNAL_VERTICES = {0}", INTERNAL_VERTICES ? "True" : "False"));

            if (SIMULATOR != 0)
            {
                Trace.WriteLine(string.Format("  SIMULATOR = {0:d}", SIMULATOR));
                Trace.WriteLine(string.Format("  SIMULATION_ITNS = {0:d}", SIMULATION_ITNS));
                Trace.WriteLine(string.Format("  CUT_COST = {0:g}", CUT_COST));
                Trace.WriteLine(string.Format("  HOP_COST = {0:g}", HOP_COST));
                Trace.WriteLine(string.Format("  BDY_COST = {0:g}", BDY_COST));
                Trace.WriteLine(string.Format("  BDY_HOP_COST = {0:g}", BDY_HOP_COST));
                Trace.WriteLine(string.Format("  STARTUP_COST = {0:g}", STARTUP_COST));
            }

            /* Now print out all the nonzero debug parameters. */
            if (DEBUG_CONNECTED)
            {
                Trace.WriteLine("  DEBUG_CONNECTED = True");
            }

            if (DEBUG_PERTURB)
            {
                Trace.WriteLine("  DEBUG_PERTURB = True");
            }

            if (DEBUG_ASSIGN)
            {
                Trace.WriteLine("  DEBUG_ASSIGN = True");
            }

            if (DEBUG_INERTIAL)
            {
                Trace.WriteLine("  DEBUG_INERTIAL = True");
            }

            if (DEBUG_OPTIMIZE)
            {
                Trace.WriteLine("  DEBUG_OPTIMIZE = True");
            }

            if (DEBUG_BPMATCH != DebugFlagBP.NoDebugging)
            {
                Trace.WriteLine("  DEBUG_BPMATCH = " + DEBUG_BPMATCH);
            }

            if (DEBUG_COARSEN)
            {
                Trace.WriteLine("  DEBUG_COARSEN = True");
            }

            if (DEBUG_EVECS != 0)
            {
                Trace.WriteLine(string.Format("  DEBUG_EVECS = {0:d}", DEBUG_EVECS));
            }

            if (DEBUG_KL != DebugFlagKL.NoDebugging)
            {
                Trace.WriteLine("  DEBUG_KL = " + DEBUG_KL);
            }

            if (DEBUG_INTERNAL)
            {
                Trace.WriteLine("  DEBUG_INTERNAL = True");
            }

            if (DEBUG_REFINE_PART)
            {
                Trace.WriteLine("  DEBUG_REFINE_PART = True");
            }

            if (DEBUG_REFINE_MAP)
            {
                Trace.WriteLine("  DEBUG_REFINE_MAP = True");
            }

            if (DEBUG_TRACE)
            {
                Trace.WriteLine(string.Format("  DEBUG_TRACE = {0:d}", DEBUG_TRACE));
            }

            if (DEBUG_MACH_PARAMS)
            {
                Trace.WriteLine("  DEBUG_MACH_PARAMS = True");
            }
        }
    }
}
