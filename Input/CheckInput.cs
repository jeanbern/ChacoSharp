#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
using System;
using System.Diagnostics;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Graph.CheckGraph;
// ReSharper disable HeuristicUnreachableCode
#pragma warning disable 162

namespace ChacoSharp.Input
{
    public static unsafe class CheckInput
    {
        /* Check graph and input options and parameters. */
        public static bool check_input(vtx_data** graph, /* linked lists of vertex data */
            int nvtxs, /* number of vertices */
            int nedges, /* number of edges */
            int igeom, /* geometric dimension for inertial method */
            float** coords, /* coordinates for inertial method */
            int* assignment, /* set numbers if read-from-file */
            double[] goal, /* desired sizes of different sets */
            int architecture, /* 0=> hypercube, d=> d-dimensional mesh */
            int ndims_tot, /* number of hypercube dimensions */
            int[] mesh_dims /*[3]*/, /* size of mesh in each dimension */
            PartitioningStrategy partitioningStrategy, /* global partitioning algorithm */
            LocalPartitioningStrategy localParitioningStrategy, /* local partitioning algorithm */
            bool rqi_flag, /* flag for RQI/symmlq eigensolver */
            int* vmax, /* smallest acceptable coarsened nvtxs */
            int ndims, /* partitioning level */
            double eigtol /* tolerance for eigen-pairs */
        )
        {
            if (architecture > 0 && mesh_dims == null)
            {
                throw new ArgumentNullException(nameof(mesh_dims));
            }

            if (DEBUG_TRACE)
            {
                Trace.WriteLine("<Entering check_input>");
            }

            /* First check for consistency in the graph. */
            bool graphError; /* does graph check out OK? */
            if (graph != null)
            {
                graphError = check_graph(graph, nvtxs, nedges);
                if (graphError)
                {
                    Trace.WriteLine("ERRORS in graph.");
                }
                else
                {
                    Trace.WriteLine("Graph check OK");
                }
            }
            else
            {
                /* Only allowed if simple or inertial w/o KL and no weights. */
                graphError = false;
                if (partitioningStrategy == PartitioningStrategy.Multilevel_KL || partitioningStrategy == PartitioningStrategy.Spectral || localParitioningStrategy == LocalPartitioningStrategy.KernighanLin)
                {
                    Trace.WriteLine("No graph input.  Only allowed for inertial or simple methods without KL.");
                    graphError = true;
                }
            }

            /* Now check the input values. */
            var errorFound = false;
            var assignmentError = false;

            if (architecture < 0 || architecture > 3)
            {
                Trace.WriteLine($"Machine architecture parameter = {architecture:d}, must be in [0,3].");
                errorFound = true;
            }
            else if (architecture == 0)
            {
                if (ndims_tot < 0)
                {
                    Trace.WriteLine($"Dimension of hypercube = {ndims_tot:d}, must be at least 1.");
                    errorFound = true;
                }
            }
            else if (architecture > 0)
            {
                if (architecture == 1 && mesh_dims[0] <= 0)
                {
                    Trace.WriteLine($"Size of 1-D mesh improperly specified, {mesh_dims[0]:d}.");
                    errorFound = true;
                }

                if (architecture == 2 && (mesh_dims[0] <= 0 || mesh_dims[1] <= 0))
                {
                    Trace.WriteLine($"Size of 2-D mesh improperly specified, {mesh_dims[0]:d}x{mesh_dims[1]:d}.");
                    errorFound = true;
                }

                if (architecture == 2 && (mesh_dims[0] <= 0 || mesh_dims[1] <= 0 || mesh_dims[2] <= 0))
                {
                    Trace.WriteLine($"Size of 3-D mesh improperly specified, {mesh_dims[0]:d}x{mesh_dims[1]:d}x{mesh_dims[2]:d}.");
                    errorFound = true;
                }
            }

            if (ndims < 1 || ndims > MAXDIMS)
            {
                Trace.WriteLine($"Partitioning at each step = {ndims:d}, should be in [1,{MAXDIMS:d}].");
                errorFound = true;
            }

            var nprocs = 0; /* number of processors partitioning for */
            if (architecture == 0)
            {
                if (!errorFound)
                {
                    nprocs = 1 << ndims_tot;
                }
            }
            else if (architecture > 0)
            {
                nprocs = mesh_dims[0] * mesh_dims[1] * mesh_dims[2];
            }

            if (1 << ndims > nprocs)
            {
                Trace.WriteLine($"Partitioning step {ndims:d} too large for {nprocs:d} processors.");
                errorFound = true;
            }

            if ((int) partitioningStrategy < 1 || (int) partitioningStrategy > 7)
            {
                Trace.WriteLine($"Global partitioning method = {(int) partitioningStrategy:d}, must be in [1,7].");
                errorFound = true;
            }

            if ((int) localParitioningStrategy < 1 || (int) localParitioningStrategy > 2)
            {
                Trace.WriteLine($"Local partitioning method = {(int) localParitioningStrategy:d}, must be in [1,2].");
                errorFound = true;
            }

            if (partitioningStrategy == PartitioningStrategy.Multilevel_KL || (partitioningStrategy == PartitioningStrategy.Spectral && rqi_flag))
            {
                var i = 2 * (1 << ndims);
                if (*vmax < i)
                {
                    Trace.WriteLine($"WARNING: Number of vertices in coarse graph ({*vmax:d}) being reset to {i:d}.");
                    *vmax = i;
                }
            }

            if ((partitioningStrategy == PartitioningStrategy.Multilevel_KL || partitioningStrategy == PartitioningStrategy.Spectral) && eigtol <= 0)
            {
                Trace.WriteLine($"Eigen tolerance ({eigtol:g}) must be positive value");
                errorFound = true;
            }

            if (partitioningStrategy == PartitioningStrategy.Inertial ||
                (MATCH_TYPE == MatchingRoutine.maxmatch5_geometric && (partitioningStrategy == PartitioningStrategy.Multilevel_KL || (partitioningStrategy == PartitioningStrategy.Spectral && rqi_flag))))
            {
                if (igeom < 1 || igeom > 3)
                {
                    Trace.WriteLine("Geometry must be 1-, 2- or 3-dimensional");
                    errorFound = true;
                }

                if (igeom > 0 && coords == null)
                {
                    Trace.WriteLine("No coordinates given");
                    errorFound = false;
                }
                else if (igeom > 0 && coords[0] == null)
                {
                    Trace.WriteLine("No X-coordinates given");
                    errorFound = true;
                }
                else if (igeom > 1 && coords[1] == null)
                {
                    Trace.WriteLine("No Y-coordinates given");
                    errorFound = true;
                }
                else if (igeom > 2 && coords[2] == null)
                {
                    Trace.WriteLine("No Z-coordinates given");
                    errorFound = true;
                }
            }

            if (partitioningStrategy == PartitioningStrategy.ReadFromFile && localParitioningStrategy == LocalPartitioningStrategy.KernighanLin)
            {
                if (nprocs > 1 << ndims)
                {
                    Trace.WriteLine("Can only use local method on single level of read-in assignment,");
                    Trace.WriteLine($"  but ndims =  {ndims:d}, while number of processors = {nprocs:d}.");
                    errorFound = true;
                }
            }

            /* Now check for consistency in the goal array. */
            double vertexWeightSum;
            if (graph != null)
            {
                vertexWeightSum = 0;
                for (var i = 1; i <= nvtxs; i++)
                {
                    vertexWeightSum += graph[i]->vwgt;
                }
            }
            else
            {
                vertexWeightSum = nvtxs;
            }

            double vertexGoalSum = 0;
            for (var i = 0; i < nprocs; i++)
            {
                if (goal[i] < 0)
                {
                    Trace.WriteLine($"goal[{i:d}] is {goal[i]:g}, but should be nonnegative.");
                    errorFound = true;
                }

                vertexGoalSum += goal[i];
            }

            if (Math.Abs(vertexWeightSum - vertexGoalSum) > 1e-5 * (vertexWeightSum + vertexGoalSum))
            {
                Trace.WriteLine($"Sum of values in goal ({vertexGoalSum:g}) not equal to sum of vertex weights ({vertexWeightSum:g}).");
                errorFound = true;
            }

            /* Check assignment file if read in. */
            if (partitioningStrategy == PartitioningStrategy.ReadFromFile && !errorFound)
            {
                assignmentError = check_assignment(assignment, nvtxs, nprocs, ndims, localParitioningStrategy);
            }

            /* Add some checks for model parameters */

            /* Finally, check the parameters. */
            var parameterError = check_params(partitioningStrategy, localParitioningStrategy, rqi_flag, ndims);

            errorFound = errorFound || graphError || assignmentError || parameterError;

            return errorFound;
        }

        static bool check_params(PartitioningStrategy partitioningStrategy, /* global partitioning algorithm */
            LocalPartitioningStrategy localPartitioningAlgorithm, /* local partitioning algorithm */
            bool rqi_flag, /* use multilevel eigensolver? */
            int ndims /* number of eigenvectors */
        )
        {
            var parameterErrorDetected = false;

            if (OUTPUT_TIME < 0 || OUTPUT_TIME > 2)
            {
                Trace.WriteLine($"WARNING: OUTPUT_TIME ({OUTPUT_TIME:d}) should be in [0,2].");
            }

            if (partitioningStrategy == PartitioningStrategy.Multilevel_KL || partitioningStrategy == PartitioningStrategy.Spectral)
            {
                if (EXPERT)
                {
                    if ((int) LANCZOS_TYPE < 1 || (int) LANCZOS_TYPE > 4)
                    {
                        Trace.WriteLine($"LANCZOS_TYPE ({(int) LANCZOS_TYPE:d}) should be in [1,4].");
                        parameterErrorDetected = true;
                    }
                }
                else
                {
                    if ((int) LANCZOS_TYPE < 1 || (int) LANCZOS_TYPE > 3)
                    {
                        Trace.WriteLine($"LANCZOS_TYPE ({(int) LANCZOS_TYPE:d}) should be in [1,3].");
                        parameterErrorDetected = true;
                    }
                }

                if (EIGEN_TOLERANCE <= 0)
                {
                    Trace.WriteLine($"EIGEN_TOLERANCE ({EIGEN_TOLERANCE:g}) should be positive.");
                    parameterErrorDetected = true;
                }

                if (LANCZOS_SO_INTERVAL <= 0)
                {
                    Trace.WriteLine($"LANCZOS_SO_INTERVAL ({LANCZOS_SO_INTERVAL:d}) should be positive.");
                    parameterErrorDetected = true;
                }

                if (LANCZOS_SO_INTERVAL == 1)
                {
                    Trace.WriteLine("WARNING: More efficient if LANCZOS_SO_INTERVAL = 2, not 1.");
                }

                if (BISECTION_SAFETY <= 0)
                {
                    Trace.WriteLine($"BISECTION_SAFETY ({BISECTION_SAFETY:g}) should be positive.");
                    parameterErrorDetected = true;
                }

                if (LANCZOS_CONVERGENCE_MODE < 0 || LANCZOS_CONVERGENCE_MODE > 1)
                {
                    Trace.WriteLine($"LANCZOS_CONVERGENCE_MODE ({LANCZOS_CONVERGENCE_MODE:d}) should be in [0,1].");
                    parameterErrorDetected = true;
                }

                if (WARNING_ORTHTOL <= 0.0d)
                {
                    Trace.WriteLine($"WARNING: WARNING_ORTHTOL ({WARNING_ORTHTOL:g}) should be positive.");
                }

                if (WARNING_MISTOL <= 0.0d)
                {
                    Trace.WriteLine($"WARNING: WARNING_MISTOL ({WARNING_MISTOL:g}) should be positive.");
                }

                if (LANCZOS_SO_PRECISION < 1 || LANCZOS_SO_PRECISION > 2)
                {
                    Trace.WriteLine($"LANCZOS_SO_PRECISION ({LANCZOS_SO_PRECISION:d}) should be in [1,2].");
                    parameterErrorDetected = true;
                }

                if (PERTURB)
                {
                    if (NPERTURB < 0)
                    {
                        Trace.WriteLine($"NPERTURB ({NPERTURB:d}) should be nonnegative.");
                        parameterErrorDetected = true;
                    }

                    if (NPERTURB > 0 && PERTURB_MAX < 0)
                    {
                        Trace.WriteLine($"PERTURB_MAX ({PERTURB_MAX:g}) should be nonnegative.");
                        parameterErrorDetected = true;
                    }
                }

                if ((int) MAPPING_TYPE < 0 || (int) MAPPING_TYPE > 3)
                {
                    Trace.WriteLine($"MAPPING_TYPE ({(int) MAPPING_TYPE:d}) should be in [0,3].");
                    parameterErrorDetected = true;
                }

                if (ndims == 3 && OPT3D_NTRIES <= 0)
                {
                    Trace.WriteLine($"OPT3D_NTRIES ({OPT3D_NTRIES:d}) should be positive.");
                    parameterErrorDetected = true;
                }

                if (partitioningStrategy == PartitioningStrategy.Spectral && rqi_flag)
                {
                    if (COARSE_NLEVEL_RQI <= 0)
                    {
                        Trace.WriteLine($"COARSE_NLEVEL_RQI ({COARSE_NLEVEL_RQI:d}) should be positive.");
                        parameterErrorDetected = true;
                    }

                    if (RQI_CONVERGENCE_MODE < 0 || RQI_CONVERGENCE_MODE > 1)
                    {
                        Trace.WriteLine($"RQI_CONVERGENCE_MODE ({RQI_CONVERGENCE_MODE:d}) should be in [0,1].");
                        parameterErrorDetected = true;
                    }

                    if (TERM_PROP)
                    {
                        Trace.WriteLine("WARNING: Using default Lanczos for extended eigenproblem, not RQI/Symmlq.");
                    }
                }

                if (partitioningStrategy == PartitioningStrategy.Multilevel_KL && COARSE_NLEVEL_KL <= 0)
                {
                    Trace.WriteLine($"COARSE_NLEVEL_KL ({COARSE_NLEVEL_KL:d}) should be positive.");
                    parameterErrorDetected = true;
                }

                if (partitioningStrategy == PartitioningStrategy.Multilevel_KL || partitioningStrategy == PartitioningStrategy.Spectral && rqi_flag)
                {
                    if (COARSEN_RATIO_MIN < .5)
                    {
                        Trace.WriteLine($"COARSEN_RATIO_MIN ({COARSEN_RATIO_MIN:g}) should be at least 1/2.");
                        parameterErrorDetected = true;
                    }

                    if ((int)MATCH_TYPE < 1 || (int)MATCH_TYPE > 9)
                    {
                        Trace.WriteLine($"MATCH_TYPE ({MATCH_TYPE:d}) should be in [1,9].");
                        parameterErrorDetected = true;
                    }
                }
            }

            if (partitioningStrategy == PartitioningStrategy.Multilevel_KL || localPartitioningAlgorithm == LocalPartitioningStrategy.KernighanLin)
            {
                if ((int)KL_METRIC < 1 || (int)KL_METRIC > 2)
                {
                    Trace.WriteLine($"KL_METRIC ({KL_METRIC:d}) should be in [1,2].");
                    parameterErrorDetected = true;
                }

                if (KL_BAD_MOVES < 0)
                {
                    Trace.WriteLine($"KL_BAD_MOVES ({KL_BAD_MOVES:d}) should be non-negative.");
                    parameterErrorDetected = true;
                }

                if (KL_NTRIES_BAD < 0)
                {
                    Trace.WriteLine($"KL_NTRIES_BAD ({KL_NTRIES_BAD:d}) should be non-negative.");
                    parameterErrorDetected = true;
                }

                if (KL_IMBALANCE < 0.0 || KL_IMBALANCE > 1.0)
                {
                    Trace.WriteLine($"KL_IMBALANCE ({KL_IMBALANCE:g}) should be in [0,1].");
                    parameterErrorDetected = true;
                }
            }

            if (SIMULATOR < 0 || SIMULATOR > 3)
            {
                Trace.WriteLine($"SIMULATOR ({SIMULATOR:d}) should be in [0,3].");
                parameterErrorDetected = true;
            }

            // ReSharper disable once InvertIf
            if (TERM_PROP)
            {
                if (CUT_TO_HOP_COST <= 0)
                {
                    Trace.WriteLine($"CUT_TO_HOP_COST ({CUT_TO_HOP_COST:g}) should be positive.");
                    parameterErrorDetected = true;
                }

                if (ndims > 1)
                {
                    Trace.WriteLine("WARNING: May ignore terminal propagation in spectral quadri/octa section");
                }
            }

            return parameterErrorDetected;
        }

        static bool check_assignment(int* assignment, /* set numbers if read-from-file */
            int nvtxs, /* number of vertices */
            int nsets_tot, /* total number of desired sets */
            int ndims, /* partitioning level */
            LocalPartitioningStrategy localPartitioningStrategy /* local partitioning algorithm */
        )
        {
            var nsets = 1 << ndims;
            var flag = false;

            for (var i = 1; i <= nvtxs && !flag; i++)
            {
                if (assignment[i] < 0)
                {
                    Trace.WriteLine($"Assignment[{i:d}] = {assignment[i]:d} less than zero.");
                    flag = true;
                }
                else if (assignment[i] >= nsets_tot)
                {
                    Trace.WriteLine($"Assignment[{i:d}] = {assignment[i]:d}, too large for {nsets_tot:d} sets.");
                    flag = true;
                }
                else if (localPartitioningStrategy == LocalPartitioningStrategy.KernighanLin && assignment[i] >= nsets)
                {
                    Trace.WriteLine("Can only use local method on single level of read-in assignment,");
                    Trace.WriteLine($"  but assignment[{i:d}] =  {assignment[i]:d}.");
                    flag = true;
                }
            }

            return flag;
        }
    }
}
