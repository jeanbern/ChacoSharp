#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
using System;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Graph.CheckGraph;

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
            string graphname, /* graph input file name */
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

            if (DEBUG_TRACE || FullTrace)
            {
                Console.WriteLine("<Entering check_input>");
            }

            /* First check for consistency in the graph. */
            bool graphError; /* does graph check out OK? */
            if (graph != null)
            {
                graphError = check_graph(graph, nvtxs, nedges);
                if (graphError)
                {
                    Console.WriteLine("ERRORS in graph.");
                }
                else
                {
                    Console.WriteLine("Graph check OK");
                }
            }
            else
            {
                /* Only allowed if simple or inertial w/o KL and no weights. */
                graphError = false;
                if (partitioningStrategy == PartitioningStrategy.Multilevel_KL || partitioningStrategy == PartitioningStrategy.Spectral || localParitioningStrategy == LocalPartitioningStrategy.KernighanLin)
                {
                    Console.WriteLine("No graph input.  Only allowed for inertial or simple methods without KL.");
                    graphError = true;
                }
            }

            /* Now check the input values. */
            var errorFound = false;
            var assignmentError = false;

            if (architecture < 0 || architecture > 3)
            {
                Console.WriteLine("Machine architecture parameter = {0:d}, must be in [0,3].", architecture);
                errorFound = true;
            }
            else if (architecture == 0)
            {
                if (ndims_tot < 0)
                {
                    Console.WriteLine("Dimension of hypercube = {0:d}, must be at least 1.", ndims_tot);
                    errorFound = true;
                }
            }
            else if (architecture > 0)
            {
                if (architecture == 1 && mesh_dims[0] <= 0)
                {
                    Console.WriteLine("Size of 1-D mesh improperly specified, {0:d}.", mesh_dims[0]);
                    errorFound = true;
                }

                if (architecture == 2 && (mesh_dims[0] <= 0 || mesh_dims[1] <= 0))
                {
                    Console.WriteLine("Size of 2-D mesh improperly specified, {0:d}x{1:d}.", mesh_dims[0], mesh_dims[1]);
                    errorFound = true;
                }

                if (architecture == 2 && (mesh_dims[0] <= 0 || mesh_dims[1] <= 0 || mesh_dims[2] <= 0))
                {
                    Console.WriteLine("Size of 3-D mesh improperly specified, {0:d}x{1:d}x{2:d}.", mesh_dims[0], mesh_dims[1], mesh_dims[2]);
                    errorFound = true;
                }
            }

            if (ndims < 1 || ndims > MAXDIMS)
            {
                Console.WriteLine("Partitioning at each step = {0:d}, should be in [1,{1:d}].", ndims, MAXDIMS);
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
                Console.WriteLine("Partitioning step {0:d} too large for {1:d} processors.", ndims, nprocs);
                errorFound = true;
            }

            if ((int) partitioningStrategy < 1 || (int) partitioningStrategy > 7)
            {
                Console.WriteLine("Global partitioning method = {0:d}, must be in [1,7].", (int) partitioningStrategy);
                errorFound = true;
            }

            if ((int) localParitioningStrategy < 1 || (int) localParitioningStrategy > 2)
            {
                Console.WriteLine("Local partitioning method = {0:d}, must be in [1,2].", (int) localParitioningStrategy);
                errorFound = true;
            }

            if (partitioningStrategy == PartitioningStrategy.Multilevel_KL || (partitioningStrategy == PartitioningStrategy.Spectral && rqi_flag))
            {
                var i = 2 * (1 << ndims);
                if (*vmax < i)
                {
                    Console.WriteLine("WARNING: Number of vertices in coarse graph ({0:d}) being reset to {1:d}.", *vmax, i);
                    *vmax = i;
                }
            }

            if ((partitioningStrategy == PartitioningStrategy.Multilevel_KL || partitioningStrategy == PartitioningStrategy.Spectral) && eigtol <= 0)
            {
                Console.WriteLine("Eigen tolerance ({0:g}) must be positive value", eigtol);
                errorFound = true;
            }

            if (partitioningStrategy == PartitioningStrategy.Inertial ||
                (MATCH_TYPE == MatchingRoutine.maxmatch5_geometric && (partitioningStrategy == PartitioningStrategy.Multilevel_KL || (partitioningStrategy == PartitioningStrategy.Spectral && rqi_flag))))
            {
                if (igeom < 1 || igeom > 3)
                {
                    Console.WriteLine("Geometry must be 1-, 2- or 3-dimensional");
                    errorFound = true;
                }

                if (igeom > 0 && coords == null)
                {
                    Console.WriteLine("No coordinates given");
                    errorFound = false;
                }
                else if (igeom > 0 && coords[0] == null)
                {
                    Console.WriteLine("No X-coordinates given");
                    errorFound = true;
                }
                else if (igeom > 1 && coords[1] == null)
                {
                    Console.WriteLine("No Y-coordinates given");
                    errorFound = true;
                }
                else if (igeom > 2 && coords[2] == null)
                {
                    Console.WriteLine("No Z-coordinates given");
                    errorFound = true;
                }
            }

            if (partitioningStrategy == PartitioningStrategy.ReadFromFile && localParitioningStrategy == LocalPartitioningStrategy.KernighanLin)
            {
                if (nprocs > 1 << ndims)
                {
                    Console.WriteLine("Can only use local method on single level of read-in assignment,");
                    Console.WriteLine("  but ndims =  {0:d}, while number of processors = {1:d}.", ndims, nprocs);
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
                    Console.WriteLine("goal[{0:d}] is {1:g}, but should be nonnegative.", i, goal[i]);
                    errorFound = true;
                }

                vertexGoalSum += goal[i];
            }

            if (Math.Abs(vertexWeightSum - vertexGoalSum) > 1e-5 * (vertexWeightSum + vertexGoalSum))
            {
                Console.WriteLine("Sum of values in goal ({0:g}) not equal to sum of vertex weights ({1:g}).", vertexGoalSum, vertexWeightSum);
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
                Console.WriteLine("WARNING: OUTPUT_TIME ({0:d}) should be in [0,2].", OUTPUT_TIME);
            }

            if (partitioningStrategy == PartitioningStrategy.Multilevel_KL || partitioningStrategy == PartitioningStrategy.Spectral)
            {
                if (EXPERT)
                {
                    if ((int) LANCZOS_TYPE < 1 || (int) LANCZOS_TYPE > 4)
                    {
                        Console.WriteLine("LANCZOS_TYPE ({0:d}) should be in [1,4].", (int) LANCZOS_TYPE);
                        parameterErrorDetected = true;
                    }
                }
                else
                {
                    if ((int) LANCZOS_TYPE < 1 || (int) LANCZOS_TYPE > 3)
                    {
                        Console.WriteLine("LANCZOS_TYPE ({0:d}) should be in [1,3].", (int) LANCZOS_TYPE);
                        parameterErrorDetected = true;
                    }
                }

                if (EIGEN_TOLERANCE <= 0)
                {
                    Console.WriteLine("EIGEN_TOLERANCE ({0:g}) should be positive.", EIGEN_TOLERANCE);
                    parameterErrorDetected = true;
                }

                if (LANCZOS_SO_INTERVAL <= 0)
                {
                    Console.WriteLine("LANCZOS_SO_INTERVAL ({0:d}) should be positive.", LANCZOS_SO_INTERVAL);
                    parameterErrorDetected = true;
                }

                if (LANCZOS_SO_INTERVAL == 1)
                {
                    Console.WriteLine("WARNING: More efficient if LANCZOS_SO_INTERVAL = 2, not 1.");
                }

                if (BISECTION_SAFETY <= 0)
                {
                    Console.WriteLine("BISECTION_SAFETY ({0:g}) should be positive.", BISECTION_SAFETY);
                    parameterErrorDetected = true;
                }

                if (LANCZOS_CONVERGENCE_MODE < 0 || LANCZOS_CONVERGENCE_MODE > 1)
                {
                    Console.WriteLine("LANCZOS_CONVERGENCE_MODE ({0:d}) should be in [0,1].", LANCZOS_CONVERGENCE_MODE);
                    parameterErrorDetected = true;
                }

                if (WARNING_ORTHTOL <= 0.0d)
                {
                    Console.WriteLine("WARNING: WARNING_ORTHTOL ({0:g}) should be positive.", WARNING_ORTHTOL);
                }

                if (WARNING_MISTOL <= 0.0d)
                {
                    Console.WriteLine("WARNING: WARNING_MISTOL ({0:g}) should be positive.", WARNING_MISTOL);
                }

                if (LANCZOS_SO_PRECISION < 1 || LANCZOS_SO_PRECISION > 2)
                {
                    Console.WriteLine("LANCZOS_SO_PRECISION ({0:d}) should be in [1,2].", LANCZOS_SO_PRECISION);
                    parameterErrorDetected = true;
                }

                if (PERTURB)
                {
                    if (NPERTURB < 0)
                    {
                        Console.WriteLine("NPERTURB ({0:d}) should be nonnegative.", NPERTURB);
                        parameterErrorDetected = true;
                    }

                    if (NPERTURB > 0 && PERTURB_MAX < 0)
                    {
                        Console.WriteLine("PERTURB_MAX ({0:g}) should be nonnegative.", PERTURB_MAX);
                        parameterErrorDetected = true;
                    }
                }

                if ((int) MAPPING_TYPE < 0 || (int) MAPPING_TYPE > 3)
                {
                    Console.WriteLine("MAPPING_TYPE ({0:d}) should be in [0,3].", (int) MAPPING_TYPE);
                    parameterErrorDetected = true;
                }

                if (ndims == 3 && OPT3D_NTRIES <= 0)
                {
                    Console.WriteLine("OPT3D_NTRIES ({0:d}) should be positive.", OPT3D_NTRIES);
                    parameterErrorDetected = true;
                }

                if (partitioningStrategy == PartitioningStrategy.Spectral && rqi_flag)
                {
                    if (COARSE_NLEVEL_RQI <= 0)
                    {
                        Console.WriteLine("COARSE_NLEVEL_RQI ({0:d}) should be positive.", COARSE_NLEVEL_RQI);
                        parameterErrorDetected = true;
                    }

                    if (RQI_CONVERGENCE_MODE < 0 || RQI_CONVERGENCE_MODE > 1)
                    {
                        Console.WriteLine("RQI_CONVERGENCE_MODE ({0:d}) should be in [0,1].", RQI_CONVERGENCE_MODE);
                        parameterErrorDetected = true;
                    }

                    if (TERM_PROP)
                    {
                        Console.WriteLine("WARNING: Using default Lanczos for extended eigenproblem, not RQI/Symmlq.");
                    }
                }

                if (partitioningStrategy == PartitioningStrategy.Multilevel_KL && COARSE_NLEVEL_KL <= 0)
                {
                    Console.WriteLine("COARSE_NLEVEL_KL ({0:d}) should be positive.", COARSE_NLEVEL_KL);
                    parameterErrorDetected = true;
                }

                if (partitioningStrategy == PartitioningStrategy.Multilevel_KL || partitioningStrategy == PartitioningStrategy.Spectral && rqi_flag)
                {
                    if (COARSEN_RATIO_MIN < .5)
                    {
                        Console.WriteLine("COARSEN_RATIO_MIN ({0:g}) should be at least 1/2.", COARSEN_RATIO_MIN);
                        parameterErrorDetected = true;
                    }

                    if ((int)MATCH_TYPE < 1 || (int)MATCH_TYPE > 9)
                    {
                        Console.WriteLine("MATCH_TYPE ({0:d}) should be in [1,9].", MATCH_TYPE);
                        parameterErrorDetected = true;
                    }
                }
            }

            if (partitioningStrategy == PartitioningStrategy.Multilevel_KL || localPartitioningAlgorithm == LocalPartitioningStrategy.KernighanLin)
            {
                if ((int)KL_METRIC < 1 || (int)KL_METRIC > 2)
                {
                    Console.WriteLine("KL_METRIC ({0:d}) should be in [1,2].", KL_METRIC);
                    parameterErrorDetected = true;
                }

                if (KL_BAD_MOVES < 0)
                {
                    Console.WriteLine("KL_BAD_MOVES ({0:d}) should be non-negative.", KL_BAD_MOVES);
                    parameterErrorDetected = true;
                }

                if (KL_NTRIES_BAD < 0)
                {
                    Console.WriteLine("KL_NTRIES_BAD ({0:d}) should be non-negative.", KL_NTRIES_BAD);
                    parameterErrorDetected = true;
                }

                if (KL_IMBALANCE < 0.0 || KL_IMBALANCE > 1.0)
                {
                    Console.WriteLine("KL_IMBALANCE ({0:g}) should be in [0,1].", KL_IMBALANCE);
                    parameterErrorDetected = true;
                }
            }

            if (SIMULATOR < 0 || SIMULATOR > 3)
            {
                Console.WriteLine("SIMULATOR ({0:d}) should be in [0,3].", SIMULATOR);
                parameterErrorDetected = true;
            }

            // ReSharper disable once InvertIf
            if (TERM_PROP)
            {
                if (CUT_TO_HOP_COST <= 0)
                {
                    Console.WriteLine("CUT_TO_HOP_COST ({0:g}) should be positive.", CUT_TO_HOP_COST);
                    parameterErrorDetected = true;
                }

                if (ndims > 1)
                {
                    Console.WriteLine("WARNING: May ignore terminal propagation in spectral quadri/octa section");
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
                    Console.WriteLine("Assignment[{0:d}] = {1:d} less than zero.", i, assignment[i]);
                    flag = true;
                }
                else if (assignment[i] >= nsets_tot)
                {
                    Console.WriteLine("Assignment[{0:d}] = {1:d}, too large for {2:d} sets.", i, assignment[i], nsets_tot);
                    flag = true;
                }
                else if (localPartitioningStrategy == LocalPartitioningStrategy.KernighanLin && assignment[i] >= nsets)
                {
                    Console.WriteLine("Can only use local method on single level of read-in assignment,");
                    Console.WriteLine("  but assignment[{0:d}] =  {1:d}.", i, assignment[i]);
                    flag = true;
                }
            }

            return flag;
        }
    }
}
