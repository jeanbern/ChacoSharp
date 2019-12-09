using System;
using System.Diagnostics;
using System.Runtime.InteropServices;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Utilities.Timer;
using static ChacoSharp.Graph.FreeGraph;
using static ChacoSharp.Graph.Reformat;
using static ChacoSharp.SubMain.SubMainHelper;

namespace ChacoSharp
{
    public static unsafe class ChacoInterface
    {
        /// <summary>
        ///
        /// </summary>
        /// <param name="nvtxs">number of vertices in full graph</param>
        /// <param name="start">start of edge list for each vertex</param>
        /// <param name="adjacency">edge list data</param>
        /// <param name="vwgts">weights for all vertices</param>
        /// <param name="ewgts">weights for all edges</param>
        /// <param name="x">coordinates for inertial method</param>
        /// <param name="y">coordinates for inertial method</param>
        /// <param name="z">coordinates for inertial method</param>
        /// <param name="assignment">set number of each vtx (length n)</param>
        /// <param name="architecture">0 => hypercube, d => d-dimensional mesh</param>
        /// <param name="ndims_tot">total number of cube dimensions to divide</param>
        /// <param name="mesh_dims">dimensions of mesh of processors</param>
        /// <param name="goal">desired set sizes for each set</param>
        /// <param name="globalMethod">global partitioning algorithm</param>
        /// <param name="localMethod">local partitioning algorithm</param>
        /// <param name="rqi_flag">should I use RQI/Symmlq eigensolver?</param>
        /// <param name="vmax">how many vertices to coarsen down to?</param>
        /// <param name="ndims">number of eigenvectors (2^d sets)</param>
        /// <param name="eigtol">tolerance on eigenvectors</param>
        /// <param name="seed">for random graph mutations</param>
        /// <returns>A flag indicating if an error occured.</returns>
        public static bool INTERFACE(int nvtxs,
            int* start,
            int* adjacency,
            int* vwgts,
            float* ewgts,
            float* x, float* y, float* z,
            int* assignment,
            int architecture,
            int ndims_tot,
            int[] mesh_dims/*[3]*/,
            double[] goal,
            PartitioningStrategy globalMethod,
            LocalPartitioningStrategy localMethod,
            bool rqi_flag,
            int vmax,
            int ndims,
            double eigtol,
            int seed
        )
        {
            vtx_data** graph = null; /* graph data structure */
            float** coords = null; /* coordinates for vertices if used */
            bool flag; /* return code from balance */
            int nedges; /* number of edges in graph */
            var totalSetsCreated = 0; /* total number of sets being created */
            int igeom; /* geometric dimension for inertial method */
            bool defaultGoal; /* using default goals? */

            if (DEBUG_TRACE)
            {
                Trace.WriteLine("<Entering INTERFACE>");
            }

            if (goal == null)
            {
                /* If not passed in, default goals have equal set sizes. */
                defaultGoal = true;
                if (architecture == 0)
                {
                    totalSetsCreated = 1 << ndims_tot;
                }
                else if (architecture == 1)
                {
                    totalSetsCreated = mesh_dims[0];
                }
                else if (architecture == 2)
                {
                    totalSetsCreated = mesh_dims[0] * mesh_dims[1];
                }
                else if (architecture > 2)
                {
                    totalSetsCreated = mesh_dims[0] * mesh_dims[1] * mesh_dims[2];
                }

                double vwgt_sum; /* sum of vertex weights */
                if (MAKE_VWGTS && start != null)
                {
                    vwgt_sum = start[nvtxs] - start[0] + nvtxs;
                }
                else if (vwgts == null)
                {
                    vwgt_sum = nvtxs;
                }
                else
                {
                    vwgt_sum = 0;
                    var vptr = vwgts; /* loops through vertex weights */
                    for (var i = nvtxs; i != 0; i--)
                    {
                        vwgt_sum += *(vptr++);
                    }
                }

                if (totalSetsCreated > 0)
                {
                    vwgt_sum /= totalSetsCreated;
                }

                goal = new double[totalSetsCreated];
                if (goal == null)
                {
                    Trace.WriteLine("\nERROR: No room to make goals.");
                    flag = true;
                    goto skip;
                }

                for (var i = 0; i < totalSetsCreated; i++)
                {
                    goal[i] = vwgt_sum;
                }
            }
            else
            {
                defaultGoal = false;
            }

            if (MAKE_VWGTS)
            {
                /* Generate vertex weights equal to degree of node. */
                if (vwgts != null)
                {
                    Trace.WriteLine("WARNING: Vertex weights being overwritten by vertex degrees.");
                }

                vwgts = (int*) Marshal.AllocHGlobal(nvtxs * sizeof(int));
                if (vwgts == null)
                {
                    Trace.WriteLine("\nERROR: No room to make vertex weights.");
                    flag = true;
                    goto skip;
                }

                if (start != null)
                {
                    for (var i = 0; i < nvtxs; i++)
                    {
                        vwgts[i] = 1 + start[i + 1] - start[i];
                    }
                }
                else
                {
                    for (var i = 0; i < nvtxs; i++)
                    {
                        vwgts[i] = 1;
                    }
                }
            }

            var using_vwgts = (vwgts != null);
            var useEdgeWeights = (ewgts != null);

            if (start != null || vwgts != null)
            {
                /* Reformat into our data structure. */
                double time = seconds(); /* timing variable */
                flag = reformat(start, adjacency, nvtxs, &nedges, vwgts, ewgts, &graph);

                if (flag)
                {
                    Trace.WriteLine("\nERROR: No room to reformat graph.");
                    goto skip;
                }

                reformat_time += seconds() - time;
            }
            else
            {
                nedges = 0;
            }

            if (FREE_GRAPH)
            {
                /* Free old graph data structures. */
                Marshal.FreeHGlobal((IntPtr) start);
                start = null;
                Marshal.FreeHGlobal((IntPtr) adjacency);
                adjacency = null;
                Marshal.FreeHGlobal((IntPtr) vwgts);
                vwgts = null;
                Marshal.FreeHGlobal((IntPtr) ewgts);
                ewgts = null;
            }

            if (globalMethod == PartitioningStrategy.Inertial ||
                (MATCH_TYPE == MatchingRoutine.maxmatch5_geometric && (globalMethod == PartitioningStrategy.Multilevel_KL || (globalMethod == PartitioningStrategy.Spectral && rqi_flag))))
            {
                if (x == null)
                {
                    igeom = 0;
                }
                else
                {
                    /* Set up coordinate data structure. */
                    coords = (float**) Marshal.AllocHGlobal(3 * sizeof(float*));
                    if (coords == null)
                    {
                        Trace.WriteLine("\nERROR: No room to make coordinate array.");
                        flag = true;
                        goto skip;
                    }

                    /* Minus 1's are to allow remainder of program to index with 1. */
                    coords[0] = x - 1;
                    igeom = 1;
                    if (y != null)
                    {
                        coords[1] = y - 1;
                        igeom = 2;
                        if (z != null)
                        {
                            coords[2] = z - 1;
                            igeom = 3;
                        }
                    }
                }
            }
            else
            {
                igeom = 0;
            }

            /* Subtract from assignment to allow code to index from 1. */
            assignment = assignment - 1;
            flag = submain(graph, nvtxs, nedges, using_vwgts, useEdgeWeights, igeom, coords,
                assignment, goal, architecture, ndims_tot, mesh_dims, globalMethod,
                localMethod, rqi_flag, vmax, ndims, eigtol, seed);

            skip:
            Marshal.FreeHGlobal((IntPtr) coords);

            if (defaultGoal)
            {
                goal = null;
            }

            if (graph != null)
            {
                free_graph(graph);
            }

            if (flag && FREE_GRAPH)
            {
                Marshal.FreeHGlobal((IntPtr) start);
                Marshal.FreeHGlobal((IntPtr) adjacency);
                Marshal.FreeHGlobal((IntPtr) vwgts);
                Marshal.FreeHGlobal((IntPtr) ewgts);
            }

            return (flag);
        }

    }
}
