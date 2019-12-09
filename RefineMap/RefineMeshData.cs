using static ChacoSharp.Graph.CheckGraph;
using static ChacoSharp.RefineMap.RefineMapHelper;

namespace ChacoSharp.RefineMap
{
    public static unsafe class RefineMeshData
    {
        /// <summary>
        ///
        /// </summary>
        /// <param name="vertex">vertex in comm_graph</param>
        /// <param name="dim">direction of edge from node</param>
        /// <param name="edata">data structure for edge preferences</param>
        /// <param name="mesh_dims">dimensions of mesh</param>
        /// <param name="vtx2node">maps comm_graph vtxs to processors</param>
        private static refine_edata* find_edge_mesh(int vertex, int dim, refine_edata* edata, int[] mesh_dims, int* vtx2node)
        {
            refine_edata* eguy; /* returned pointer to edge info */
            int dir; /* higher or lower direction? */
            int my_node; /* processor vertex assigned to */
            int[] my_loc = new int[3]; /* location of my processor */
            int index = 0; /* computed index into edata */

            if (dim < 0)
            {
                dir = -1;
                dim = -(dim + 1);
            }
            else
            {
                dir = 1;
            }

            my_node = vtx2node[vertex];
            my_loc[0] = my_node % mesh_dims[0];
            my_loc[1] = (my_node / mesh_dims[0]) % mesh_dims[1];
            my_loc[2] = my_node / (mesh_dims[0] * mesh_dims[1]);

            if ((my_loc[dim] == 0 && dir == -1) || (my_loc[dim] == mesh_dims[dim] - 1 && dir == 1))
            {
                eguy = null;
            }

            else
            {
                /* Figure out where edge is in data structure. */
                /* Note: indexing must match with that in init_mesh_edata. */
                if (dir < 0)
                {
                    --my_loc[dim];
                }

                if (dim == 0)
                {
                    /* Edge in x-direction. */
                    index = (mesh_dims[0] - 1) * mesh_dims[1] * my_loc[2] + (mesh_dims[0] - 1) * my_loc[1] +
                            my_loc[0];
                }

                else if (dim == 1)
                {
                    /* Edge in y-direction. */
                    index = (mesh_dims[0] - 1) * mesh_dims[1] * mesh_dims[2] +
                            mesh_dims[0] * (mesh_dims[1] - 1) * my_loc[2] + mesh_dims[0] * my_loc[1] + my_loc[0];
                }
                else if (dim == 2)
                {
                    /* Edge in z-direction. */
                    index = (mesh_dims[0] - 1) * mesh_dims[1] * mesh_dims[2] +
                            mesh_dims[0] * (mesh_dims[1] - 1) * mesh_dims[2] +
                            mesh_dims[0] * mesh_dims[1] * my_loc[2] + mesh_dims[0] * my_loc[1] + my_loc[0];
                }

                eguy = &(edata[index]);
            }

            return (eguy);
        }

        /// <summary>
        ///
        /// </summary>
        /// <param name="vertex">graph vertex being worked on</param>
        /// <param name="dim">mesh dimension to be adjusted</param>
        /// <param name="edata">data structure for edge preferences</param>
        /// <param name="vdata">data structure for vertex preferences</param>
        /// <param name="comm_graph">communication graph</param>
        /// <param name="mesh_dims">extent of mesh</param>
        /// <param name="node2vtx">maps mesh nodes to comm_graph vtxs</param>
        /// <param name="vtx2node">maps mesh nodes to comm_graph vtxs</param>
        /// <param name="best_desire">best desire seen</param>
        /// <param name="imax">offset in desire_ptr array</param>
        /// <param name="desire_ptr">buckets for desire values</param>
        public static void update_mesh_edata(int vertex,int dim, refine_edata* edata, refine_vdata* vdata, vtx_data** comm_graph, int[] mesh_dims /*[3]*/, int* node2vtx, int* vtx2node, double* best_desire, int imax, refine_edata** desire_ptr)
        {
            int i; /* loop counter */

            for (i = 0; i < 2; i++)
            {
                /* Have to adjust two edges. */
                dim = -(dim + 1);
                var eguy = find_edge_mesh(vertex, dim, edata, mesh_dims, vtx2node); /* data for desired edge */
                if (eguy == null)
                {
                    continue;
                }

                var old_desire = eguy->swap_desire; /* original desire for edge to flip */
                var new_desire = (float) compute_mesh_edata(eguy, vdata, mesh_dims, comm_graph, node2vtx); /* new desire for edge to flip */

                if (new_desire == old_desire)
                {
                    continue;
                }

                /* Update linked list if necessary. */
                eguy->swap_desire = new_desire;

                if (new_desire > *best_desire)
                {
                    *best_desire = new_desire;
                }

                /* Remove eguy from it's current place in list. */
                int k; /* loop counter */
                if (eguy->prev == null)
                {
                    /* Round up for index into desire_ptr. */
                    if (old_desire >= 0)
                    {
                        k = (int) old_desire;
                        if (k != old_desire)
                        {
                            k++;
                        }
                    }
                    else
                    {
                        k = (int) -old_desire;
                        if (k != -old_desire)
                        {
                            k++;
                        }

                        k = -k;
                    }

                    k += imax;
                    desire_ptr[k] = eguy->next;
                }
                else
                {
                    eguy->prev->next = eguy->next;
                }

                if (eguy->next != null)
                {
                    eguy->next->prev = eguy->prev;
                }

                /* Now add eguy to it's new desire bucket. */
                if (new_desire >= 0)
                {
                    k = (int) new_desire;
                    if (k != new_desire)
                    {
                        k++;
                    }
                }
                else
                {
                    k = (int) -new_desire;
                    if (k != -new_desire)
                    {
                        k++;
                    }

                    k = -k;
                }

                k += imax;

                eguy->prev = null;
                eguy->next = desire_ptr[k];
                if (desire_ptr[k] != null)
                {
                    desire_ptr[k]->prev = eguy;
                }

                desire_ptr[k] = eguy;
            }
        }

        public static void update_mesh_vdata(int old_loc, /* previous node for moved vertex in moved dimension */
            int new_loc, /* new node for moved vertex in moved dimension */
            int dim, /* dimension that was changed */
            double ewgt, /* weight of edge */
            refine_vdata* vdata, /* array of vertex data */
            int[] mesh_dims /*[3]*/, /* size of processor mesh */
            int neighbor, /* vertex impacted by flip */
            int* vtx2node /* mapping from comm_graph vtxs to processors */
        )
        {
            refine_vdata* vptr = null; /* correct element in vdata */
            int offset = 0; /* index into vdata array */
            int my_loc = 0; /* my location in relevant dimension */
            int neighbor_node = 0; /* processor neighbor assigned to */

            neighbor_node = vtx2node[neighbor];

            if (dim == 0)
            {
                offset = 0;
                my_loc = neighbor_node % mesh_dims[0];
            }
            else if (dim == 1)
            {
                offset = mesh_dims[0] * mesh_dims[1] * mesh_dims[2];
                my_loc = (neighbor_node / mesh_dims[0]) % mesh_dims[1];
            }
            else if (dim == 2)
            {
                offset = 2 * mesh_dims[0] * mesh_dims[1] * mesh_dims[2];
                my_loc = neighbor_node / (mesh_dims[0] * mesh_dims[1]);
            }

            vptr = &vdata[offset + neighbor];

            /* If I'm far away from the flipped edge, I'm not effected. */
            if (!((old_loc < my_loc && new_loc < my_loc) || (old_loc > my_loc && new_loc > my_loc)))
            {
                if (old_loc < my_loc)
                {
                    /* Old moves to the right to line up with me. */
                    vptr->same += (float) ewgt;
                    vptr->below -= (float) ewgt;
                }

                else if (old_loc > my_loc)
                {
                    /* Old moves to the left to line up with me. */
                    vptr->same += (float) ewgt;
                    vptr->above -= (float) ewgt;
                }

                else if (new_loc < my_loc)
                {
                    /* Old moves to the left to pass me. */
                    vptr->same -= (float) ewgt;
                    vptr->below += (float) ewgt;
                }

                else if (new_loc > my_loc)
                {
                    /* Old moves to the right to pass me. */
                    vptr->same -= (float) ewgt;
                    vptr->above += (float) ewgt;
                }
            }
        }


        /* Initialize the mapping of sets to endpoints of wires in the mesh. */
        public static void init_mesh_edata(refine_edata* edata, /* desire data for all edges */
            int[] mesh_dims /*[3]*/ /* dimensions of mesh */
        )
        {
            int wire; /* loops through wires */
            int i, j, k; /* loop counters */

            wire = 0;
            /* First do all the x-axis wires. */
            for (k = 0; k < mesh_dims[2]; k++)
            {
                for (j = 0; j < mesh_dims[1]; j++)
                {
                    for (i = 0; i < mesh_dims[0] - 1; i++)
                    {
                        edata[wire].node1 = (short) (i + mesh_dims[0] * (j + k * mesh_dims[1]));
                        edata[wire].node2 = (short) (i + 1 + mesh_dims[0] * (j + k * mesh_dims[1]));
                        edata[wire].dim = 0;
                        wire++;
                    }
                }
            }

            /* Now do all the y-axis wires. */
            for (k = 0; k < mesh_dims[2]; k++)
            {
                for (j = 0; j < mesh_dims[1] - 1; j++)
                {
                    for (i = 0; i < mesh_dims[0]; i++)
                    {
                        edata[wire].node1 = (short) (i + mesh_dims[0] * (j + k * mesh_dims[1]));
                        edata[wire].node2 = (short) (i + mesh_dims[0] * (j + 1 + k * mesh_dims[1]));
                        edata[wire].dim = 1;
                        wire++;
                    }
                }
            }

            /* Finally, do all the z-axis wires. */
            for (k = 0; k < mesh_dims[2] - 1; k++)
            {
                for (j = 0; j < mesh_dims[1]; j++)
                {
                    for (i = 0; i < mesh_dims[0]; i++)
                    {
                        edata[wire].node1 = (short) (i + mesh_dims[0] * (j + k * mesh_dims[1]));
                        edata[wire].node2 = (short) (i + mesh_dims[0] * (j + (k + 1) * mesh_dims[1]));
                        edata[wire].dim = 2;
                        wire++;
                    }
                }
            }
        }

        public static double compute_mesh_edata(refine_edata* edata, /* desire data for current edge */
            refine_vdata* vdata, /* data for all vertices */
            int[] mesh_dims /*[3]*/, /* dimensions of processor mesh */
            vtx_data** comm_graph, /* communication graph */
            int* node2vtx /* maps mesh nodes to graph vertices */
        )
        {
            double desire; /* edge's interest in flipping */
            float ewgt; /* edge weight */
            int vtx1, vtx2; /* vertices on either side of wire */
            int off; /* index into vdata */

            vtx1 = node2vtx[edata->node1];
            vtx2 = node2vtx[edata->node2];

            off = edata->dim * mesh_dims[0] * mesh_dims[1] * mesh_dims[2];

            desire = (vdata[off + vtx1].above - vdata[off + vtx1].below - vdata[off + vtx1].same) +
                     (vdata[off + vtx2].below - vdata[off + vtx2].above - vdata[off + vtx2].same);

            /* Subtract off potential doubly counted edge. */
            if (is_an_edge(comm_graph[vtx1], vtx2, &ewgt))
            {
                desire -= 2 * ewgt;
            }

            return (desire);
        }

        public static void compute_mesh_vdata(refine_vdata* vdata, /* preference data for a vertex */
            vtx_data** comm_graph, /* communication graph data structure */
            int vtx, /* current vertex */
            int* vtx2node, /* maps graph vtxs to mesh nodes */
            int[] mesh_dims /*[3]*/, /* size of mesh */
            int dim /* dimension we are currently working in */
        )
        {
            float above; /* my preference to move up in each dimension */
            float below; /* my preference to move down in each dimension */
            float same; /* my preference to stay where I am */
            int my_loc; /* my location in mesh */
            int neighb_loc; /* neighbor's location in mesh */
            float ewgt; /* weight of an edge */
            int node; /* set vertex is assigned to */
            int neighbor; /* neighboring vtx in comm_graph */
            int j; /* loop counter */

            node = vtx2node[vtx];

            neighb_loc = 0;
            my_loc = 0;

            if (dim == 0)
            {
                my_loc = node % mesh_dims[0];
            }
            else if (dim == 1)
            {
                my_loc = (node / mesh_dims[0]) % mesh_dims[1];
            }
            else if (dim == 2)
            {
                my_loc = node / (mesh_dims[0] * mesh_dims[1]);
            }

            below = above = same = 0;
            for (j = 1; j < comm_graph[vtx]->nedges; j++)
            {
                neighbor = comm_graph[vtx]->edges[j];
                ewgt = comm_graph[vtx]->ewgts[j];
                node = vtx2node[neighbor];

                if (dim == 0)
                {
                    neighb_loc = node % mesh_dims[0];
                }
                else if (dim == 1)
                {
                    neighb_loc = (node / mesh_dims[0]) % mesh_dims[1];
                }
                else if (dim == 2)
                {
                    neighb_loc = node / (mesh_dims[0] * mesh_dims[1]);
                }

                if (neighb_loc < my_loc)
                {
                    below += ewgt;
                }
                else if (neighb_loc > my_loc)
                {
                    above += ewgt;
                }
                else
                {
                    same += ewgt;
                }
            }

            vdata->below = below;
            vdata->above = above;
            vdata->same = same;
        }
    }
}
