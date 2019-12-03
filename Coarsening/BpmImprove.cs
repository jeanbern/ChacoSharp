#pragma warning disable HAA0101 // Array allocation for params parameter
#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
using System;
using System.Runtime.InteropServices;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Coarsening.MatchingHelper;
using static ChacoSharp.Coarsening.Flow;
using static ChacoSharp.Coarsening.MakeBpGraph;

namespace ChacoSharp.Coarsening
{
    public static unsafe class BpmImprove
    {
        /* Refine a vertex separator by finding a maximum bipartite matching. */
        public static void bpm_improve(vtx_data** graph, /* list of graph info for each vertex */
            int* sets, /* local partitioning of vtxs */
            double[] goal, /* desired set sizes */
            int max_dev, /* largest deviation from balance allowed */
            int** bndy_list, /* list of vertices on boundary (0 ends) */
            double[] weights, /* vertex weights in each set */
            bool using_vwgts /* invoke weighted cover routines? */
        )
        {
            var separatorSize = 0;
            while ((*bndy_list)[separatorSize] != 0)
            {
                separatorSize++;
            }

            var separatorWeight = 0;
            if (using_vwgts)
            {
                for (var i = 0; i < separatorSize; i++)
                {
                    separatorWeight += graph[(*bndy_list)[i]]->vwgt;
                }
            }
            else
            {
                separatorWeight = separatorSize;
            }

            if (DEBUG_COVER > 1)
            {
                Console.WriteLine("Before first matching, sep_size = {0:d}, sep_weight = {1:d},  Sizes = {2:g}/{3:g}", separatorSize,
                    separatorWeight, weights[0], weights[1]);
            }

            var ratio = (weights[0] + weights[1]) / (goal[0] + goal[1]);
            var deltaplus = Math.Abs(weights[0] - goal[0] * ratio);
            var deltaminus = Math.Abs(weights[1] - goal[1] * ratio);
            var imbalance = deltaplus + deltaminus;
            var oldCost = sep_cost(weights[0], weights[1], (double) separatorWeight, (double) max_dev);

            var isImprovedmentFound = true;
            while (isImprovedmentFound)
            {
                /* First match towards the larger side, then the smaller. */
                int setBig; /* side of graph I'm matching against */
                int setSmall; /* side of graph I'm not matching against */
                if (goal[0] - weights[0] >= goal[1] - weights[1])
                {
                    setBig = 1;
                    setSmall = 0;
                }
                else
                {
                    setBig = 0;
                    setSmall = 1;
                }

                isImprovedmentFound = bpm_improve1(graph, sets, bndy_list, weights, setBig, setSmall, goal, max_dev,
                    &imbalance, &separatorSize, &separatorWeight, using_vwgts, &oldCost);

                if (DEBUG_COVER != 0)
                {
                    Console.WriteLine("After big matching, sep_size = {0:d}, sep_weight = {1:d},  Sizes = {2:g}/{3:g}", separatorSize,
                        separatorWeight, weights[0], weights[1]);
                }

                if (VERTEX_COVER)
                {
                    break;
                }

                if (!isImprovedmentFound)
                {
                    /* If balanced, try the other direction. */
                    if (imbalance < max_dev)
                    {
                        isImprovedmentFound = bpm_improve1(graph, sets, bndy_list, weights, setSmall, setBig, goal, max_dev,
                            &imbalance, &separatorSize, &separatorWeight, using_vwgts, &oldCost);

                        if (DEBUG_COVER != 0)
                        {
                            Console.WriteLine("After small matching, sep_size = {0:d},  Sizes = {1:g}/{2:g}", separatorSize, weights[0], weights[1]);
                        }
                    }
                }
            }

            if (DEBUG_COVER != 0)
            {
#pragma warning disable HAA0101 // Array allocation for params parameter
#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
                Console.WriteLine("After all matchings, sep_size = {0:d}, sep_weight = {1:d},  Sizes = {2:g}/{3:g}\n", separatorSize, separatorWeight, weights[0], weights[1]);
            }
        }

        static bool bpm_improve1(vtx_data** graph, /* list of graph info for each vertex */
            int* sets, /* local partitioning of vtxs */
            int** pbndy_list, /* list of vertices on boundary (0 ends) */
            double[] weights, /* vertex weights in each set */
            int set_match, /* side of graph I'm matching against */
            int set_other, /* side of graph I'm not matching against */
            double[] goal, /* desired set sizes */
            int max_dev, /* largest deviation from balance allowed */
            double* pimbalance, /* imbalance of current partition */
            int* sep_size, /* separator size */
            int* sep_weight, /* weight of separator */
            bool using_vwgts, /* use weighted model? */
            double* pcost /* cost of current separator */
        )
        {
            double[] new_weights = new double[2]; /* weights associated with new separator */
            double ratio; /* fraction of non-separator vertices */
            double deltaplus; /* amount set is too big */
            double deltaminus; /* amount set is too small */
            double new_imbalance; /* imbalance of new partition */
            double new_cost; /* cost of new separator */
            int* pointers; /* start/stop indices into adjacencies */
            int* indices; /* adjacencies for each bipartite vertex */
            int* vweight; /* vertex weights if needed */
            int* loc2glob; /* mapping from bp graph to original */
            int* new_bndy_list; /* new list of boundary vertices */
            int old_sep_size; /* previous separator size */
            int old_sep_weight; /* previous separator weight */
            int vtx; /* vertex in graph */
            bool change; /* does this routine alter separator? */
            int nleft, nright; /* # vtxs in two sides on bp graph */
            int i, j; /* loop counter */

            make_bpgraph(graph, sets, *pbndy_list, *sep_size, set_match, &pointers, &indices, &vweight,
                &loc2glob, &nleft, &nright, using_vwgts);

            old_sep_size = *sep_size;
            old_sep_weight = *sep_weight;
            if (!using_vwgts)
            {
                new_bndy_list = (int*) Marshal.AllocHGlobal((*sep_size + 1) * sizeof(int));
                new_bndy_list[0] = nleft + nright;
                bpcover(nleft, nright, pointers, indices, sep_size, new_bndy_list);
                *sep_weight = *sep_size;
            }
            else
            {
                wbpcover(nleft, nright, pointers, indices, vweight, sep_size, sep_weight, &new_bndy_list);
            }

            /* Update weights. */
            new_weights[0] = weights[0];
            new_weights[1] = weights[1];
            for (j = 0; j < new_bndy_list[0]; j++)
            {
                /* First handle nodes numbered less than separator nodes. */
                vtx = loc2glob[j];
                if (sets[vtx] == 2)
                {
                    new_weights[set_other] += graph[vtx]->vwgt;
                }
            }

            for (i = 0; i < *sep_size; i++)
            {
                vtx = loc2glob[new_bndy_list[i]];
                if (sets[vtx] == set_match)
                {
                    new_weights[set_match] -= graph[vtx]->vwgt;
                }

                if (i != 0)
                {
                    for (j = new_bndy_list[i - 1] + 1; j < new_bndy_list[i]; j++)
                    {
                        vtx = loc2glob[j];
                        if (sets[vtx] == 2)
                        {
                            new_weights[set_other] += graph[vtx]->vwgt;
                        }
                    }
                }
            }

            if (*sep_size != 0)
            {
                i = new_bndy_list[*sep_size - 1] + 1;
            }
            else
            {
                i = 0;
            }

            for (j = i; j < nleft + nright; j++)
            {
                vtx = loc2glob[j];
                if (sets[vtx] == 2)
                {
                    new_weights[set_other] += graph[vtx]->vwgt;
                }
            }

            /* Check to see if new partition is acceptably balanced. */
            ratio = (new_weights[0] + new_weights[1]) / (goal[0] + goal[1]);
            deltaplus = Math.Abs(new_weights[0] - goal[0] * ratio);
            deltaminus = Math.Abs(new_weights[1] - goal[1] * ratio);
            new_imbalance = deltaplus + deltaminus;

            new_cost = sep_cost(weights[0], weights[1], (double) *sep_weight, (double) max_dev);

            if (DEBUG_COVER > 1)
            {
                Console.WriteLine("Sides {0:f}, {1:f}: sep {2:d} total {3:f} {4:f}", new_weights[0], new_weights[1], *sep_size, new_weights[0] + new_weights[1], new_weights[0] + new_weights[1] + *sep_size);
            }

            /* if (new_cost < *pcost) { */
            if ((new_cost < *pcost && new_imbalance <= max_dev) ||
                (new_cost <= *pcost && new_imbalance < *pimbalance))
            {
                /* Update set values. */
                change = true;
                *pcost = new_cost;
                for (j = 0; j < new_bndy_list[0]; j++)
                {
                    /* First handle nodes numbered  less than separator nodes. */
                    vtx = loc2glob[j];
                    if (sets[vtx] == 2)
                    {
                        sets[vtx] = set_other;
                    }
                }

                for (i = 0; i < *sep_size; i++)
                {
                    vtx = loc2glob[new_bndy_list[i]];
                    if (sets[vtx] == set_match)
                    {
                        sets[vtx] = 2;
                    }

                    if (i != 0)
                    {
                        for (j = new_bndy_list[i - 1] + 1; j < new_bndy_list[i]; j++)
                        {
                            vtx = loc2glob[j];
                            if (sets[vtx] == 2)
                            {
                                sets[vtx] = set_other;
                            }
                        }
                    }
                }

                if (*sep_size != 0)
                {
                    i = new_bndy_list[*sep_size - 1] + 1;
                }
                else
                {
                    i = 0;
                }

                for (j = i; j < nleft + nright; j++)
                {
                    vtx = loc2glob[j];
                    if (sets[vtx] == 2)
                    {
                        sets[vtx] = set_other;
                    }
                }

                /* Restore bndy_list to global numbering. */
                for (i = 0; i < *sep_size; i++)
                {
                    new_bndy_list[i] = loc2glob[new_bndy_list[i]];
                }

                new_bndy_list[*sep_size] = 0;

                Marshal.FreeHGlobal((IntPtr) (*pbndy_list));

                *pbndy_list = new_bndy_list;
                *pimbalance = new_imbalance;

                weights[0] = new_weights[0];
                weights[1] = new_weights[1];
            }

            else
            {
                change = false;
                Marshal.FreeHGlobal((IntPtr) new_bndy_list);
                *sep_size = old_sep_size;
                *sep_weight = old_sep_weight;
            }

            Marshal.FreeHGlobal((IntPtr) vweight);
            Marshal.FreeHGlobal((IntPtr) loc2glob);
            Marshal.FreeHGlobal((IntPtr) indices);
            Marshal.FreeHGlobal((IntPtr) pointers);

            return (change);
        }

        private static double sep_cost(double v1, double v2, double v3, double max_dev)
        {
            return 0.0d;
        }

        /* Routine that can be modified to allow different cost functions. */

        static double sep_cost(double size_sep /* maximum allowed imbalance */ )
        {
            return (size_sep);
        }
    }
}
