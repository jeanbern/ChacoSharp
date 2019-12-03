#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
#pragma warning disable HAA0101 // Array allocation for params parameter
using System;
using System.Runtime.InteropServices;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Coarsening.BiListHelper;
using static ChacoSharp.Utilities.Timer;
using static ChacoSharp.Utilities.Randomize;
using static ChacoSharp.Coarsening.BucketSort;

namespace ChacoSharp.Coarsening
{
    public static unsafe class NWayKlv
    {
        /*
   Keep guys moved in and guys moving out of separator.
   To restore, move all undesirable guys back.
     (1) guys moved out of separator get put back in.
     (2) guys moved into separator (bspace) get put back.
         Note: Should be done in this order.
         Note: No neighbors need be considered.
     (3) To clear dvals, I should compute touch all guys that
         were ever in separator (bspace) and their neighbors.
*/

        public static bool nway_klv(vtx_data** graph, /* data structure for graph */
            int nvtxs, /* number of vtxs in graph */
            bilist** lbuckets, /* array of lists for bucket sort */
            bilist** rbuckets, /* array of lists for bucket sort */
            bilist* llistspace, /* list data structure for each vertex */
            bilist* rlistspace, /* list data structure for each vertex */
            int* ldvals, /* d-values for each transition */
            int* rdvals, /* d-values for each transition */
            int* sets, /* processor each vertex is assigned to */
            int maxdval, /* maximum d-value for a vertex */
            double[] goal, /* desired set sizes */
            int max_dev, /* largest allowed deviation from balance */
            int** bndy_list, /* list of vertices on boundary (0 ends) */
            double[] weightsum /* sum of vweights in each set (in and out) */
        )
        {
            bilist** to_buckets; /* buckets I'm moving to */
            bilist** from_buckets; /* buckets I'm moving from */
            bilist* to_listspace; /* list structure I'm moving to */
            bilist* from_listspace; /* list structure I'm moving from */
            bilist* out_list; /* list of vtxs moved out of separator */
            int* to_dvals; /* d-values I'm moving to */
            int* from_dvals; /* d-values I'm moving from */
            int* bspace; /* list of active vertices for bucketsort */
            int* edges; /* edge list for a vertex */
            int* edges2; /* edge list for a vertex */
            int* bdy_ptr; /* loops through bndy_list */
            double total_weight; /* weight of all vertices */
            double partial_weight; /* weight of vertices not in separator */
            double ratio; /* fraction of weight not in separator */
            double time; /* timing parameter */
            double delta0; /* largest negative deviation from goal size */
            double delta1; /* largest negative deviation from goal size */
            double left_imbalance = 0.0; /* imbalance if I move to the left */
            double right_imbalance; /* imbalance if I move to the right */
            double balance_val; /* how imbalanced is it */
            double balance_best; /* best balance yet if trying hard */
            bool flag; /* condition indicator */
            int to = -1, from; /* sets moving into / out of */
            int rtop, ltop; /* top of each set of buckets */
            int* to_top; /* ptr to top of set moving to */
            int lvtx, rvtx; /* next vertex to move left/right */
            int lweight, rweight; /* weights of moving vertices */
            int weightfrom = 0; /* weight moving out of a set */
            int list_length; /* how long is list of vertices to bucketsort? */
            bool balanced; /* is partition balanced? */
            bool temp_balanced; /* is intermediate partition balanced? */
            bool ever_balanced; /* has any partition been balanced? */
            int bestvtx = -1; /* best vertex to move */
            int bestval = -1; /* best change in value for a vtx move */
            int vweight; /* weight of best vertex */
            int gtotal; /* sum of changes from moving */
            int improved; /* total improvement from KL */
            double bestg; /* maximum gtotal found in KL loop */
            double bestg_min; /* smaller than any possible bestg */
            int beststep; /* step where maximum value occurred */
            int bestlength; /* step where maximum value occurred */
            int neighbor; /* neighbor of a vertex */
            int step_cutoff; /* number of negative steps in a row allowed */
            int cost_cutoff; /* amount of negative d-values allowed */
            int neg_steps; /* number of negative steps in a row */
            int neg_cost; /* decrease in sum of d-values */
            int vtx; /* vertex number */
            int dval; /* dval of a vertex */
            int group, group2; /* set that a vertex is assigned to */
            bool left_too_big; /* is left set too large? */
            bool right_too_big; /* is right set too large? */
            int vwgt; /* weight of a vertex */
            int gain; /* reduction in separator due to a move */
            int neighbor2; /* neighbor of a vertex */
            int step; /* loops through movements of vertices */
            bool parity; /* sort forwards or backwards? */
            bool done; /* has termination criteria been achieved? */
            int nbad; /* number of unhelpful passes in a row */
            int npass; /* total number of passes */
            int nbadtries; /* number of unhelpful passes before quitting */
            bool enforce_balance; /* force a balanced partition? */
            bool enforce_balance_hard; /* really force a balanced partition? */
            int i, j, k; /* loop counters */

            nbadtries = KL_NTRIES_BAD;

            enforce_balance = false;
            enforce_balance_hard = false;

            total_weight = goal[0] + goal[1];

            bspace = (int*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));

            if (bspace == null)
            {
                return true;
            }

            bdy_ptr = *bndy_list;
            list_length = 0;
            while (*bdy_ptr != 0)
            {
                bspace[list_length++] = *bdy_ptr++;
            }

            Marshal.FreeHGlobal((IntPtr) (*bndy_list));

            clear_dvals(graph, nvtxs, ldvals, rdvals, bspace, list_length);

            step_cutoff = KL_BAD_MOVES;
            cost_cutoff = maxdval * step_cutoff / 7;
            if (cost_cutoff < step_cutoff)
            {
                cost_cutoff = step_cutoff;
            }

            partial_weight = weightsum[0] + weightsum[1];
            ratio = partial_weight / total_weight;
            delta0 = Math.Abs(weightsum[0] - goal[0] * ratio);
            delta1 = Math.Abs(weightsum[1] - goal[1] * ratio);
            balanced =
                (delta0 + delta1 <= max_dev) && weightsum[0] != total_weight && weightsum[1] != total_weight;

            bestg_min = -2.0 * nvtxs * maxdval;
            parity = false;
            nbad = 0;
            npass = 0;
            improved = 0;
            done = false;
            while (!done)
            {
                npass++;
                ever_balanced = false;
                balance_best = delta0 + delta1;

                /* Initialize various quantities. */
                ltop = rtop = 2 * maxdval;

                gtotal = 0;
                bestg = bestg_min;
                beststep = -1;
                bestlength = list_length;
                out_list = null;

                neg_steps = 0;

                /* Compute the initial d-values, and bucket-sort them. */
                time = seconds();
                bucketsortsv(graph, nvtxs, lbuckets, rbuckets, llistspace, rlistspace, ldvals, rdvals, sets,
                    maxdval, parity, bspace, list_length);
                parity = !parity;
                kl_bucket_time += seconds() - time;

                if (DEBUG_KL == DebugFlagKL.PrintBucket)
                {
                    Console.WriteLine("After sorting, left buckets:");
                    p1bucket(lbuckets, llistspace, maxdval);
                    Console.WriteLine("              right buckets:");
                    p1bucket(rbuckets, rlistspace, maxdval);
                }

                /* Now determine the set of vertex moves. */

                for (step = 1;; step++)
                {

                    /* Find the highest d-value in each set. */
                    /* But only consider moves from large to small sets, or moves */
                    /* in which balance is preserved. */
                    /* Break ties in some nonarbitrary manner. */
                    bestval = -maxdval - 1;

                    partial_weight = weightsum[0] + weightsum[1];
                    ratio = partial_weight / total_weight;
                    left_too_big = (weightsum[0] > (goal[0] + .5 * max_dev) * ratio);
                    right_too_big = (weightsum[1] > (goal[1] + .5 * max_dev) * ratio);

                    while (ltop >= 0 && lbuckets[ltop] == null)
                    {
                        --ltop;
                    }

                    if (ltop >= 0 && !left_too_big)
                    {
                        lvtx = (int) (((long) lbuckets[ltop] - (long) llistspace) / sizeof(bilist));
                        lweight = graph[lvtx]->vwgt;
                        rweight = lweight - (ltop - maxdval);
                        weightfrom = rweight;
                        to = 0;
                        bestvtx = lvtx;
                        bestval = ltop - maxdval;
                        partial_weight = weightsum[0] + lweight + weightsum[1] - rweight;
                        ratio = partial_weight / total_weight;
                        left_imbalance = Math.Max(Math.Abs(weightsum[0] + lweight - goal[0] * ratio),
                            Math.Abs(weightsum[1] - rweight - goal[1] * ratio));
                    }

                    while (rtop >= 0 && rbuckets[rtop] == null)
                    {
                        --rtop;
                    }

                    if (rtop >= 0 && !right_too_big)
                    {
                        rvtx = (int) (((long) rbuckets[rtop] - (long) rlistspace) / sizeof(bilist));
                        rweight = graph[rvtx]->vwgt;
                        lweight = rweight - (rtop - maxdval);
                        partial_weight = weightsum[0] - lweight + weightsum[1] + rweight;
                        ratio = partial_weight / total_weight;
                        right_imbalance = Math.Max(Math.Abs(weightsum[0] - lweight - goal[0] * ratio),
                            Math.Abs(weightsum[1] + rweight - goal[1] * ratio));
                        if (rtop - maxdval > bestval || (rtop - maxdval == bestval &&
                                                         (right_imbalance < left_imbalance ||
                                                          (right_imbalance == left_imbalance && drandom() < .5))))
                        {
                            to = 1;
                            weightfrom = lweight;
                            bestvtx = rvtx;
                            bestval = rtop - maxdval;
                        }
                    }

                    if (bestval == -maxdval - 1)
                    {
                        /* No allowed moves */
                        if (DEBUG_KL != DebugFlagKL.NoDebugging)
                        {
                            Console.WriteLine("No KLV moves at step {0:d}.  bestg = {1:g} at step {2:d}.", step, bestg, beststep);
                        }

                        break;
                    }

                    if (to == 0)
                    {
                        from = 1;
                        to_listspace = llistspace;
                        from_listspace = rlistspace;
                        to_dvals = ldvals;
                        from_dvals = rdvals;
                        to_buckets = lbuckets;
                        from_buckets = rbuckets;
                        to_top = &ltop;
                    }
                    else
                    {
                        from = 0;
                        to_listspace = rlistspace;
                        from_listspace = llistspace;
                        to_dvals = rdvals;
                        from_dvals = ldvals;
                        to_buckets = rbuckets;
                        from_buckets = lbuckets;
                        to_top = &rtop;
                    }

                    vweight = graph[bestvtx]->vwgt;

                    weightsum[to] += vweight;
                    weightsum[from] -= weightfrom;

                    /* Check if this partition is balanced. */
                    partial_weight = weightsum[0] + weightsum[1];
                    ratio = partial_weight / total_weight;
                    delta0 = Math.Abs(weightsum[0] - goal[0] * ratio);
                    delta1 = Math.Abs(weightsum[1] - goal[1] * ratio);
                    temp_balanced = (delta0 + delta1 <= max_dev) && weightsum[0] != total_weight &&
                                    weightsum[1] != total_weight;
                    ever_balanced = (ever_balanced || temp_balanced);
                    balance_val = delta0 + delta1;

                    gtotal += bestval;

                    if ((gtotal > bestg && temp_balanced) ||
                        (enforce_balance_hard && balance_val < balance_best))
                    {
                        bestg = gtotal;
                        beststep = step;
                        if (balance_val < balance_best)
                        {
                            balance_best = balance_val;
                        }

                        if (temp_balanced)
                        {
                            enforce_balance_hard = false;
                        }
                    }

                    if (DEBUG_KL == DebugFlagKL.MoreInfo || DEBUG_KL == DebugFlagKL.PrintBucket)
                    {
                        Console.WriteLine("At KLV step {0:d}, bestvtx={1:d}, bestval={2:d} (2->{3:d}), wt0 = {4:g,} wt1 = {5:g}", step, bestvtx, bestval, to, weightsum[0], weightsum[1]);
                    }

                    /* Monitor the stopping criteria. */
                    if (bestval < 0)
                    {
                        if (!enforce_balance || ever_balanced)
                        {
                            neg_steps++;
                        }

                        if (bestg != bestg_min)
                        {
                            neg_cost = (int) bestg - gtotal;
                        }
                        else
                        {
                            neg_cost = -maxdval - 1;
                        }

                        if ((neg_steps > step_cutoff || neg_cost > cost_cutoff) &&
                            !(enforce_balance && bestg == bestg_min))
                        {
                            if (DEBUG_KL != DebugFlagKL.NoDebugging)
                            {
                                if (neg_steps > step_cutoff || neg_cost > cost_cutoff)
                                {
                                    Console.WriteLine("KLV step cutoff at step {0:d}.  bestg = {1:g} at step {2:d}.", step, bestg, beststep);
                                }
                            }

                            weightsum[to] -= vweight;
                            weightsum[from] += weightfrom;
                            break;
                        }
                    }
                    else if (bestval > 0)
                    {
                        neg_steps = 0;
                    }

                    /* Remove vertex from its buckets, and flag it as finished. */
                    sets[bestvtx] = to;
                    removebilist(&to_listspace[bestvtx], &to_buckets[bestval + maxdval]);
                    /*
                                printf("After to removebilist\n");
                                p1bucket(to_buckets, to_listspace, maxdval);
                    */

                    if (from_dvals[bestvtx] != -maxdval - 1)
                    {
                        removebilist(&from_listspace[bestvtx], &from_buckets[from_dvals[bestvtx] + maxdval]);
                        /*
                                        printf("After from removebilist\n");
                                        p1bucket(from_buckets, from_listspace, maxdval);
                        */
                    }

                    from_dvals[bestvtx] = -maxdval - 1;

                    /* Now keep track of vertices moved out of separator so */
                    /* I can restore them as needed. */
                    llistspace[bestvtx].next = out_list;
                    out_list = &(llistspace[bestvtx]);

                    /* Now update the d-values of all the neighbors */
                    /* And neighbors of neighbors ... */

                    /* If left move:
                       1. Separator neighbors right gain => infinity
                       2. Left neighbors unaffected.
                       3. Right neighbors move into separator.
                          A. Right gain = infinity.
                          B. Left gain = computed.
                          C. For any of their neighbors in separator increase left gain.
                    */

                    edges = graph[bestvtx]->edges;
                    for (j = graph[bestvtx]->nedges - 1; j != 0; j--)
                    {
                        neighbor = *(++edges);

                        group = sets[neighbor];

                        if (group == 2)
                        {
                            /* In separator. */
                            gain = from_dvals[neighbor] + maxdval;
                            /* Gain in the from direction => -infinity */
                            if (gain >= 0)
                            {
                                removebilist(&from_listspace[neighbor], &from_buckets[gain]);
                                /*
                                                        printf("\n  After removing %d\n", neighbor);
                                                        p1bucket(from_buckets, from_listspace, maxdval);
                                */
                                from_dvals[neighbor] = -maxdval - 1;
                            }
                        }
                        else if (group == from)
                        {
                            /* Gain in the from direction => -infinity */
                            sets[neighbor] = 2;
                            from_dvals[neighbor] = -maxdval - 1;

                            if (to == 0)
                            {
                                bspace[list_length++] = -neighbor;
                            }
                            else
                            {
                                bspace[list_length++] = neighbor;
                            }

                            edges2 = graph[neighbor]->edges;
                            vwgt = graph[neighbor]->vwgt;
                            gain = graph[neighbor]->vwgt;
                            flag = false;
                            for (k = graph[neighbor]->nedges - 1; k != 0; k--)
                            {
                                neighbor2 = *(++edges2);
                                group2 = sets[neighbor2];
                                if (group2 == 2)
                                {
                                    dval = to_dvals[neighbor2] + maxdval;
                                    if (dval >= 0)
                                    {
                                        movebilist(&to_listspace[neighbor2], &to_buckets[dval], &to_buckets[dval + vwgt]);
                                        /*
                                                                        printf("\n  After moving %d from bucket %d to bucket
                                           %d\n", neighbor2, dval, dval + vwgt);
                                                                        p1bucket(to_buckets, to_listspace, maxdval);
                                        */
                                        to_dvals[neighbor2] += vwgt;
                                        dval += vwgt;
                                        if (dval > *to_top)
                                        {
                                            *to_top = dval;
                                        }
                                    }
                                }
                                else if (group2 == from)
                                {
                                    gain -= graph[neighbor2]->vwgt;
                                    if (to_dvals[neighbor2] + maxdval < 0)
                                    {
                                        flag = true;
                                    }
                                }
                            }

                            if (flag)
                            {
                                /* Not allowed to move further. */
                                to_dvals[neighbor] = -maxdval - 1;
                            }
                            else
                            {
                                to_dvals[neighbor] = gain;
                                /* place in appropriate bucket */

                                gain += maxdval;
                                add2bilist(&to_listspace[neighbor], &to_buckets[gain]);
                                /*
                                                        printf("\nAfter adding %d to bucket %d\n", neighbor, gain -
                                   maxdval);
                                                        p1bucket(to_buckets, to_listspace, maxdval);
                                */

                                if (gain > *to_top)
                                {
                                    *to_top = gain;
                                }
                            }
                        }
                    }

                    if (beststep == step)
                    {
                        bestlength = list_length;
                    }

                    if (DEBUG_KL == DebugFlagKL.PrintBucket)
                    {
                        Console.WriteLine("\n-- After step, left buckets:");
                        p1bucket(lbuckets, llistspace, maxdval);
                        Console.WriteLine("             right buckets:");
                        p1bucket(rbuckets, rlistspace, maxdval);
                    }
                }

                /* Done with a pass; should we actually perform any swaps? */
                if (bestg > 0 || (bestg != bestg_min && !balanced && enforce_balance))
                {
                    improved += (int) bestg;
                }
                else
                {
                    if (enforce_balance_hard)
                    {
                        /* I've done the best I can, give up. */
                        done = true;
                    }

                    if (enforce_balance)
                    {
                        enforce_balance_hard = true;
                    }

                    enforce_balance = true;
                    nbad++;
                }

                /* Work backwards, undoing all the undesirable moves. */

                /* First reset vertices moved out of the separator. */
                if (out_list != null)
                {
                    if (beststep < 0)
                    {
                        beststep = 0;
                    }

                    for (i = step - 1; i > beststep; i--)
                    {
                        vtx = (int) (((long) out_list - (long) llistspace) / sizeof(bilist));
                        if (sets[vtx] != 2)
                        {
                            weightsum[sets[vtx]] -= graph[vtx]->vwgt;
                        }

                        sets[vtx] = 2;
                        out_list = out_list->next;
                    }
                }

                for (i = list_length - 1; i >= bestlength; i--)
                {
                    vtx = bspace[i];
                    if (vtx < 0)
                    {
                        if (sets[-vtx] == 2)
                        {
                            weightsum[1] += graph[-vtx]->vwgt;
                        }

                        sets[-vtx] = 1;
                    }
                    else
                    {
                        if (sets[vtx] == 2)
                        {
                            weightsum[0] += graph[vtx]->vwgt;
                        }

                        sets[vtx] = 0;
                    }
                }

                partial_weight = weightsum[0] + weightsum[1];
                ratio = partial_weight / total_weight;
                delta0 = Math.Abs(weightsum[0] - goal[0] * ratio);
                delta1 = Math.Abs(weightsum[1] - goal[1] * ratio);
                balanced = (delta0 + delta1 <= max_dev) && weightsum[0] != total_weight &&
                           weightsum[1] != total_weight;

                done = done || (nbad >= nbadtries && balanced);
                if (KL_MAX_PASS > 0)
                {
                    done = done || (npass == KL_MAX_PASS && balanced);
                }

                if (!done)
                {
                    /* Rezero dval values. */
                    clear_dvals(graph, nvtxs, ldvals, rdvals, bspace, list_length);
                }

                /* Construct list of separator vertices to pass to buckets or return */
                list_length = make_sep_list(bspace, list_length, sets);

                if (done)
                {
                    bspace[list_length] = 0;
                    bspace = (int*) Marshal.ReAllocHGlobal((IntPtr) bspace, new IntPtr((list_length + 1) * sizeof(int)));
                    *bndy_list = bspace;
                }

                gain = 0;
                j = k = 0;
                for (i = 1; i <= nvtxs; i++)
                {
                    if (sets[i] == 0)
                    {
                        j += graph[i]->vwgt;
                    }
                    else if (sets[i] == 1)
                    {
                        k += graph[i]->vwgt;
                    }
                    else if (sets[i] == 2)
                    {
                        gain += graph[i]->vwgt;
                    }
                }

                /*
                        printf("\nAfter pass of KLV: sets = %d/%d, sep = %d  (bestg = %g)\n\n\n", j, k, gain, bestg);
                */
            }

            if (DEBUG_KL != DebugFlagKL.NoDebugging)
            {
                Console.WriteLine("   KLV required {0:d} passes to improve by {1:d}.", npass, improved);
            }

            return false;
        }


        private static void p1bucket(bilist** a, bilist* b, int c)
        {
            throw new NotImplementedException("I couldn't find the implementation for this in Chaco source.");
        }

        private static void clear_dvals(vtx_data** graph, /* data structure for graph */
            int nvtxs, /* number of vtxs in graph */
            int* ldvals, /* d-values for each transition */
            int* rdvals, /* d-values for each transition */
            int* bspace, /* list of activated vertices */
            int list_length /* number of activated vertices */
        )
        {
            int* edges; /* loops through edge list */
            int vtx; /* vertex in bspace */
            int neighbor; /* neighbor of vtx */
            int i, j; /* loop counters */

            if (list_length > .05 * nvtxs)
            {
                /* Do it directly. */
                for (i = 1; i <= nvtxs; i++)
                {
                    ldvals[i] = rdvals[i] = 0;
                }
            }
            else
            {
                /* Do it more carefully */
                for (i = 0; i < list_length; i++)
                {
                    vtx = bspace[i];
                    if (vtx < 0)
                    {
                        vtx = -vtx;
                    }

                    ldvals[vtx] = rdvals[vtx] = 0;
                    edges = graph[vtx]->edges;
                    for (j = graph[vtx]->nedges - 1; j != 0; j--)
                    {
                        neighbor = *(++edges);
                        ldvals[neighbor] = rdvals[neighbor] = 0;
                    }
                }
            }
        }

        private static int make_sep_list(int* bspace, /* list of vtxs to be moved */
            int list_length, /* current length of bspace */
            int* sets /* processor each vertex is assigned to */
        )
        {
            int vtx; /* vertex in list */
            int i, k; /* loop counters */

            /* Compress out the actual boundary vertices. */
            k = 0;
            for (i = 0; i < list_length; i++)
            {
                vtx = bspace[i];
                if (vtx < 0)
                {
                    vtx = -vtx;
                }

                if (sets[vtx] == 2)
                {
                    bspace[k++] = vtx;
                }
            }

            bspace[k] = 0;
            return (k);
        }
    }
}
