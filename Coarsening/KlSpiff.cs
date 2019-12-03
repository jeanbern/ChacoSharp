using System;
using System.Runtime.InteropServices;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Utilities.Timer;
using static ChacoSharp.Utilities.CountHelper;
using static ChacoSharp.Coarsening.NWayKl;
using static ChacoSharp.Utilities.ArrayAlloc2d;

namespace ChacoSharp.Coarsening
{
    public static unsafe class KlSpiff
    {
        public static void klspiff(vtx_data** graph, /* list of graph info for each vertex */
            int nvtxs, /* number of vertices in graph */
            int* sets, /* local partitioning of vtxs */
            int nsets, /* number of sets at each level */
            int[][] hops, /* hop cost between sets */
            double[] goal, /* desired set sizes */
            float*[] term_wgts, /* weights for terminal propagation */
            int max_dev, /* largest deviation from balance allowed */
            double maxdeg, /* largest weighted vertex degree */
            bool useEdgeWeights, /* are edge weights being used? */
            int** bndy_list, /* list of vertices on boundary (0 ends) */
            double[] weights /* vertex weights in each set */
        )
        {
            bilist**** buckets; /* space for bucket sorts */
            bilist** listspace; /* space for all bidirectional elements */
            float* twptr; /* loops through term_wgts */
            int** dvals; /* change in penalty for each possible move */
            int** tops; /* starting dval for each type of move */
            double time, time1; /* timing variables */
            float maxterm; /* largest terminal propagation preference */
            int maxhop; /* maximum hops between sets */
            int maxdval; /* largest transition cost for a vertex */
            double cut_cost; /* relative importance of cuts/hops */
            double hop_cost; /* relative importance of hops/cuts */
            bool error; /* out of space? */
            int i, j; /* loop counters */

            time = seconds();

            if (DEBUG_TRACE)
            {
                Console.WriteLine("<Entering klspiff, nvtxs = {0:d}>", nvtxs);
            }

            /* Find the largest hop value. */
            maxhop = 0;
            for (i = 0; i < nsets; i++)
            {
                for (j = 0; j < nsets; j++)
                {
                    if (hops[i][j] > maxhop)
                    {
                        maxhop = hops[i][j];
                    }
                }
            }

            maxterm = 0;
            cut_cost = hop_cost = 1;
            if (term_wgts[1] != null)
            {
                for (j = 1; j < nsets; j++)
                {
                    twptr = term_wgts[j];
                    for (i = nvtxs; i != 0; i--)
                    {
                        ++twptr;
                        if (*twptr > maxterm)
                        {
                            maxterm = *twptr;
                        }
                        else if (-*twptr > maxterm)
                        {
                            maxterm = -*twptr;
                        }
                    }
                }

                if (CUT_TO_HOP_COST > 1)
                {
                    cut_cost = CUT_TO_HOP_COST;
                }
                else
                {
                    hop_cost = 1.0 / CUT_TO_HOP_COST;
                }
            }

            maxdval = (int) ((2 * maxterm * hop_cost + .5) + (maxdeg * cut_cost + .5) * maxhop);

            /* Allocate a bunch of space for KL. */
            time1 = seconds();
            error = kl_init(&buckets, &listspace, &dvals, &tops, nvtxs, nsets, maxdval);
            kl_init_time += seconds() - time1;

            if (!error)
            {
                if (DEBUG_KL != DebugFlagKL.NoDebugging)
                {
                    Console.WriteLine(" Before KL: ");
                    count(graph, nvtxs, sets, nsets, hops, false, useEdgeWeights);
                }

                time1 = seconds();
                error = nway_kl(graph, nvtxs, buckets, listspace, tops, dvals, sets, maxdval, nsets, goal,
                    term_wgts, hops, max_dev, useEdgeWeights, bndy_list, weights);
                nway_kl_time += seconds() - time1;

                if (DEBUG_KL == DebugFlagKL.MoreInfo || DEBUG_KL == DebugFlagKL.PrintBucket)
                {
                    Console.WriteLine(" After KL:");
                    count(graph, nvtxs, sets, nsets, hops, false, useEdgeWeights);
                }
            }

            if (error)
            {
                Console.WriteLine("\nWARNING: No space to perform KL on graph with {0:d} vertices.", nvtxs);
                Console.WriteLine("         NO LOCAL REFINEMENT PERFORMED.\n");
            }

            free_kl(buckets, listspace, dvals, tops);

            kl_total_time += seconds() - time;
        }

        static void free_kl(
            /* Free everything malloc'd for KL. */
            bilist**** buckets, /* space for bucket sorts */
            bilist** listspace, /* space for all bidirectional elements */
            int** dvals, /* change in penalty for each possible move */
            int** tops /* starting dval for each type of move */
        )
        {

            Marshal.FreeHGlobal((IntPtr) dvals);
            Marshal.FreeHGlobal((IntPtr) tops);

            Marshal.FreeHGlobal((IntPtr) listspace[0]);
            Marshal.FreeHGlobal((IntPtr) buckets[0][1]);
            Marshal.FreeHGlobal((IntPtr) listspace);
            Marshal.FreeHGlobal((IntPtr) buckets);
        }

        public static void make_bndy_list(vtx_data** graph, /* data structure for graph */
            bilist* movelist, /* list of vtxs to be moved */
            bilist**** buckets, /* array of lists for bucket sort */
            bilist** listspace, /* list data structure for each vertex */
            int* sets, /* processor each vertex is assigned to */
            int nsets, /* number of sets divided into */
            int* bspace, /* list of active vertices for bucketsort */
            int** tops, /* top of each set of buckets */
            int** bndy_list /* list of boundary vertices returned */
        )
        {
            bilist* bptr; /* loops through bspace */
            int vtx; /* vertex that was moved */
            int set; /* set a vertex is in */
            int list_length; /* number of active vertices */
            int bndy_length; /* number of vertices actually on boundary */
            int size; /* array spacing */
            int i, j, k; /* loop counters */

            /* First push all the moved vertices onto list, so they can be flagged. */
            /* They've already been removed from buckets, so want to avoid them. */
            size = (int) (&(listspace[0][1]) - &(listspace[0][0]));
            bptr = movelist;
            list_length = 0;
            while (bptr != null)
            {
                vtx = ((int) (bptr - listspace[0])) / size;
                bspace[list_length++] = vtx;
                bptr = bptr->next;
            }

            /* Now get all the vertices still in the bucket lists. */
            for (k = tops[0][1]; k >= 0; k--)
            {
                bptr = buckets[0][1][k];
                while (bptr != null)
                {
                    vtx = ((int) (bptr - listspace[0])) / size;
                    bspace[list_length++] = vtx;
                    bptr = bptr->next;
                }
            }

            for (i = 1; i < nsets; i++)
            {
                for (k = tops[i][0]; k >= 0; k--)
                {
                    bptr = buckets[i][0][k];
                    while (bptr != null)
                    {
                        vtx = ((int) (bptr - listspace[0])) / size;
                        bspace[list_length++] = vtx;
                        bptr = bptr->next;
                    }
                }
            }

            /* Now that list is constructed, go reconstruct all the set numbers. */
            bndy_length = 0;
            for (i = 0; i < list_length; i++)
            {
                vtx = bspace[i];
                set = sets[vtx];
                for (j = 1; j < graph[vtx]->nedges; j++)
                {
                    if (sets[graph[vtx]->edges[j]] != set)
                    {
                        bspace[bndy_length++] = vtx;
                        break;
                    }
                }
            }

            /* Finally, copy boundary vertices into boundary list. */
            *bndy_list = (int*) Marshal.AllocHGlobal((bndy_length + 1) * sizeof(int));
            for (i = 0; i < bndy_length; i++)
            {
                (*bndy_list)[i] = bspace[i];
            }

            (*bndy_list)[bndy_length] = 0;
        }

        private static bool kl_init(bilist***** bucket_ptrs, /* space for multiple bucket sorts */
            bilist*** listspace, /* space for all elements of linked lists */
            int*** dvals, /* change in cross edges for each move */
            int*** tops, /* top dval for each type of move */
            int nvtxs, /* number of vertices in the graph */
            int nsets, /* number of sets created at each step */
            int maxchange /* maximum change by moving a vertex */
        )
        {
            bilist* spacel; /* space for all listspace entries */
            bilist** spaceb; /* space for all buckets entries */
            int sizeb; /* size of set of buckets */
            int sizel; /* size of set of pointers for all vertices */
            int i, j; /* loop counters */

            /* Allocate appropriate data structures for buckets, and listspace. */

            *bucket_ptrs = (bilist****) array_alloc_2D_ret<IntPtr>(nsets, nsets, sizeof(bilist*));

            *dvals = (int**) array_alloc_2D_ret<int>(nvtxs + 1, nsets - 1, sizeof(int));

            *tops = (int**) array_alloc_2D_ret<int>(nsets, nsets, sizeof(int));

            /* By using '-1' in the next line, I save space, but I need to */
            /* be careful to get the right element in listspace each time. */
            *listspace = (bilist**) Marshal.AllocHGlobal((nsets - 1) * sizeof(bilist*));

            sizeb = (2 * maxchange + 1) * sizeof(bilist*);
            sizel = (nvtxs + 1) * sizeof(bilist);
            spacel = (bilist*) Marshal.AllocHGlobal((nsets - 1) * sizel);
            spaceb = (bilist**) Marshal.AllocHGlobal(nsets * (nsets - 1) * sizeb);

            if (*bucket_ptrs == null || *dvals == null || *tops == null || *listspace == null ||
                spacel == null || spaceb == null)
            {
                Marshal.FreeHGlobal((IntPtr) spacel);
                Marshal.FreeHGlobal((IntPtr) spaceb);
                return true;
            }

            for (i = 0; i < nsets; i++)
            {
                if (i != nsets - 1)
                {
                    (*listspace)[i] = spacel;
                    spacel += nvtxs + 1;
                }

                for (j = 0; j < nsets; j++)
                {
                    if (i != j)
                    {
                        (*bucket_ptrs)[i][j] = spaceb;
                        spaceb += 2 * maxchange + 1;
                    }
                }
            }

            return false;
        }


        static float *old_ewgts; /* space for old edge weights */

public static void compress_ewgts(vtx_data **graph,      /* list of graph info for each vertex */
                    int               nvtxs,      /* number of vertices in graph */
                    int               nedges,     /* number of edges in graph */
                    double            ewgt_max,   /* largest edge weight */
                    bool               useEdgeWeights /* are edge weights being used? */
)
{
  float *       old_ewptr;      /* loops through old edge weights */
  float *       new_ewptr;      /* loops through old edge weights */
  float *       self_ptr;       /* points to self edge in new_ewgts */
  float *       new_ewgts;      /* space for new edge weights */
  int           ewgt;           /* new edge weight value */
  double        sum;            /* sum of all the new edge weights */
  double        ratio;          /* amount to reduce edge weights */
  int           i, j;           /* loop counter */

  /* Check easy cases first. */
  if (!useEdgeWeights) {
    old_ewgts = null;
  }
  else if (ewgt_max < EWGT_RATIO_MAX * nvtxs) {
    /* If not too heavy, leave it alone. */
    old_ewgts = null;
    Console.WriteLine("In compress_ewgts, but not too heavy, ewgt_max = {0:g}, nvtxs = {1:d}", ewgt_max, nvtxs);
  }

  else { /* Otherwise, compress edge weights. */
    /* Allocate space for new edge weights */
    old_ewgts = graph[1]->ewgts;
    new_ewgts = (float*)Marshal.AllocHGlobal((2 * nedges + nvtxs) * sizeof(float));
    ratio     = (EWGT_RATIO_MAX * nvtxs) / ewgt_max;
    Console.WriteLine("In compress_ewgts, ewgt_max = {0:g}, nvtxs = {1:d}, ratio = {2:e}", ewgt_max, nvtxs, ratio);
    old_ewptr = old_ewgts;
    new_ewptr = new_ewgts;
    for (i = 1; i <= nvtxs; i++) {
      self_ptr = new_ewptr++;
      old_ewptr++;
      sum = 0;
      for (j = graph[i]->nedges - 1; j != 0; j--) {
        ewgt           = (int)(*(old_ewptr++) * ratio + 1.0);
        *(new_ewptr++) = (float)ewgt;
        sum += (float)ewgt;
      }
      *self_ptr       = (float)-sum;
      graph[i]->ewgts = self_ptr;
    }
  }
}

public static void restore_ewgts(vtx_data **graph, /* list of graph info for each vertex */
                   int               nvtxs  /* number of vertices in graph */
)
{
  int i; /* loop counter */

  /* Check easy case first. */
  if (old_ewgts == null) {
    return;
  }
  else { /* otherwise, compress edge weights. */
    Marshal.FreeHGlobal((IntPtr)(graph[1]->ewgts));
    for (i = 1; i <= nvtxs; i++) {
      graph[i]->ewgts = old_ewgts;
      old_ewgts += graph[i]->nedges;
    }
    old_ewgts = null;
  }
}
    }
}
