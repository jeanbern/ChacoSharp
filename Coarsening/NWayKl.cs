#pragma warning disable HAA0101 // Array allocation for params parameter
#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
using System;
using System.Runtime.InteropServices;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Utilities.Timer;
using static ChacoSharp.Utilities.Randomize;
using static ChacoSharp.Coarsening.KlSpiff;
using static ChacoSharp.Coarsening.BucketSort;

namespace ChacoSharp.Coarsening
{
    public static unsafe class NWayKl
    {
        /* Idea:
   'buckets[i][j]' is a set of buckets to sort moves from i to j.
   listspace[i] is space for lists in buckets[i][j].
   Loop through all nonequal pairs [i][j], taking the first element
   in each list.  Compare them all to find the largest allowed move.
   Make that move, and save it in movelist.
*/


public static bool nway_kl(vtx_data **graph,       /* data structure for graph */
            int               nvtxs,       /* number of vtxs in graph */
            bilist ****buckets,     /* array of lists for bucket sort */
            bilist **  listspace,   /* list data structure for each vertex */
            int **            tops,        /* 2-D array of top of each set of buckets */
            int **            dvals,       /* d-values for each transition */
            int *             sets,        /* processor each vertex is assigned to */
            int               maxdval,     /* maximum d-value for a vertex */
            int               nsets,       /* number of sets divided into */
            double []          goal,        /* desired set sizes */
            float *[]           term_wgts, /* weights for terminal propagation */
            int [][] hops,          /* cost of set transitions */
            int     max_dev,               /* largest allowed deviation from balance */
            bool     useEdgeWeights,           /* are edge weights being used? */
            int **  bndy_list,             /* list of vertices on boundary (0 ends) */
            double [] startweight            /* sum of vweights in each set (in and out) */
)

/* Suaris and Kedem algorithm for quadrisection, generalized to an */
/* arbitrary number of sets, with intra-set cost function specified by hops. */
/* Note: this is for a single divide step. */
/* Also, sets contains an initial (possibly crummy) partitioning. */

{
  bilist * movelist;                   /* list of vtxs to be moved */
  bilist **endlist;                    /* end of movelists */
  bilist * bestptr;                    /* best vertex in linked list */
  bilist * bptr;                       /* loops through bucket list */
  float *         ewptr     = null;           /* loops through edge weights */
  double *        locked    = null;           /* weight of vertices locked in a set */
  double *        loose     = null;           /* weight of vtxs that can move from a set */
  int *           bspace    = null;           /* list of active vertices for bucketsort */
  double *        weightsum = null;           /* sum of vweights for each partition */
  int *           edges     = null;           /* edge list for a vertex */
  int *           bdy_ptr   = null;           /* loops through bndy_list */
  double          time;                       /* timing parameter */
  double          delta;                      /* desire of sets to change size */
  double          bestdelta = -1;             /* strongest delta value */
  double          deltaplus;                  /* largest negative deviation from goal size */
  double          deltaminus;                 /* largest negative deviation from goal size */
  int             list_length;                /* how long is list of vertices to bucketsort? */
  bool             balanced;                   /* is partition balanced? */
  bool             temp_balanced;              /* is intermediate partition balanced? */
  bool             ever_balanced;              /* has any partition been balanced? */
  int             bestvtx  = -1;              /* best vertex to move */
  int             bestval  = -1;              /* best change in value for a vtx move */
  int             bestfrom = -1, bestto = -1; /* sets best vertex moves between */
  int             vweight;                    /* weight of best vertex */
  int             gtotal;                     /* sum of changes from moving */
  int             improved;                   /* total improvement from KL */
  double          balance_val = 0.0;          /* how imbalanced is it */
  double          balance_best;               /* best balance yet if trying hard */
  double          bestg;                      /* maximum gtotal found in KL loop */
  double          bestg_min;                  /* smaller than any possible bestg */
  int             beststep;                   /* step where maximum value occurred */
  int             neighbor;                   /* neighbor of a vertex */
  int             step_cutoff;                /* number of negative steps in a row allowed */
  int             cost_cutoff;                /* amount of negative d-values allowed */
  int             neg_steps;                  /* number of negative steps in a row */
  int             neg_cost;                   /* decrease in sum of d-values */
  int             vtx;                        /* vertex number */
  int             dval;                       /* dval of a vertex */
  int             group;                      /* set that a vertex is assigned to */
  double          cut_cost;                   /* if term_prop; relative cut/hop importance */
  int             diff;                       /* change in a d-value */
  int             stuck1st, stuck2nd;         /* how soon will moves be disallowed? */
  int             beststuck1 = -1, beststuck2 = -1; /* best stuck values for tie-breaking */
  int             eweight;                          /* a particular edge weight */
  bool             worth_undoing;                    /* is it worth undoing list? */
  float           undo_frac;            /* fraction of vtxs indicating worth of undoing */
  int             step;                 /* loops through movements of vertices */
  bool             parity;               /* sort forwards or backwards? */
  bool             done;                 /* has termination criteria been achieved? */
  int             nbad;                 /* number of unhelpful passes in a row */
  int             npass;                /* total number of passes */
  int             nbadtries;            /* number of unhelpful passes before quitting */
  bool             enforce_balance;      /* force a balanced partition? */
  bool             enforce_balance_hard; /* really force a balanced partition? */
  bool             balance_trouble;      /* even balance_hard isn't working */
  int             size;                 /* array spacing */
  int             i, j, k, l;           /* loop counters */


  nbadtries = KL_NTRIES_BAD;

  enforce_balance      = false;
  temp_balanced        = false;
  enforce_balance_hard = false;
  balance_trouble      = false;

  size = (int)(&(listspace[0][1]) - &(listspace[0][0]));

  undo_frac = 0.3f;

  cut_cost = 1;
  if (term_wgts[1] != null) {
    if (CUT_TO_HOP_COST > 1) {
      cut_cost = CUT_TO_HOP_COST;
    }
  }

  bspace    = (int*)Marshal.AllocHGlobal(nvtxs * sizeof(int));
  weightsum = (double*)Marshal.AllocHGlobal(nsets * sizeof(double));
  locked    = (double*)Marshal.AllocHGlobal(nsets * sizeof(double));
  loose     = (double*)Marshal.AllocHGlobal(nsets * sizeof(double));

  if (bspace == null || weightsum == null || locked == null || loose == null) {
    Marshal.FreeHGlobal((IntPtr)loose);
    Marshal.FreeHGlobal((IntPtr)locked);
    Marshal.FreeHGlobal((IntPtr)weightsum);
    Marshal.FreeHGlobal((IntPtr)bspace);
    return true;
  }

  if (*bndy_list != null) {
    bdy_ptr     = *bndy_list;
    list_length = 0;
    while (*bdy_ptr != 0) {
      bspace[list_length++] = *bdy_ptr++;
    }
    Marshal.FreeHGlobal((IntPtr)(*bndy_list));

    if (list_length == 0) { /* No boundary -> make everybody bndy. */
      for (i = 0; i < nvtxs; i++) {
        bspace[i] = i + 1;
      }
      list_length = nvtxs;
    }
    /* Set dvals to flag uninitialized vertices. */
    for (i = 1; i <= nvtxs; i++) {
      dvals[i][0] = 3 * maxdval;
    }
  }
  else {
    list_length = nvtxs;
  }

  step_cutoff = KL_BAD_MOVES;
  cost_cutoff = maxdval * step_cutoff / 7;
  if (cost_cutoff < step_cutoff) {
    cost_cutoff = step_cutoff;
  }

  deltaminus = deltaplus = 0;
  for (i = 0; i < nsets; i++) {
    if (startweight[i] - goal[i] > deltaplus) {
      deltaplus = startweight[i] - goal[i];
    }
    else if (goal[i] - startweight[i] > deltaminus) {
      deltaminus = goal[i] - startweight[i];
    }
  }
  balanced = (deltaplus + deltaminus <= max_dev);

  bestg_min = -2.0 * nvtxs * maxdval;
  parity    = false;
  eweight   = (int)(cut_cost + .5);
  nbad      = 0;
  npass     = 0;
  improved  = 0;
  done      = false;
  while (!done) {
    npass++;
    ever_balanced = false;

    /* Initialize various quantities. */
    balance_best = 0;
    for (i = 0; i < nsets; i++) {
      for (j = 0; j < nsets; j++) {
        tops[i][j] = 2 * maxdval;
      }
      weightsum[i] = startweight[i];
      loose[i]     = weightsum[i];
      locked[i]    = 0;
      balance_best += goal[i];
    }

    gtotal   = 0;
    bestg    = bestg_min;
    beststep = -1;

    movelist = null;
    endlist  = &movelist;

    neg_steps = 0;

    /* Compute the initial d-values, and bucket-sort them. */
    time = seconds();
    if (nsets == 2) {
      bucketsorts_bi(graph, nvtxs, buckets, listspace, dvals, sets, term_wgts, maxdval, nsets,
                     parity, hops, bspace, list_length, npass, useEdgeWeights);
    }
    else {
      bucketsorts(graph, nvtxs, buckets, listspace, dvals, sets, term_wgts, maxdval, nsets, parity,
                  hops, bspace, list_length, npass, useEdgeWeights);
    }
    parity = !parity;
    kl_bucket_time += seconds() - time;

    if (DEBUG_KL == DebugFlagKL.PrintBucket) {
      pbuckets(buckets, listspace, maxdval, nsets);
    }

    /* Now determine the set of K-L moves. */

    for (step = 1;; step++) {

      /* Find the highest d-value in each set. */
      /* But only consider moves from large to small sets, or moves */
      /* in which balance is preserved. */
      /* Break ties in some nonarbitrary manner. */
      bestval = -maxdval - 1;
      for (i = 0; i < nsets; i++) {
        for (j = 0; j < nsets; j++) {
          /* Only allow moves from large sets to small sets, or */
          /* moves which preserve balance. */
          if (i != j) {
            /* Find the best move from i to j. */
            for (k = tops[i][j]; k >= 0 && buckets[i][j][k] == null; k--) {
              ;
            }
            tops[i][j] = k;

            if (k >= 0) {
              l       = (j > i) ? j - 1 : j;
              vtx     = ((int)(buckets[i][j][k] - listspace[l])) / size;
              vweight = graph[vtx]->vwgt;

              if ((enforce_balance_hard && weightsum[i] >= goal[i] && weightsum[j] <= goal[j] &&
                   weightsum[i] - goal[i] - (weightsum[j] - goal[j]) > max_dev) ||
                  (!enforce_balance_hard && weightsum[i] >= goal[i] && weightsum[j] <= goal[j]) ||
                  (!enforce_balance_hard &&
                   weightsum[i] - vweight - goal[i] > -(double)((max_dev + 1) / 2) &&
                   weightsum[j] + vweight - goal[j] < (double)((max_dev + 1) / 2))) {

                /* Is it the best move seen so far? */
                if (k - maxdval > bestval) {
                  bestval = k - maxdval;
                  bestvtx = vtx;
                  bestto  = j;
                  /* DO I NEED ALL THIS DATA?  Just to break ties. */
                  bestdelta = Math.Abs(weightsum[i] - vweight - goal[i]) +
                              Math.Abs(weightsum[j] + vweight - goal[j]);
                  beststuck1 = (int)Math.Min(loose[i], goal[j] - locked[j]);
                  beststuck2 = (int)Math.Max(loose[i], goal[j] - locked[j]);
                }

                else if (k - maxdval == bestval) {
                  /* Tied.  Is better balanced than current best? */
                  /* If tied, move among sets with most freedom. */
                  stuck1st = (int)Math.Min(loose[i], goal[j] - locked[j]);
                  stuck2nd = (int)Math.Max(loose[i], goal[j] - locked[j]);
                  delta    = Math.Abs(weightsum[i] - vweight - goal[i]) +
                             Math.Abs(weightsum[j] + vweight - goal[j]);

                  /* NOTE: Randomization in this check isn't ideal */
                  /* if more than two guys are tied. */
                  if (delta < bestdelta ||
                      (delta == bestdelta &&
                       (stuck1st > beststuck1 ||
                        (stuck1st == beststuck1 &&
                         (stuck2nd > beststuck2 ||
                          (stuck2nd == beststuck2 && (KL_RANDOM && drandom() < .5))))))) {
                    bestval    = k - maxdval;
                    bestvtx    = vtx;
                    bestto     = j;
                    bestdelta  = delta;
                    beststuck1 = stuck1st;
                    beststuck2 = stuck2nd;
                  }
                }
              }
            }
          }
        }
      }
      if (bestval == -maxdval - 1) { /* No allowed moves */
        if (DEBUG_KL != DebugFlagKL.NoDebugging) {
          Console.WriteLine("No KL moves at step {0:d}.  bestg = {1:g} at step {2:d}.\n", step, bestg, beststep);
        }
        break;
      }

      bestptr  = &(listspace[0][bestvtx]);
      bestfrom = sets[bestvtx];

      vweight = graph[bestvtx]->vwgt;
      weightsum[bestto] += vweight;
      weightsum[bestfrom] -= vweight;
      loose[bestfrom] -= vweight;
      locked[bestto] += vweight;

      if (enforce_balance) { /* Check if this partition is balanced. */
        deltaminus = deltaplus = 0;
        for (i = 0; i < nsets; i++) {
          if (weightsum[i] - goal[i] > deltaplus) {
            deltaplus = weightsum[i] - goal[i];
          }
          else if (goal[i] - weightsum[i] > deltaminus) {
            deltaminus = goal[i] - weightsum[i];
          }
        }
        balance_val   = deltaminus + deltaplus;
        temp_balanced = (balance_val <= max_dev);
        ever_balanced = (ever_balanced || temp_balanced);
      }

      gtotal += bestval;
      if (((gtotal > bestg && (!enforce_balance || temp_balanced)) ||
           (enforce_balance_hard && balance_val < balance_best)) &&
          step != nvtxs) {
        bestg    = gtotal;
        beststep = step;
        if (enforce_balance_hard) {
          balance_best = balance_val;
        }
        if (temp_balanced) {
          enforce_balance_hard = false;
        }
      }

      if (DEBUG_KL == DebugFlagKL.MoreInfo || DEBUG_KL == DebugFlagKL.PrintBucket)
      {
          Console.WriteLine("At KL step {0:d}, bestvtx={1:d}, bestval={2:d} ({3:d}-> {4:d})", step, bestvtx, bestval, bestfrom, bestto);
      }

      /* Monitor the stopping criteria. */
      if (bestval < 0) {
        if (!enforce_balance || ever_balanced) {
          neg_steps++;
        }
        if (bestg != bestg_min) {
          neg_cost = (int)(bestg - gtotal);
        }
        else {
          neg_cost = -maxdval - 1;
        }
        if ((neg_steps > step_cutoff || neg_cost > cost_cutoff) &&
            !(enforce_balance && bestg == bestg_min) && (beststep != step)) {
          if (DEBUG_KL != DebugFlagKL.NoDebugging) {
            if (neg_steps > step_cutoff) {
                Console.WriteLine("KL step cutoff at step {0:d}.  bestg = {1:g} at step {2:d}.", step, bestg, beststep);
            }
            else if (neg_cost > cost_cutoff) {
                Console.WriteLine("KL cost cutoff at step {0:d}.  bestg = {1:g} at step {2:d}.", step, bestg, beststep);
            }
          }
          break;
        }
      }
      else if (bestval > 0) {
        neg_steps = 0;
      }

      /* Remove vertex from its buckets, and flag it as finished. */
      l = 0;
      for (k = 0; k < nsets; k++) {
        if (k != bestfrom) {
          dval = dvals[bestvtx][l] + maxdval;
          removebilist(&listspace[l][bestvtx], &buckets[bestfrom][k][dval]);
          l++;
        }
      }

      /* Is there a better way to do this? */
      sets[bestvtx] = -sets[bestvtx] - 1;

      /* Set up the linked list of moved vertices. */
      bestptr->next = null;
      bestptr->prev = (bilist *)(ulong)bestto;
      *endlist      = bestptr;
      endlist       = &(bestptr->next);

      /* Now update the d-values of all the neighbors */
      edges = graph[bestvtx]->edges;
      if (useEdgeWeights) {
        ewptr = graph[bestvtx]->ewgts;
      }
      for (j = graph[bestvtx]->nedges - 1; j != 0; j--) {
        neighbor = *(++edges);
        if (useEdgeWeights) {
          eweight = (int)(*(++ewptr) * cut_cost + .5);
        }

        /* First make sure neighbor is alive. */
        if (sets[neighbor] >= 0) {
          group = sets[neighbor];

          if (dvals[neighbor][0] >= 3 * maxdval) {
            /* New vertex, not yet in buckets. */
            /* Can't be neighbor of moved vtx, so compute */
            /* initial dvals and buckets, then update. */
            bucketsort1(graph, neighbor, buckets, listspace, dvals, sets, term_wgts, maxdval, nsets,
                        hops, useEdgeWeights);
          }

          l = 0;
          for (k = 0; k < nsets; k++) {
            if (k != group) {
              diff = eweight * (hops[k][bestfrom] - hops[group][bestfrom] + hops[group][bestto] -
                                hops[k][bestto]);
              dval = dvals[neighbor][l] + maxdval;
              movebilist(&listspace[l][neighbor], &buckets[group][k][dval],
                         &buckets[group][k][dval + diff]);
              dvals[neighbor][l] += diff;
              dval += diff;
              if (dval > tops[group][k]) {
                tops[group][k] = dval;
              }
              l++;
            }
          }
        }
      }
      if (DEBUG_KL  == DebugFlagKL.PrintBucket) {
        pbuckets(buckets, listspace, maxdval, nsets);
      }
    }

    /* Done with a pass; should we actually perform any swaps? */
    bptr = movelist;
    if (bestg > 0 || (bestg != bestg_min && !balanced && enforce_balance) ||
        (bestg != bestg_min && balance_trouble)) {
      improved += (int)bestg;
      for (i = 1; i <= beststep; i++) {
        vtx    = ((int)(bptr - listspace[0])) / size;
        bestto = (int)(ulong)bptr->prev;
        startweight[bestto] += graph[vtx]->vwgt;
        startweight[-sets[vtx] - 1] -= graph[vtx]->vwgt;
        sets[vtx] = bestto;
        bptr      = bptr->next;
      }

      deltaminus = deltaplus = 0;
      for (i = 0; i < nsets; i++) {
        if (startweight[i] - goal[i] > deltaplus) {
          deltaplus = startweight[i] - goal[i];
        }
        else if (goal[i] - startweight[i] > deltaminus) {
          deltaminus = goal[i] - startweight[i];
        }
      }
      /*
      printf(" deltaplus = %f, deltaminus = %f, max_dev = %d\n", deltaplus, deltaminus, max_dev);
      */
      balanced = (deltaplus + deltaminus <= max_dev);
    }
    else {
      nbad++;
    }

    if (!balanced || bptr == movelist) {
      if (enforce_balance) {
        if (enforce_balance_hard) {
          balance_trouble = true;
        }
        enforce_balance_hard = true;
      }
      enforce_balance = true;
      nbad++;
    }

    worth_undoing = (step < undo_frac * nvtxs);
    done          = (nbad >= nbadtries && balanced);
    if (KL_MAX_PASS > 0) {
      done = done || (npass == KL_MAX_PASS && balanced);
    }
    if (!done) { /* Prepare for next pass. */
      if (KL_UNDO_LIST && worth_undoing && !balance_trouble) {
        /* Make a list of modified vertices for next bucketsort. */
        /* Also, ensure these vertices are removed from their buckets. */
        list_length =
            make_kl_list(graph, movelist, buckets, listspace, sets, nsets, bspace, dvals, maxdval);
      }
    }
    if (done || !(KL_UNDO_LIST && worth_undoing && !balance_trouble)) {
      /* Restore set numbers of remaining, altered vertices. */
      while (bptr != null) {
        vtx       = ((int)(bptr - listspace[0])) / size;
        sets[vtx] = -sets[vtx] - 1;
        bptr      = bptr->next;
      }
      list_length = nvtxs;
    }

    if (done && *bndy_list != null) {
      make_bndy_list(graph, movelist, buckets, listspace, sets, nsets, bspace, tops, bndy_list);
    }
  }

  if (DEBUG_KL != DebugFlagKL.NoDebugging) {
    Console.WriteLine("   KL required {0:d} passes to improve by {1:d}.", npass, improved);
  }

  Marshal.FreeHGlobal((IntPtr)loose);
  Marshal.FreeHGlobal((IntPtr)locked);
  Marshal.FreeHGlobal((IntPtr)weightsum);
  Marshal.FreeHGlobal((IntPtr)bspace);
  return false;
}

private static int make_kl_list(vtx_data **graph,     /* data structure for graph */
                 bilist *   movelist,  /* list of vtxs to be moved */
                 bilist ****buckets,   /* array of lists for bucket sort */
                 bilist **  listspace, /* list data structure for each vertex */
                 int *             sets,      /* processor each vertex is assigned to */
                 int               nsets,     /* number of sets divided into */
                 int *             bspace,    /* list of active vertices for bucketsort */
                 int **            dvals,     /* d-values for each transition */
                 int               maxdval    /* maximum d-value for a vertex */
)
{
  bilist **list;        /* bucket to erase element from */
  bilist * vptr;        /* loops through movelist */
  int *           bptr;        /* loops through bspace */
  int *           iptr;        /* loops through edge list */
  int             vtx;         /* vertex that was moved */
  int             neighbor;    /* neighbor of a vertex */
  int             myset;       /* set a vertex is in */
  int             newset;      /* loops through other sets */
  int             list_length; /* number of values put in bspace */
  int             size;        /* array spacing */
  int             i, l;        /* loop counter */

  /* First push all the moved vertices onto list, so they can be flagged. */
  /* They've already been removed from buckets, so want to avoid them. */
  size        = (int)(&(listspace[0][1]) - &(listspace[0][0]));
  vptr        = movelist;
  bptr        = bspace;
  list_length = 0;
  while (vptr != null) {
    vtx     = ((int)(vptr - listspace[0])) / size;
    *bptr++ = vtx;
    if (sets[vtx] >= 0) {
      sets[vtx] = -sets[vtx] - 1;
    }
    ++list_length;
    vptr = vptr->next;
  }

  /* Now look at all the neighbors of moved vertices. */
  vptr = movelist;
  while (vptr != null) {
    vtx = ((int)(vptr - listspace[0])) / size;

    iptr = graph[vtx]->edges;
    for (i = graph[vtx]->nedges - 1; i != 0; i--) {
      neighbor = *(++iptr);
      if (sets[neighbor] >= 0) {
        *bptr++ = neighbor;
        ++list_length;
        myset          = sets[neighbor];
        sets[neighbor] = -sets[neighbor] - 1;

        /* Remove neighbor entry from all his buckets. */
        /* Note: vertices in movelist already removed from buckets. */
        l = 0;
        for (newset = 0; newset < nsets; newset++) {
          if (newset != myset) {
            list = &buckets[myset][newset][dvals[neighbor][l] + maxdval];
            removebilist(&listspace[l][neighbor], list);
            l++;
          }
        }
      }
    }
    vptr = vptr->next;
  }

  /* Now that list is constructed, go reconstruct all the set numbers. */
  bptr = bspace;
  for (i = list_length; i!= 0; i--) {
    vtx       = *bptr++;
    sets[vtx] = -sets[vtx] - 1;
  }

  return (list_length);
}

/*static void p1bucket();*/

private static void pbuckets(bilist ****buckets,   /* pointers to bucket lists */
bilist **  listspace, /* elements within buckets */
int               maxdeg,    /* maximum degree of a vertex */
int               nsets      /* number of sets being divided into */
)
{
    bilist *lptr; /* points to correct listspace */
    int            i, j; /* loop counter */

    Console.WriteLine();
    for (i = 0; i < nsets; i++) {
        for (j = 0; j < nsets; j++) {
            if (i != j) {
                Console.WriteLine("For transition {0:d} -> {1:d}", i, j);
                if (j > i) {
                    lptr = listspace[j - 1];
                }
                else {
                    lptr = listspace[j];
                }
                p1bucket(buckets[i][j], lptr, maxdeg);
                Console.WriteLine();
            }
        }
    }
    Console.WriteLine();
}

private static void p1bucket(bilist **bucket, /* buckets holding bucket list */
bilist * lptr,   /* elements within bucket */
int             maxdeg  /* maximum degree of a vertex */
)
{
    bilist *bptr; /* loops through list at a bucket */
    int            val;  /* element in a bucket */
    int            size; /* array spacing */
    int            i;    /* loop counter */

    size = (int)(&(lptr[1]) - &(lptr[0]));
    for (i = 2 * maxdeg; i >= 0; i--) {
        if (bucket[i] != null) {
            Console.Write("  Bucket {0:d}:", i - maxdeg);
            for (bptr = bucket[i]; bptr != null; bptr = bptr->next) {
                val = ((int)(bptr - lptr)) / size;
                Console.Write(" {0:d}", val);
            }

            Console.WriteLine();
        }
    }
}

/* Note: bi-directional lists aren't assumed to be sorted. */
private static void add2bilist(                      /* add val to unsorted list */
bilist * lptr, /* element to add */
bilist **list  /* list added to */
)
{
    lptr->next = *list;
    if (*list != null) {
        (*list)->prev = lptr;
    }
    lptr->prev = null;
    *list      = lptr;
}

private static void removebilist(bilist * lptr, /* ptr to element to remove */
bilist **list  /* head of list to remove it from */
)

/* Remove an element from a bidirectional list. */
{
    if (lptr->next != null) {
        lptr->next->prev = lptr->prev;
    }
    if (lptr->prev != null) {
        lptr->prev->next = lptr->next;
    }
    else {
        *list = lptr->next;
    }
}

private static void movebilist(bilist * lptr,    /* ptr to element to move */
bilist **oldlist, /* head of list to remove it from */
bilist **newlist  /* head of list to add it to */
)

/* Move an element from a old bidirectional list to new one. */
{
    removebilist(lptr, oldlist);

    add2bilist(lptr, newlist);
}
    }
}
