using static ChacoSharp.StaticConstants;
using static ChacoSharp.Utilities.Randomize;
using static ChacoSharp.Coarsening.BiListHelper;

namespace ChacoSharp.Coarsening
{
    public static unsafe class BucketSort
    {
        /* Idea:
   'buckets[i][j]' is a set of buckets to sort moves from i to j.
   listspace[i] is space for lists in buckets[i][j].
   Loop through all nonequal pairs [i][j], taking the first element
   in each list.  Compare them all to find the largest allowed move.
   Make that move, and save it in movelist.
*/

public static void bucketsorts(vtx_data **graph,       /* graph data structure */
                 int               nvtxs,       /* number of vertices */
                 bilist ****buckets,     /* array of lists for bucket sort */
                 bilist **  listspace,   /* list data structure for each vertex */
                 int **            dvals,       /* d-values for each vertex for removing */
                 int *             sets,        /* processor each vertex is assigned to */
                 float *[]           term_wgts, /* weights for terminal prapogation */
                 int               maxdval,     /* maximum possible dvalue for a vertex */
                 int               nsets,       /* number of sets being divided into */
                 bool               parity,      /* work in forward or backward direction? */
                 int [][] hops,          /* hop cost between sets */
                 int *bspace,                   /* indices for randomly ordering vtxs */
                 int  list_length,              /* number of values in bspace to work with */
                 int  npass,                    /* which pass through KL is this? */
                 bool  useEdgeWeights               /* are edge weights being used? */
)
{
  bilist **bptr  = null;    /* loops through set of buckets */
  bilist * lptr  = null;    /* pointer to an element in listspace */
  float *         ewptr = null;    /* loops through edge weights */
  int *           bsptr = null;    /* loops through bspace */
  int *           edges = null;    /* edge list for a vertex */
  int             myset;           /* set that current vertex belongs to */
  int             newset;          /* set current vertex could move to */
  int             set;             /* set that neighboring vertex belongs to */
  int             weight;          /* edge weight for a particular edge */
  int             vtx;             /* vertex in graph */
  float           tval;            /* terminal propagation value */
  int             val;             /* terminal propagation rounded value */
  double          cut_cost;        /* relative cut/hop importance */
  double          hop_cost;        /* relative hop/cut importance */
  int             myhop;           /* hops associated with current vertex */
  int             i, j, l;         /* loop counters */

  /* For each vertex, compute d-values for each possible transition. */
  /* Then store them in each appropriate bucket. */

  if (npass == 1 || !KL_UNDO_LIST || list_length == nvtxs) {
    /* Empty all the buckets. */
    /* Last clause catches case where lists weren't undone. */
    bptr = buckets[0][1];
    for (i = nsets * (nsets - 1) * (2 * maxdval + 1); i!= 0; i--) {
      *bptr++ = null;
    }
  }

  /* Randomize the order of the vertices */

  if (list_length == nvtxs || !KL_UNDO_LIST) {
    bsptr       = bspace;
    list_length = nvtxs;
    if (parity) {
      for (i = 1; i <= nvtxs; i++) {
        *bsptr++ = i;
      }
    }
    else {
      for (i = nvtxs; i != 0; i--) {
        *bsptr++ = i;
      }
    }
  }
  if (KL_RANDOM) {
    randomize(bspace - 1, list_length);
  }

  /* Now compute d-vals by seeing which sets neighbors belong to. */
  cut_cost = hop_cost = 1;
  if (term_wgts[1] != null) {
    if (CUT_TO_HOP_COST > 1) {
      cut_cost = CUT_TO_HOP_COST;
    }
    else {
      hop_cost = 1.0 / CUT_TO_HOP_COST;
    }
  }
  weight = (int)(cut_cost + .5);
  bsptr  = bspace;
  for (i = 0; i < list_length; i++) { /* Loop through vertices. */
    vtx   = *bsptr++;
    myset = sets[vtx];

    /* Initialize all the preference values. */
    if (term_wgts[1] != null) {
      /* Using terminal propagation. */
      if (myset == 0) { /* No terminal value. */
        for (newset = 1; newset < nsets; newset++) {
          tval = (term_wgts[newset])[vtx];
          if (tval < 0) {
            val = (int)(-tval * hop_cost + .5);
            val = -val;
          }
          else {
            val = (int)(tval * hop_cost + .5);
          }
          dvals[vtx][newset - 1] = val;
        }
      }
      else {
        tval = -(term_wgts[myset])[vtx];
        if (tval < 0) {
          val = (int)(-tval * hop_cost + .5);
          val = -val;
        }
        else {
          val = (int)(tval * hop_cost + .5);
        }
        dvals[vtx][0] = val;
        l             = 1;
        for (newset = 1; newset < nsets; newset++) {
          if (newset != myset) {
            tval = (term_wgts[newset])[vtx] - (term_wgts[myset])[vtx];
            if (tval < 0) {
              val = (int)(-tval * hop_cost + .5);
              val = -val;
            }
            else {
              val = (int)(tval * hop_cost + .5);
            }
            dvals[vtx][l] = val;
            l++;
          }
        }
      }
    }
    else {
      for (j = 0; j < nsets - 1; j++) {
        dvals[vtx][j] = 0;
      }
    }

    /* First count the neighbors in each set. */
    edges = graph[vtx]->edges;
    if (useEdgeWeights) {
      ewptr = graph[vtx]->ewgts;
    }
    for (j = graph[vtx]->nedges - 1; j != 0; j--) {
      set = sets[*(++edges)];
      if (set < 0) {
        set = -set - 1;
      }
      if (useEdgeWeights) {
        weight = (int)(*(++ewptr) * cut_cost + .5);
      }
      myhop = hops[myset][set];

      l = 0;
      for (newset = 0; newset < nsets; newset++) {
        if (newset != myset) {
          dvals[vtx][l] += weight * (myhop - hops[newset][set]);
          l++;
        }
      }
    }

    /* Now add to appropriate buckets. */
    l = 0;
    for (newset = 0; newset < nsets; newset++) {
      if (newset != myset) {
        lptr = listspace[l];
        add2bilist(&lptr[vtx], &buckets[myset][newset][dvals[vtx][l] + maxdval]);
        ++l;
      }
    }
  }
}

        public static void bucketsortsv(vtx_data** graph, /* graph data structure */
            int nvtxs, /* number of vertices */
            bilist** lbuckets, /* array of lists for bucket sort */
            bilist** rbuckets, /* array of lists for bucket sort */
            bilist* llistspace, /* list data structure for each vertex */
            bilist* rlistspace, /* list data structure for each vertex */
            int* ldvals, /* d-values for each vertex for removing */
            int* rdvals, /* d-values for each vertex for removing */
            int* sets, /* processor each vertex is assigned to */
            int maxdval, /* maximum possible dvalue for a vertex */
            bool parity, /* work in forward or backward direction? */
            int* bspace, /* indices for randomly ordering vtxs */
            int list_length /* number of values in bspace to work with */
        )
        {
            bilist** lbptr; /* loops through set of buckets */
            bilist** rbptr; /* loops through set of buckets */
            int* bsptr; /* loops through bspace */
            int* edges; /* edge list for a vertex */
            int left_weight; /* my neighbors in 0 set */
            int right_weight; /* my neighbors in 1 set */
            int vtx; /* vertex in graph */
            int neighbor; /* neighbor of vertex */
            int set; /* set that neighboring vertex belongs to */
            int i, j; /* loop counters */

            /* For each vertex, compute d-values and store in buckets. */

            /* Empty all the buckets. */
            rbptr = lbuckets;
            lbptr = rbuckets;
            for (i = (2 * maxdval + 1); i != 0; i--)
            {
                *lbptr++ = null;
                *rbptr++ = null;
            }

            /* Randomize the order of the vertices */

            if ((KL_UNDO_LIST && list_length == nvtxs) || (!KL_UNDO_LIST && !KL_RANDOM) || list_length == 0)
            {
                /* Don't need to reoder if about to randomize. */
                list_length = nvtxs;
                bsptr = bspace;
                if (parity)
                {
                    for (i = 1; i <= nvtxs; i++)
                    {
                        *bsptr++ = i;
                    }
                }
                else
                {
                    for (i = nvtxs; i != 0; i--)
                    {
                        *bsptr++ = i;
                    }
                }
            }

            if (KL_RANDOM)
            {
                randomize(bspace - 1, list_length);
            }

            /* Now compute d-vals by seeing which sets neighbors belong to. */

            bsptr = bspace;
            for (i = 0; i < list_length; i++)
            {
                /* Loop through vertices. */
                vtx = *bsptr++;
                if (sets[vtx] == 2)
                {
                    /* Only worry about separator vertices. */

                    /* Initialize all the preference values. */
                    left_weight = right_weight = 0;

                    /* First count the neighbors in each set. */
                    edges = graph[vtx]->edges;
                    for (j = graph[vtx]->nedges - 1; j != 0; j--)
                    {
                        neighbor = *(++edges);
                        set = sets[neighbor];
                        if (set < 0)
                        {
                            set = -set - 1;
                        }

                        if (set == 0)
                        {
                            left_weight += graph[neighbor]->vwgt;
                        }
                        else if (set == 1)
                        {
                            right_weight += graph[neighbor]->vwgt;
                        }
                    }

                    ldvals[vtx] = graph[vtx]->vwgt - right_weight;
                    rdvals[vtx] = graph[vtx]->vwgt - left_weight;

                    /* Now add to appropriate buckets. */
                    add2bilist(&llistspace[vtx], &lbuckets[ldvals[vtx] + maxdval]);
                    add2bilist(&rlistspace[vtx], &rbuckets[rdvals[vtx] + maxdval]);
                }
            }
        }

        /* Idea:
   'buckets[i][j]' is a set of buckets to sort moves from i to j.
   listspace[i] is space for lists in buckets[i][j].
   Loop through all nonequal pairs [i][j], taking the first element
   in each list.  Compare them all to find the largest allowed move.
   Make that move, and save it in movelist.
*/

/* Routine slightly streamlined for case with only two sets. */

public static void bucketsorts_bi(vtx_data **graph,       /* graph data structure */
                    int               nvtxs,       /* number of vertices */
                    bilist ****buckets,     /* array of lists for bucket sort */
                    bilist **  listspace,   /* list data structure for each vertex */
                    int **            dvals,       /* d-values for each vertex for removing */
                    int *             sets,        /* processor each vertex is assigned to */
                    float *[]           term_wgts, /* weights for terminal propagation */
                    int               maxdval,     /* maximum possible dvalue for a vertex */
                    int               nsets,       /* number of sets being divided into */
                    bool               parity,      /* work in forward or backward direction? */
                    int [][] hops,          /* hop cost between sets */
                    int *bspace,                   /* indices for randomly ordering vtxs */
                    int  list_length,              /* number of values in bspace to work with */
                    int  npass,                    /* which pass through KL is this? */
                    bool  useEdgeWeights               /* are edge weights being used? */
)
{
  bilist **bptr  = null;    /* loops through set of buckets */
  bilist * lptr  = null;    /* pointer to an element in listspace */
  float *         ewptr = null;    /* loops through edge weights */
  float *         twptr = null;    /* weights for terminal propagation */
  int *           bsptr = null;    /* loops through bspace */
  int *           edges = null;    /* edge list for a vertex */
  int             myset;           /* set current vertex belongs to */
  int             other_set;       /* set current vertex doesn't belong to */
  int             set;             /* set that neighboring vertex belongs to */
  int             weight;          /* edge weight for a particular edge */
  int             vtx;             /* vertex in graph */
  int             val;             /* terminal propagation rounded value */
  double          cut_cost;        /* relative cut/hop importance */
  double          hop_cost;        /* relative hop/cut importance */
  int             myhop;           /* hops associated with current vertex */
  int             i, j;            /* loop counters */

  /* For each vertex, compute d-values for each possible transition. */
  /* Then store them in each appropriate bucket. */

  if (npass == 1 || !KL_UNDO_LIST || list_length == nvtxs) {
    /* Empty all the buckets. */
    /* Last clause catches case where lists weren't undone. */
    bptr = buckets[0][1];
    for (i = nsets * (nsets - 1) * (2 * maxdval + 1); i != 0; i--) {
      *bptr++ = null;
    }
  }

  /* Randomize the order of the vertices */

  if ((KL_UNDO_LIST && list_length == nvtxs) || !KL_UNDO_LIST) {
    list_length = nvtxs;
    bsptr       = bspace;
    if (parity) {
      for (i = 1; i <= nvtxs; i++) {
        *bsptr++ = i;
      }
    }
    else {
      for (i = nvtxs; i!= 0; i--) {
        *bsptr++ = i;
      }
    }
  }
  if (KL_RANDOM) {
    randomize(bspace - 1, list_length);
  }

  /* Now compute d-vals by seeing which sets neighbors belong to. */
  cut_cost = hop_cost = 1;
  if (term_wgts[1] != null) {
    if (CUT_TO_HOP_COST > 1) {
      cut_cost = CUT_TO_HOP_COST;
    }
    else {
      hop_cost = 1.0 / CUT_TO_HOP_COST;
    }
  }

  weight = (int)(cut_cost + .5);

  bsptr = bspace;
  twptr = term_wgts[1];
  for (i = 0; i < list_length; i++) { /* Loop through vertices. */
    vtx       = *bsptr++;
    myset     = sets[vtx];
    other_set = myset == 0 ? 1 : 0;

    /* Initialize all the preference values. */
    if (twptr != null) {
      /* Using terminal propagation.  Round to integer value. */
      if (twptr[vtx] < 0) {
        val = (int)(-twptr[vtx] * hop_cost + .5);
        val = -val;
      }
      else {
        val = (int) (twptr[vtx] * hop_cost + .5);
      }
      if (myset == 0) {
        dvals[vtx][0] = val;
      }
      else {
        dvals[vtx][0] = -val;
      }
    }
    else {
      dvals[vtx][0] = 0;
    }

    /* First count the neighbors in each set. */
    edges = graph[vtx]->edges;
    if (useEdgeWeights) {
      ewptr = graph[vtx]->ewgts;
    }
    for (j = graph[vtx]->nedges - 1; j != 0; j--) {
      set = sets[*(++edges)];
      if (set < 0) {
        set = -set - 1;
      }
      if (useEdgeWeights) {
        weight = (int)((*(++ewptr)) * cut_cost + .5);
      }
      myhop = hops[myset][set];

      dvals[vtx][0] += weight * (myhop - hops[other_set][set]);
    }

    /* Now add to appropriate buckets. */
    lptr = listspace[0];
    add2bilist(&lptr[vtx], &buckets[myset][other_set][dvals[vtx][0] + maxdval]);
  }
}

/* Idea:
   'buckets[i][j]' is a set of buckets to sort moves from i to j.
   listspace[i] is space for lists in buckets[i][j].
   Loop through all nonequal pairs [i][j], taking the first element
   in each list.  Compare them all to find the largest allowed move.
   Make that move, and save it in movelist.
*/

public static void bucketsort1(vtx_data **graph,       /* graph data structure */
                 int               vtx,         /* vertex being added to lists */
                 bilist ****buckets,     /* array of lists for bucket sort */
                 bilist **  listspace,   /* list data structure for each vertex */
                 int **            dvals,       /* d-values for each vertex for removing */
                 int *             sets,        /* processor each vertex is assigned to */
                 float *[]           term_wgts, /* weights for terminal propagation */
                 int               maxdval,     /* maximum possible dvalue for a vertex */
                 int               nsets,       /* number of sets being divided into */
                 int [][] hops,          /* hop cost between sets */
                 bool useEdgeWeights                /* are edge weights being used? */
)
{
  bilist *lptr  = null;    /* pointer to an element in listspace */
  float *        ewptr = null;    /* loops through edge weights */
  int *          edges = null;    /* edge list for a vertex */
  int            myset;           /* set that current vertex belongs to */
  int            newset;          /* set current vertex could move to */
  int            set;             /* set that neighboring vertex belongs to */
  int            weight;          /* edge weight for a particular edge */
  float          tval;            /* terminal propagation value */
  int            val;             /* terminal propagation rounded value */
  double         cut_cost;        /* relative cut/hop importance */
  double         hop_cost;        /* relative hop/cut importance */
  int            myhop;           /* hops associated with current vertex */
  int            j, l;            /* loop counters */

  /* Compute d-vals by seeing which sets neighbors belong to. */
  cut_cost = hop_cost = 1;
  if (term_wgts[1] != null) {
    if (CUT_TO_HOP_COST > 1) {
      cut_cost = CUT_TO_HOP_COST;
    }
    else {
      hop_cost = 1.0 / CUT_TO_HOP_COST;
    }
  }

  weight = (int)(cut_cost + .5);
  myset  = sets[vtx];

  /* Initialize all the preference values. */
  if (term_wgts[1] != null) {
    /* Using terminal propagation. */
    if (myset == 0) { /* No terminal value. */
      for (newset = 1; newset < nsets; newset++) {
        tval = (term_wgts[newset])[vtx];
        if (tval < 0) {
          val = (int)(-tval * hop_cost + .5);
          val = -val;
        }
        else {
          val = (int)(tval * hop_cost + .5);
        }
        dvals[vtx][newset - 1] = val;
      }
    }
    else {
      tval = -(term_wgts[myset])[vtx];
      if (tval < 0) {
        val = (int)(-tval * hop_cost + .5);
        val = -val;
      }
      else {
        val = (int)(tval * hop_cost + .5);
      }
      dvals[vtx][0] = val;
      l             = 1;
      for (newset = 1; newset < nsets; newset++) {
        if (newset != myset) {
          tval = (term_wgts[newset])[vtx] - (term_wgts[myset])[vtx];
          if (tval < 0) {
            val = (int)(-tval * hop_cost + .5);
            val = -val;
          }
          else {
            val = (int)(tval * hop_cost + .5);
          }
          dvals[vtx][l] = val;
          l++;
        }
      }
    }
  }
  else {
    for (j = 0; j < nsets - 1; j++) {
      dvals[vtx][j] = 0;
    }
  }

  /* First count the neighbors in each set. */
  edges = graph[vtx]->edges;
  if (useEdgeWeights) {
    ewptr = graph[vtx]->ewgts;
  }
  for (j = graph[vtx]->nedges - 1; j != 0; j--) {
    set = sets[*(++edges)];
    if (set < 0) {
      set = -set - 1;
    }
    if (useEdgeWeights) {
      weight = (int)(*(++ewptr) * cut_cost + .5);
    }
    myhop = hops[myset][set];

    l = 0;
    for (newset = 0; newset < nsets; newset++) {
      if (newset != myset) {
        dvals[vtx][l] += weight * (myhop - hops[newset][set]);
        l++;
      }
    }
  }

  /* Now add to appropriate buckets. */
  l = 0;
  for (newset = 0; newset < nsets; newset++) {
    if (newset != myset) {
      lptr = listspace[l];
      add2bilist(&lptr[vtx], &buckets[myset][newset][dvals[vtx][l] + maxdval]);
      ++l;
    }
  }
}
    }
}
