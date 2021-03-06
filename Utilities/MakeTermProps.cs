﻿using System;
using System.Runtime.InteropServices;
using static ChacoSharp.StaticConstants;

namespace ChacoSharp.Utilities
{
    public static unsafe class MakeTermProps
    {
        /* Compute the terminal constraints for next partition. */
public static void make_term_props(vtx_data **graph,        /* data structure for graph */
                     int               sub_nvtxs,    /* number of vtxs in subgraph */
                     int *             loc2glob,     /* mapping from subgraph to graph */
                     int *             assignment,   /* set for each vertex */
                     int               architecture, /* 0 => hypercube, 1 => mesh */
                     int               ndims_tot,    /* total hypercube dimensions */
                     int               ndims,        /* number of dimensions at this step */
                     set_info * set_info,     /* data about all the sets */
                     int               setnum,       /* number of set being divided */
                     int               nsets,        /* number of subsets being created */
                     int               set_max,      /* largest set created so far */
                     int []             subsets,      /* subsets being created */
                     float *[]           term_wgts,  /* set of terminal weights for each vertex */
                     bool               useEdgeWeights   /* are edge weights being used? */
)
{
  double[] term_wgt = new double[MAXSETS]; /* terminal weights */
  float *twptr;             /* one of the term_wgts vectors */
  float *[] dists = new float*[MAXSETS];    /* distances from my subsets to other sets */
  float *dist;              /* one of the dists arrays */
  float  edge_wgt;          /* weight of an edge */
  int    vtx;               /* vertex number */
  int    neighbor;          /* neighboring vertex number */
  int    neighbor_setnum;   /* set neighboring vertex is in */
  int    i, j, k;           /* loop counters */

  /* First compute average distance between my subsets and all other sets. */
  dist = (float*)Marshal.AllocHGlobal(nsets * (set_max + 1) * sizeof(float));
  for (i = 0; i < nsets; i++) {
    dists[i] = dist;
    dist += set_max + 1;
  }

  for (k = 0; k < MAXSETS; k++) {
    term_wgt[k] = 0;
  }

  if (architecture == 0) {
    avg_dists_cube(ndims_tot, ndims, set_info, nsets, set_max, subsets, dists);
  }
  else if (architecture > 0) {
    avg_dists_mesh(architecture, set_info, nsets, set_max, subsets, dists);
  }

  edge_wgt = 1;
  for (i = 1; i <= sub_nvtxs; i++) {
    for (k = 1; k < nsets; k++) {
      term_wgt[k] = 0;
    }

    vtx = loc2glob[i];

    for (j = 1; j < graph[vtx]->nedges; j++) {
      neighbor        = graph[vtx]->edges[j];
      neighbor_setnum = assignment[neighbor];
      if (neighbor_setnum != setnum) {
        if (useEdgeWeights) {
          edge_wgt = graph[vtx]->ewgts[j];
        }
        for (k = 1; k < nsets; k++) {
          dist = dists[k];
          term_wgt[k] += edge_wgt * dist[neighbor_setnum];
        }
      }
    }

    for (k = 1; k < nsets; k++) {
      twptr    = term_wgts[k];
      twptr[i] = (float)term_wgt[k];
    }
  }

  Marshal.FreeHGlobal((IntPtr)dists[0]);
}

static void avg_dists_cube(int              ndims_tot, /* total number of hypercube dimensions */
                           int              ndims,     /* number of dimensions created this step */
                           set_info *set_info,  /* data about all the sets */
                           int              nsets,     /* number of subsets being created */
                           int              set_max,   /* largest set created so far */
                           int []            subsets,   /* subsets being created */
                           float *[] dists/*[MAXSETS]*/       /* distances from my subsets to other sets */
)
{
  float *dist0;      /* first of dists vectors */
  float *dist;       /* one of dists vectors */
  int    ndims_old;  /* hypercube dimensions not relevant */
  int    ndims_left; /* hypercube dimensions left to do */
  int    myset;      /* subset being analyzed */
  int    start;      /* bit difference between two sets */
  int    val;        /* number of differing bits */
  int    set;        /* loops through all other sets */
  int    i;          /* loop counter */

  /* First compute distances for subset 0. */
  myset      = subsets[0];
  dist0      = dists[0];
  ndims_left = set_info[myset].ndims;
  ndims_old  = ndims_tot - ndims_left - ndims;
  for (set = 0; set < set_max; set++) {
    if (set_info[set].ndims >= 0) {
      val = 0;
      if (ndims_left == set_info[set].ndims) {
        start = (myset ^ set) >> ndims_old;
        while (start != 0) {
          if ((start & 1) != 0) {
            val++;
          }
          start >>= 1;
        }
      }
      dist0[set] = val;
    }
  }

  /* Now compute all distances relative to subset 0. */
  for (i = 1; i < nsets; i++) {
    myset = subsets[i];
    dist  = dists[i];

    for (set = 0; set < set_max; set++) {
      if (set_info[set].ndims >= 0) {
        val = 0;
        if (ndims_left == set_info[set].ndims) {
          start = (myset ^ set) >> ndims_old;
          while (start != 0) {
            if ((start & 1) != 0) {
              val++;
            }
            start >>= 1;
          }
        }
        /* Note: this is net preference for set over set 0. */
        dist[set] = dist0[set] - val;
      }
    }
  }
}

static void avg_dists_mesh(int              architecture, /* dimensions of mesh */
                           set_info *set_info,     /* data about all the sets */
                           int              nsets,        /* number of subsets being created */
                           int              set_max,      /* largest set created so far */
                           int []            subsets,      /* subsets being created */
                           float *[] dists/*[MAXSETS]*/ /* distances from my subsets to other sets */
)
{
  float *dist0; /* first of dists vectors */
  float *dist;  /* one of dists vectors */
  double val;   /* distance from subset to set */
  double sep;   /* distance between two subsets */
  int    set;   /* loops through all other sets */
  int    i;     /* loop counter */

  /* First compute distances for subset 0. */
  dist0 = dists[0];

  for (set = 0; set < set_max; set++) {
    if (set_info[set].span[0] >= 0) {
      val        = avg_dist_mesh(&set_info[subsets[0]], &set_info[set], architecture);
      dist0[set] = (float)val;
    }
  }

  /* Now compute all distances relative to subset 0. */
  for (i = 1; i < nsets; i++) {
    dist = dists[i];
    sep  = avg_dist_mesh(&set_info[subsets[i]], &set_info[subsets[0]], architecture);

    for (set = 0; set < set_max; set++) {
      if (set_info[set].span[0] >= 0) {
        val = avg_dist_mesh(&set_info[subsets[i]], &set_info[set], architecture);
        /* Note: this is net preference for set over 0. */
        dist[set] = (float)((dist0[set] - val) / sep);
      }
    }
  }
}

/* Compute the average distance between two subsets of mesh processors. */
private static double avg_dist_mesh(set_info *set1,        /* data about all first set */
                            set_info *set2,        /* data about all second set */
                            int              architecture /* dimension of mesh */
)
{
  double val; /* distance returned */
  int    i;   /* loop counter */

  val = 0;

  for (i = 0; i < architecture; i++) {
    val += avg_dist_interval(set1->low[i], set1->span[i], set2->low[i], set2->span[i]);
  }

  return (val);
}

/* Compute the average distance between two intervals */
private static double avg_dist_interval(int set1_low,  /* lowest point for first interval */
                                int set1_span, /* highest point for first interval */
                                int set2_low,  /* lowest point for second interval */
                                int set2_span  /* highest point for second interval */
)
{
  double set1_high; /* length of first interval */
  double set1_avg;  /* average value in first interval */
  double set2_high; /* length of second interval */
  double set2_avg;  /* average value in second interval */
  double val;       /* average distance between intervals */

  val       = 0;
  set1_high = set1_low + set1_span - 1;
  set1_avg  = .5 * (set1_high + set1_low);
  set2_high = set2_low + set2_span - 1;
  set2_avg  = .5 * (set2_high + set2_low);

  if (set1_low > set2_high || set2_low > set1_high) {
    val = Math.Abs(set1_avg - set2_avg);
  }

  else {
    if (set1_high > set2_high) {
      val += .5 * (set2_high - set2_low + 1) * (set1_high - set2_high) * (set1_high - set2_low + 1);
      set1_high = set2_high;
    }
    else if (set2_high > set1_high) {
      val += .5 * (set1_high - set1_low + 1) * (set2_high - set1_high) * (set2_high - set1_low + 1);
      set2_high = set1_high;
    }
    if (set1_low < set2_low) {
      val += .5 * (set2_high - set2_low + 1) * (set2_low - set1_low) * (set2_high - set1_low + 1);
      set1_low = set2_low;
    }
    else if (set2_low < set1_low) {
      val += .5 * (set1_high - set1_low + 1) * (set1_low - set2_low) * (set1_high - set2_low + 1);
    }
    val += (set1_high - set1_low) * (set1_high - set1_low + 1) * (set1_high - set1_low + 2) / 3.0;

    val /= set1_span * set2_span;
  }

  return (val);
}
    }
}
