using System;
using System.Diagnostics;
using System.Runtime.InteropServices;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Graph.CountWeights;
using static ChacoSharp.Graph.SimplePartition;
using static ChacoSharp.Graph.FreeGraph;
using static ChacoSharp.Graph.FindMaxDeg;
using static ChacoSharp.Coarsening.Coarsen1;
using static ChacoSharp.Coarsening.KlSpiff;
using static ChacoSharp.Coarsening.KlvSpiff;
using static ChacoSharp.Coarsening.BpmImprove;


namespace ChacoSharp.Coarsening
{
    public static unsafe class CoarsenKlv
    {

/* Once, or iteratively? */

public static void coarsen_klv(
    /* Coarsen until nvtxs < vmax, compute and uncoarsen. */
    vtx_data **graph,        /* array of vtx data for graph */
    int               nvtxs,        /* number of vertices in graph */
    int               nedges,       /* number of edges in graph */
    bool               using_vwgts,  /* are vertices weights being used? */
    bool               useEdgeWeights,  /* are edge weights being used? */
    float *[]           term_wgts,  /* weights for terminal propagation */
    int               igeom,        /* dimension for geometric information */
    float **          coords,       /* coordinates for vertices */
    int               vwgt_max,     /* largest vertex weight */
    int *             assignment,   /* processor each vertex gets assigned to */
    double []          goal,         /* desired set sizes */
    int               architecture, /* 0 => hypercube, d => d-dimensional mesh */
    int [][] hops,           /* cost of edge between sets */
    LanczosType     solver_flag,            /* which eigensolver to use */
    int     ndims,                  /* number of eigenvectors to calculate */
    int     nsets,                  /* number of sets being divided into */
    int     vmax,                   /* largest subgraph to stop coarsening */
    MappingType     mediantype,             /* flag for different assignment strategies */
    bool     mkconnected,            /* enforce connectivity before eigensolve? */
    double  eigtol,                 /* tolerance in eigen calculation */
    int     nstep,                  /* number of coarsenings between RQI steps */
    int     step,                   /* current step number */
    int **  pbndy_list,             /* pointer to returned boundary list */
    double[] weights,                /* weights of vertices in each set */
    bool    give_up                 /* has coarsening bogged down? */
)
{
  vtx_data **cgraph = null;                 /* array of vtx data for coarsened graph */
  double[]            new_goal = new double[MAXSETS];      /* new goal if not using vertex weights */
  double []          real_goal;              /* chooses between goal and new_goal */
  double            total_weight;           /* total weight of vertices */
  double            goal_weight;            /* total weight of vertices in goal */
  float *[]           cterm_wgts = new float*[MAXSETS];    /* terminal weights for coarse graph */
  float *[]           new_term_wgts = new float*[MAXSETS]; /* modified for Bui's method */
  float *[]          real_term_wgts;         /* which of previous two to use */
  float             ewgt_max;               /* largest edge weight in graph */
  float *           twptr      = null;      /* loops through term_wgts */
  float *           twptr_save = null;      /* copy of twptr */
  float *           ctwptr;                 /* loops through cterm_wgts */
  float **          ccoords;                /* coarse graph coordinates */
  int *             v2cv = null;                   /* mapping from vtxs to coarse vtxs */
  bool *             flag;                   /* scatter array for coarse bndy vtxs */
  int *             bndy_list;              /* list of vertices on boundary */
  int *             cbndy_list = null;             /* list of vertices of coarse graph on boundary */
  int *             cassignment;            /* set assignments for coarsened vertices */
  bool               flattened;              /* was this graph flattened? */
  int               list_length;            /* length of boundary vtx list */
  int               cnvtxs = 0;                 /* number of vertices in coarsened graph */
  int               cnedges = 0;                /* number of edges in coarsened graph */
  int               cvwgt_max;              /* largest vertex weight in coarsened graph */
  int               nextstep;               /* next step in RQI test */
  int               max_dev;                /* largest allowed deviation from balance */
  int               i, j;                   /* loop counters */

  if (DEBUG_COARSEN || DEBUG_TRACE) {
      Trace.WriteLine($"<Entering coarsen_kl, {nameof(step)}={step:d}, {nameof(nvtxs)}={nvtxs:d}, {nameof(nedges)}={nedges:d}, {nameof(vmax)}={vmax:d}>\n");
  }

  /* Is problem small enough to solve? */
  if (nvtxs <= vmax || give_up) {
    real_goal = goal;

    simple_part(graph, nvtxs, assignment, nsets, MappingType.MinCost, real_goal);
    list_length = find_bndy(graph, nvtxs, assignment, 2, &bndy_list);

    count_weights(graph, nvtxs, assignment, nsets + 1, weights, true);

    max_dev      = (step == 0) ? vwgt_max : 5 * vwgt_max;
    total_weight = 0;
    for (i = 0; i < nsets; i++) {
      total_weight += real_goal[i];
    }

    if (max_dev > total_weight) {
      max_dev = vwgt_max;
    }
    goal_weight = total_weight * KL_IMBALANCE / nsets;
    if (goal_weight > max_dev) {
      max_dev = (int)goal_weight;
    }

    if (COARSE_KLV) {
      klvspiff(graph, nvtxs, assignment, real_goal, max_dev, &bndy_list, weights);
    }
    if (COARSE_BPM) {
      bpm_improve(graph, assignment, real_goal, max_dev, &bndy_list, weights, using_vwgts);
    }
    *pbndy_list = bndy_list;
    return;
  }

  /* Otherwise I have to coarsen. */
  flattened = false;
  if (coords != null) {
    ccoords = (float**)Marshal.AllocHGlobal(igeom * sizeof(float *));
  }
  else {
    ccoords = null;
  }
  if (FLATTEN && step == 0) {
    flattened = flatten(graph, nvtxs, nedges, &cgraph, &cnvtxs, &cnedges, &v2cv,
                        useEdgeWeights && COARSEN_EWGTS, igeom, coords, ccoords);
  }
  if (!flattened) {
    coarsen1(graph, nvtxs, nedges, &cgraph, &cnvtxs, &cnedges, &v2cv, igeom, coords, ccoords,
             useEdgeWeights && COARSEN_EWGTS);
  }

  if (term_wgts[1] != null) {
    twptr      = (float*)Marshal.AllocHGlobal((cnvtxs + 1) * (nsets - 1) * sizeof(float));
    twptr_save = twptr;
    for (i = (cnvtxs + 1) * (nsets - 1); i != 0; i--) {
      *twptr++ = 0;
    }
    twptr = twptr_save;
    for (j = 1; j < nsets; j++) {
      cterm_wgts[j] = twptr;
      twptr += cnvtxs + 1;
    }
    for (j = 1; j < nsets; j++) {
      ctwptr = cterm_wgts[j];
      twptr  = term_wgts[j];
      for (i = 1; i < nvtxs; i++) {
        ctwptr[v2cv[i]] += twptr[i];
      }
    }
  }

  else {
    cterm_wgts[1] = null;
  }

  /* If coarsening isn't working very well, give up and partition. */
  give_up = false;
  if (nvtxs * COARSEN_RATIO_MIN < cnvtxs && cnvtxs > vmax && !flattened) {
      Trace.WriteLine($"WARNING: Coarsening not making enough progress, {nameof(nvtxs)} = {nvtxs:d}, {nameof(cnvtxs)} = {cnvtxs:d}.");
    Trace.WriteLine("         Recursive coarsening being stopped prematurely.");
    give_up = true;
  }

  /* Now recurse on coarse subgraph. */
  if (COARSEN_VWGTS) {
    cvwgt_max = 0;
    for (i = 1; i <= cnvtxs; i++) {
      if (cgraph[i]->vwgt > cvwgt_max) {
        cvwgt_max = cgraph[i]->vwgt;
      }
    }
  }

  else {
    cvwgt_max = 1;
  }

  cassignment = (int*)Marshal.AllocHGlobal((cnvtxs + 1) * sizeof(int));
  if (flattened) {
    nextstep = step;
  }
  else {
    nextstep = step + 1;
  }
  coarsen_klv(cgraph, cnvtxs, cnedges, COARSEN_VWGTS, COARSEN_EWGTS, cterm_wgts, igeom, ccoords,
              cvwgt_max, cassignment, goal, architecture, hops, solver_flag, ndims, nsets, vmax,
              mediantype, mkconnected, eigtol, nstep, nextstep, &cbndy_list, weights, give_up);

  /* Interpolate assignment back to fine graph. */
  for (i = 1; i <= nvtxs; i++) {
    assignment[i] = cassignment[v2cv[i]];
  }

  /* Construct boundary list from coarse boundary list. */
  flag = (bool*)Marshal.AllocHGlobal((cnvtxs + 1) * sizeof(bool));
  for (i = 1; i <= cnvtxs; i++) {
    flag[i] = false;
  }
  for (i = 0; cbndy_list[i] != 0; i++) {
    flag[cbndy_list[i]] = true;
  }

  list_length = 0;
  for (i = 1; i <= nvtxs; i++) {
    if (flag[v2cv[i]]) {
      ++list_length;
    }
  }

  bndy_list = (int*)Marshal.AllocHGlobal((list_length + 1) * sizeof(int));

  list_length = 0;
  for (i = 1; i <= nvtxs; i++) {
    if (flag[v2cv[i]]) {
      bndy_list[list_length++] = i;
    }
  }
  bndy_list[list_length] = 0;

  Marshal.FreeHGlobal((IntPtr)flag);
  Marshal.FreeHGlobal((IntPtr)cbndy_list);

  /* Free the space that was allocated. */
  Marshal.FreeHGlobal((IntPtr)cassignment);
  if (twptr_save != null) {
      Marshal.FreeHGlobal((IntPtr)twptr_save);
    twptr_save = null;
  }
  free_graph(cgraph);
  Marshal.FreeHGlobal((IntPtr)v2cv);

  /* Smooth using KL or BPM every nstep steps. */
  if ((step % nstep) == 0 && !flattened) {
    if (!COARSEN_VWGTS && step != 0) { /* Construct new goal */
      goal_weight = 0;
      for (i = 0; i < nsets; i++) {
        goal_weight += goal[i];
      }
      for (i = 0; i < nsets; i++) {
        new_goal[i] = goal[i] * (nvtxs / goal_weight);
      }
      real_goal = new_goal;
    }
    else {
      real_goal = goal;
    }

    if (LIMIT_KL_EWGTS) {
      find_maxdeg(graph, nvtxs, useEdgeWeights, &ewgt_max);
      compress_ewgts(graph, nvtxs, nedges, ewgt_max, useEdgeWeights);
    }

    /* If not coarsening ewgts, then need care with term_wgts. */
    if (!useEdgeWeights && term_wgts[1] != null && step != 0) {
      twptr      = (float*)Marshal.AllocHGlobal((nvtxs + 1) * (nsets - 1) * sizeof(float));
      twptr_save = twptr;
      for (j = 1; j < nsets; j++) {
        new_term_wgts[j] = twptr;
        twptr += nvtxs + 1;
      }

      for (j = 1; j < nsets; j++) {
        twptr  = term_wgts[j];
        ctwptr = new_term_wgts[j];
        for (i = 1; i <= nvtxs; i++) {
          if (twptr[i] > .5) {
            ctwptr[i] = 1;
          }
          else if (twptr[i] < -.5) {
            ctwptr[i] = -1;
          }
          else {
            ctwptr[i] = 0;
          }
        }
      }
      real_term_wgts = new_term_wgts;
    }
    else {
      real_term_wgts   = term_wgts;
      new_term_wgts[1] = null;
    }

    max_dev      = (step == 0) ? vwgt_max : 5 * vwgt_max;
    total_weight = 0;
    for (i = 0; i < nsets; i++) {
      total_weight += real_goal[i];
    }
    if (max_dev > total_weight) {
      max_dev = vwgt_max;
    }
    goal_weight = total_weight * KL_IMBALANCE / nsets;
    if (goal_weight > max_dev) {
      max_dev = (int)goal_weight;
    }

    if (!COARSEN_VWGTS) {
      count_weights(graph, nvtxs, assignment, nsets + 1, weights, (vwgt_max != 1));
    }

    if (COARSE_KLV) {
      klvspiff(graph, nvtxs, assignment, real_goal, max_dev, &bndy_list, weights);
    }
    if (COARSE_BPM) {
      bpm_improve(graph, assignment, real_goal, max_dev, &bndy_list, weights, using_vwgts);
    }

    if (real_term_wgts != term_wgts && new_term_wgts[1] != null) {
        Marshal.FreeHGlobal((IntPtr)real_term_wgts[1]);
    }

    if (LIMIT_KL_EWGTS) {
      restore_ewgts(graph, nvtxs);
    }
  }

  *pbndy_list = bndy_list;

  if (twptr_save != null) {
      Marshal.FreeHGlobal((IntPtr)twptr_save);
  }

  /* Free the space that was allocated. */
  if (ccoords != null) {
    for (i = 0; i < igeom; i++) {
        Marshal.FreeHGlobal((IntPtr)ccoords[i]);
    }
    Marshal.FreeHGlobal((IntPtr)ccoords);
  }

  if (DEBUG_COARSEN)
  {
      Trace.WriteLine($" Leaving coarsen_klv, {nameof(step)}={step:d}");
  }
}

public static void print_sep_size(int *list, vtx_data **graph /* array of vtx data for graph */
)
{
  int sep_size, sep_weight;
  int i;

  sep_size = sep_weight = 0;
  for (i = 0; list[i] != 0; i++) {
    sep_size++;
    sep_weight += graph[list[i]]->vwgt;
  }
  Trace.WriteLine($" {nameof(sep_size)} = {sep_size:d}, {nameof(sep_weight)} = {sep_weight:d}");
}
    }
}
