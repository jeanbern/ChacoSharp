using System;
using System.Runtime.InteropServices;
using static ChacoSharp.Graph.Reformat;

namespace ChacoSharp.RefineMap
{
    public static unsafe class MakeCommGraph
    {
        /* Construct a weighted quotient graph representing the inter-set communication. */
public static bool make_comm_graph(vtx_data ***pcomm_graph, /* graph for communication requirements */
                    vtx_data ** graph,       /* graph data structure */
                    int                nvtxs,       /* number of vertices in graph */
                    bool                useEdgeWeights, /* are edge weights being used? */
                    int *              assign,      /* current assignment */
                    int                nsets_tot    /* total number of sets */
)
{
  float  ewgt;               /* edge weight in graph */
  int ** edges_list  = null; /* lists of edges */
  int ** ewgts_list  = null; /* lists of edge weights */
  int *  edges       = null; /* edges in communication graph */
  int *  ewgts       = null; /* edge weights in communication graph */
  float *float_ewgts = null; /* edge weights in floating point */
  int *  adj_sets    = null; /* weights connecting sets */
  int *  order       = null; /* ordering of vertices by set */
  int *  sizes       = null; /* sizes of different sets */
  int *  start       = null; /* pointers into adjacency data */
  int *  adjacency   = null; /* array with all the edge info */
  int *  eptr        = null; /* loops through edges in graph */
  int *  ewptr       = null; /* loop through edge weights */
  int    set, set2;          /* sets two vertices belong to */
  int    vertex;             /* vertex in graph */
  int    ncomm_edges;        /* number of edges in communication graph */
  bool    error;              /* out of space? */
  int    i, j;               /* loop counters */

  error        = true;
  *pcomm_graph = null;

  /* First construct some mappings to ease later manipulations. */
  sizes = (int*)Marshal.AllocHGlobal(nsets_tot * sizeof(int));
  if (sizes == null) {
    goto skip;
  }

  for (i = 0; i < nsets_tot; i++) {
    sizes[i] = 0;
  }
  for (i = 1; i <= nvtxs; i++) {
    ++(sizes[assign[i]]);
  }

  /* Now make sizes reflect the start index for each set. */
  for (i = 1; i < nsets_tot - 1; i++) {
    sizes[i] += sizes[i - 1];
  }
  for (i = nsets_tot - 1; i != 0; i--) {
    sizes[i] = sizes[i - 1];
  }
  sizes[0] = 0;

  /* Now construct list of all vertices in set 0, all in set 1, etc. */
  order = (int*)Marshal.AllocHGlobal(nvtxs * sizeof(int));
  if (order == null) {
    goto skip;
  }
  for (i = 1; i <= nvtxs; i++) {
    set               = assign[i];
    order[sizes[set]] = i;
    ++sizes[set];
  }

  /* For each set, find total weight to all neighbors. */
  adj_sets   = (int*)Marshal.AllocHGlobal(nsets_tot * sizeof(int));
  edges_list = (int**)Marshal.AllocHGlobal(nsets_tot * sizeof(int *));
  ewgts_list = (int**)Marshal.AllocHGlobal(nsets_tot * sizeof(int *));
  start      = (int*)Marshal.AllocHGlobal((nsets_tot + 1) * sizeof(int));
  if (adj_sets == null || edges_list == null || ewgts_list == null || start == null) {
    goto skip;
  }

  start[0]    = 0;
  ewgt        = 1;
  ncomm_edges = 0;

  for (set = 0; set < nsets_tot; set++) {
    edges_list[set] = null;
    ewgts_list[set] = null;
  }

  for (set = 0; set < nsets_tot; set++) {
    for (i = 0; i < nsets_tot; i++) {
      adj_sets[i] = 0;
    }
    for (i = (set != 0 ? sizes[set - 1] : 0); i < sizes[set]; i++) {
      vertex = order[i];
      for (j = 1; j < graph[vertex]->nedges; j++) {
        set2 = assign[graph[vertex]->edges[j]];
        if (set2 != set) {
          if (useEdgeWeights) {
            ewgt = graph[vertex]->ewgts[j];
          }
          adj_sets[set2] += (int)ewgt;
        }
      }
    }

    /* Now save adj_sets data to later construct graph. */
    j = 0;
    for (i = 0; i < nsets_tot; i++) {
      if (adj_sets[i] != 0) {
        j++;
      }
    }
    ncomm_edges += j;
    start[set + 1] = ncomm_edges;
    if (j != 0) {
      edges_list[set] = edges = (int*)Marshal.AllocHGlobal(j * sizeof(int));
      ewgts_list[set] = ewgts = (int*)Marshal.AllocHGlobal(j * sizeof(int));
      if (edges == null || ewgts == null) {
        goto skip;
      }
    }
    j = 0;
    for (i = 0; i < nsets_tot; i++) {
      if (adj_sets[i] != 0) {
        edges[j] = i + 1;
        ewgts[j] = adj_sets[i];
        j++;
      }
    }
  }

  Marshal.FreeHGlobal((IntPtr)adj_sets);
  Marshal.FreeHGlobal((IntPtr)order);
  Marshal.FreeHGlobal((IntPtr)sizes);
  adj_sets = order = sizes = null;

  /* I now need to pack the edge and weight data into single arrays. */
  adjacency   = (int*)Marshal.AllocHGlobal((ncomm_edges + 1) * sizeof(int));
  float_ewgts = (float*)Marshal.AllocHGlobal((ncomm_edges + 1) * sizeof(float));
  if (adjacency == null || float_ewgts == null) {
    goto skip;
  }

  for (set = 0; set < nsets_tot; set++) {
    j     = start[set];
    eptr  = edges_list[set];
    ewptr = ewgts_list[set];
    for (i = start[set]; i < start[set + 1]; i++) {
      adjacency[i]   = eptr[i - j];
      float_ewgts[i] = ewptr[i - j];
    }
    if (start[set] != start[set + 1]) {
        Marshal.FreeHGlobal((IntPtr)edges_list[set]);
        Marshal.FreeHGlobal((IntPtr)ewgts_list[set]);
    }
  }
  Marshal.FreeHGlobal((IntPtr)edges_list);
  Marshal.FreeHGlobal((IntPtr)ewgts_list);
  edges_list = ewgts_list = null;

  error =
      reformat(start, adjacency, nsets_tot, &ncomm_edges, (int *)null, float_ewgts, pcomm_graph);

skip:
Marshal.FreeHGlobal((IntPtr)adj_sets);
Marshal.FreeHGlobal((IntPtr)order);
Marshal.FreeHGlobal((IntPtr)sizes);
  if (edges_list != null) {
    for (set = nsets_tot - 1; set >= 0; set--) {
      if (edges_list[set] != null) {
          Marshal.FreeHGlobal((IntPtr)edges_list[set]);
      }
    }
    Marshal.FreeHGlobal((IntPtr)edges_list);
  }

  if (ewgts_list != null) {
    for (set = nsets_tot - 1; set >= 0; set--) {
      if (ewgts_list[set] != null) {
          Marshal.FreeHGlobal((IntPtr)ewgts_list[set]);
      }
    }
    Marshal.FreeHGlobal((IntPtr)ewgts_list);
  }

  Marshal.FreeHGlobal((IntPtr)float_ewgts);
  Marshal.FreeHGlobal((IntPtr)adjacency);
  Marshal.FreeHGlobal((IntPtr)start);

  return (error);
}

    }
}
