using System;
using System.Runtime.InteropServices;
using static ChacoSharp.Connect.FindComps;
using static ChacoSharp.Connect.AddEdges;

namespace ChacoSharp.Connect
{
    public static unsafe class Connected
    {
        public static void make_connected(
    /* Add edges to make graph connected. */
    vtx_data **    graph,      /* graph data structure */
    int                   nvtxs,      /* number of vertices in graph */
    int *                 nedges,     /* number of edges in graph */
    int *                 mark,       /* space for nvtxs+1 ints */
    int *                 vtxlist,    /* space for nvtxs ints */
    connect_data **cdata,      /* space for connectivity data */
    bool                   useEdgeWeights /* are edges of graph weighted? */
)
{
  edgeslist *new_edges; /* list of edges connecting graph */
  edgeslist *prev_edge; /* pointer for manipulating edge list */
  edgeslist *curr_edge; /* pointer for manipulating edge list */
  edgeslist *next_edge; /* pointer for manipulating edge list */
  int               nadded;    /* number of edges being added */

  /* First find edges needed to make graph connected. */
  nadded = find_edges(graph, nvtxs, mark, vtxlist, &new_edges);

  /* Now add these needed edges to graph data structure if needed. */
  if (nadded == 0) {
    *cdata = null;
  }
  else {
    *cdata              = (connect_data*)Marshal.AllocHGlobal(sizeof(connect_data));
    (*cdata)->old_edges = null;
    (*cdata)->old_ewgts = null;
    add_edges(graph, new_edges, &(*cdata)->old_edges, &(*cdata)->old_ewgts, useEdgeWeights);
    *nedges += nadded;

    /* Now, reverse the order of the new_edges list for consistency with */
    /* the removal order. */
    curr_edge       = new_edges->next;
    new_edges->next = null;
    prev_edge       = new_edges;
    while (curr_edge != null) {
      next_edge       = curr_edge->next;
      curr_edge->next = prev_edge;
      prev_edge       = curr_edge;
      curr_edge       = next_edge;
    }
    (*cdata)->new_edges = prev_edge;
  }
}

public static void make_unconnected(
    /* Restore graph to its pristine state and free space for connectivity. */
    vtx_data **    graph,      /* graph data structure */
    int *                 nedges,     /* number of edges in graph */
    connect_data **cdata,      /* space for connectivity data */
    bool                   useEdgeWeights /* are edges of graph weighted? */
)
{
  ilists *   old_edges = null; /* edges overwritten for connecting */
  flists *   old_ewgts = null; /* weights of edges overwritten */
  edgeslist *new_edges;        /* list of edges connecting graph */
  ilists *   tempi;            /* used for freeing space */
  flists *   tempf;            /* used for freeing space */
  edgeslist *tempe;            /* used for freeing edgelist space */
  edgeslist *edges;            /* loops through new edges */
  int               vtx;              /* vertex in an added edge */
  int               j;                /* loop counters */

  if (*cdata == null) {
    return;
  }

  old_edges = (*cdata)->old_edges;
  old_ewgts = (*cdata)->old_ewgts;
  new_edges = (*cdata)->new_edges;
  Marshal.FreeHGlobal((IntPtr)(*cdata));
  *cdata = null;

  edges = new_edges;
  while (edges != null) {
    /* Restore edges and weights to original status. */
    (*nedges)--;
    for (j = 0; j < 2; j++) {
      if (j == 0) {
        vtx = edges->vtx2;
      }
      else {
        vtx = edges->vtx1;
      }

      Marshal.FreeHGlobal((IntPtr)(graph[vtx]->edges));
      graph[vtx]->edges = old_edges->list;
      graph[vtx]->nedges--;
      tempi     = old_edges;
      old_edges = old_edges->next;
      Marshal.FreeHGlobal((IntPtr)(tempi));

      if (useEdgeWeights) {
          Marshal.FreeHGlobal((IntPtr)(graph[vtx]->ewgts));
        graph[vtx]->ewgts = old_ewgts->list;
        tempf             = old_ewgts;
        old_ewgts         = old_ewgts->next;
        Marshal.FreeHGlobal((IntPtr)(tempf));
      }
    }
    tempe = edges;
    edges = edges->next;
    Marshal.FreeHGlobal((IntPtr)(tempe));
  }
}

/* Print out the added edges. */
public static void print_connected(connect_data *cdata /* space for connectivity data */
)
{
  edgeslist *edges; /* loops through new edges */

  if (cdata == null) {
      Console.WriteLine("No phantom edges\n");
  }
  else {
    Console.WriteLine("Phantom edges: ");
    edges = cdata->new_edges;
    while (edges != null) {
        Console.WriteLine("({0:d},{1:d}) ", edges->vtx1, edges->vtx2);
      edges = edges->next;
    }
    Console.WriteLine("\n");
  }
}

/* Free the edge list created by find_edges. */
public static void free_edgeslist(edgeslist *edge_list /* list to be freed */
)
{
  edgeslist *next_list; /* next guy in list */

  while (edge_list != null) {
    next_list = edge_list->next;
    Marshal.FreeHGlobal((IntPtr)(edge_list));
    edge_list = next_list;
  }
}
    }
}
