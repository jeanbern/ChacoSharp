using System.Runtime.InteropServices;

namespace ChacoSharp.Connect
{
    public static unsafe class AddEdges
    {
        public static void add_edges(vtx_data **graph,      /* graph data structure */
               edgeslist *new_edges,  /* list of edges connecting graph */
               ilists **  old_edges,  /* edges data overwritten for connecting */
               flists **  old_ewgts,  /* weights of edges overwritten */
               bool               useEdgeWeights /* are edge weights being used? */
)
{
  ilists *   save_list;  /* space to save old edge list */
  flists *   save_ewgts; /* space to save old edge weights */
  edgeslist *edges;      /* loops through new edges */
  float *           new_ewgts;  /* new edge weights */
  int *             new_list;   /* new edge list */
  int               nedges;     /* number of edges a vertex has */
  int               vtx, vtx2;  /* two vertices in edge to be added */
  int               i, j;       /* loop counter */

  *old_edges = null;
  *old_ewgts = null;
  edges      = new_edges;
  while (edges != null) {
    for (j = 0; j < 2; j++) {
      if (j == 0) {
        vtx  = edges->vtx1;
        vtx2 = edges->vtx2;
      }
      else {
        vtx  = edges->vtx2;
        vtx2 = edges->vtx1;
      }

      /* Copy old edge list to new edge list. */
      nedges   = graph[vtx]->nedges;
      new_list = (int*)Marshal.AllocHGlobal((nedges + 1) * sizeof(int));
      for (i = 0; i < nedges; i++) {
        new_list[i] = graph[vtx]->edges[i];
      }
      new_list[nedges] = vtx2;

      /* Save old edges. */
      save_list       = (ilists*)Marshal.AllocHGlobal(sizeof(ilists));
      save_list->list = graph[vtx]->edges;

      /* Add new list at FRONT of linked list to facilitate uncoarsening. */
      save_list->next = *old_edges;
      *old_edges      = save_list;

      /* Now modify graph to have new edges list. */
      graph[vtx]->nedges++;
      graph[vtx]->edges = new_list;

      /* If using edge weights, I have to modify those too. */
      if (useEdgeWeights) {
        new_ewgts = (float*)Marshal.AllocHGlobal((nedges + 1) * sizeof(float));
        for (i = 1; i < nedges; i++) {
          new_ewgts[i] = graph[vtx]->ewgts[i];
        }
        new_ewgts[nedges] = 1;
        new_ewgts[0]      = graph[vtx]->ewgts[0] - new_ewgts[nedges];

        /* Save old edge weights. */
        save_ewgts       = (flists*)Marshal.AllocHGlobal(sizeof(flists));
        save_ewgts->list = graph[vtx]->ewgts;

        save_ewgts->next = *old_ewgts;
        *old_ewgts       = save_ewgts;

        /* Finally, modify graph to have new edge weights. */
        graph[vtx]->ewgts = new_ewgts;
      }
    }
    edges = edges->next;
  }
}
    }
}
