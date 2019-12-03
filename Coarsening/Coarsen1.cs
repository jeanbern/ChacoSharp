using System;
using System.Runtime.InteropServices;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Coarsening.MakeV2Cv;
using static ChacoSharp.Coarsening.MaxMatch;
using static ChacoSharp.Coarsening.MakeFGraph;
using static ChacoSharp.Utilities.Timer;

namespace ChacoSharp.Coarsening
{
    public static unsafe class Coarsen1
    {
        public static void coarsen1(vtx_data ** graph,      /* array of vtx data for graph */
              int                nvtxs,      /* number of vertices in graph */
              int                nedges,     /* number of edges in graph */
              vtx_data ***pcgraph,    /* coarsened version of graph */
              int *              pcnvtxs,    /* number of vtxs in coarsened graph */
              int *              pcnedges,   /* number of edges in coarsened graph */
              int **             pv2cv,      /* pointer to v2cv */
              int                igeom,      /* dimension for geometric information */
              float **           coords,     /* coordinates for vertices */
              float **           ccoords,    /* coordinates for coarsened vertices */
              bool                useEdgeWeights /* are edge weights being used? */
)
{
  double        time;    /* time routine is entered */
  int *         v2cv;    /* maps from vtxs to cvtxs */
  int *         mflag;   /* flag indicating vtx matched or not */
  int           cnvtxs;  /* number of vtxs in coarse graph */
  int           nmerged; /* number of edges contracted */

  time = seconds();

  /* Allocate and initialize space. */
  v2cv  = (int*)Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));
  mflag = (int*)Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));

  /* Find a maximal matching in the graph. */
  nmerged = maxmatch(graph, nvtxs, nedges, mflag, useEdgeWeights, igeom, coords);
  match_time += seconds() - time;

  /* Now construct coarser graph by contracting along matching edges. */
  /* Pairs of values in mflag array indicate matched vertices. */
  /* A zero value indicates that vertex is unmatched. */

  /*
      makecgraph(graph, nvtxs, pcgraph, pcnvtxs, pcnedges, mflag, *pv2cv, nmerged, useEdgeWeights, igeom, coords, ccoords);
      makecgraph2(graph, nvtxs, nedges, pcgraph, pcnvtxs, pcnedges, mflag, *pv2cv, nmerged, useEdgeWeights, igeom, coords, ccoords);
  */

  makev2cv(mflag, nvtxs, v2cv);

  Marshal.FreeHGlobal((IntPtr)mflag);

  cnvtxs = nvtxs - nmerged;
  makefgraph(graph, nvtxs, nedges, pcgraph, cnvtxs, pcnedges, v2cv, useEdgeWeights, igeom, coords,
             ccoords);

  *pcnvtxs = cnvtxs;
  *pv2cv   = v2cv;
  coarsen_time += seconds() - time;
}
    }
}
