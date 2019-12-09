using System;
using System.Diagnostics;
using System.Runtime.InteropServices;
using static ChacoSharp.Graph.FindMaxDeg;
using static ChacoSharp.Graph.FreeGraph;
using static ChacoSharp.RefineMap.RefineMesh;
using static ChacoSharp.RefineMap.RefineCube;
using static ChacoSharp.RefineMap.MakeCommGraph;

namespace ChacoSharp.RefineMap
{
    public static unsafe class RefineMapHelper
    {
        /* Given a partition, refine the mapping in a locally greedy fashion. */

public static void refine_map(vtx_data **graph,        /* graph data structure */
                int               nvtxs,        /* number of vertices in graph */
                bool               useEdgeWeights,  /* are edge weights being used? */
                int *             assign,       /* current assignment */
                int               cube_or_mesh, /* 0 => hypercube, d => d-dimensional mesh */
                int               ndims_tot,    /* if hypercube, number of dimensions */
                int[]               mesh_dims/*[3]*/  /* if mesh, dimensions of mesh */
)
{
  vtx_data **comm_graph;       /* graph for communication requirements */
  int               nsets_tot = 0;    /* total number of sets */
  int *             vtx2node  = null; /* mapping of comm_graph vtxs to processors */
  int *             node2vtx  = null; /* mapping of sets to comm_graph vtxs */
  double            maxdesire = 0.0;  /* largest possible desire to flip an edge */
  bool               error     = false;    /* out of space? */
  int               i;                /* loop counter */

  if (cube_or_mesh == 0) {
    nsets_tot = 1 << ndims_tot;
  }
  else if (cube_or_mesh == 1) {
    nsets_tot = mesh_dims[0];
  }
  else if (cube_or_mesh == 2) {
    nsets_tot = mesh_dims[0] * mesh_dims[1];
  }
  else if (cube_or_mesh == 3) {
    nsets_tot = mesh_dims[0] * mesh_dims[1] * mesh_dims[2];
  }

  node2vtx = vtx2node = null;

  /* Construct the weighted quotient graph representing communication. */
  error = make_comm_graph(&comm_graph, graph, nvtxs, useEdgeWeights, assign, nsets_tot);

  if (!error) {
    maxdesire = 2 * find_maxdeg(comm_graph, nsets_tot, true, (float *)null);

    vtx2node = (int*)Marshal.AllocHGlobal((nsets_tot + 1) * sizeof(int));
    node2vtx = (int*)Marshal.AllocHGlobal(nsets_tot * sizeof(int));
    if (node2vtx == null || vtx2node == null) {
      error = true;
      goto skip;
    }

    for (i = 1; i <= nsets_tot; i++) {
      vtx2node[i]     = i - 1;
      node2vtx[i - 1] = i;
    }

    if (cube_or_mesh > 0) {
      error = refine_mesh(comm_graph, cube_or_mesh, mesh_dims, maxdesire, vtx2node, node2vtx);
    }

    else if (cube_or_mesh == 0) {
      error = refine_cube(comm_graph, ndims_tot, maxdesire, vtx2node, node2vtx);
    }

    if (!error) {
      for (i = 1; i <= nvtxs; i++) {
        assign[i] = vtx2node[assign[i] + 1];
      }
    }
  }

skip:

  if (error) {
    Trace.WriteLine("\nWARNING: No space to refine mapping to processors.");
    Trace.WriteLine("         NO MAPPING REFINEMENT PERFORMED.\n");
  }

  Marshal.FreeHGlobal((IntPtr)node2vtx);
  Marshal.FreeHGlobal((IntPtr)vtx2node);
  free_graph(comm_graph);
}

        public struct refine_vdata
        {
            public float above; /* sum of edge weights pulling me higher */
            public float below; /* sum of edge weights pulling me lower */
            public float same;  /* sum of edge weights keeping me here */
        };

        public struct refine_edata
        {
            public short                node1, node2; /* nodes in mesh connected by this edge */
            public short                dim;          /* which dimension of mesh does wire span? */
            public float                swap_desire;  /* reduction in hops if edge is flipped */
            public refine_edata *prev;         /* pointer to previous guy in list */
            public refine_edata *next;         /* pointer to next guy in list */
        };
    }
}
