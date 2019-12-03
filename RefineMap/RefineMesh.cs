using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using static ChacoSharp.Utilities.MergeSort;
using static ChacoSharp.Graph.FindMaxDeg;
using static ChacoSharp.RefineMap.RefineMeshData;
using static ChacoSharp.RefineMap.RefineMapHelper;

namespace ChacoSharp.RefineMap
{
    public static unsafe class RefineMesh
    {
        /* Use a greedy strategy to swap assignments to reduce hops. */
/* Note that because of our graph data structure, set assignments in the graph */
/* begin at 1 instead of at 0. */
public static bool refine_mesh(vtx_data **comm_graph,   /* graph for communication requirements */
                int               cube_or_mesh, /* number of dimensions in mesh */
                int[]               mesh_dims/*[3]*/, /* dimensions of mesh */
                double            maxdesire,    /* largest possible desire to flip an edge */
                int *             vtx2node,     /* mapping from comm_graph vtxs to mesh nodes */
                int *             node2vtx      /* mapping from mesh nodes to comm_graph vtxs */
)
{
  refine_vdata * vdata = null;      /* desire data for all vertices */
  refine_vdata * vptr;              /* loops through vdata */
  refine_edata * edata = null;      /* desire data for all edges */
  refine_edata * eguy;              /* one element in edata array */
  refine_edata **desire_ptr = null; /* array of desire buckets */
  double *              desires    = null; /* each edge's inclination to flip */
  int *                 indices    = null; /* sorted list of desire values */
  int *                 space      = null; /* used for sorting disire values */
  double                best_desire;       /* highest desire of edge to flip */
  int                   imax;              /* maxdesire rounded up */
  int                   nsets_tot;         /* total number of sets/processors */
  int                   neighbor;          /* neighboring vertex */
  int                   dim;               /* loops over mesh dimensions */
  int                   nwires;            /* number of wires in processor mesh */
  int                   wire;              /* loops through all wires */
  int                   node1, node2;      /* processors joined by a wire */
  int                   vtx1, vtx2;        /* corresponding vertices in comm_graph */
  int                   loc1, loc2;        /* location of vtxs in flipping dimension */
  bool                   error;             /* out of space? */
  int                   i, j, k;           /* loop counter */

  error = true;

  nsets_tot = mesh_dims[0] * mesh_dims[1] * mesh_dims[2];

  imax = (int)maxdesire;
  if (imax != maxdesire) {
    imax++;
  }

  vdata = (refine_vdata*) Marshal.AllocHGlobal((cube_or_mesh * nsets_tot + 1) * sizeof(refine_vdata));
  if (vdata == null) {
    goto skip;
  }

  /* Compute each node's desires to move or stay put. */
  vptr = vdata;
  for (dim = 0; dim < cube_or_mesh; dim++) {
    for (i = 1; i <= nsets_tot; i++) {
      compute_mesh_vdata(++vptr, comm_graph, i, vtx2node, mesh_dims, dim);
    }
  }

  nwires = (mesh_dims[0] - 1) * mesh_dims[1] * mesh_dims[2] +
           mesh_dims[0] * (mesh_dims[1] - 1) * mesh_dims[2] +
           mesh_dims[0] * mesh_dims[1] * (mesh_dims[2] - 1);

  edata   = (refine_edata*)Marshal.AllocHGlobal((nwires + 1) * sizeof(refine_edata));
  desires = (double*)Marshal.AllocHGlobal(nwires * sizeof(double));
  if (desires == null) {
    goto skip;
  }

  /* Initialize all the edge values. */
  init_mesh_edata(edata, mesh_dims);
  for (wire = 0; wire < nwires; wire++) {
    desires[wire] = edata[wire].swap_desire =
        (float)compute_mesh_edata(&(edata[wire]), vdata, mesh_dims, comm_graph, node2vtx);
  }

  /* Set special value for end pointer. */
  edata[nwires].swap_desire = 2.0f * (float)find_maxdeg(comm_graph, nsets_tot, true, (float *)null);

  /* I now need to sort all the wire preference values */
  indices = (int*)Marshal.AllocHGlobal(nwires * sizeof(int));
  space   = (int*)Marshal.AllocHGlobal(nwires * sizeof(int));
  if (indices == null || space == null) {
    goto skip;
  }

  ch_mergesort(desires, nwires, indices, space);

  Marshal.FreeHGlobal((IntPtr)space);
  Marshal.FreeHGlobal((IntPtr)desires);
  space   = null;
  desires = null;

  best_desire = (edata[indices[nwires - 1]]).swap_desire;

  /* Now construct a buckets of linked lists with desire values. */

  if (best_desire > 0) {
    desire_ptr =
        (refine_edata **)Marshal.AllocHGlobal((2 * imax + 1) * sizeof(refine_edata *));
    if (desire_ptr == null) {
      goto skip;
    }

    for (i = 2 * imax; i >= 0; i--) {
      desire_ptr[i] = null;
    }

    for (i = nwires - 1; i >= 0; i--) {
      eguy = &(edata[indices[i]]);
      /* Round the swap desire up. */
      if (eguy->swap_desire >= 0) {
        k = (int)eguy->swap_desire;
        if (k != eguy->swap_desire) {
          k++;
        }
      }
      else {
        k = (int)-eguy->swap_desire;
        if (k != -eguy->swap_desire) {
          k++;
        }
        k = -k;
      }

      k += imax;

      eguy->prev = null;
      eguy->next = desire_ptr[k];
      if (desire_ptr[k] != null) {
        desire_ptr[k]->prev = eguy;
      }
      desire_ptr[k] = eguy;
    }
  }
  else {
    desire_ptr = null;
  }

  Marshal.FreeHGlobal((IntPtr)indices);
  indices = null;

  loc1 = 0;
  loc2 = 0;

  /* Everything is now set up.  Swap sets across wires until no more improvement. */

  while (best_desire > 0) {
    k = (int)(best_desire + 1 + imax);
    if (k > 2 * imax) {
      k = 2 * imax;
    }
    while (k > imax && desire_ptr[k] == null) {
      k--;
    }
    eguy = desire_ptr[k];

    dim   = eguy->dim;
    node1 = eguy->node1;
    node2 = eguy->node2;
    vtx1  = node2vtx[node1];
    vtx2  = node2vtx[node2];
    if (dim == 0) {
      loc1 = node1 % mesh_dims[0];
      loc2 = node2 % mesh_dims[0];
    }
    else if (dim == 1) {
      loc1 = (node1 / mesh_dims[0]) % mesh_dims[1];
      loc2 = (node2 / mesh_dims[0]) % mesh_dims[1];
    }
    else if (dim == 2) {
      loc1 = node1 / (mesh_dims[0] * mesh_dims[1]);
      loc2 = node2 / (mesh_dims[0] * mesh_dims[1]);
    }

    /* Now swap the vertices. */
    node2vtx[node1] = vtx2;
    node2vtx[node2] = vtx1;
    vtx2node[vtx1]  = node2;
    vtx2node[vtx2]  = node1;

    /* First update all the vdata fields for vertices effected by this flip. */
    for (j = 1; j < comm_graph[vtx1]->nedges; j++) {
      neighbor = comm_graph[vtx1]->edges[j];
      if (neighbor != vtx2) {
        update_mesh_vdata(loc1, loc2, dim, comm_graph[vtx1]->ewgts[j], vdata, mesh_dims, neighbor,
                          vtx2node);
      }
    }

    for (j = 1; j < comm_graph[vtx2]->nedges; j++) {
      neighbor = comm_graph[vtx2]->edges[j];
      if (neighbor != vtx1) {
        update_mesh_vdata(loc2, loc1, dim, comm_graph[vtx2]->ewgts[j], vdata, mesh_dims, neighbor,
                          vtx2node);
      }
    }

    /* Now recompute all preferences for vertices that were moved. */
    for (j = 0; j < cube_or_mesh; j++) {
      compute_mesh_vdata(&(vdata[j * nsets_tot + vtx1]), comm_graph, vtx1, vtx2node, mesh_dims, j);
      compute_mesh_vdata(&(vdata[j * nsets_tot + vtx2]), comm_graph, vtx2, vtx2node, mesh_dims, j);
    }

    /* Now I can update the values of all the edges associated with all the
       effected vertices.  Note that these include mesh neighbors of node1 and
       node2 in addition to the dim-edges of graph neighbors of vtx1 and vtx2. */

    /* For each neighbor vtx, look at -1 and +1 edge.  If desire hasn't changed,
       return.  Otherwise, pick him up and move him. Similarly for all
       directional neighbors of node1 and node2. */

    for (j = 1; j < comm_graph[vtx1]->nedges; j++) {
      neighbor = comm_graph[vtx1]->edges[j];
      if (neighbor != vtx2) {
        update_mesh_edata(neighbor, dim, edata, vdata, comm_graph, mesh_dims, node2vtx, vtx2node,
                          &best_desire, imax, desire_ptr);
      }
    }

    for (j = 1; j < comm_graph[vtx2]->nedges; j++) {
      neighbor = comm_graph[vtx2]->edges[j];
      if (neighbor != vtx1) {
        update_mesh_edata(neighbor, dim, edata, vdata, comm_graph, mesh_dims, node2vtx, vtx2node,
                          &best_desire, imax, desire_ptr);
      }
    }
    for (j = 0; j < cube_or_mesh; j++) {
      update_mesh_edata(vtx1, j, edata, vdata, comm_graph, mesh_dims, node2vtx, vtx2node,
                        &best_desire, imax, desire_ptr);
      update_mesh_edata(vtx2, j, edata, vdata, comm_graph, mesh_dims, node2vtx, vtx2node,
                        &best_desire, imax, desire_ptr);
    }

    k = (int)(best_desire + 1 + imax);
    if (k > 2 * imax) {
      k = 2 * imax;
    }
    while (k > imax && desire_ptr[k] == null) {
      k--;
    }
    best_desire = k - imax;
  }
  error = false;

skip:
  Marshal.FreeHGlobal((IntPtr)indices);
  Marshal.FreeHGlobal((IntPtr)space);
  Marshal.FreeHGlobal((IntPtr)desires);
  Marshal.FreeHGlobal((IntPtr)desire_ptr);
  Marshal.FreeHGlobal((IntPtr)vdata);
  Marshal.FreeHGlobal((IntPtr)edata);

  return (error);
}
    }
}
