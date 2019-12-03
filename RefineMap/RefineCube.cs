using System;
using System.Runtime.InteropServices;
using static ChacoSharp.RefineMap.RefineMapHelper;
using static ChacoSharp.Graph.CheckGraph;
using static ChacoSharp.Graph.FindMaxDeg;
using static ChacoSharp.Utilities.MergeSort;

namespace ChacoSharp.RefineMap
{
    public static unsafe class RefineCube
    {
        /* Use a greedy strategy to swap assignments to reduce hops. */
/* Note that because of our graph data structure, set assignments in the graph */
/* begin at 1 instead of at 0. */
public static bool refine_cube(vtx_data **comm_graph, /* graph for communication requirements */
                int               ndims_tot,  /* dimensionality of hypercube */
                double            maxdesire,  /* largest possible desire to flip an edge */
                int *             vtx2node,   /* mapping from comm_graph vtxs to processors */
                int *             node2vtx    /* mapping from processors to comm_graph vtxs */
)
{
  refine_vdata * vdata = null;      /* desire data for vertices */
  refine_vdata * vptr;              /* loops through vdata */
  refine_edata * edata = null;      /* desire data for edges */
  refine_edata * eptr;              /* loops through edata */
  refine_edata * eguy;              /* one element in edata array */
  refine_edata **desire_ptr = null; /* array of desire buckets */
  double *              desires    = null; /* each edge's inclination to flip */
  double *              dptr;              /* loops through desire */
  int *                 indices = null;    /* sorted list of desire values */
  int *                 space   = null;    /* used for sorting disire values */
  double                best_desire;       /* desire of max edge to flip */
  int                   imax;              /* maxdesire rounded up */
  int                   nsets_tot;         /* total number of sets/processors */
  int                   neighbor;          /* neighboring vertex */
  int                   dim;               /* loops over cube dimensions */
  int                   mask;              /* bit set for current dimension */
  int                   side;              /* side of hypercube node is on */
  int                   nwires;            /* number of wires in dimension of hypercube */
  int                   nwires_tot;        /* total number of wires in hypercube */
  int                   wire;              /* loops through all wires */
  int                   node1, node2;      /* processors joined by a wire */
  int                   vtx1, vtx2;        /* corresponding vertices in comm_graph */
  bool                   error;             /* out of space? */
  int                   i, j, k;           /* loop counter */

  nsets_tot = 1 << ndims_tot;
  error     = true;

  imax = (int)maxdesire;
  if (imax != maxdesire) {
    imax++;
  }

  /* This is really just ndims_tot different 1-D problems. */

  /* Allocate space for and inititalize the vertex data. */
  vdata =
      (refine_vdata *)Marshal.AllocHGlobal((ndims_tot * nsets_tot + 1) * sizeof(refine_vdata));
  if (vdata == null) {
    goto skip;
  }

  /* Compute each node's desires to move or stay put in each direction. */
  vptr = vdata;
  for (dim = 0; dim < ndims_tot; dim++) {
    mask = 1 << dim;
    for (i = 1; i <= nsets_tot; i++) {
      compute_cube_vdata(++vptr, comm_graph, i, mask, vtx2node);
    }
  }

  /* Now allocate space for and initialize the wire data. */
  nwires     = nsets_tot / 2;
  nwires_tot = nwires * ndims_tot;

  edata = (refine_edata *)Marshal.AllocHGlobal((nwires_tot + 1) * sizeof(refine_edata));
  if (edata == null) {
    goto skip;
  }

  desires = (double*)Marshal.AllocHGlobal(nwires_tot * sizeof(double));
  if (desires == null) {
    goto skip;
  }

  /* Initialize all the wire swap_desire values. */
  eptr = edata;
  dptr = desires;
  i    = 0;
  for (dim = 0; dim < ndims_tot; dim++) {
    mask = 1 << dim;
    for (wire = 0; 2 * wire < nsets_tot; wire++) {
      /* Insert zero bit at position dim. */
      j = (wire >> dim) << (dim + 1);
      i = (wire << 1) ^ j;
      j ^= (i >> 1);
      init_cube_edata(eptr, j, dim, mask);
      *dptr++ = eptr->swap_desire =
          (float) compute_cube_edata(eptr, vdata, nsets_tot, comm_graph, node2vtx);
      eptr++;
    }
  }

  /* Set value for end pointer larger than all others. */
  edata[nwires_tot].swap_desire = (float)(2 * find_maxdeg(comm_graph, nsets_tot, true, (float *)null));

  /* I now need to sort all the wire preference values */
  indices = (int*)Marshal.AllocHGlobal(nwires_tot * sizeof(int));
  space   = (int*)Marshal.AllocHGlobal(nwires_tot * sizeof(int));
  if (indices == null || space == null) {
    goto skip;
  }

  ch_mergesort(desires, nwires_tot, indices, space);

  Marshal.FreeHGlobal((IntPtr)space);
  Marshal.FreeHGlobal((IntPtr)desires);
  space   = null;
  desires = null;

  best_desire = (edata[indices[nwires_tot - 1]]).swap_desire;

  /* Now construct buckets of linked lists with desire values. */

  if (best_desire > 0) {
    desire_ptr =
        (refine_edata **)Marshal.AllocHGlobal((2 * imax + 1) * sizeof(refine_edata *));
    if (desire_ptr == null) {
      goto skip;
    }
    for (i = 2 * imax; i >= 0; i--) {
      desire_ptr[i] = null;
    }

    for (i = nwires_tot - 1; i >= 0; i--) {
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
    mask  = 1 << dim;
    node1 = eguy->node1;
    node2 = eguy->node2;
    vtx1  = node2vtx[node1];
    vtx2  = node2vtx[node2];

    /* Swap the sets. */
    node2vtx[node1] = vtx2;
    node2vtx[node2] = vtx1;
    vtx2node[vtx1]  = node2;
    vtx2node[vtx2]  = node1;

    /* Update all the vdata fields for vertices effected by this flip. */
    /* First do the vertices adjacent to swapped guys, in swapped dimension. */
    side = node1 & mask;
    for (j = 1; j < comm_graph[vtx1]->nedges; j++) {
      neighbor = comm_graph[vtx1]->edges[j];
      if (neighbor != vtx2) {
        update_cube_vdata(side, mask, vtx2node[neighbor], comm_graph[vtx1]->ewgts[j],
                          &(vdata[dim * nsets_tot + neighbor]));
      }
    }

    side = node2 & mask;
    for (j = 1; j < comm_graph[vtx2]->nedges; j++) {
      neighbor = comm_graph[vtx2]->edges[j];
      if (neighbor != vtx1) {
        update_cube_vdata(side, mask, vtx2node[neighbor], comm_graph[vtx2]->ewgts[j],
                          &(vdata[dim * nsets_tot + neighbor]));
      }
    }

    /* Now recompute all preferences for vertices that were moved. */
    for (j = 0; j < ndims_tot; j++) {
      k = 1 << j;
      compute_cube_vdata(&(vdata[j * nsets_tot + vtx1]), comm_graph, vtx1, k, vtx2node);
      compute_cube_vdata(&(vdata[j * nsets_tot + vtx2]), comm_graph, vtx2, k, vtx2node);
    }

    /* Now I can update the values of all the edges associated with all the
       effected vertices.  Note that these include cube neighbors of node1 and
       node2 in addition to the dim-edges of graph neighbors of vtx1 and vtx2. */

    /* For each neighbor vtx, look at wire in this direction.  If desire hasn't changed,
       return.  Otherwise, pick him up and move him in desire list. Similarly for all
       directional neighbors of node1 and node2. */

    for (j = 1; j < comm_graph[vtx1]->nedges; j++) {
      neighbor = comm_graph[vtx1]->edges[j];
      if (neighbor != vtx2) {
        update_cube_edata(neighbor, dim, edata, vdata, comm_graph, node2vtx, vtx2node, nsets_tot,
                          &best_desire, imax, desire_ptr);
      }
    }

    for (j = 1; j < comm_graph[vtx2]->nedges; j++) {
      neighbor = comm_graph[vtx2]->edges[j];
      if (neighbor != vtx1) {
        update_cube_edata(neighbor, dim, edata, vdata, comm_graph, node2vtx, vtx2node, nsets_tot,
                          &best_desire, imax, desire_ptr);
      }
    }
    for (j = 0; j < ndims_tot; j++) {
      update_cube_edata(vtx1, j, edata, vdata, comm_graph, node2vtx, vtx2node, nsets_tot,
                        &best_desire, imax, desire_ptr);
      update_cube_edata(vtx2, j, edata, vdata, comm_graph, node2vtx, vtx2node, nsets_tot,
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
  Marshal.FreeHGlobal((IntPtr)space);
  Marshal.FreeHGlobal((IntPtr)desires);
  Marshal.FreeHGlobal((IntPtr)indices);
  Marshal.FreeHGlobal((IntPtr)desire_ptr);
  Marshal.FreeHGlobal((IntPtr)vdata);
  Marshal.FreeHGlobal((IntPtr)edata);

  return (error);
}

private static void init_cube_edata(refine_edata *edata, /* desire data for current edge */
int                  node1, /* processor incident to current wire */
int                  dim,   /* direction of wire */
int                  mask   /* bit set in wire dimension */
)
{

    edata->node1 = (short)node1;
    edata->node2 = (short)(node1 ^ mask);
    edata->dim   = (short)dim;
}

private static void compute_cube_vdata(refine_vdata *vdata,      /* preference data for a vertex */
vtx_data **   comm_graph, /* communication graph data structure */
int                  vtx,        /* current vertex */
int                  mask,    /* bit set in current hypercube dimension */
int *                vtx2node /* maps graph vtxs to mesh nodes */
)
{
    float same;        /* my preference to stay where I am */
    float change;      /* my preference to change this bit */
    float ewgt;        /* weight of an edge */
    int   neighbor;    /* neighboring vtx in comm_graph */
    int   my_side;     /* which side of hypercube I'm on */
    int   neighb_side; /* which side of hypercube neighbor's on */
    int   j;           /* loop counter */

    my_side = (vtx2node[vtx] & mask);

    change = 0;
    same   = 0;
    for (j = 1; j < comm_graph[vtx]->nedges; j++) {
        neighbor = comm_graph[vtx]->edges[j];
        ewgt     = comm_graph[vtx]->ewgts[j];

        neighb_side = (vtx2node[neighbor] & mask);

        if (neighb_side != my_side) {
            change += ewgt;
        }
        else {
            same += ewgt;
        }
    }
    vdata->same  = same;
    vdata->above = change;
}

private static double compute_cube_edata(refine_edata *edata,      /* desire data for current edge */
refine_vdata *vdata,      /* data for all vertices */
int                  nsets_tot,  /* total number of processors */
vtx_data **   comm_graph, /* communication graph */
int *                node2vtx    /* maps mesh nodes to graph vertices */
)
{
    double desire;     /* edge's interest in flipping */
    float  ewgt;       /* edge weight */
    int    offset;     /* offset into vdata array */
    int    vtx1, vtx2; /* vertices on either side of wire */

    vtx1   = node2vtx[edata->node1];
    vtx2   = node2vtx[edata->node2];
    offset = nsets_tot * edata->dim;

    desire = (vdata[offset + vtx1].above - vdata[offset + vtx1].same) +
             (vdata[offset + vtx2].above - vdata[offset + vtx2].same);

    /* Subtract off potential doubly counted edge. */
    if (is_an_edge(comm_graph[vtx1], vtx2, &ewgt)) {
        desire -= 2 * ewgt;
    }

    return (desire);
}

private static refine_edata *
find_edge_cube(int                  node,     /* processor node */
    int                  dim,      /* direction of edge from node */
refine_edata *edata,    /* data structure for edge preferences */
int                  nsets_tot /* total number of processors */
)
{
    refine_edata *eguy;  /* returned pointer to edge info */
    int                  index; /* computed index into edata */

    /* Squeeze out bit dim from node number. */
    index = node ^ ((node >> dim) << dim);
    index ^= ((node >> (dim + 1)) << dim);
    index += dim * nsets_tot / 2;

    eguy = &(edata[index]);

    return (eguy);
}

private static void update_cube_edata(int                   vertex,     /* graph vertex being worked on */
                       int                   dim,        /* mesh dimension to be adjusted */
                       refine_edata * edata,      /* data structure for edge preferences */
                       refine_vdata * vdata,      /* data structure for vertex preferences */
                       vtx_data **    comm_graph, /* communication graph */
                       int *                 node2vtx,   /* maps processors to comm_graph vtxs */
                       int *                 vtx2node,   /* maps comm_graph vtxs to processors */
                       int                   nsets_tot,  /* total number of processors */
                       double *              best_desire, /* best desire seen */
                       int                   imax,        /* offset in desire_ptr array */
                       refine_edata **desire_ptr   /* buckets for desire values */
)
{
  refine_edata *eguy;       /* data for desired edge */
  float                old_desire; /* original desire for edge to flip */
  float                new_desire; /* new desire for edge to flip */
  int                  node;       /* node number of vertex */
  int                  k;          /* index into desire_ptr array */

  node = vtx2node[vertex];
  eguy = find_edge_cube(node, dim, edata, nsets_tot);

  old_desire = eguy->swap_desire;
  new_desire = (float)compute_cube_edata(eguy, vdata, nsets_tot, comm_graph, node2vtx);

  if (new_desire != old_desire) { /* Update linked list if necessary. */
    eguy->swap_desire = new_desire;
    if (new_desire > *best_desire) {
      *best_desire = new_desire;
    }

    /* Remove eguy from it's current place in list. */
    if (eguy->prev == null) {
      /* Round up for index into desire_ptr. */
      if (old_desire >= 0) {
        k = (int)old_desire;
        if (k != old_desire) {
          k++;
        }
      }
      else {
        k = (int)-old_desire;
        if (k != -old_desire) {
          k++;
        }
        k = -k;
      }
      k += imax;
      desire_ptr[k] = eguy->next;
    }
    else {
      eguy->prev->next = eguy->next;
    }
    if (eguy->next != null) {
      eguy->next->prev = eguy->prev;
    }

    /* Now add eguy to it's new desire bucket. */
    if (new_desire >= 0) {
      k = (int)new_desire;
      if (k != new_desire) {
        k++;
      }
    }
    else {
      k = (int)-new_desire;
      if (k != -new_desire) {
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

private static void update_cube_vdata(int    old_side,      /* previous side for moved vertex in moved dimension */
    int    mask,          /* bit set in current dimension */
    int    neighbor_node, /* node neighbor vertex assigned to */
    double ewgt,          /* weight of edge */
refine_vdata *vdata /* neighbor connected by that edge */
)
{
    int neighbor_side; /* side of cube neighbor is on */

    neighbor_side = (neighbor_node & mask);

    if (neighbor_side == old_side) {
        vdata->above += (float)ewgt;
        vdata->same -= (float)ewgt;
    }
    else {
        vdata->above -= (float)ewgt;
        vdata->same += (float)ewgt;
    }
}

    }
}
