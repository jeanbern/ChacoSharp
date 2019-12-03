using System;
using System.Runtime.InteropServices;
using static ChacoSharp.Assignment.RecMedian;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Utilities.Timer;
using static ChacoSharp.Inertial.EigenVec;

namespace ChacoSharp.Inertial
{
    public static unsafe class InertialHelper
    {

public static void inertial(vtx_data **graph,        /* graph data structure */
              int               nvtxs,        /* number of vtxs in graph */
              bool               cube_or_mesh, /* 0 => hypercube, d => d-dimensional mesh */
              int               nsets,        /* number of sets to cut into */
              int               igeom,        /* 1, 2 or 3 dimensional geometry? */
              float **          coords,       /* x, y and z coordinates of vertices */
              int *             sets,         /* set each vertex gets assigned to */
              double []          goal,         /* desired set sizes */
              bool               using_vwgts   /* are vertex weights being used? */
)
{
  double        time;            /* timing parameter */
  float *[]       inert_coords = new float*[3]; /* coord arrays passed down */
  int           i, j;            /* loop counters */

  time = seconds();

  if (DEBUG_TRACE) {
    Console.WriteLine("<Entering inertial, nvtxs = {0:d}>", nvtxs);
  }

  if (PROJECTION_AXIS == 0) {
    for (i = 0; i < igeom; i++) {
      inert_coords[i] = coords[i];
    }
  }

  else { /* project out an axis to get long regions */
    j = 0;
    for (i = 0; i < igeom; i++) {
      if (PROJECTION_AXIS != i + 1) {
        inert_coords[j] = coords[i];
        j++;
      }
    }
    --igeom;
  }

  if (igeom == 1) {
    inertial1d(graph, nvtxs, cube_or_mesh, nsets, inert_coords[0], sets, goal, using_vwgts);
  }
  else if (igeom == 2) {
    inertial2d(graph, nvtxs, cube_or_mesh, nsets, inert_coords[0], inert_coords[1], sets, goal,
               using_vwgts);
  }
  else if (igeom == 3) {
    inertial3d(graph, nvtxs, cube_or_mesh, nsets, inert_coords[0], inert_coords[1], inert_coords[2],
               sets, goal, using_vwgts);
  }
  inertial_time += seconds() - time;
}
private static void inertial1d(vtx_data **graph,        /* graph data structure */
int               nvtxs,        /* number of vtxs in graph */
bool               cube_or_mesh, /* 0 => hypercube, d => d-dimensional mesh */
int               nsets,        /* number of sets to divide into */
float *           x,            /* x coordinates of vertices */
int *             sets,         /* set each vertex gets assigned to */
double []          goal,         /* desired set sizes */
bool               using_vwgts   /* are vertex weights being used? */
)
{
double *      value;       /* values passed to median routine */
double        time;        /* timing variables */
int *         space;       /* space required by median routine */
int           i;           /* loop counter */

value = (double*)Marshal.AllocHGlobal((nvtxs + 1) * sizeof(double));

    /* Copy values into double precision array. */
    for (i = 1; i <= nvtxs; i++) {
    value[i] = x[i];
}

/* Now find the median value and partition based upon it. */
space = (int*)Marshal.AllocHGlobal(nvtxs * sizeof(int));

time = seconds();
rec_median_1(graph, value, nvtxs, space, cube_or_mesh, nsets, goal, using_vwgts, sets, true);
median_time += seconds() - time;

Marshal.FreeHGlobal((IntPtr)space);
Marshal.FreeHGlobal((IntPtr)value);
    }


private static void inertial2d(vtx_data **graph,        /* graph data structure for weights */
                int               nvtxs,        /* number of vtxs in graph */
                bool               cube_or_mesh, /* 0 => hypercube, d => d-dimensional mesh */
                int               nsets,        /* number of sets to divide into */
                float *x, float *y,             /* x and y coordinates of vertices */
                int *   sets,                   /* set each vertex gets assigned to */
                double[] goal,                   /* desired set sizes */
                bool     using_vwgts             /* are vertex weights being used? */
)
{
  double[][]        tensor = new double[2][];       /* inertial tensor */
  tensor[0] = new double[2];
  tensor[1] = new double[2];
  double[]        evec = new double[2];            /* eigenvector of tensor */
  double *      value;              /* values along selected direction to sort */
  double        xcm, ycm;           /* center of mass in each direction */
  double        xx, yy, xy;         /* elements of inertial tensor */
  double        xdif, ydif;         /* deviation from center of mass */
  double        eval, res;          /* eigenvalue and error in eval calculation */
  double        vwgt_sum;           /* sum of all the vertex weights */
  double        time;               /* timing parameters */
  int *         space;              /* space required by median routine */
  int           i;                  /* loop counter */

  /* Compute center of mass and total mass. */

  time = seconds();
  xcm = ycm = 0.0;
  if (using_vwgts) {
    vwgt_sum = 0.0;
    for (i = 1; i <= nvtxs; i++) {
      vwgt_sum += graph[i]->vwgt;
      xcm += graph[i]->vwgt * x[i];
      ycm += graph[i]->vwgt * y[i];
    }
  }
  else {
    vwgt_sum = nvtxs;
    for (i = 1; i <= nvtxs; i++) {
      xcm += x[i];
      ycm += y[i];
    }
  }

  xcm /= vwgt_sum;
  ycm /= vwgt_sum;

  /* Generate 3 elements of Inertial tensor. */
  xx = yy = xy = 0.0;
  if (using_vwgts) {
    for (i = 1; i <= nvtxs; i++) {
      xdif = x[i] - xcm;
      ydif = y[i] - ycm;
      xx += graph[i]->vwgt * xdif * xdif;
      yy += graph[i]->vwgt * ydif * ydif;
      xy += graph[i]->vwgt * xdif * ydif;
    }
  }
  else {
    for (i = 1; i <= nvtxs; i++) {
      xdif = x[i] - xcm;
      ydif = y[i] - ycm;
      xx += xdif * xdif;
      yy += ydif * ydif;
      xy += xdif * ydif;
    }
  }

  /* Compute eigenvector with maximum eigenvalue. */

  tensor[0][0] = xx;
  tensor[1][1] = yy;
  tensor[1][0] = tensor[0][1] = xy;
  evals2(tensor, &res, &eval);
  eigenvec2(tensor, eval, evec, &res);

  inertial_axis_time += seconds() - time;

  if (DEBUG_INERTIAL) {
    Console.WriteLine("Principle Axis = ({0:g}, {1:g}), Eval={2:g}, Residual={3:e}", evec[0], evec[1], eval, res);
  }

  /* Allocate space for value array. */

  value = (double*)Marshal.AllocHGlobal((nvtxs + 1) * sizeof(double));

  /* Calculate value to sort/split on for each cell. */
  /* This is inner product with eigenvector. */
  for (i = 1; i <= nvtxs; i++) {
    value[i] = (x[i] - xcm) * evec[0] + (y[i] - ycm) * evec[1];
  }

  /* Now find the median value and partition based upon it. */
  space = (int*)Marshal.AllocHGlobal(nvtxs * sizeof(int));
  time  = seconds();
  rec_median_1(graph, value, nvtxs, space, cube_or_mesh, nsets, goal, using_vwgts, sets, true);
  median_time += seconds() - time;

  Marshal.FreeHGlobal((IntPtr)space);
  Marshal.FreeHGlobal((IntPtr)value);
}

private static void inertial3d(vtx_data **graph,        /* graph data structure */
                int               nvtxs,        /* number of vtxs in graph */
                bool               cube_or_mesh, /* 0 => hypercube, d => d-dimensional mesh */
                int               nsets,        /* number of sets to divide into */
                float *x, float *y, float *z,   /* x, y and z coordinates of vertices */
                int *   sets,                   /* set each vertex gets assigned to */
                double[] goal,                   /* desired set sizes */
                bool     using_vwgts             /* are vertex weights being used? */
)
{
  double[][]        tensor = new double[3][];       /* inertia tensor */
  tensor[0] = new double[3];
  tensor[1] = new double[3];
  tensor[2] = new double[3];
  double[]        evec = new double[3];            /* eigenvector */
  double *      value;              /* values along selected direction to sort */
  double        xcm, ycm, zcm;      /* center of mass in each direction */
  double        xx, yy, zz;         /* elements of inertial tensor */
  double        xy, xz, yz;         /* elements of inertial tensor */
  double        xdif, ydif;         /* deviation from center of mass */
  double        zdif;               /* deviation from center of mass */
  double        eval, res;          /* eigenvalue and error in eval calculation */
  double        vwgt_sum;           /* sum of all the vertex weights */
  double        time;               /* timing parameter */
  int *         space;              /* space required by median routine */
  int           i;                  /* loop counter */

  /* Compute center of mass and total mass. */
  time = seconds();
  xcm = ycm = zcm = 0.0;
  if (using_vwgts) {
    vwgt_sum = 0;
    for (i = 1; i <= nvtxs; i++) {
      vwgt_sum += graph[i]->vwgt;
      xcm += graph[i]->vwgt * x[i];
      ycm += graph[i]->vwgt * y[i];
      zcm += graph[i]->vwgt * z[i];
    }
  }
  else {
    vwgt_sum = nvtxs;
    for (i = 1; i <= nvtxs; i++) {
      xcm += x[i];
      ycm += y[i];
      zcm += z[i];
    }
  }

  xcm /= vwgt_sum;
  ycm /= vwgt_sum;
  zcm /= vwgt_sum;

  /* Generate 6 elements of Inertial tensor. */
  xx = yy = zz = xy = xz = yz = 0.0;
  if (using_vwgts) {
    for (i = 1; i <= nvtxs; i++) {
      xdif = x[i] - xcm;
      ydif = y[i] - ycm;
      zdif = z[i] - zcm;
      xx += graph[i]->vwgt * xdif * xdif;
      yy += graph[i]->vwgt * ydif * ydif;
      zz += graph[i]->vwgt * zdif * zdif;
      xy += graph[i]->vwgt * xdif * ydif;
      xz += graph[i]->vwgt * xdif * zdif;
      yz += graph[i]->vwgt * ydif * zdif;
    }
  }
  else {
    for (i = 1; i <= nvtxs; i++) {
      xdif = x[i] - xcm;
      ydif = y[i] - ycm;
      zdif = z[i] - zcm;
      xx += xdif * xdif;
      yy += ydif * ydif;
      zz += zdif * zdif;
      xy += xdif * ydif;
      xz += xdif * zdif;
      yz += ydif * zdif;
    }
  }

  /* Compute eigenvector with maximum eigenvalue. */

  tensor[0][0] = xx;
  tensor[1][1] = yy;
  tensor[2][2] = zz;
  tensor[0][1] = tensor[1][0] = xy;
  tensor[0][2] = tensor[2][0] = xz;
  tensor[1][2] = tensor[2][1] = yz;
  ch_evals3(tensor, &res, &res, &eval);
  ch_eigenvec3(tensor, eval, evec, &res);

  inertial_axis_time += seconds() - time;

  if (DEBUG_INERTIAL) {
    Console.WriteLine("Principle Axis = ({0:g}, {1:g}, {2:g})\n  Eval={3:g}, Residual={4:e}", evec[0], evec[1], evec[2], eval, res);
  }

  /* Allocate space for value array. */

  value = (double*)Marshal.AllocHGlobal((nvtxs + 1) * sizeof(double));

  /* Calculate value to sort/split on for each cell. */
  /* This is inner product with eigenvector. */
  for (i = 1; i <= nvtxs; i++) {
    value[i] = (x[i] - xcm) * evec[0] + (y[i] - ycm) * evec[1] + (z[i] - zcm) * evec[2];
  }

  /* Now find the median value and partition based upon it. */
  space = (int*)Marshal.AllocHGlobal(nvtxs * sizeof(int));
  time  = seconds();
  rec_median_1(graph, value, nvtxs, space, cube_or_mesh, nsets, goal, using_vwgts, sets, true);
  median_time += seconds() - time;

  Marshal.FreeHGlobal((IntPtr)space);
  Marshal.FreeHGlobal((IntPtr)value);
}

    }
}
