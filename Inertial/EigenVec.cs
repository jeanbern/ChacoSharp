using System;
using static ChacoSharp.Optimize.Determinant;

namespace ChacoSharp.Inertial
{
    public static unsafe class EigenVec
    {
        private static int sign(double value)
        {
            return value < 0 ? -1 : 1;
        }

        /* Find eigenvalues of 3x3 symmetric system by solving cubic. */
public static void ch_evals3(double[][]  H/*[3][3]*/, /* 3x3 sym matrix for lowest eigenvalue */
               double *eval1,   /* smallest eigenvalue */
               double *eval2,   /* middle eigenvalue */
               double *eval3    /* largest eigenvalue */
)
{
  double[][] mat = {new []{0.0d, 0.0d, 0.0d}, new []{0.0d, 0.0d, 0.0d}, new[]{0.0d, 0.0d, 0.0d}}; /* scaled version of H */
  double a1, a2, a3;          /* coefficients of cubic equation */
  double q, r;                /* intermediate terms */
  double q3, r2;              /* powers of q and r */
  double theta;               /* angular parameter */
  double root1, root2, root3; /* three roots of cubic */
  double tol = 1.0e-6;        /* allowed deviation */
  double xmax;                /* largest matrix element for scaling */
  int    i, j;                /* loop indices */

  /* This first requires solving a cubic equation. */
  /* Normalize to avoid any numerical problems. */
  xmax = 0.0;
  for (i = 0; i < 3; i++) {
    for (j = i; j < 3; j++) {
      if (Math.Abs(H[i][j]) > xmax) {
        xmax = Math.Abs(H[i][j]);
      }
    }
  }
  if (xmax != 0) {
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        mat[i][j] = H[i][j] / xmax;
      }
    }
  }

  a1 = -(mat[0][0] + mat[1][1] + mat[2][2]);
  a2 = (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]) +
       (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) +
       (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]);
  a3 = -determinant(mat, 3);

  if (a3 == 0) {
    root1 = 0;
    /* Solve quadratic. */
    q     = -.5 * (a1 + sign(a1) * Math.Sqrt(a1 * a1 - 4 * a2));
    root2 = q;
    root3 = a2 / q;
  }

  else { /* solve cubic */
    q  = (a1 * a1 - 3 * a2) / 9;
    r  = (2 * a1 * a1 * a1 - 9 * a1 * a2 + 27 * a3) / 54;
    q3 = q * q * q;
    r2 = r * r;

    /* To avoid missing a root, check for roundoff. */
    if ((q3 < r2) && Math.Abs(q3 - r2) < tol * (Math.Abs(q3) + Math.Abs(r2))) {
      q3 = r2;
    }

    if (q3 >= r2) { /* Three real roots. */
      if (r == 0) {
        theta = Math.PI/2;
      }
      else {
        q3 = Math.Sqrt(q3);
        if (q3 < Math.Abs(r)) {
          q3 = Math.Abs(r);
        }
        theta = Math.Acos(r / q3);
      }
      q = -2 * Math.Sqrt(q);

      root1 = q * Math.Cos(theta / 3) - a1 / 3;
      root2 = q * Math.Cos((theta + Math.PI*2) / 3) - a1 / 3;
      root3 = q * Math.Cos((theta + 2 * Math.PI*2) / 3) - a1 / 3;
    }
    else { /* Only one real root. */
      theta = Math.Sqrt(r2 - q3) + Math.Abs(r);
      theta = Math.Pow(theta, 1.0 / 3.0);

      root1 = root2 = root3 = -sign(r) * (theta + q / theta) - a1 / 3;
    }
  }
  root1 *= xmax;
  root2 *= xmax;
  root3 *= xmax;
  *eval1 = Math.Min(root1, root2);
  *eval1 = Math.Min(*eval1, root3);
  *eval3 = Math.Max(root1, root2);
  *eval3 = Math.Max(*eval3, root3);
  if (root1 != *eval1 && root1 != *eval3) {
    *eval2 = root1;
  }
  else if (root2 != *eval1 && root2 != *eval3) {
    *eval2 = root2;
  }
  else {
    *eval2 = root3;
  }
}

private static void kramer3(                                         /* Use Kramer's rule to solve 3x3 */
             double[][] A/*[3][3]*/, double[] b/*[3]*/, double[] x/*[3]*/ /* Solve Ax=b */
)
{
  double det; /* determinant of system */

  det = 1.0 / determinant(A, 3);

  x[0] = (b[0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
          b[1] * (A[0][1] * A[2][2] - A[0][2] * A[2][1]) +
          b[2] * (A[0][1] * A[1][2] - A[0][2] * A[1][1])) *
         det;

  x[1] = -(b[0] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) -
           b[1] * (A[0][0] * A[2][2] - A[0][2] * A[2][0]) +
           b[2] * (A[0][0] * A[1][2] - A[0][2] * A[1][0])) *
         det;

  x[2] = (b[0] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]) -
          b[1] * (A[0][0] * A[2][1] - A[0][1] * A[2][0]) +
          b[2] * (A[0][0] * A[1][1] - A[0][1] * A[1][0])) *
         det;
}

/* Find the eigenvector of symmetric 3x3 matrix w/ given eigenvalue. */
public static void ch_eigenvec3(double[][]  A/*[3][3]*/, /* matrix to find eigenvector of */
                  double  eval,    /* eigenvalue */
                  double[]  evec/*[3]*/, /* eigenvector returned */
                  double *res      /* returned error estimate */
)
{
  double[][] mat = {new []{0.0d, 0.0d, 0.0d}, new []{0.0d, 0.0d, 0.0d}, new[]{0.0d, 0.0d, 0.0d}};            /* copy of A to write over */
  int[]    ind = new int[3];               /* permutation indices */
  double ex, ey, ez;           /* elements of eigenvector returned */
  double xmax;                 /* maximum value in matrix */
  double tmp;                  /* intermediate values */
  double norm;                 /* norm of eigenvector */
  double res1, res2, res3;     /* elements of residual vector */
  double tol  = 1.0e-6;        /* smaller value assumed to be zero */
  int    imax = -1, jmax = -1; /* indices of Math.Max value in matrix */
  int    i, j;                 /* loop counters */

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      mat[i][j] = A[i][j];
    }
  }
  for (i = 0; i < 3; i++) {
    mat[i][i] -= eval;
  }

  ind[0] = 0;
  ind[1] = 1;
  ind[2] = 2;

  /* Find the largest element in the matrix. */
  xmax = 0.0;
  for (i = 0; i < 3; i++) {
    for (j = i; j < 3; j++) {
      if (Math.Abs(mat[i][j]) > xmax) {
        imax = i;
        jmax = j;
        xmax = Math.Abs(mat[i][j]);
      }
    }
  }

  if (xmax == 0.0) { /* Handle completely degenerate case first. */
    evec[0] = 1.0;
    evec[1] = evec[2] = 0.0;
  }

  else {
    /* Scale the matrix so largest value is 1.0 */
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        mat[i][j] /= xmax;
      }
    }

    /* Swap rows if necessary to move Math.Max element to first row. */
    if (imax != 0) {
      for (j = 0; j < 3; j++) {
        tmp          = mat[0][j];
        mat[0][j]    = mat[imax][j];
        mat[imax][j] = tmp;
      }
    }

    /* Swap columns if necessary to move Math.Max element to first position. */
    if (jmax != 0) {
      for (i = 0; i < 3; i++) {
        tmp          = mat[i][0];
        mat[i][0]    = mat[i][jmax];
        mat[i][jmax] = tmp;
      }
      ind[0]    = jmax;
      ind[jmax] = 0;
    }

    /* Reduce matrix to 2x2 by subtracting off first row. */
    for (i = 1; i < 3; i++) {
      for (j = 1; j < 3; j++) {
        mat[i][j] = mat[0][0] * mat[i][j] - mat[i][0] * mat[0][j];
      }
    }

    /* Find maximum element in reduced 2x2 matrix. */
    xmax = 0.0;
    for (i = 1; i < 3; i++) {
      for (j = i; j < 3; j++) {
        if (Math.Abs(mat[i][j]) > xmax) {
          imax = i;
          jmax = j;
          xmax = Math.Abs(mat[i][j]);
        }
      }
    }

    if (xmax < tol) { /* Handle 2-fold degenerate case - skip to end. */
      ey = 1.0;
      ex = ez = 0;
    }

    else {
      /* Swap rows 2 and 3 to move Math.Max element to 2nd row. */
      if (imax != 1) {
        for (j = 0; j < 3; j++) {
          tmp          = mat[1][j];
          mat[1][j]    = mat[imax][j];
          mat[imax][j] = tmp;
        }
      }

      /* Swap columns to move Math.Max element to (1,1) position. */
      if (jmax != 1) {
        for (i = 0; i < 3; i++) {
          tmp       = mat[i][1];
          mat[i][1] = mat[i][2];
          mat[i][2] = tmp;
        }
        i      = ind[1];
        ind[1] = ind[2];
        ind[2] = i;
      }

      /* Compute eigenvector from numerically stabilized matrix. */
      ez = mat[0][0] * mat[1][1];
      ey = -mat[1][2] * mat[0][0];
      ex = mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1];
    }
    /* Reorder the e-vector to undo pivoting - end of 3-D case. */
    evec[ind[0]] = ex;
    evec[ind[1]] = ey;
    evec[ind[2]] = ez;
  }

  /* Normalize eigenvector and calculate a normalized eigen-residual. */
  norm = Math.Sqrt(evec[0] * evec[0] + evec[1] * evec[1] + evec[2] * evec[2]);
  for (i = 0; i < 3; i++) {
    evec[i] /= norm;
  }
  res1 = (A[0][0] - eval) * evec[0] + A[0][1] * evec[1] + A[0][2] * evec[2];
  res2 = A[1][0] * evec[0] + (A[1][1] - eval) * evec[1] + A[1][2] * evec[2];
  res3 = A[2][0] * evec[0] + A[2][1] * evec[1] + (A[2][2] - eval) * evec[2];
  *res = Math.Sqrt(res1 * res1 + res2 * res2 + res3 * res3);

  /* Now normalize the residual */
  res1 = Math.Abs(A[0][0]) + Math.Abs(A[0][1]) + Math.Abs(A[0][2]);
  res2 = Math.Abs(A[1][0]) + Math.Abs(A[1][1]) + Math.Abs(A[1][2]);
  res3 = Math.Abs(A[2][0]) + Math.Abs(A[2][1]) + Math.Abs(A[2][2]);
  res2 = Math.Max(res2, res3);
  *res /= Math.Max(res1, res2);
}

/* Find eigenvalues of 2x2 symmetric system by solving quadratic. */
public static void evals2(double[][]  H/*[2][2]*/, /* symmetric matrix for eigenvalues */
            double *eval1,   /* smallest eigenvalue */
            double *eval2    /* middle eigenvalue */
)
{
  double[][] M = new []{new []{0.0d, 0.0d}, new []{0.0d, 0.0d} };      /* normalized version of matrix */
  double b, c;         /* coefficients of cubic equation */
  double root1, root2; /* roots of quadratic */
  double xmax;         /* largest matrix element */
  int    i, j;         /* loop counters */

  M[0][0] = M[1][0] = M[0][1] = M[1][1] = 0.0;

  xmax = 0.0;
  for (i = 0; i < 2; i++) {
    for (j = i; j < 2; j++) {
      if (Math.Abs(H[i][j]) > xmax) {
        xmax = Math.Abs(H[i][j]);
      }
    }
  }
  if (xmax != 0) {
    for (i = 0; i < 2; i++) {
      for (j = 0; j < 2; j++) {
        M[i][j] = H[i][j] / xmax;
      }
    }
  }

  b     = -M[0][0] - M[1][1];
  c     = M[0][0] * M[1][1] - M[1][0] * M[1][0];
  root1 = -.5 * (b + sign(b) * Math.Sqrt(b * b - 4 * c));
  root2 = c / root1;

  root1 *= xmax;
  root2 *= xmax;
  *eval1 = Math.Min(root1, root2);
  *eval2 = Math.Max(root1, root2);
}

/* Solve for eigenvector of SPD 2x2 matrix, with given eigenvalue. */
public static void eigenvec2(double[][]  A/*[2][2]*/, /* matrix */
               double  eval,    /* eigenvalue */
               double[]  evec/*[2]*/, /* eigenvector returned */
               double *res      /* normalized residual */
)
{
  double norm;       /* norm of eigenvector */
  double res1, res2; /* components of residual vector */
  int    i;          /* loop counter */

  if (Math.Abs(A[0][0] - eval) > Math.Abs(A[1][1] - eval)) {
    evec[0] = -A[1][0];
    evec[1] = A[0][0] - eval;
  }
  else {
    evec[0] = A[1][1] - eval;
    evec[1] = -A[1][0];
  }

  /* Normalize eigenvector and calculate a normalized eigen-residual. */
  norm = Math.Sqrt(evec[0] * evec[0] + evec[1] * evec[1]);
  if (norm == 0) {
    evec[0] = 1;
    evec[1] = 0;
    norm    = 1;
  }
  for (i = 0; i < 2; i++) {
    evec[i] /= norm;
  }
  res1 = (A[0][0] - eval) * evec[0] + A[1][0] * evec[1];
  res2 = A[1][0] * evec[0] + (A[1][1] - eval) * evec[1];
  *res = Math.Sqrt(res1 * res1 + res2 * res2);

  res1 = Math.Abs(A[0][0]) + Math.Abs(A[1][0]);
  res2 = Math.Abs(A[1][1]) + Math.Abs(A[1][0]);
  *res /= Math.Max(res1, res2);
}
    }
}
