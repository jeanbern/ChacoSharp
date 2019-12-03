using System;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Eigen.SturmCount;
using static ChacoSharp.Eigen.CkSturmCount;

namespace ChacoSharp.Eigen
{
    public static unsafe class Bisect
    {
        /* Finds selected eigenvalues of T using Sturm sequence bisection. Based
    on Wilkinson's algorithm, AEP, p.302. Returns 1 if sturmcnt() fails and
    2 if it hasn't converged in max_steps of bisection. If neither of these
    errors is detected the return value is 0. */
public static int bisect(double *alpha,        /* vector of Lanczos scalars */
           double *beta,         /* vector of Lanczos scalars */
           int     j,            /* number of Lanczos iterations taken */
           double  Anorm,        /* Gershgorin estimate */
           double *workj,        /* work vector for Sturm sequence */
           double *ritz,         /* array holding evals */
           int     nevals_left,  /* number of evals on right to find */
           int     nevals_right, /* number of evals on left to find */
           double  tol,          /* tolerance on bracket width */
           double *ritz_sav,     /* space to copy ritzvals for debugging */
           int     max_steps     /* maximum number of bisection steps allowed */
)
{
  int           index;        /* index of sturm polynomial */
  int           i;            /* loop index */
  double *      pntr;         /* pntr to double array */
  double        x1, x2;       /* the bracketing interval */
  int           x1cnt, x2cnt; /* Sturm counts at x1 and x2 */
  double        x;            /* the inserted point */
  int           xcnt;         /* the Sturm count at x */
  int           steps;        /* number of bisection steps for a Ritzval */
  int           tot_steps;    /* number of bisection steps for all Ritzvals */
  int           numbracketed; /* number of evals between x1 and x2 */
  int           x1ck;         /* debugging check on x1cnt */
  int           x2ck;         /* debugging check on x2cnt */
  int           numck;        /* debugging check on numbracketed */
  double        diff;         /* debugging register */
  int           ii = 0;           /* debugging loop counter */

  /* If space has been allocated for a copy of the ritz values, assume
     we are to check the Sturm sequence counts directly using ql(). */
  if (ritz_sav != null) {
    Console.WriteLine("\nAnorm {0:g} j {1:d} nevals_left {2:d}", Anorm, j, nevals_left);
    Console.WriteLine("step              x1                 x2         x1cnt  ck  x2cnt  ck  brack   ck   x2-x1");
  }

  /* Initialize portion of ritz we will use (use max double so scanmin will work
     properly when called later on) */
  pntr = &ritz[1];
  for (i = j; i != 0; i--) {
    *pntr++ = double.MaxValue;
  }

  tot_steps = 0;

  /* find evals on left in decreasing index order */
  x2           = Anorm;
  x2cnt        = j;
  numbracketed = j;
  for (index = nevals_left; index >= 1; index--) {
    x1    = 0;
    x1cnt = 0;
    steps = 1; /* ... since started with Anorm bracketing j roots */
    while ((x2 - x1) > tol || numbracketed > 1) {
      x    = 0.5 * (x1 + x2);
      xcnt = sturmcnt(alpha, beta, j, x, workj);
      if (xcnt == -1) {
        return (1);
        /* ... sturmcnt() failed; bail out with error code */
      }
      if (xcnt >= index) {
        x2    = x;
        x2cnt = xcnt;
      }
      else {
        x1    = x;
        x1cnt = xcnt;
      }
      numbracketed = x2cnt - x1cnt;
      steps++;
      if (steps == max_steps) {
        return (2);
        /* ... not converging; bail out with error code */
      }

      if (ritz_sav != null) {
        diff = x2 - x1;
        cksturmcnt(ritz_sav, 1, j, x1, x2, &x1ck, &x2ck, &numck);
        Console.WriteLine("{0:d} {1:f} {2:f}   {3:d}   {4:d}  {5:d}   {6:d}   {7:d}   {8:d}   {9:g}", ii++, x1, x2, x1cnt,x1ck, x2cnt, x2ck, numbracketed, numck, diff);
        if (x1cnt != x1ck || x2cnt != x2ck || numbracketed != numck) {
            Console.WriteLine("**");
        }
        else {
            Console.WriteLine();
        }
      }
    }
    ritz[index] = 0.5 * (x1 + x2);
    if (ritz_sav != null) {
        Console.WriteLine("Ritzval #{0:d}:", index);
      Console.WriteLine("            bisection {0:f}", ritz[index]);
      Console.WriteLine("                   ql {0:f}", ritz_sav[index]);
      Console.WriteLine("           difference {0:f}", ritz[index] - ritz_sav[index]);
      Console.WriteLine("---------------------------------------------------");
    }
    if (DEBUG_EVECS > 2) {
      Console.WriteLine("    index {0:d}, bisection steps {1:d}, root {2:f}", index, steps, ritz[index]);
    }
    tot_steps += steps;
  }

  /* find evals on right in increasing index order */
  x1    = 0;
  x1cnt = 0;
  for (index = j - nevals_right + 1; index <= j; index++) {
    x2    = Anorm;
    x2cnt = j;
    steps = 1; /* ... since started with Anorm bracketing j roots */
    while ((x2 - x1) > tol || numbracketed > 1) {
      x    = 0.5 * (x1 + x2);
      xcnt = sturmcnt(alpha, beta, j, x, workj);
      if (xcnt == -1) {
        return (1);
        /* ... sturmcnt() failed; bail out with error code */
      }
      if (xcnt >= index) {
        x2    = x;
        x2cnt = xcnt;
      }
      else {
        x1    = x;
        x1cnt = xcnt;
      }
      numbracketed = x2cnt - x1cnt;
      steps++;
      if (steps == max_steps) {
        return (2);
        /* ... not converging; bail out with error code */
      }

      if (ritz_sav != null) {
        diff = x2 - x1;
        cksturmcnt(ritz_sav, 1, j, x1, x2, &x1ck, &x2ck, &numck);
        Console.WriteLine("{0:d} {1:f} {2:f}   {3:d}   {4:d}  {5:d}   {6:d}   {7:d}   {8:d}   {9:g}", ii++, x1, x2, x1cnt,x1ck, x2cnt, x2ck, numbracketed, numck, diff);
        if (x1cnt != x1ck || x2cnt != x2ck || numbracketed != numck) {
            Console.WriteLine("**");
        }
        else {
            Console.WriteLine("");
        }
      }
    }
    ritz[index] = 0.5 * (x1 + x2);
    if (ritz_sav != null) {
        Console.WriteLine("Ritzval #{0:d}:", index);
        Console.WriteLine("            bisection {0:f}", ritz[index]);
        Console.WriteLine("                   ql {0:f}", ritz_sav[index]);
        Console.WriteLine("           difference {0:f}", ritz[index] - ritz_sav[index]);
        Console.WriteLine("---------------------------------------------------");
    }
    if (DEBUG_EVECS > 2) {
        Console.WriteLine("    index {0:d}, bisection steps {1:d}, root {2:f}", index, steps, ritz[index]);
    }
    tot_steps += steps;
  }
  if (DEBUG_EVECS > 2) {
      Console.WriteLine("  Total number of bisection steps {0:d}.", tot_steps);
  }

  return (0); /* ... things seem ok. */
}
    }
}
