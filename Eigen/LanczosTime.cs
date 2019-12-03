using System;
using static ChacoSharp.Utilities.Timer;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Eigen.Orthogonalization;
using static ChacoSharp.Utilities.Dot;

namespace ChacoSharp.Eigen
{
    public static unsafe class LanczosTime
    {
        public static double lanc_seconds()
        {
            if (LANCZOS_TIME)
            {
                return seconds();
            }

            return 0;
        }

        /* Determine whether to pause in Lanczos */
        //TODO: remove this from LanczosTime, it has nothing to do with time.
public static bool lanpause(int      j,         /* current step */
             int      lastpause, /* when last paused */
             int      interval,  /* interval between pauses */
             double **q,         /* the Lanczos vectors */
             int      n,         /* length of Lanczos vectors */
             int *    pausemode, /* which pausing criterion to use */
             int      version,   /* which version of sel. orth. we are using */
             double   beta       /* current off-diagonal value */
)
{
  double        paige_dot;      /* q[j]^T q[1] */
  double        paigetol;       /* pause if paigedot > paigetol */

  /* Check orthogonality of last Lanczos vector against previous ones */
  if (DEBUG_EVECS > 3) {
    checkorth(q, n, j);
  }

  /* periodic reorthogonalization */
  if (version == 1 || version == 2) {
    if ((j - lastpause) == interval || beta < 1000 * DOUBLE_EPSILON) {
      return (true);
    }

    return (false);
  }

  /* Run until orthogonality with first Lanczos vector deteriorates, then switch
     switch to periodic reorthog. */
  if (version == 3) {
    paigetol = 1.0e-3;
    if (*pausemode == 1) {
      paige_dot = Math.Abs(dot(q[1], 1, n, q[j]));
      if ((paige_dot > paigetol && j > 1) || beta < 1000 * DOUBLE_EPSILON) {
        if (DEBUG_EVECS > 1) {
          Console.WriteLine("  Pausing on step {0:d} with Paige prod. = {1:g}", j, paige_dot);
        }
        *pausemode = 2;
        return (true);
      }

      return (false);
    }
    if (*pausemode == 2) {
      if ((j - lastpause) == interval || beta < 1000 * DOUBLE_EPSILON) {
        return (true);
      }

      return (false);
    }
  }

  /* shouldn't ever get this far, but alint really wants a return value */
  return (false);
}

public static bool lanpause_float(int     j,         /* current step */
                   int     lastpause, /* when last paused */
                   int     interval,  /* interval between pauses */
                   float **q,         /* the Lanczos vectors */
                   int     n,         /* length of Lanczos vectors */
                   int *   pausemode, /* which pausing criterion to use */
                   int     version,   /* which version of sel. orth. we are using */
                   double  beta       /* current off-diagonal value */
)
{
    double        paige_dot;      /* q[j]^T q[1] */
    double        paigetol;       /* pause if paigedot > paigetol */

  /* Check orthogonality of last Lanczos vector against previous ones */
  if (DEBUG_EVECS > 3) {
    checkorth_float(q, n, j);
  }

  /* periodic reorthogonalization */
  if (version == 1 || version == 2) {
    if ((j - lastpause) == interval || beta < 1000 * DOUBLE_EPSILON) {
      return (true);
    }

    return (false);
  }

  /* Run until orthogonality with first Lanczos vector deteriorates, then switch
     switch to periodic reorthog. */
  if (version == 3) {
    paigetol = 1.0e-3;
    if (*pausemode == 1) {
      paige_dot = Math.Abs(dot_float(q[1], 1, n, q[j]));
      if ((paige_dot > paigetol && j > 1) || beta < 1000 * DOUBLE_EPSILON) {
        if (DEBUG_EVECS > 1) {
          Console.WriteLine("  Pausing on step {0:d} with Paige prod. = {1:g}", j, paige_dot);
        }
        *pausemode = 2;
        return (true);
      }

      return (false);
    }
    if (*pausemode == 2) {
      if ((j - lastpause) == interval || beta < 1000 * DOUBLE_EPSILON) {
        return (true);
      }

      return (false);
    }
  }

  /* shouldn't ever get this far, but alint really wants a return value */
  return (false);
}

    }
}
