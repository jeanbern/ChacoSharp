using System;
using static ChacoSharp.Utilities.Norm;

namespace ChacoSharp.Eigen
{
    public static unsafe class BiDir
    {
        /* NOTE: This should only be called if j >= 2. It returns a residual and
         the point where the forward and backward recurrences met. */

/* NOTE: The paper has an error in the residual calculation. And the
         heuristic of increasing the cut-off for the local maximum
         by a factor of 10 when the residual tolerance is not met
         appears to work poorly - letting the recursion run further
         to find a bigger local max often results in <less> accurate
         results. Static cut-off of 1 (recall we set the first element
         to 1) worked better. */

/* Finds eigenvector s of T and returns residual norm. */
        public static double bidir(double* alpha, /* vector of Lanczos scalars */
            double* beta, /* vector of Lanczos scalars */
            int j, /* number of Lanczos iterations taken */
            double ritz, /* approximate eigenvalue  of T */
            double* s, /* approximate eigenvector of T */
            double hurdle /* hurdle for local maximum in recurrence */
        )

        {
            int i; /* index */
            int k = 0; /* meeting point for bidirect. recurrence */
            double sav_val = 0.0; /* value at meeting point */
            bool iterate; /* keep iterating? */
            double scale; /* scale factor for bidirectional recurrence */
            double residual; /* eigen residual */

            s[2] = -(alpha[1] - ritz) / beta[2];
            i = 3;
            iterate = true;
            while (i <= j && iterate)
            {
                /* follow the forward recurrence until a local maximum > 1.0 */
                s[i] = -((alpha[i - 1] - ritz) * s[i - 1] + beta[i - 1] * s[i - 2]) / beta[i];
                if (Math.Abs(s[i - 1]) > hurdle && i > 2)
                {
                    if (Math.Abs(s[i]) < Math.Abs(s[i - 1]) && Math.Abs(s[i - 1]) > Math.Abs(s[i - 2]))
                    {
                        iterate = false;
                        k = i - 1;
                        sav_val = s[k];
                    }
                }

                i++;
            }

            if (!iterate)
            {
                /* we stopped at an intermediate point */
                s[j] = 1.0;
                s[j - 1] = -(alpha[j] - ritz) / beta[j];
                for (i = j; i >= k + 2; i--)
                {
                    s[i - 2] = -((alpha[i - 1] - ritz) * s[i - 1] + beta[i] * s[i]) / beta[i - 1];
                }

                scale = s[k] / sav_val;
                for (i = 1; i <= k - 1; i++)
                {
                    s[i] = s[i] * scale;
                }

                /* paper gets this expression wrong */
                residual = beta[k] * s[k - 1] + (alpha[k] - ritz) * s[k] + beta[k + 1] * s[k + 1];
            }
            else
            {
                /* we went all the way forward */
                residual = (alpha[j] - ritz) * s[j] + beta[j] * s[j - 1];
            }

            residual = Math.Abs(residual) / sign_normalize(s, 1, j, j);
            return (residual);
        }
    }
}
