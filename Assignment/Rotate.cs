using System;

namespace ChacoSharp.Assignment
{
    public static unsafe class Rotate
    {
        public static void Rotate2d(
            double*[] yvecs, /* ptr to list of y-vectors (rotated) */
            int nmyvtxs, /* length of yvecs */
            double theta /* angle to rotate by */
        )
        {
            var sinTheta = Math.Sin(theta);
            var cosTheta = Math.Cos(theta);

            for (var i = 1; i <= nmyvtxs; i++)
            {
                var temp1 = yvecs[1][i];
                yvecs[1][i] = cosTheta * temp1 + sinTheta * yvecs[2][i];
                yvecs[2][i] = -sinTheta * temp1 + cosTheta * yvecs[2][i];
            }
        }

        public static void Rotate3d(
            double*[] yvecs, /* ptr to list of y-vectors (to be rotated) */
            int nmyvtxs, /* length of yvecs */
            double theta, double phi, double gamma2 /* rotational parameters */
        )
        {
            var stheta = Math.Sin(theta); /* sine of theta */
            var ctheta = Math.Cos(theta); /* cosine of theta */
            var sphi = Math.Sin(phi); /* sine of phi */
            var cphi = Math.Cos(phi); /* cosine of phi */
            var sgamma = Math.Sin(gamma2); /* sine of gamma */
            var cgamma = Math.Cos(gamma2); /* cosine of gamma */

            var onemcg = 1.0 - cgamma; /* 1.0 - cosine(gamma) */

            /* rotation matrix entries */
            var a1 = cgamma + cphi * ctheta * onemcg * cphi * ctheta;
            var a2 = sgamma * sphi + cphi * stheta * onemcg * cphi * ctheta;
            var a3 = -sgamma * cphi * stheta + sphi * onemcg * cphi * ctheta;

            /* rotation matrix entries */
            var b1 = -sgamma * sphi + cphi * ctheta * onemcg * cphi * stheta;
            var b2 = cgamma + cphi * stheta * onemcg * cphi * stheta;
            var b3 = sgamma * cphi * ctheta + sphi * onemcg * cphi * stheta;

            /* rotation matrix entries */
            var c1 = sgamma * cphi * stheta + cphi * ctheta * onemcg * sphi;
            var c2 = -sgamma * cphi * ctheta + cphi * stheta * onemcg * sphi;
            var c3 = cgamma + sphi * onemcg * sphi;

            for (var i = 1; i <= nmyvtxs; i++)
            {
                var temp1 = yvecs[1][i];
                var temp2 = yvecs[2][i];

                yvecs[1][i] = a1 * temp1 + b1 * temp2 + c1 * yvecs[3][i];
                yvecs[2][i] = a2 * temp1 + b2 * temp2 + c2 * yvecs[3][i];
                yvecs[3][i] = a3 * temp1 + b3 * temp2 + c3 * yvecs[3][i];
            }
        }
    }
}
