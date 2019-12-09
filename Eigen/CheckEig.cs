using System.Diagnostics;
using static ChacoSharp.Eigen.Splarax;
using static ChacoSharp.Utilities.Scadd;
using static ChacoSharp.Utilities.Norm;
using static ChacoSharp.Utilities.CpVec;
using static ChacoSharp.Utilities.VecScale;
using static ChacoSharp.StaticConstants;

namespace ChacoSharp.Eigen
{
    public static unsafe class CheckEig
    {
        /* Check an eigenpair of A by direct multiplication.  */
        public static double checkeig(double* err, vtx_data** A, double* y, int n, double lambda, double* vwsqrt, double* work)
        {
            splarax(err, A, n, y, vwsqrt, work);
            scadd(err, 1, n, -lambda, y);
            return ch_norm(err, 1, n) / ch_norm(y, 1, n);
        }

        /* Check an extended eigenpair of A by direct multiplication. Uses
   the Ay = extval*y + Dg form of the problem for convenience. */

        public static double checkeig_ext(double* err, double* work, /* work vector of length n */
            vtx_data** A, double* y, int n, double extval, double* vwsqrt,
            double* gvec, double eigtol,
            bool warnings /* don't want to see warning messages in one of the
                                    contexts this is called */
        )
        {
            double resid; /* the extended eigen residual */

            splarax(err, A, n, y, vwsqrt, work);
            scadd(err, 1, n, -extval, y);
            cpvec(work, 1, n, gvec); /* only need if going to re-use gvec */
            scale_diag(work, 1, n, vwsqrt);
            scadd(err, 1, n, -1.0, work);
            resid = ch_norm(err, 1, n);

            if (DEBUG_EVECS > 0)
            {
                {
                    Trace.WriteLine($"  extended residual: {resid:g}");
                }
            }

            if (warnings && WARNING_EVECS > 0 && resid > eigtol)
            {
                Trace.WriteLine($"WARNING: Extended residual ({resid:g}) greater than tolerance ({eigtol:g}).");
            }

            return (resid);
        }
    }
}
