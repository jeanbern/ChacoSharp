#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
using System;
using System.Diagnostics;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Utilities.MkVec;
using static ChacoSharp.Utilities.Dot;
using static ChacoSharp.Utilities.Scadd;
using static ChacoSharp.Utilities.Perturb;
using static ChacoSharp.Eigen.Splarax;
using static ChacoSharp.Utilities.VecRan;
using static ChacoSharp.Utilities.Norm;
using static ChacoSharp.Utilities.Update;

namespace ChacoSharp.Utilities
{
    public static unsafe class Timer
    {
        public static int seconds()
        {
            return (int) (DateTime.Now.Ticks / 10000000);
        }

        /// <summary>
        /// Benchmark certain kernel operations
        /// </summary>
        /// <param name="A">matrix/graph being analyzed</param>
        /// <param name="n">number of rows/columns in matrix</param>
        /// <param name="vwsqrt">square roots of vertex weights</param>
        public static void time_kernels(vtx_data** A, int n, double* vwsqrt)
        {
            const double minTime = 0.5;
            const double targetTime = 1.0;

            int i;
            float* vwsqrtFloat;
            double time;

            if (DEBUG_TRACE)
            {
                Trace.WriteLine("<Entering time_kernels>");
            }

            const int beg = 1;
            var end = n;

            var dvec1 = mkvec(beg, end);
            var dvec2 = mkvec(beg, end);
            var dvec3 = mkvec(beg - 1, end);
            var svec1 = mkvec_float(beg, end);
            var svec2 = mkvec_float(beg, end);
            var svec3 = mkvec_float(beg - 1, end);

            if (vwsqrt == null)
            {
                vwsqrtFloat = null;
            }
            else
            {
                vwsqrtFloat = mkvec_float(beg - 1, end);
                for (i = beg - 1; i <= end; i++)
                {
                    vwsqrtFloat[i] = (float) vwsqrt[i];
                }
            }

            vecran(dvec1, beg, end);
            vecran(dvec2, beg, end);
            vecran(dvec3, beg, end);
            for (i = beg; i <= end; i++)
            {
                svec1[i] = (float) dvec1[i];
                svec2[i] = (float) dvec2[i];
                svec3[i] = (float) dvec3[i];
            }

            // Set number of loops so that ch_norm() takes about one second.
            // This should insulate against inaccurate timings on faster machines.

            var normDvec = 0.0d;
            var loops = 1;
            double timeDvec = 0;
            while (timeDvec < minTime)
            {
                time = seconds();
                for (i = loops; i != 0; i--)
                {
                    normDvec = ch_norm(dvec1, beg, end);
                }

                timeDvec = seconds() - time;
                if (timeDvec < minTime)
                {
                    loops = 10 * loops;
                }
            }

            loops = (int) ((targetTime / timeDvec) * loops);
            if (loops < 1)
            {
                loops = 1;
            }

            Trace.WriteLine("                Kernel benchmarking");
            Trace.WriteLine($"Time (in seconds) for {loops:d} loops of each operation:\n");

            Trace.WriteLine("Routine      Double     Float      Discrepancy      Description");
            Trace.WriteLine("-------      ------     -----      -----------      -----------");

            /* Norm operation */
            time = seconds();
            for (i = loops; i != 0; i--)
            {
                normDvec = ch_norm(dvec1, beg, end);
            }

            timeDvec = seconds() - time;

            time = seconds();
            var normSvec = 0.0d;
            for (i = loops; i != 0; i--)
            {
                normSvec = norm_float(svec1, beg, end);
            }

            var timeSvec = seconds() - time;

            var diff = normDvec - normSvec;
            Trace.Write($"norm        {timeDvec:f}    {timeSvec:f}    {diff:e}");
            Trace.WriteLine("      2 norm");

            /* Dot operation */
            time = seconds();
            var dotDvec = 0.0d;
            for (i = loops; i != 0; i--)
            {
                dotDvec = dot(dvec1, beg, end, dvec2);
            }

            timeDvec = seconds() - time;

            time = seconds();
            var dotSvec = 0.0d;
            for (i = loops; i != 0; i--)
            {
                dotSvec = dot_float(svec1, beg, end, svec2);
            }

            timeSvec = seconds() - time;

            diff = dotDvec - dotSvec;
            Trace.Write($"dot         {timeDvec:f}    {timeSvec:f}    {diff:e}");
            Trace.WriteLine("      scalar product");

            /* Scadd operation */
            const double factor = 1.01d;
            const float factorFloat = (float) factor;

            var fac = factor;
            time = seconds();
            for (i = loops; i != 0; i--)
            {
                scadd(dvec1, beg, end, fac, dvec2);
                fac = -fac; /* to keep things in scale */
            }

            timeDvec = seconds() - time;

            var facFloat = factorFloat;
            time = seconds();
            for (i = loops; i != 0; i--)
            {
                scadd_float(svec1, beg, end, facFloat, svec2);
                facFloat = -facFloat; /* to keep things in scale */
            }

            timeSvec = seconds() - time;

            diff = checkvec(dvec1, beg, end, svec1);
            Trace.Write($"scadd       {timeDvec:f}    {timeSvec:f}    {diff:e}");
            Trace.WriteLine("      vec1 <- vec1 + alpha*vec2");

            /* Update operation */
            time = seconds();
            for (i = loops; i != 0; i--)
            {
                update(dvec1, beg, end, dvec2, factor, dvec3);
            }

            timeDvec = seconds() - time;

            time = seconds();
            for (i = loops; i != 0; i--)
            {
                update_float(svec1, beg, end, svec2, factorFloat, svec3);
            }

            timeSvec = seconds() - time;

            diff = checkvec(dvec1, beg, end, svec1);
            Trace.Write($"update      {timeDvec:f}    {timeSvec:f}    {diff:g}");
            Trace.WriteLine("      vec1 <- vec2 + alpha*vec3");

            /* splarax operation */
            if (PERTURB)
            {
                if (NPERTURB > 0 && PERTURB_MAX > 0.0)
                {
                    perturb_init(n);
                    if (DEBUG_PERTURB)
                    {
                        Trace.WriteLine($"Matrix being perturbed with scale {PERTURB_MAX:e}");
                    }
                }
                else if (DEBUG_PERTURB)
                {
                    Trace.WriteLine("Matrix not being perturbed");
                }
            }

            time = seconds();
            for (i = loops; i != 0; i--)
            {
                splarax(dvec1, A, n, dvec2, vwsqrt, dvec3);
            }

            timeDvec = seconds() - time;

            time = seconds();
            for (i = loops; i != 0; i--)
            {
                splarax_float(svec1, A, n, svec2, vwsqrtFloat, svec3);
            }

            timeSvec = seconds() - time;

            diff = checkvec(dvec1, beg, end, svec1);
            Trace.Write($"splarax     {timeDvec:f}    {timeSvec:f}    {diff:e}");
            Trace.WriteLine("      sparse matrix vector multiply");

            if (PERTURB && NPERTURB > 0 && PERTURB_MAX > 0.0)
            {
                perturb_clear();
            }

            Trace.WriteLine("");

            /* Free memory */
            frvec(dvec1, 1);
            frvec(dvec2, 1);
            frvec(dvec3, 0);
            frvec_float(svec1, 1);
            frvec_float(svec2, 1);
            frvec_float(svec3, 0);
            if (vwsqrtFloat != null)
            {
                frvec_float(vwsqrtFloat, beg - 1);
            }
        }

        /// <summary>
        /// Compute norm of difference between a double and float vector.
        /// </summary>
        static double checkvec(double* dvec, int beg, int end, float* svec)
        {
            int i;

            double sum = 0;
            for (i = beg; i <= end; i++)
            {
                var diff = dvec[i] - svec[i];
                sum += diff * diff;
            }

            return Math.Sqrt(sum);
        }
    }
}
