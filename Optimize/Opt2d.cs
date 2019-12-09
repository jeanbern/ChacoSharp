#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
using System;
using System.Diagnostics;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Optimize.Func2d;

namespace ChacoSharp.Optimize
{
    public static unsafe class Opt2d
    {
        private static double Sign(double x)
        {
            return x < 0.0 ? -1.0 : 1.0;
        }

        /* Compute rotation angle to minimize distance to discrete points. */
        public static double opt2d(
            vtx_data** graph, /* data structure with vertex weights */
            double*[] yvecs, /* eigenvectors */
            int nvtxs, /* total number of vertices */
            int nmyvtxs /* number of vertices I own */
        )
        {
            double* aptr; /* loop through yvecs */
            double* bptr; /* loop through yvecs */
            double[] coeffs = new double[5]; /* various products of yvecs */
            double func = 0.0d; /* value of function to be minimized */
            double grad, hess; /* first and 2nd derivatives of function */
            double grad_min; /* acceptably small gradient */
            double theta; /* angle being optimized */
            double step; /* change in angle */
            double step_max; /* maximum allowed step */
            double step_min; /* minimum step => convergence */
            double hess_min; /* value for hessian is < 0 */
            double hfact; /* scaling for min tolerated hessian */
            double w; /* vertex weight squared */
            double pdtol; /* allowed error in hessian pd-ness */
            bool pdflag; /* is hessian positive semi-definite? */

            /* Set parameters. */
            step_max = Math.PI / 8;
            step_min = 2.0e-5;
            grad_min = 1.0e-7;
            pdtol = 1.0e-8;
            hfact = 2;

            for (var i = 0; i < 5; i++)
            {
                coeffs[i] = 0;
            }

            aptr = yvecs[1] + 1;
            bptr = yvecs[2] + 1;
            for (var i = 1; i <= nmyvtxs; i++)
            {
                var a = *aptr++; /* temporary values */
                var b = *bptr++; /* temporary values */
                w = graph[i]->vwgt;
                if (w == 1)
                {
                    coeffs[0] += a * a * a * a;
                    coeffs[1] += a * a * a * b;
                    coeffs[2] += a * a * b * b;
                    coeffs[3] += a * b * b * b;
                    coeffs[4] += b * b * b * b;
                }
                else
                {
                    w = 1 / (w * w);
                    coeffs[0] += a * a * a * a * w;
                    coeffs[1] += a * a * a * b * w;
                    coeffs[2] += a * a * b * b * w;
                    coeffs[3] += a * b * b * b * w;
                    coeffs[4] += b * b * b * b * w;
                }
            }

            /* Adjust for normalization of eigenvectors. */
            /* This should make tolerances independent of vector length */
            for (var i = 0; i < 5; i++)
            {
                coeffs[i] *= nvtxs;
            }

            var passes = 0;
            theta = 0.0;
            step = step_max;
            pdflag = false;
            grad = 0;
            while (Math.Abs(step) >= step_min && (!pdflag || Math.Abs(grad) > grad_min))
            {
                func = func2d(coeffs, theta);
                grad = grad2d(coeffs, theta);
                hess = hess2d(coeffs);

                if (hess < -pdtol)
                {
                    pdflag = false;
                }
                else
                {
                    pdflag = true;
                }

                hess_min = hfact * Math.Abs(grad) / step_max;
                if (hess < hess_min)
                {
                    hess = hess_min;
                }

                if (Math.Abs(grad) > Math.Abs(hess * step_max))
                {
                    step = -step_max * Sign(grad);
                }
                else
                {
                    step = -grad / hess;
                }

                theta += step;
                if (Math.Abs(step) < step_min && !pdflag)
                {
                    /* Convergence to non-min. */
                    step = step_min;
                    theta += step;
                }

                passes++;
            }

            if (DEBUG_OPTIMIZE)
            {
                Trace.WriteLine($"After {passes:d} passes, func={func:e}, theta = {theta:f}");
            }

            return theta;
        }
    }
}
