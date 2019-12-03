using System;
using System.Runtime.InteropServices;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Utilities.Randomize;

namespace ChacoSharp.Utilities
{
    public static unsafe class Perturb
    {
        private static ipairs* pedges; /* perturbed edges */
        private static double* pvals; /* perturbed values */

        /* Inititialize the perturbation */
        public static void perturb_init(
            int n /* graph size at this level */
        )
        {
            int i, j; /* loop counter */

            /* Initialize the diagonal perturbation weights */
            pedges = (ipairs*) Marshal.AllocHGlobal(NPERTURB * sizeof(ipairs));
            pvals = (double*) Marshal.AllocHGlobal(NPERTURB * sizeof(double));

            if (n <= 1)
            {
                for (i = 0; i < NPERTURB; i++)
                {
                    pedges[i].val1 = pedges[i].val2 = 0;
                    pvals[i] = 0;
                }

                return;
            }

            for (i = 0; i < NPERTURB; i++)
            {
                pedges[i].val1 = 1 + (int)(n * drandom());

                /* Find another vertex to define an edge. */
                j = 1 + (int)(n * drandom());
                while (j == i)
                {
                    j = 1 + (int)(n * drandom());
                }

                pedges[i].val2 = 1 + (int)(n * drandom());

                pvals[i] = PERTURB_MAX * drandom();
            }
        }

        public static void perturb_clear()
        {
            Marshal.FreeHGlobal((IntPtr)pedges);
            Marshal.FreeHGlobal((IntPtr)pvals);
            pedges = null;
            pvals = null;
        }

        /* Modify the result of splarax to break any graph symmetry */
        public static void perturb(double* result, /* result of matrix-vector multiply */
            double* vec /* vector matrix multiplies */
        )
        {
            int i; /* loop counter */

            for (i = 0; i < NPERTURB; i++)
            {
                result[pedges[i].val1] += pvals[i] * (vec[pedges[i].val2] - vec[pedges[i].val1]);
                result[pedges[i].val2] += pvals[i] * (vec[pedges[i].val1] - vec[pedges[i].val2]);
            }
        }

        /* Modify the result of splarax to break any graph symmetry, float version */
        public static void perturb_float(float* result, /* result of matrix-vector multiply */
            float* vec /* vector matrix multiplies */
        )
        {
            int i; /* loop counter */

            for (i = 0; i < NPERTURB; i++)
            {
                result[pedges[i].val1] += (float) (pvals[i] * (vec[pedges[i].val2] - vec[pedges[i].val1]));
                result[pedges[i].val2] += (float) (pvals[i] * (vec[pedges[i].val1] - vec[pedges[i].val2]));
            }
        }

    }
}
