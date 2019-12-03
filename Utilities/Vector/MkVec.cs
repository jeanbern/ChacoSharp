using System;
using System.Runtime.InteropServices;

namespace ChacoSharp.Utilities
{
    public static unsafe class MkVec
    {
        /* Allocates a double vector with range [nl..nh]. Dies. */
        public static double* mkvec(int nl, int nh)
        {
            var v = (double*) Marshal.AllocHGlobal((nh - nl + 1) * sizeof(double));
            return v - nl;
        }

        /* Allocates a double vector with range [nl..nh]. Returns error code. */
        public static double* mkvec_ret(int nl, int nh)
        {
            var v = (double*) Marshal.AllocHGlobal((nh - nl + 1) * sizeof(double));
            if (v == null)
            {
                throw new InvalidOperationException();
            }

            return v - nl;
        }

        /* Free a double vector with range [nl..nh]. */
        public static void frvec(double* v, int nl)
        {
            Marshal.FreeHGlobal((IntPtr) (v + nl));
        }

        /* Allocates a float vector with range [nl..nh]. Dies. */
        public static float* mkvec_float(int nl, int nh)
        {
            var v = (float*) Marshal.AllocHGlobal((nh - nl + 1) * sizeof(float));
            return v - nl;
        }

        /* Allocates a float vector with range [nl..nh]. Returns error code. */
        public static float* mkvec_ret_float(int nl, int nh)
        {
            var v = (float*) Marshal.AllocHGlobal((nh - nl + 1) * sizeof(float));
            if (v == null)
            {
                throw new InvalidOperationException();
            }

            return v - nl;
        }

        /* Free a float vector with range [nl..nh]. */
        public static void frvec_float(float* v, int nl)
        {
            Marshal.FreeHGlobal((IntPtr) (v + nl));
        }
    }
}
