using System;
using static ChacoSharp.Utilities.Dot;

namespace ChacoSharp.Utilities
{
    public static unsafe class Norm
    {
        /// <summary>
        /// Print vertically range of double vector.
        /// </summary>
        public static void vecout(double* vec, int beg, int end, string tag, string file_name)
        {
            int i;

            Console.WriteLine(tag);
            for (i = beg; i <= end; i++)
            {
                if (Math.Abs(vec[i]) >= 1.0e-16)
                {
                    Console.WriteLine("{0:d}.   {1:f}", i, vec[i]);
                }
                else
                {
                    Console.WriteLine("{0:d}.         {1:g}", i, vec[i]);
                }
            }
        }

        /// <summary>
        /// Scale the eigenvector such that the first element is non-negative.
        /// This is an attempt to reduce machine-to-machine variance which can result from calculated eigenvectors being mirror-images of each other due to small roundoff..
        /// </summary>
        public static void vecnorm(double* vec, int beg, int end)
        {
            // ReSharper disable once InvertIf
            if (vec[beg] < 0.0)
            {
                for (var i = beg; i <= end; i++)
                {
                    vec[i] *= -1.0;
                }
            }
        }

        /// <summary>
        /// Normalizes a double n-vector over range.
        /// </summary>
        public static double ch_normalize(double* vec, int beg, int end)
        {
            var scale = ch_norm(vec, beg, end);
            vec += beg;
            for (var i = end - beg + 1; i != 0; i--)
            {
                *vec /= scale;
                vec++;
            }

            return scale;
        }

        /// <summary>
        /// Normalizes such that element k is positive
        /// </summary>
        public static double sign_normalize(double* vec, int beg, int end, int k)
        {
            double scale2;
            var scale = ch_norm(vec, beg, end);
            if (vec[k] < 0)
            {
                scale2 = -scale;
            }
            else
            {
                scale2 = scale;
            }

            vec += beg;
            for (var i = end - beg + 1; i != 0; i--)
            {
                *vec /= scale2;
                vec++;
            }

            return scale;
        }

        /// <summary>
        /// Normalizes a float n-vector over range.
        /// </summary>
        public static double normalize_float(float* vec, int beg, int end)
        {
            var scale = (float) norm_float(vec, beg, end);
            vec += beg;
            for (var i = end - beg + 1; i != 0; i--)
            {
                *vec /= scale;
                vec++;
            }

            return scale;
        }

        /// <summary>
        /// Returns 2-norm of a double n-vector over range.
        /// </summary>
        public static double ch_norm(double* vec, int beg, int end)
        {
            return Math.Sqrt(dot(vec, beg, end, vec));
        }

        /// <summary>
        /// Returns 2-norm of a float n-vector over range.
        /// </summary>
        public static double norm_float(float* vec, int beg, int end)
        {
            return Math.Sqrt(dot_float(vec, beg, end, vec));
        }
    }
}
