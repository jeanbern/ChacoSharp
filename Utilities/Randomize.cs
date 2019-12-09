using System;
using System.Diagnostics;

namespace ChacoSharp.Utilities
{
    public static unsafe class Randomize
    {
        /// <summary>
        /// Randomly permute elements of an array.
        /// </summary>
        /// <param name="array">array of integer values</param>
        /// <param name="n">number of values</param>
        public static void randomize(int* array, int n)
        {
            for (var i = 1; i <= n; i++)
            {
                var value = drandom();
                var index = (int)(n * value) + 1;
                var temp = array[i];
                array[i] = array[index];
                array[index] = temp;
            }
        }

        private static Random _rand = new Random();

        public static double drandom()
        {
            return _rand.NextDouble();
        }

        public static void setrandom(int seed)
        {
            Trace.WriteLine($"Setting random seed to: {seed}");
            _rand = new Random(seed);
        }
    }
}
