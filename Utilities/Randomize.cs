using System;

namespace ChacoSharp.Utilities
{
    public static unsafe class Randomize
    {
        /* Randomly permute elements of an array. */
        public static void randomize(int* array, int n)
/* array of integer values */
/* number of values */
        {
            int i; /* loop counter */

            for (i = 1; i <= n; i++)
            {
                double value = drandom();
                int index = (int)(n * value) + 1;
                int temp = array[i];
                array[i] = array[index];
                array[index] = temp;
            }
        }

        private static Random rand = new Random();

        public static double drandom()
        {
            return rand.NextDouble();
        }

        public static void setrandom(int seed)
        {
            Console.WriteLine("Setting random to: " + seed);
            rand = new Random(seed);
        }
    }
}
