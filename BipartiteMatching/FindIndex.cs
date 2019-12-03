namespace ChacoSharp.BipartiteMatching
{
    public static unsafe class FindIndex
    {
        public static int findindex(int* indices, /* indices sorting values */
            double* vals, /* values sorted by indices */
            double target, /* target value */
            int nvals /* number of values */
        )
        {
            if (target <= vals[indices[0]])
            {
                return (0);
            }

            if (target >= vals[indices[nvals - 1]])
            {
                return nvals - 1;
            }

            var low = 0;/* range left to search */
            var high = nvals - 1;/* range left to search */

            while (high - low > 1)
            {
                var vlow = vals[indices[low]]; /* values at limits of search range */
                var vhigh = vals[indices[high]]; /* values at limits of search range */
                if (vlow == vhigh)
                {
                    return (int) ((vlow + vhigh) / 2);
                }

                var ratio = (target - vlow) / (vhigh - vlow); /* interpolation parameter */
                var newLimit = (int) (low + ratio * (high - low)); /* new index limit */
                if (newLimit == low)
                {
                    ++newLimit;
                }
                else if (newLimit == high)
                {
                    --newLimit;
                }

                if (vals[indices[newLimit]] < target)
                {
                    low = newLimit;
                }
                else
                {
                    high = newLimit;
                }
            }

            return target == vals[indices[high]] ? high : low;
        }
    }
}
