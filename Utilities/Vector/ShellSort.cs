namespace ChacoSharp.Utilities
{
    public static unsafe class ShellSort
    {
        public static void sort_double(int count, double* ra)
        {
            int start, end;

            /* heapify */
            for (start = (count - 2) / 2; start >= 0; start--)
            {
                siftDown(ra, start, count);
            }

            for (end = count - 1; end > 0; end--)
            {
                var temp = ra[end];
                ra[end] = ra[0];
                ra[0] = temp;
                siftDown(ra, 0, end);
            }
        }

        static void siftDown(double* a, int start, int end)
        {
            int root = start;

            while (root * 2 + 1 < end)
            {
                int child = 2 * root + 1;
                if ((child + 1 < end) && (a[child] < a[child + 1]))
                {
                    child += 1;
                }

                if (a[root] < a[child])
                {
                    var temp = a[child];
                    a[child] = a[root];
                    a[root] = temp;
                    root = child;
                }
                else
                {
                    return;
                }
            }
        }

        public static void shell_sort(int n, double* arr)
        {
            sort_double(n, arr);
        }
    }
}