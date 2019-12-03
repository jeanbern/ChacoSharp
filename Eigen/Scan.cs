using System.Runtime.InteropServices;

namespace ChacoSharp.Eigen
{
    public static unsafe class Scan
    {
        public static scanlink* mkscanlist(int depth)
        {
            scanlink* prevlnk;
            scanlink* newlnk;
            int i;

            prevlnk = (scanlink*) Marshal.AllocHGlobal(sizeof(scanlink));
            prevlnk->pntr = null;
            newlnk = prevlnk; /* in case the list is one long */
            for (i = 1; i <= (depth - 1); i++)
            {
                newlnk = (scanlink*) Marshal.AllocHGlobal(sizeof(scanlink));
                newlnk->pntr = prevlnk;
                newlnk->indx = -1;
                prevlnk = newlnk;
            }

            return (newlnk);
        }

        /* Return minimum of vector over range */
        public static void scanmin(double* vec, /* vector to scan */
            int beg, int end, /* index range */
            scanlink** scanlist /* pntr to list holding results of scan */
        )
        {
            scanlink* top;
            scanlink* curlnk;
            scanlink* prevlnk;
            double val;
            int i;

            curlnk = *scanlist;
            while (curlnk != null)
            {
                curlnk->indx = 0;
                curlnk->val = double.MaxValue;
                curlnk = curlnk->pntr;
            }

            /* Note: Uses current top link (which would need to be deleted anyway) each time
                     an insertion to the list is required. */

            for (i = beg; i <= end; i++)
            {
                /* consider each element for insertion */
                top = *scanlist;
                val = vec[i];
                if (val < top->val)
                {
                    if (top->pntr == null)
                    {
                        /* the list is only one long, so just replace */
                        top->val = val;
                        top->indx = i;
                    }
                    else
                    {
                        /* beats top element; scan for insertion point */
                        if (val < (top->pntr)->val)
                        {
                            /* 2nd link becomes list pntr; otherwise stays same */
                            *scanlist = top->pntr;
                        }

                        prevlnk = curlnk = top;
                        while ((val < curlnk->val) && (curlnk->pntr != null))
                        {
                            prevlnk = curlnk;
                            curlnk = curlnk->pntr;
                        }

                        if (val < curlnk->val)
                        {
                            /* got to end of list; add top to bottom */
                            curlnk->pntr = top;
                            top->val = val;
                            top->indx = i;
                            top->pntr = null;
                        }
                        else
                        {
                            /* stopped within list; insert top here */
                            prevlnk->pntr = top;
                            top->val = val;
                            top->indx = i;
                            top->pntr = curlnk;
                        }
                    }
                }
            }
        }


        /* Return maximum entries of vector over range */
        public static void scanmax(double* vec, /* vector to scan */
            int beg, int end, /* index range */
            scanlink** scanlist /* pntr to list holding results of scan */
        )
        {
            scanlink* top;
            scanlink* curlnk;
            scanlink* prevlnk;
            double val;
            int i;

            curlnk = *scanlist;
            while (curlnk != null)
            {
                curlnk->indx = 0;
                curlnk->val = double.MinValue;
                curlnk = curlnk->pntr;
            }

            /* Note: Uses current top link (which would need to be deleted anyway) each time
                     an insertion to the list is required. */

            for (i = beg; i <= end; i++)
            {
                /* consider each element for insertion */
                top = *scanlist;
                val = vec[i];
                if (val > top->val)
                {
                    if (top->pntr == null)
                    {
                        /* the list is only one long, so just replace */
                        top->val = val;
                        top->indx = i;
                    }
                    else
                    {
                        /* beats top element; scan for insertion point */
                        if (val > (top->pntr)->val)
                        {
                            /* 2nd link becomes list pntr; otherwise stays same */
                            *scanlist = top->pntr;
                        }

                        prevlnk = curlnk = top;
                        while ((val > curlnk->val) && (curlnk->pntr != null))
                        {
                            prevlnk = curlnk;
                            curlnk = curlnk->pntr;
                        }

                        if (val > curlnk->val)
                        {
                            /* got to end of list; add top to bottom */
                            curlnk->pntr = top;
                            top->val = val;
                            top->indx = i;
                            top->pntr = null;
                        }
                        else
                        {
                            /* stopped within list; insert top here */
                            prevlnk->pntr = top;
                            top->val = val;
                            top->indx = i;
                            top->pntr = curlnk;
                        }
                    }
                }
            }
        }
    }
}
