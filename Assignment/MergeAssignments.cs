using System;
using System.Collections.Generic;
using System.Text;

namespace ChacoSharp.Assignment
{
    public static unsafe class MergeAssignments
    {
        /* Combine the old assignment value with the new partition. */
        public static void merge_assignments(int* assignment, /* assignment list for graph */
            int* subassign, /* subgraph assignment list */
            int[] subsets, /* mapping from local to global sets */
            int subnvtxs, /* number of vtxs in subgraph */
            int* loc2glob /* subgraph -> graph numbering map */
        )
        {
            int i; /* loop counter */

            for (i = 1; i <= subnvtxs; i++)
            {
                assignment[loc2glob[i]] = subsets[subassign[i]];
            }
        }
    }
}
