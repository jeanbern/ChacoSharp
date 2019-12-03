#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
using System;
using static ChacoSharp.StaticConstants;

namespace ChacoSharp.Utilities
{
    public static unsafe class MergeSort
    {
        /// <summary>
        /// Merge sort values in vals, returning sorted indices.
        /// </summary>
        /// <param name="vals">values to be sorted</param>
        /// <param name="nvals">number of values</param>
        /// <param name="indices">indices of sorted list</param>
        /// <param name="space">space for nvals integers</param>
        public static void ch_mergesort(double* vals, int nvals, int* indices, int* space)
        {
            for (var i = 0; i < nvals; i++)
            {
                indices[i] = i;
            }

            RecurseSort(vals, nvals, indices, space);

            if (DEBUG_BPMATCH == DebugFlagBP.NoDebugging)
            {
                return;
            }

            var errorFlag = false; /* has sorting screwed up? */
            for (var i = 1; i < nvals; i++)
            {
                if (vals[indices[i - 1]] > vals[indices[i]])
                {
                    errorFlag = true;
                }
            }

            if (!errorFlag)
            {
                return;
            }

            Console.WriteLine("List improperly sorted in mergesort");
            if (DEBUG_BPMATCH == DebugFlagBP.Logging)
            {
                return;
            }

            for (var i = 1; i < nvals; i++)
            {
                Console.WriteLine("{0:d}  {1:f}", indices[i], vals[indices[i]]);
            }
        }

        /// <summary>
        /// Recursive implementation of mergesort.
        /// </summary>
        /// <param name="vals">values to be sorted</param>
        /// <param name="nvals">number of values</param>
        /// <param name="indices">indices of sorted list</param>
        /// <param name="space">space for nvals integers</param>
        private static void RecurseSort(double* vals, int nvals, int* indices, int* space)
        {
            /* First consider base cases */
            if (nvals <= 1)
            {
                return;
            }

            if (nvals == 2)
            {
                if (vals[indices[0]] <= vals[indices[1]])
                {
                    return;
                }

                var temp = indices[0];
                indices[0] = indices[1];
                indices[1] = temp;
                return;
            }

            var sublist1Length = nvals / 2;
            var sublist2Length = nvals - sublist1Length;
            RecurseSort(vals, sublist1Length, indices, space);
            RecurseSort(vals, sublist2Length, indices + sublist1Length, space);
            Merge(vals, indices, sublist1Length, sublist2Length, space);
        }

        /// <summary>
        /// Merge two sorted lists to create longer sorted list.
        /// </summary>
        /// <param name="vals">values to be sorted</param>
        /// <param name="indices">start of first sorted list</param>
        /// <param name="length1">number of values in first list</param>
        /// <param name="length2">number of values in second list</param>
        /// <param name="space">sorted answer</param>
        private static void Merge(double* vals, int* indices, int length1, int length2, int* space)
        {
            int list1Index = 0, list2Index = 0;
            var spaceptr = space;
            var index1 = indices;
            var index2 = indices + length1;

            while (list1Index < length1 && list2Index < length2)
            {
                if (vals[*index1] <= vals[*index2])
                {
                    *spaceptr++ = *index1++;
                    ++list1Index;
                }
                else
                {
                    *spaceptr++ = *index2++;
                    ++list2Index;
                }
            }

            /* Now add on remaining elements */
            while (list1Index < length1)
            {
                *spaceptr++ = *index1++;
                ++list1Index;
            }

            while (list2Index < length2)
            {
                *spaceptr++ = *index2++;
                ++list2Index;
            }

            spaceptr = space;
            index1 = indices;
            for (list1Index = length1 + length2; list1Index != 0; list1Index--)
            {
                *index1++ = *spaceptr++;
            }
        }
    }
}
