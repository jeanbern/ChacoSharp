using System;

namespace ChacoSharp.Connect
{
    public static unsafe class HeapHelper
    {
        private static int left(int i)
        {
            return 2 * i;
        }
        private static int right(int i)
        {
            return (2 * i) + 1;
        }

        private static int parent(int i)
        {
            return i/2;
        }
/* NOTE: heap assumes indices are 1-based. */

/* Note further: All the map arguments and manipulations allow me to */
/* quickly find a particular tag value.  This assumes that tags */
/* are integers between 0 and nvals. */

public static void heapify(heap *heap,  /* array of vals/tag to make into heap */
             int          index, /* root of subtree to heapify */
             int          nvals, /* number of values in array */
             int *        map    /* maps from tag values to heap indices */
)
{
  double swap_val; /* temporary storage for swapping values */
  double swap_tag; /* temporary storage for swapping values */
  int    l, r;     /* indices of left & right children */
  int    largest;  /* index of largest value */

  l = left(index);
  r = right(index);

  if (l <= nvals && heap[l].val > heap[index].val) {
    largest = l;
  }
  else {
    largest = index;
  }

  if (r <= nvals && heap[r].val > heap[largest].val) {
    largest = r;
  }

  if (largest != index) { /* swap index with largest and recurse */
    swap_val = heap[index].val;
    swap_tag = heap[index].tag;

    heap[index].val = heap[largest].val;
    heap[index].tag = heap[largest].tag;

    heap[largest].val = swap_val;
    heap[largest].tag = (int)swap_tag;

    if (map != null) { /* update pointers to heap tags */
      map[heap[index].tag]   = index;
      map[heap[largest].tag] = largest;
    }

    heapify(heap, largest, nvals, map);
  }
}

/* Construct a heap from an unordered set of values. */
        public static void heap_build(heap *heap,  /* array of vals/tag to make into heap */
                int          nvals, /* number of values in array */
                int *        map    /* maps from tag values to heap indices */
)
{
  int i; /* loop counter */

  for (i = nvals / 2; i != 0; i--) {
    heapify(heap, i, nvals, (int *)null);
  }

  if (map != null) {
    for (i = 1; i <= nvals; i++) {
      map[heap[i].tag] = i;
    }
  }
}

public static double heap_extract_max(heap *heap,  /* array of vals/tag in a heap */
                        int          nvals, /* number of values in array */
                        int *        ptag,  /* tag associated with return value */
                        int *        map    /* maps from tag values to heap indices */
)
{
  double maxval; /* return value */

  if (nvals < 1) {
    Console.WriteLine("Heap underflow");
    throw new InvalidOperationException("Heap underflow");
  }

  if (map != null) { /* turn off map value for extracted tag */
    map[heap[1].tag] = 0;
  }

  maxval = heap[1].val;
  *ptag  = heap[1].tag;

  heap[1].val = heap[nvals].val;
  heap[1].tag = heap[nvals].tag;

  if (map != null) { /* update map value for root */
    map[heap[1].tag] = 1;
  }

  heapify(heap, 1, nvals - 1, map);

  return (maxval);
}

public static void heap_update_val(heap *heap,   /* array of vals/tag in a heap */
                     int          index,  /* index of value to update */
                     double       newval, /* new value to insert */
                     int *        map     /* maps from tag values to heap indices */
)
{
  int tag; /* tag value associated with updated val */
  int dad; /* parent of a tree node */

  tag = heap[index].tag;

  dad = parent(index);
  while (index > 1 && heap[dad].val < newval) {
    heap[index].val = heap[dad].val;
    heap[index].tag = heap[dad].tag;
    if (map != null) {
      map[heap[index].tag] = index;
    }
    index = dad;
    dad   = parent(index);
  }

  heap[index].val = newval;
  heap[index].tag = tag;
  if (map != null) {
    map[tag] = index;
  }
}
    }
}
