using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using static ChacoSharp.Utilities.MergeSort;
using static ChacoSharp.StaticConstants;

namespace ChacoSharp.Internal
{
    public static unsafe class ForceInternal
    {
        public struct bidint
        {
            public int val;
            public bidint* prev;
            public bidint* next;
        };

        /* Greedily increase the number of internal vtxs in each set. */
        public static void force_internal(vtx_data** graph, /* graph data structure */
            int nvtxs, /* number of vertices in graph */
            bool useEdgeWeights, /* are edge weights being used? */
            int* assign, /* current assignment */
            double[] goal, /* desired set sizes */
            int nsets_tot, /* total number of sets */
            int npasses_max /* number of passes to make */
        )
        {
            bidint* prev; /* back pointer for setting up lists */
            bidint* int_list = null; /* internal vwgt in each set */
            bidint* vtx_elems = null; /* linked lists of vtxs in each set */
            bidint* set_list = null; /* headers for vtx_elems lists */
            double* internal_vwgt = null; /* total internal vwgt in each set */
            int* total_vwgt = null; /* total vertex weight in each set */
            int* indices = null; /* orders sets by internal vwgt */
            bool* locked = null; /* is vertex allowed to switch sets? */
            bool isInternalVertex; /* is a vertex internal or not? */
            int* space = null; /* space for mergesort */
            int npasses; /* number of callse to improve_internal */
            int nlocked; /* number of vertices that can't move */
            int set, set2; /* sets two vertices belong to */
            bool any_change; /* did pass improve # internal vtxs? */
            int niter; /* counts calls to improve_internal */
            int vwgt_max; /* largest vertex weight in graph */
            bool progress; /* am I improving # internal vertices? */
            bool error; /* out of space? */
            int size; /* array spacing */
            int i, j; /* loop counters */

            error = true;

            /* For each set, compute the total weight of internal vertices. */

            if (DEBUG_TRACE)
            {
                Console.WriteLine("<Entering force_internal>");
            }

            indices = (int*) Marshal.AllocHGlobal(nsets_tot * sizeof(int));
            internal_vwgt = (double*) Marshal.AllocHGlobal(nsets_tot * sizeof(double));
            total_vwgt = (int*) Marshal.AllocHGlobal(nsets_tot * sizeof(int));
            if (indices == null || internal_vwgt == null || total_vwgt == null)
            {
                goto skip;
            }

            for (set = 0; set < nsets_tot; set++)
            {
                internal_vwgt[set] = 0.0d;
                total_vwgt[set] = 0;
                indices[set] = set;
            }

            vwgt_max = 0;
            for (i = 1; i <= nvtxs; i++)
            {
                isInternalVertex = true;
                set = assign[i];
                for (j = 1; j < graph[i]->nedges && isInternalVertex; j++)
                {
                    set2 = assign[graph[i]->edges[j]];
                    isInternalVertex = (set2 == set);
                }

                total_vwgt[set] += graph[i]->vwgt;
                if (isInternalVertex)
                {
                    internal_vwgt[set] += graph[i]->vwgt;
                }

                if (graph[i]->vwgt > vwgt_max)
                {
                    vwgt_max = graph[i]->vwgt;
                }
            }

            /* Now sort all the internal_vwgt values. */
            space = (int*) Marshal.AllocHGlobal(nsets_tot * sizeof(int));
            if (space == null)
            {
                goto skip;
            }

            ch_mergesort(internal_vwgt, nsets_tot, indices, space);
            Marshal.FreeHGlobal((IntPtr) space);
            space = null;

            /* Now construct a doubly linked list of sorted, internal_vwgt values. */
            int_list = (bidint*) Marshal.AllocHGlobal((nsets_tot + 1) * sizeof(bidint));
            if (int_list == null)
            {
                goto skip;
            }

            prev = &(int_list[nsets_tot]);
            prev->prev = null;
            for (i = 0; i < nsets_tot; i++)
            {
                set = indices[i];
                int_list[set].prev = prev;
                int_list[set].val = (int) internal_vwgt[set];
                prev->next = &(int_list[set]);
                prev = &(int_list[set]);
            }

            prev->next = null;
            int_list[nsets_tot].val = -1;

            Marshal.FreeHGlobal((IntPtr) internal_vwgt);
            Marshal.FreeHGlobal((IntPtr) indices);
            internal_vwgt = null;
            indices = null;

            /* Set up convenient data structure for navigating through sets. */
            set_list = (bidint*) Marshal.AllocHGlobal(nsets_tot * sizeof(bidint));
            vtx_elems = (bidint*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(bidint));
            if (set_list == null || vtx_elems == null)
            {
                goto skip;
            }

            for (i = 0; i < nsets_tot; i++)
            {
                set_list[i].next = null;
            }

            for (i = 1; i <= nvtxs; i++)
            {
                set = assign[i];
                vtx_elems[i].next = set_list[set].next;
                if (vtx_elems[i].next != null)
                {
                    vtx_elems[i].next->prev = &(vtx_elems[i]);
                }

                vtx_elems[i].prev = &(set_list[set]);
                set_list[set].next = &(vtx_elems[i]);
            }

            locked = (bool*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(bool));
            if (locked == null)
            {
                goto skip;
            }

            nlocked = 0;
            size = (int) (&(int_list[1]) - &(int_list[0]));

            any_change = true;
            npasses = 1;
            while (any_change && npasses <= npasses_max)
            {
                for (i = 1; i <= nvtxs; i++)
                {
                    locked[i] = false;
                }

                /* Now select top guy off the list and improve him. */
                any_change = false;
                progress = true;
                niter = 1;
                while (progress)
                {
                    prev = int_list[nsets_tot].next;
                    set = ((int) (prev - int_list)) / size;

                    if (DEBUG_INTERNAL)
                    {
                        Console.WriteLine("Before iteration {0:d}, nlocked = {1:d}, int[{2:d}] = {3:d}", niter, nlocked, set, prev->val);
                    }

                    if (DEBUG_INTERNAL)
                    {
                        check_internal(graph, nvtxs, int_list, set_list, vtx_elems, total_vwgt, assign, nsets_tot);
                    }

                    progress = improve_internal(graph, nvtxs, assign, goal, int_list, set_list, vtx_elems, set,
                        locked, &nlocked, useEdgeWeights, vwgt_max, total_vwgt);
                    if (progress)
                    {
                        any_change = true;
                    }

                    niter++;
                }

                npasses++;
            }

            error = false;

            skip:

            if (error)
            {
                Console.WriteLine("\nWARNING: No space to increase internal vertices.");
                Console.WriteLine("         NO INTERNAL VERTEX INCREASE PERFORMED.");
            }

            Marshal.FreeHGlobal((IntPtr) internal_vwgt);
            Marshal.FreeHGlobal((IntPtr) indices);
            Marshal.FreeHGlobal((IntPtr) locked);
            Marshal.FreeHGlobal((IntPtr) total_vwgt);
            Marshal.FreeHGlobal((IntPtr) vtx_elems);
            Marshal.FreeHGlobal((IntPtr) int_list);
            Marshal.FreeHGlobal((IntPtr) set_list);
        }

        private static void check_internal(vtx_data** graph, /* graph data structure */
            int nvtxs, /* number of vertices in graph */
            bidint* int_list, /* sorted list of internal weights per set */
            bidint* set_list, /* head of vtx_elems lists */
            bidint* vtx_elems, /* start of vertex data */
            int* total_vwgt, /* total weight in each set */
            int* assign, /* current assignment */
            int nsets_tot /* total number of sets */
        )
        {
            bidint* ptr;
            bidint* ptr2; /* elements in int_list */
            bidint* old_ptr;
            bidint* old_ptr2; /* elements in set_list */
            int vwgt_sum; /* sum of vertex weights */
            int set, set2; /* sets two vertices are in */
            int sum; /* sum of internal weights */
            int nseen; /* number of vertices found in set_lists */
            int old_val, val; /* consecutive values in int_list */
            int vtx; /* vertex in set list */
            bool isVectorInternal; /* is a vertex internal or not? */
            int size; /* array spacing */
            int j, k; /* loop counters */

            k = 0;
            size = (int) (&(int_list[1]) - &(int_list[0]));
            nseen = 0;
            old_val = -1;
            old_ptr = &(int_list[nsets_tot]);
            for (ptr = int_list[nsets_tot].next; ptr != null; ptr = ptr->next)
            {
                set = ((int) (ptr - int_list)) / size;
                val = ptr->val;
                if (val < old_val)
                {
                    Console.WriteLine("int_list out of order, k={0:d}, set = {1:d}, old_val={2:d}, val = {3:d}", k, set, old_val, val);
                }

                if (ptr->prev != old_ptr)
                {
                    Console.WriteLine(" int_list back link screwed up, set={0:d}, k={1:d}, old_ptr={2:ld}, ptr->prev = {3:ld}", set, k, (long) old_ptr, (long) ptr->prev);
                }

                old_ptr = ptr;
                old_val = val;

                vwgt_sum = 0;
                sum = 0;
                old_ptr2 = &(set_list[set]);
                for (ptr2 = set_list[set].next; ptr2 != null; ptr2 = ptr2->next)
                {
                    vtx = ((int) (ptr2 - vtx_elems)) / size;
                    vwgt_sum += graph[vtx]->vwgt;
                    if (ptr2->prev != old_ptr2)
                    {
                        Console.WriteLine(" set_list back link screwed up, set={0:d}, k={1:d}, old_ptr2={2:ld}, ptr2->prev = {3:ld}",
                            set, k, (long) old_ptr2, (long) ptr2->prev);
                    }

                    old_ptr2 = ptr2;

                    ++nseen;
                    if (assign[vtx] != set)
                    {
                        Console.WriteLine("assign[{0:d}] = {1:d}, but in set_list[{2:d}]", vtx, assign[vtx], set);
                    }

                    isVectorInternal = true;
                    for (j = 1; j < graph[vtx]->nedges && isVectorInternal; j++)
                    {
                        set2 = assign[graph[vtx]->edges[j]];
                        isVectorInternal = (set2 == set);
                    }

                    if (isVectorInternal)
                    {
                        sum += graph[vtx]->vwgt;
                    }
                }

                if (sum != val)
                {
                    Console.WriteLine("set = {0:d}, val = {1:d}, but I compute internal = {2:d}", set, val, sum);
                }

                if (vwgt_sum != total_vwgt[set])
                {
                    Console.WriteLine(" vwgt_sum = {0:d}, but total_vwgt[{1:d}] = {2:d}", vwgt_sum, set, total_vwgt[set]);
                }

                k++;
            }

            if (k != nsets_tot)
            {
                Console.WriteLine(" Only {0:d} sets in int_sets list, but nsets_tot = {1:d}", k, nsets_tot);
            }

            if (nseen != nvtxs)
            {
                Console.WriteLine(" Only {0:d} vertices found in int_sets lists, but nvtxs = {1:d}", nseen, nvtxs);
            }
        }

        public static bool improve_internal(vtx_data **graph,       /* graph data structure */
                     int               nvtxs,       /* number of vertices in graph */
                     int *             assign,      /* current assignment */
                     double[]           goal,        /* desired set sizes */
                     bidint *   int_list,    /* sorted list of internal vtx values */
                     bidint *   set_list,    /* headers of vtx_elems lists */
                     bidint *   vtx_elems,   /* lists of vtxs in each set */
                     int               set1,        /* set to try to improve */
                     bool *             locked,      /* indicates vertices not allowed to move */
                     int *             nlocked,     /* number of vertices that can't move */
                     bool               useEdgeWeights, /* are edge weights being used? */
                     int               vwgt_max,    /* largest vertex weight */
                     int *             total_vwgt   /* total vertex weight in each set */
)
{
  bidint *move_list;          /* list of vertices changing sets */
  bidint* ptr;bidint *ptr2;         /* loop through bidints */
  bidint *changed_sets;       /* list of sets that were modified */
  double         vwgt_avg   = 0.0;   /* average vertex weight in current set */
  double         degree_avg = 0.0;   /* average vertex degree in current set */
  double         frac       = .4;    /* fraction of neighbors acceptable to move. */
  double         cost, min_cost;     /* cost of making a vertex internal */
  double         min_cost_start;     /* larger than any possible cost */
  double         cost_limit;         /* acceptable cost of internalization */
  double         ratio;              /* possible wgt / desired wgt */
  float          ewgt;               /* weight of an edge */
  int            set2, set3;         /* sets of two vertices */
  int            vtx, best_vtx = -1; /* vertex to make internal */
  int            move_vtx = -1;      /* vertex to move between sets */
  int            neighbor;           /* neighbor of a vertex */
  int            nguys = 0;          /* number of vertices in current set */
  bool            isInternalVertex;           /* is a vertex internal or not? */
  bool            balanced;           /* are two sets balanced? */
  bool            flag;               /* did I improve things: return code */
  int            size;               /* array spacing */
  int            i, j;               /* loop counters */

  /* First find best candidate vertex to internalize. */
  /* This is vertex which is already most nearly internal. */
  min_cost_start = 2.0 * vwgt_max * nvtxs;
  min_cost       = min_cost_start;
  size           = (int)(&(vtx_elems[1]) - &(vtx_elems[0]));
  for (ptr = set_list[set1].next; ptr != null; ptr = ptr->next) {
    ++nguys;
    vtx = ((int)(ptr - vtx_elems)) / size;
    vwgt_avg += graph[vtx]->vwgt;
    degree_avg += (graph[vtx]->nedges - 1);
    cost = 0;
    for (i = 1; i < graph[vtx]->nedges; i++) {
      neighbor = graph[vtx]->edges[i];
      set2     = assign[neighbor];
      if (set2 != set1) {
        if (locked[neighbor]) {
          cost = min_cost_start;
        }
        else {
          cost += graph[neighbor]->vwgt;
        }
      }
    }
    if (cost == 0) { /* Lock vertex and all it's neighbors. */
      for (i = 1; i < graph[vtx]->nedges; i++) {
        neighbor = graph[vtx]->edges[i];
        if (!locked[neighbor]) {
          locked[neighbor] = true;
          ++(*nlocked);
        }
      }
    }

    if (cost < min_cost && cost != 0) {
      min_cost = cost;
      best_vtx = vtx;
    }
  }

  if (nguys > 0) {
    vwgt_avg /= nguys;
    degree_avg /= nguys;
  }
  cost_limit = frac * vwgt_avg * degree_avg;

  if (min_cost > cost_limit) {
    return false;
  }

  /* Lock the candidate vertex in current set */
  if (!locked[best_vtx]) {
    locked[best_vtx] = true;
    ++(*nlocked);
  }

  /* Also lock all his neighbors in set. */
  for (i = 1; i < graph[best_vtx]->nedges; i++) {
    neighbor = graph[best_vtx]->edges[i];
    set2     = assign[neighbor];
    if (set1 == set2 && !locked[neighbor]) {
      locked[neighbor] = true;
      ++(*nlocked);
    }
    vtx_elems[neighbor].val = set1;
  }

  ewgt      = 1;
  move_list = null;

  /* Now move neighbors of best_vtx to set1. */
  for (i = 1; i < graph[best_vtx]->nedges; i++) {
    neighbor = graph[best_vtx]->edges[i];
    set2     = assign[neighbor];
    if (set2 != set1) {
      /* Add vertex to list of guys to move to set1. */
      /* Don't move it yet in case I get stuck later. */
      /* But change his assignment so that swapping vertex has current info. */
      /* Note: This will require me to undo changes if I fail. */

      locked[neighbor] = true;
      ++(*nlocked);

      /* Remove him from his set list. */
      if (vtx_elems[neighbor].next != null) {
        vtx_elems[neighbor].next->prev = vtx_elems[neighbor].prev;
      }
      if (vtx_elems[neighbor].prev != null) {
        vtx_elems[neighbor].prev->next = vtx_elems[neighbor].next;
      }

      /* Put him in list of moved vertices */
      vtx_elems[neighbor].next = move_list;
      vtx_elems[neighbor].val  = set2;
      move_list                = &(vtx_elems[neighbor]);
      assign[neighbor]         = set1;

      total_vwgt[set2] -= graph[neighbor]->vwgt;
      total_vwgt[set1] += graph[neighbor]->vwgt;
    }
  }

  /* Now check if vertices need to be handed back to restore balance. */
  flag = true;
  for (i = 1; i < graph[best_vtx]->nedges && flag; i++) {
    neighbor = graph[best_vtx]->edges[i];
    set2     = vtx_elems[neighbor].val;
    if (set2 != set1) {
      ratio    = (total_vwgt[set1] + total_vwgt[set2]) / (goal[set1] + goal[set2]);
      balanced = (total_vwgt[set1] - goal[set1] * ratio + goal[set2] * ratio - total_vwgt[set2]) <=
                 vwgt_max;
      while (!balanced && flag) {
        /* Find a vertex to move back to set2. Use a KL metric. */
        min_cost = min_cost_start;

        for (ptr = set_list[set1].next; ptr != null; ptr = ptr->next) {
          vtx = ((int)(ptr - vtx_elems)) / size;
          if (!locked[vtx]) {
            cost = 0;
            for (j = 1; j < graph[vtx]->nedges; j++) {
              neighbor = graph[vtx]->edges[j];
              if (useEdgeWeights) {
                ewgt = graph[vtx]->ewgts[j];
              }
              set3 = assign[neighbor];
              if (set3 == set1) {
                cost += ewgt;
              }
              else if (set3 == set2) {
                cost -= ewgt;
              }
            }
            if (cost < min_cost) {
              min_cost = cost;
              move_vtx = vtx;
            }
          }
        }
        if (min_cost >= min_cost_start) {
          flag = false;
        }
        else {
          /* Add move_vtx to list of guys to move to set2. */
          /* Don't move it yet in case I get stuck later. */
          /* But change assign so later decisions have up-to-date info. */
          if (vtx_elems[move_vtx].next != null) {
            vtx_elems[move_vtx].next->prev = vtx_elems[move_vtx].prev;
          }
          if (vtx_elems[move_vtx].prev != null) {
            vtx_elems[move_vtx].prev->next = vtx_elems[move_vtx].next;
          }
          vtx_elems[move_vtx].next = move_list;
          vtx_elems[move_vtx].val  = -(set2 + 1);
          move_list                = &(vtx_elems[move_vtx]);
          assign[move_vtx]         = set2;

          total_vwgt[set2] += graph[move_vtx]->vwgt;
          total_vwgt[set1] -= graph[move_vtx]->vwgt;
        }
        balanced = total_vwgt[set1] - goal[set1] + goal[set2] - total_vwgt[set2] <= vwgt_max;
      }
    }
  }

  if (!flag) {
    /* Can't rebalance sets.  Give up, but first restore the data structures. */
    /* These include vtx_lists, total_vwgts and assign. */

    for (ptr = move_list; ptr != null;) {
      ptr2 = ptr->next;
      vtx  = ((int)(ptr - vtx_elems)) / size;
      if (ptr->val >= 0) { /* Almost moved from set2 to set1. */
        set2        = ptr->val;
        assign[vtx] = set2;
        total_vwgt[set2] += graph[vtx]->vwgt;
        total_vwgt[set1] -= graph[vtx]->vwgt;
        locked[vtx] = false;
        --(*nlocked);
      }
      else { /* Almost moved from set1 to set2. */
        set2        = -(ptr->val + 1);
        assign[vtx] = set1;
        total_vwgt[set2] -= graph[vtx]->vwgt;
        total_vwgt[set1] += graph[vtx]->vwgt;
        set2 = set1;
      }

      /* Now add vertex back into its old vtx_list (now indicated by set2) */
      ptr->next = set_list[set2].next;
      if (ptr->next != null) {
        ptr->next->prev = ptr;
      }
      ptr->prev           = &(set_list[set2]);
      set_list[set2].next = ptr;

      ptr = ptr2;
    }
    return false;
  }

  /* Now perform actual moves. */
  /* First, update assignment and place vertices into their new sets. */
  changed_sets = null;
  for (ptr = move_list; ptr != null;) {
    ptr2 = ptr->next;
    vtx  = ((int)(ptr - vtx_elems)) / size;
    if (ptr->val >= 0) {
      set2 = set1;
    }
    else {
      set2 = -(ptr->val + 1);
    }

    ptr->next = set_list[set2].next;
    if (ptr->next != null) {
      ptr->next->prev = ptr;
    }
    ptr->prev           = &(set_list[set2]);
    set_list[set2].next = ptr;

    /* Pull int_list[set2] out of its list to be used later. */
    if (ptr->val >= 0) {
      set2 = ptr->val;
    }
    if (int_list[set2].val >= 0) {
      int_list[set2].val = -(int_list[set2].val + 1);
      if (int_list[set2].next != null) {
        int_list[set2].next->prev = int_list[set2].prev;
      }
      if (int_list[set2].prev != null) {
        int_list[set2].prev->next = int_list[set2].next;
      }

      int_list[set2].next = changed_sets;
      changed_sets        = &(int_list[set2]);
    }
    ptr = ptr2;
  }
  if (int_list[set1].val >= 0) {
    if (int_list[set1].next != null) {
      int_list[set1].next->prev = int_list[set1].prev;
    }
    if (int_list[set1].prev != null) {
      int_list[set1].prev->next = int_list[set1].next;
    }

    int_list[set1].next = changed_sets;
    changed_sets        = &(int_list[set1]);
  }

  /* Finally, update internal node calculations for all modified sets. */
  while (changed_sets != null) {
    set2         = ((int)(changed_sets - int_list)) / size;
    changed_sets = changed_sets->next;

    /* Next line uses fact that list has dummy header so prev isn't null. */
    int_list[set2].next = int_list[set2].prev->next;
    int_list[set2].val  = 0;
    /* Recompute internal nodes for this set */
    for (ptr = set_list[set2].next; ptr != null; ptr = ptr->next) {
      vtx      = ((int)(ptr - vtx_elems)) / size;
      isInternalVertex = true;
      for (j = 1; j < graph[vtx]->nedges && isInternalVertex; j++) {
        set3     = assign[graph[vtx]->edges[j]];
        isInternalVertex = (set3 == set2);
      }
      if (isInternalVertex) {
        int_list[set2].val += graph[vtx]->vwgt;
      }
    }

    /* Now move internal value in doubly linked list. */
    /* Move higher in list? */
    while (int_list[set2].next != null && int_list[set2].val >= int_list[set2].next->val) {
      int_list[set2].prev = int_list[set2].next;
      int_list[set2].next = int_list[set2].next->next;
    }
    /* Move lower in list? */
    while (int_list[set2].prev != null && int_list[set2].val < int_list[set2].prev->val) {
      int_list[set2].next = int_list[set2].prev;
      int_list[set2].prev = int_list[set2].prev->prev;
    }

    if (int_list[set2].next != null) {
      int_list[set2].next->prev = &(int_list[set2]);
    }
    if (int_list[set2].prev != null) {
      int_list[set2].prev->next = &(int_list[set2]);
    }
  }
  return true;
}

    }
}
