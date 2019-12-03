using System;
using System.Runtime.InteropServices;
using static ChacoSharp.Utilities.MakeSetLists;
using static ChacoSharp.Utilities.MakeMaps;
using static ChacoSharp.Graph.Subgraph;
using static ChacoSharp.Connect.FindComps;
using static ChacoSharp.Connect.HeapHelper;

namespace ChacoSharp.Connect
{
    public static unsafe class ConnectEnforce
    {
        /* Move vertices between domains to ensure each domain is connected. */
/* Note: This will likely result in load imbalance. */

public static void connect_enforce(vtx_data **graph,       /* data structure for graph */
                     int               nvtxs,       /* number of vertices in full graph */
                     bool               useEdgeWeights, /* are edge weights being used? */
                     int *             assignment,  /* set number of each vtx (length n) */
                     double []          goal,        /* desired sizes for each set */
                     int               nsets_tot,   /* total number sets created */
                     int *             total_move,  /* total number of vertices moved */
                     int *             max_move     /* largest connected component moved */
)
{
  vtx_data **subgraph;           /* data structure for domain graph */
  int               subnvtxs;           /* number of vertices in a domain */
  int               subnedges;          /* number of edges in a domain */
  heap *     heap;               /* data structure for sorting set sizes */
  int *             heap_map;           /* pointers from sets to heap locations */
  int *             list_ptrs;          /* header of vtx list for each domain */
  int *             setlists;           /* linked list of vtxs for each domain */
  int *             vtxlist;            /* space for breadth first search list */
  int *             comp_lists;         /* list of vtxs in each connected comp */
  int *             clist_ptrs;         /* pointers to heads of comp_lists */
  int []             subsets;            /* list of active domains (all of them) */
  int []             subsets2;           /* list of active domains (all of them) */
  double *          set_size;           /* weighted sizes of different partitions */
  double            size;               /* size of subset being moved to new domain */
  int *             bndy_list;          /* list of domains adjacent to component */
  double *          bndy_size;          /* size of these boundaries */
  double *          comp_size;          /* sizes of different connected components */
  double            comp_size_max;      /* size of largest connected component */
  int               comp_max_index = 0; /* which component is largest? */
  int *             glob2loc;           /* global to domain renumbering */
  int *             loc2glob;           /* domain to global renumbering */
  int *             degree;             /* number of neighbors of a vertex */
  int *             comp_flag;          /* component number for each vtx */
  double            ewgt;               /* edge weight */
  int               nbndy;              /* number of sets adjecent to component */
  int               domain;             /* which subdomain I'm working on */
  int               new_domain;         /* subdomain to move some vertices to */
  double            max_bndy;           /* max connectivity to other domain */
  int               ncomps;             /* number of connected components */
  int               change;             /* how many vertices change set? */
  int               max_change;         /* largest subset moved together */
  int               vtx;                /* vertex in a connected component */
  int               set;                /* set a neighboring vertex is in */
  int               i, j, k, l;         /* loop counters */

  change     = 0;
  max_change = 0;

  /* Allocate space & initialize some values. */

  set_size  = (double*)Marshal.AllocHGlobal(nsets_tot * sizeof(double));
  bndy_size = (double*)Marshal.AllocHGlobal(nsets_tot * sizeof(double));
  bndy_list = (int*)Marshal.AllocHGlobal(nsets_tot * sizeof(int));

  setlists  = (int*)Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));
  list_ptrs = (int*)Marshal.AllocHGlobal(nsets_tot * sizeof(int));

  glob2loc = (int*)Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));
  loc2glob = (int*)Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));
  subsets = new int[nsets_tot];
  heap     = (heap *)Marshal.AllocHGlobal((nsets_tot + 1) * sizeof(heap));
  heap_map = (int*)Marshal.AllocHGlobal(nsets_tot * sizeof(int));

  for (i = 0; i < nsets_tot; i++) {
    set_size[i]  = 0;
    bndy_list[i] = 0;
    bndy_size[i] = 0;
    subsets[i]   = i;
  }

  for (i = 1; i <= nvtxs; i++) {
    set_size[assignment[i]] += graph[i]->vwgt;
  }

  for (i = 0; i < nsets_tot; i++) {
    heap[i + 1].tag = i;
    heap[i + 1].val = set_size[i] - goal[i];
  }

  make_setlists(setlists, list_ptrs, nsets_tot, subsets, assignment, null, nvtxs, true);

  heap_build(heap, nsets_tot, heap_map);

  for (i = 0; i < nsets_tot; i++) {
    /* Find largest remaining set to work on next */
    size = heap_extract_max(heap, nsets_tot - i, &domain, heap_map);

    /* Construct subdomain graph. */
    subnvtxs = make_maps(setlists, list_ptrs, domain, glob2loc, loc2glob);
    if (subnvtxs > 1) {

      subgraph = (vtx_data **)Marshal.AllocHGlobal((subnvtxs + 1) * sizeof(vtx_data *));
      degree   = (int*)Marshal.AllocHGlobal((subnvtxs + 1) * sizeof(int));

      make_subgraph(graph, subgraph, subnvtxs, &subnedges, assignment, domain, glob2loc, loc2glob,
                    degree, useEdgeWeights);

      /* Find connected components. */
      comp_flag = (int*)Marshal.AllocHGlobal((subnvtxs + 1) * sizeof(int));
      vtxlist   = (int*)Marshal.AllocHGlobal(subnvtxs * sizeof(int));
      ncomps    = find_comps(subgraph, subnvtxs, comp_flag, vtxlist);
      Marshal.FreeHGlobal((IntPtr)vtxlist);

      /* Restore original graph */
      remake_graph(subgraph, subnvtxs, loc2glob, degree, useEdgeWeights);
      Marshal.FreeHGlobal((IntPtr)degree);
      Marshal.FreeHGlobal((IntPtr)subgraph);

      if (ncomps > 1) {

        /* Figure out sizes of components */
        comp_size = (double*)Marshal.AllocHGlobal(ncomps * sizeof(double));
        for (j = 0; j < ncomps; j++) {
          comp_size[j] = 0;
        }
        for (j = 1; j <= subnvtxs; j++) {
          comp_size[comp_flag[j]] += graph[loc2glob[j]]->vwgt;
        }
        comp_size_max = 0;
        for (j = 0; j < ncomps; j++) {
          if (comp_size[j] > comp_size_max) {
            comp_size_max  = comp_size[j];
            comp_max_index = j;
          }
        }
        for (j = 0; j < ncomps; j++) {
          if (j != comp_max_index) {
            change += (int)comp_size[j];
            if (comp_size[j] > max_change) {
              max_change = (int)comp_size[j];
            }
          }
        }
        Marshal.FreeHGlobal((IntPtr)comp_size);

        /* Make data structures for traversing components */
        comp_lists = (int*)Marshal.AllocHGlobal((subnvtxs + 1) * sizeof(int));
        clist_ptrs = (int*)Marshal.AllocHGlobal(ncomps * sizeof(int));
        if (ncomps > nsets_tot) {
          subsets2 = new int[ncomps];
          for (j = 0; j < ncomps; j++) {
            subsets2[j] = j;
          }
        }
        else {
          subsets2 = subsets;
        }
        make_setlists(comp_lists, clist_ptrs, ncomps, subsets2, comp_flag, null, subnvtxs, true);
        if (ncomps > nsets_tot) {
            subsets2 = null;
        }

        /* Move all but the largest component. */
        ewgt = 1;
        for (j = 0; j < ncomps; j++) {
          if (j != comp_max_index) {

            /* Figure out to which other domain it is most connected. */
            nbndy = 0;
            k     = clist_ptrs[j];
            while (k != 0) {
              vtx = loc2glob[k];
              for (l = 1; l <= graph[vtx]->nedges; l++) {
                set = assignment[graph[vtx]->edges[l]];
                if (set != domain) {
                  if (bndy_size[set] == 0) {
                    bndy_list[nbndy++] = set;
                  }
                  if (useEdgeWeights) {
                    ewgt = graph[vtx]->ewgts[l];
                  }
                  bndy_size[set] += ewgt;
                }
              }

              k = comp_lists[k];
            }

            /* Select a new domain. */
            /* Instead of just big boundary, penalize too-large sets. */
            /* Could be more aggressive to improve balance. */
            max_bndy   = 0;
            new_domain = -1;
            for (k = 0; k < nbndy; k++) {
              l = bndy_list[k];
              if (bndy_size[l] * goal[l] / (set_size[l] + 1) > max_bndy) {
                new_domain = bndy_list[k];
                max_bndy   = bndy_size[l] * goal[l] / (set_size[l] + 1);
              }
            }
            if (new_domain == -1) {
              Console.WriteLine("Error in connect_enforce: new_domain = -1.  Disconnected graph?");
              new_domain = domain;
            }

            /* Clear bndy_size array */
            for (k = 0; k < nbndy; k++) {
              bndy_size[bndy_list[k]] = 0;
            }

            k = clist_ptrs[j];

            size = 0;

            while (k != 0) {
              vtx             = loc2glob[k];
              assignment[vtx] = new_domain;
              size += graph[vtx]->vwgt;

              /* Finally, update setlists and list_ptrs */
              /* Note: current domain setlist now bad, but not used
                 again */
              setlists[vtx]         = list_ptrs[new_domain];
              list_ptrs[new_domain] = vtx;

              k = comp_lists[k];
            }
            /*
            printf("Updating set %d (from %d) to size %g\n", new_domain, domain,
            set_size[new_domain] + size - goal[new_domain]);
            */
            if (heap_map[new_domain] > 0) {
              set_size[new_domain] += size;
              heap_update_val(heap, heap_map[new_domain], set_size[new_domain] - goal[new_domain],
                              heap_map);
            }
          }
        }

        Marshal.FreeHGlobal((IntPtr)clist_ptrs);
        Marshal.FreeHGlobal((IntPtr)comp_lists);
      }
      Marshal.FreeHGlobal((IntPtr)comp_flag);
    }
  }

  Marshal.FreeHGlobal((IntPtr)heap_map);
  Marshal.FreeHGlobal((IntPtr)heap);
  Marshal.FreeHGlobal((IntPtr)loc2glob);
  Marshal.FreeHGlobal((IntPtr)glob2loc);
  Marshal.FreeHGlobal((IntPtr)list_ptrs);
  Marshal.FreeHGlobal((IntPtr)setlists);
  Marshal.FreeHGlobal((IntPtr)bndy_list);
  Marshal.FreeHGlobal((IntPtr)bndy_size);
  Marshal.FreeHGlobal((IntPtr)set_size);

  *total_move = change;
  *max_move   = max_change;
}
    }
}
