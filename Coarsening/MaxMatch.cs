#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
using System;
using System.Runtime.InteropServices;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Utilities.Randomize;

namespace ChacoSharp.Coarsening
{
    public static unsafe class MaxMatch
    {
        /// <summary>
        /// Find a maximal matching in a graph using one of several algorithms.
        /// </summary>
        public static int maxmatch(vtx_data** graph, /* array of vtx data for graph */
            int nvtxs, /* number of vertices in graph */
            int nedges, /* number of edges in graph */
            int* mflag, /* flag indicating vtx selected or not */
            bool useEdgeWeights, /* are edge weights being used? */
            int igeom, /* geometric dimensionality */
            float** coords /* coordinates for each vertex */
        )
        {
            int nmerged = 0; /* number of matching edges found */

            if (MATCH_TYPE == MatchingRoutine.maxmatch1 || (MATCH_TYPE == MatchingRoutine.maxmatch5_geometric && coords == null))
            {
                /* Dumb, fast routine. */
                nmerged = maxmatch1(graph, nvtxs, mflag, useEdgeWeights);
            }

            else if (MATCH_TYPE == MatchingRoutine.maxmatch2)
            {
                /* More random but somewhat slower. */
                nmerged = maxmatch2(graph, nvtxs, mflag, useEdgeWeights);
            }

            else if (MATCH_TYPE == MatchingRoutine.maxmatch3)
            {
                /* Much more random but slower still. */
                nmerged = maxmatch3(graph, nvtxs, mflag, useEdgeWeights);
            }

            else if (MATCH_TYPE == MatchingRoutine.maxmatch4_Luby)
            {
                /* Truly random but very slow. */
                nmerged = maxmatch4(graph, nvtxs, nedges, mflag, useEdgeWeights);
            }
            else if (MATCH_TYPE == MatchingRoutine.maxmatch5_geometric && coords != null)
            {
                /* Geometric nearness. */
                nmerged = maxmatch5(graph, nvtxs, mflag, igeom, coords);
            }

            else if (MATCH_TYPE == MatchingRoutine.maxmatch9_minimumVertexDegree)
            {
                /* Minimum degree of merged vertex */
                nmerged = maxmatch9(graph, nvtxs, mflag, useEdgeWeights);
            }

            if (DEBUG_COARSEN)
            {
                Console.WriteLine("Number of matching edges = {0:D}", nmerged);
            }

            return (nmerged);
        }

        /// <summary>
        /// Find a maximal matching in a graph using simple greedy algorithm.
        /// Randomly permute vertices, and then have each select an unmatched neighbor.
        /// </summary>
        public static int maxmatch1(vtx_data** graph, /* array of vtx data for graph */
            int nvtxs, /* number of vertices in graph */
            int* mflag, /* flag indicating vtx selected or not */
            bool useEdgeWeights /* are edge weights being used? */
        )
        {
            float ewgt_max; /* largest edge weight seen so far */
            int* jptr; /* loops through integer arrays */
            int vtx; /* vertex to process next */
            int neighbor; /* neighbor of a vertex */
            int nmerged; /* number of edges in matching */
            bool matched; /* is a vertex matched yet? */
            int jsave; /* best matching edge found so far */
            int i, j; /* loop counters */

            /* Initialize mflag array. */
            jptr = mflag;
            for (i = nvtxs; i != 0; i--)
            {
                *(++jptr) = 0;
            }

            nmerged = 0;

            /* Select random starting point in list of vertices. */
            vtx = 1 + (int) (drandom() * nvtxs);

            if (!useEdgeWeights || !HEAVY_MATCH)
            {
                /* Choose first neighbor */
                for (i = nvtxs; i != 0; i--)
                {
                    if (mflag[vtx] == 0)
                    {
                        /* Not already matched. */
                        /* Select first free edge. */
                        matched = false;
                        for (j = 1; !matched && j < graph[vtx]->nedges; j++)
                        {
                            neighbor = graph[vtx]->edges[j];
                            if (mflag[neighbor] == 0)
                            {
                                mflag[vtx] = neighbor;
                                mflag[neighbor] = vtx;
                                matched = true;
                                nmerged++;
                            }
                        }
                    }

                    if (++vtx > nvtxs)
                    {
                        vtx = 1;
                    }
                }
            }

            else
            {
                /* Choose heavy edge neighbor */
                for (i = nvtxs; i != 0; i--)
                {
                    if (mflag[vtx] == 0)
                    {
                        /* Not already matched. */
                        /* Select heaviest free edge. */
                        jsave = 0;
                        ewgt_max = 0;
                        for (j = 1; j < graph[vtx]->nedges; j++)
                        {
                            neighbor = graph[vtx]->edges[j];
                            if (mflag[neighbor] == 0 && graph[vtx]->ewgts[j] > ewgt_max)
                            {
                                ewgt_max = graph[vtx]->ewgts[j];
                                jsave = j;
                            }
                        }

                        if (jsave > 0)
                        {
                            neighbor = graph[vtx]->edges[jsave];
                            mflag[vtx] = neighbor;
                            mflag[neighbor] = vtx;
                            nmerged++;
                        }
                    }

                    if (++vtx > nvtxs)
                    {
                        vtx = 1;
                    }
                }
            }

            return (nmerged);
        }

        /// <summary>
        /// Find a maximal matching in a graph using simple greedy algorithm.
        /// Randomly permute vertices, and then have each select first unmatched neighbor.
        /// </summary>
        public static int maxmatch2(vtx_data** graph, /* array of vtx data for graph */
            int nvtxs, /* number of vertices in graph */
            int* mflag, /* flag indicating vtx selected or not */
            bool useEdgeWeights /* are edge weights being used? */
        )
        {
            float ewgt_max; /* heaviest edge seen so far */
            int* order; /* random ordering of vertices */
            int* iptr; /* loops through integer arrays */
            int* jptr; /* loops through integer arrays */
            bool matched; /* has vertex been matched? */
            int vtx; /* vertex to process next */
            int neighbor; /* neighbor of a vertex */
            int nmerged; /* number of edges in matching */
            int jsave; /* best edge so far */
            int i, j; /* loop counters */

            /* First, randomly permute the vertices. */
            iptr = order = (int*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));
            jptr = mflag;
            for (i = 1; i <= nvtxs; i++)
            {
                *(++iptr) = i;
                *(++jptr) = 0;
            }

            randomize(order, nvtxs);

            nmerged = 0;

            if (!useEdgeWeights || !HEAVY_MATCH)
            {
                /* Simple greedy approach. */
                for (i = 1; i <= nvtxs; i++)
                {
                    vtx = order[i];
                    if (mflag[vtx] == 0)
                    {
                        /* Not already matched. */
                        /* Find first unmatched neighbor of vtx. */
                        matched = false;
                        for (j = 1; j < graph[vtx]->nedges && !matched; j++)
                        {
                            neighbor = graph[vtx]->edges[j];
                            if (mflag[neighbor] == 0)
                            {
                                /* Found one! */
                                mflag[vtx] = neighbor;
                                mflag[neighbor] = vtx;
                                matched = true;
                                nmerged++;
                            }
                        }
                    }
                }
            }

            else
            {
                /* Find heavy edge to match */
                for (i = 1; i <= nvtxs; i++)
                {
                    vtx = order[i];
                    if (mflag[vtx] == 0)
                    {
                        /* Not already matched. */
                        /* Find heaviest unmatched neighbor of vtx. */
                        jsave = 0;
                        ewgt_max = 0;
                        for (j = 1; j < graph[vtx]->nedges; j++)
                        {
                            neighbor = graph[vtx]->edges[j];
                            if (mflag[neighbor] == 0 && graph[vtx]->ewgts[j] > ewgt_max)
                            {
                                ewgt_max = graph[vtx]->ewgts[j];
                                jsave = j;
                            }
                        }

                        if (jsave > 0)
                        {
                            neighbor = graph[vtx]->edges[jsave];
                            mflag[vtx] = neighbor;
                            mflag[neighbor] = vtx;
                            nmerged++;
                        }
                    }
                }
            }

            Marshal.FreeHGlobal((IntPtr) order);
            return nmerged;
        }

        /// <summary>
        /// Find a maximal matching in a graph using simple greedy algorithm.
        /// Randomly permute vertices, and then have each select an unmatched neighbor.
        /// </summary>
        public static int maxmatch3(vtx_data** graph, /* array of vtx data for graph */
            int nvtxs, /* number of vertices in graph */
            int* mflag, /* flag indicating vtx selected or not */
            bool useEdgeWeights /* are edge weights being used? */
        )
        {
            int* order; /* random ordering of vertices */
            int* iptr; /* loops through integer arrays */
            int* jptr; /* loops through integer arrays */
            double prob_sum; /* sum of probabilities to select from */
            double val; /* random value for selecting neighbor */
            float ewgt; /* edge weight */
            int save; /* neighbor vertex if only one active */
            int vtx; /* vertex to process next */
            int neighbor; /* neighbor of a vertex */
            int nmerged; /* number of edges in matching */
            int i, j; /* loop counters */

            /* First, randomly permute the vertices. */
            iptr = order = (int*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));
            jptr = mflag;
            for (i = 1; i <= nvtxs; i++)
            {
                *(++iptr) = i;
                *(++jptr) = 0;
            }

            randomize(order, nvtxs);

            nmerged = 0;
            if (!useEdgeWeights || !HEAVY_MATCH)
            {
                /* All edges equal. */
                for (i = 1; i <= nvtxs; i++)
                {
                    vtx = order[i];
                    if (mflag[vtx] == 0)
                    {
                        /* Not already matched. */
                        /* Add up sum of edge weights of neighbors. */
                        prob_sum = 0;
                        save = 0;
                        for (j = 1; j < graph[vtx]->nedges; j++)
                        {
                            neighbor = graph[vtx]->edges[j];
                            if (mflag[neighbor] == 0)
                            {
                                /* Set flag for single possible neighbor. */
                                if (prob_sum == 0)
                                {
                                    save = neighbor;
                                }
                                else
                                {
                                    save = 0;
                                }

                                prob_sum += 1.0;
                            }
                        }

                        if (prob_sum != 0)
                        {
                            /* Does vertex have contractible edges? */
                            nmerged++;
                            if (save != 0)
                            {
                                /* Only one neighbor, special case. */
                                mflag[vtx] = save;
                                mflag[save] = vtx;
                            }
                            else
                            {
                                /* Pick randomly neighbor. */
                                val = drandom() * prob_sum * .999999;
                                prob_sum = 0;
                                for (j = 1; mflag[vtx] == 0; j++)
                                {
                                    neighbor = graph[vtx]->edges[j];
                                    if (mflag[neighbor] == 0)
                                    {
                                        prob_sum += 1.0;
                                        if (prob_sum >= val)
                                        {
                                            mflag[vtx] = neighbor;
                                            mflag[neighbor] = vtx;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                /* Choose heavy edges preferentially. */
                for (i = 1; i <= nvtxs; i++)
                {
                    vtx = order[i];
                    if (mflag[vtx] == 0)
                    {
                        /* Not already matched. */
                        /* Add up sum of edge weights of neighbors. */
                        prob_sum = 0;
                        save = 0;
                        for (j = 1; j < graph[vtx]->nedges; j++)
                        {
                            neighbor = graph[vtx]->edges[j];
                            if (mflag[neighbor] == 0)
                            {
                                /* Set flag for single possible neighbor. */
                                if (prob_sum == 0)
                                {
                                    save = neighbor;
                                }
                                else
                                {
                                    save = 0;
                                }

                                ewgt = graph[vtx]->ewgts[j];
                                prob_sum += ewgt;
                            }
                        }

                        if (prob_sum != 0.0)
                        {
                            /* Does vertex have contractible edges? */
                            nmerged++;
                            if (save != 0)
                            {
                                /* Only one neighbor, special case. */
                                mflag[vtx] = save;
                                mflag[save] = vtx;
                            }
                            else
                            {
                                /* Pick randomly neighbor, skewed by edge weights. */
                                val = drandom() * prob_sum * .999999;
                                prob_sum = 0;
                                for (j = 1; mflag[vtx] == 0; j++)
                                {
                                    neighbor = graph[vtx]->edges[j];
                                    if (mflag[neighbor] == 0)
                                    {
                                        ewgt = graph[vtx]->ewgts[j];
                                        prob_sum += ewgt;
                                        if (prob_sum >= val)
                                        {
                                            mflag[vtx] = neighbor;
                                            mflag[neighbor] = vtx;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            Marshal.FreeHGlobal((IntPtr) order);
            return (nmerged);
        }

        /// <summary>
        /// Find a maximal matching in a graph using Luby's algorithm.  Assign a random value to each edge and include an edge in the matching if it has a higher value than all neighboring edges that aren't disallowed.
        /// (Use float instead of double values to save space.)
        /// </summary>
        /// <remarks>
        /// THIS ROUTINE IS SLOWER THAN MAXMATCH3, AND SO PROBABLY OBSOLETE.
        /// </remarks>
        public static int maxmatch4(vtx_data** graph, /* array of vtx data for graph */
            int nvtxs, /* number of vertices in graph */
            int nedges, /* number of edges in graph */
            int* mflag, /* flag indicating vtx selected or not */
            bool useEdgeWeights)
        {
            int* iptr; /* loops through integer arrays */
            float* edgevals; /* random values for all edges */
            float* evptr; /* loops through edgevals */
            double maxval; /* largest edge value for a vertex */
            int neighbor; /* neighbor of a vertex */
            int nmerged; /* number of edges in matching */
            bool change; /* any new edges in matching? */
            int* start; /* start of edgevals list for each vertex */
            int i, j, k; /* loop counters */

            /* Allocate and initialize space. */
            evptr = edgevals = (float*) Marshal.AllocHGlobal(2 * nedges * sizeof(float));

            start = (int*) Marshal.AllocHGlobal((nvtxs + 2) * sizeof(int));
            start[1] = 0;
            for (i = 1; i <= nvtxs; i++)
            {
                start[i + 1] = start[i] + graph[i]->nedges - 1;
            }

            /* Assign a random value to each edge. */
            if (!useEdgeWeights || !HEAVY_MATCH)
            {
                /* All edges are equal. */
                for (i = 1; i <= nvtxs; i++)
                {
                    for (j = 1; j < graph[i]->nedges; j++)
                    {
                        neighbor = graph[i]->edges[j];
                        if (neighbor > i)
                        {
                            *evptr = (float) drandom();
                        }
                        else
                        {
                            /* Look up already-generated value. */
                            for (k = 1; graph[neighbor]->edges[k] != i; k++)
                            { }

                            *evptr = edgevals[start[neighbor] + k - 1];
                        }

                        evptr++;
                    }
                }
            }
            else
            {
                /* Prefer heavy weight edges. */
                for (i = 1; i <= nvtxs; i++)
                {
                    for (j = 1; j < graph[i]->nedges; j++)
                    {
                        neighbor = graph[i]->edges[j];
                        if (neighbor > i)
                        {
                            *evptr = (float) (graph[i]->ewgts[j] * drandom());
                        }
                        else
                        {
                            /* Look up already-generated value. */
                            for (k = 1; graph[neighbor]->edges[k] != i; k++)
                            { }

                            *evptr = edgevals[start[neighbor] + k - 1];
                        }

                        evptr++;
                    }
                }
            }

            for (iptr = mflag, i = nvtxs; i != 0; i--)
            {
                *(++iptr) = -(nvtxs + 1);
            }

            nmerged = 0;
            change = true;
            while (change)
            {
                change = false;

                for (i = 1; i <= nvtxs; i++)
                {
                    /* Find largest valued edge of each vtx */
                    if (mflag[i] < 0)
                    {
                        maxval = 0.0;
                        k = -1;
                        evptr = &(edgevals[start[i]]);
                        for (j = 1; j < graph[i]->nedges; j++)
                        {
                            if (*evptr > maxval && mflag[graph[i]->edges[j]] < 0)
                            {
                                maxval = *evptr;
                                k = j;
                            }

                            evptr++;
                        }

                        if (k == -1)
                        {
                            mflag[i] = 0; /* No neighbors are alive. */
                        }
                        else
                        {
                            mflag[i] = -graph[i]->edges[k];
                        }
                    }
                }

                /* If vtxs agree on largest valued edge, add to independent set. */
                for (i = 1; i <= nvtxs; i++)
                {
                    if (-mflag[i] <= i || -mflag[-mflag[i]] != i)
                    {
                        continue;
                    }

                    /* Add edge to independent set. */
                    nmerged++;
                    mflag[i] = -mflag[i];
                    mflag[mflag[i]] = i;
                    change = true;
                }
            }

            /* Maximal independent set is indicated by corresponding pairs of positive values in the mflag array. */
            for (i = 1; i <= nvtxs; i++)
            {
                if (mflag[i] < 0)
                {
                    mflag[i] = 0;
                }
            }

            Marshal.FreeHGlobal((IntPtr) start);
            Marshal.FreeHGlobal((IntPtr) edgevals);
            return nmerged;
        }

        /// <summary>
        /// Find a maximal matching in a graph by geometrically near neighbors.
        /// </summary>
        public static int maxmatch5(vtx_data** graph, /* array of vtx data for graph */
            int nvtxs, /* number of vertices in graph */
            int* mflag, /* flag indicating vtx selected or not */
            int igeom, /* geometric dimensionality */
            float** coords /* coordinates of each vertex */
        )
        {
            double dist; /* distance to free neighbor */
            double min_dist; /* smallest distance to free neighbor */
            int* jptr; /* loops through integer arrays */
            int vtx; /* vertex to process next */
            int neighbor; /* neighbor of a vertex */
            int nmerged; /* number of edges in matching */
            int jsave; /* best edge so far */
            int i, j; /* loop counters */

            /* Initialize mflag array. */
            jptr = mflag;
            for (i = 1; i <= nvtxs; i++)
            {
                *(++jptr) = 0;
            }

            nmerged = 0;

            /* Select random starting point in list of vertices. */
            vtx = 1 + (int) (drandom() * nvtxs);

            if (igeom == 1)
            {
                for (i = nvtxs; i != 0; i--)
                {
                    /* Choose geometrically nearest neighbor */
                    if (mflag[vtx] == 0)
                    {
                        /* Not already matched. */
                        /* Select nearest free edge. */
                        jsave = 0;
                        min_dist = double.MaxValue;
                        for (j = 1; j < graph[vtx]->nedges; j++)
                        {
                            neighbor = graph[vtx]->edges[j];
                            if (mflag[neighbor] == 0)
                            {
                                dist = (coords[0][vtx] - coords[0][neighbor]) * (coords[0][vtx] - coords[0][neighbor]);
                                if (dist < min_dist)
                                {
                                    jsave = j;
                                    min_dist = dist;
                                }
                            }
                        }

                        if (jsave > 0)
                        {
                            neighbor = graph[vtx]->edges[jsave];
                            mflag[vtx] = neighbor;
                            mflag[neighbor] = vtx;
                            nmerged++;
                        }
                    }

                    if (++vtx > nvtxs)
                    {
                        vtx = 1;
                    }
                }
            }

            else if (igeom == 2)
            {
                for (i = nvtxs; i != 0; i--)
                {
                    /* Choose geometrically nearest neighbor */
                    if (mflag[vtx] == 0)
                    {
                        /* Not already matched. */
                        /* Select nearest free edge. */
                        jsave = 0;
                        min_dist = double.MaxValue;
                        for (j = 1; j < graph[vtx]->nedges; j++)
                        {
                            neighbor = graph[vtx]->edges[j];
                            if (mflag[neighbor] == 0)
                            {
                                dist = (coords[0][vtx] - coords[0][neighbor]) * (coords[0][vtx] - coords[0][neighbor]);
                                if (dist < min_dist)
                                {
                                    dist +=
                                        (coords[1][vtx] - coords[1][neighbor]) * (coords[1][vtx] - coords[1][neighbor]);
                                    if (dist < min_dist)
                                    {
                                        jsave = j;
                                        min_dist = dist;
                                    }
                                }
                            }
                        }

                        if (jsave > 0)
                        {
                            neighbor = graph[vtx]->edges[jsave];
                            mflag[vtx] = neighbor;
                            mflag[neighbor] = vtx;
                            nmerged++;
                        }
                    }

                    if (++vtx > nvtxs)
                    {
                        vtx = 1;
                    }
                }
            }

            else if (igeom >= 2)
            {
                for (i = nvtxs; i != 0; i--)
                {
                    /* Choose geometrically nearest neighbor */
                    if (mflag[vtx] == 0)
                    {
                        /* Not already matched. */
                        /* Select nearest free edge. */
                        jsave = 0;
                        min_dist = double.MaxValue;
                        for (j = 1; j < graph[vtx]->nedges; j++)
                        {
                            neighbor = graph[vtx]->edges[j];
                            if (mflag[neighbor] == 0)
                            {
                                dist = (coords[0][vtx] - coords[0][neighbor]) * (coords[0][vtx] - coords[0][neighbor]);
                                if (dist < min_dist)
                                {
                                    dist +=
                                        (coords[1][vtx] - coords[1][neighbor]) * (coords[1][vtx] - coords[1][neighbor]);
                                    if (dist < min_dist)
                                    {
                                        dist +=
                                            (coords[2][vtx] - coords[2][neighbor]) * (coords[2][vtx] - coords[2][neighbor]);
                                        if (dist < min_dist)
                                        {
                                            jsave = j;
                                            min_dist = dist;
                                        }
                                    }
                                }
                            }
                        }

                        if (jsave > 0)
                        {
                            neighbor = graph[vtx]->edges[jsave];
                            mflag[vtx] = neighbor;
                            mflag[neighbor] = vtx;
                            nmerged++;
                        }
                    }

                    if (++vtx > nvtxs)
                    {
                        vtx = 1;
                    }
                }
            }

            return (nmerged);
        }

        /// <summary>
        /// Find a maximal matching in a graph using simple greedy algorithm.
        /// Randomly permute vertices, and then have each select an unmatched neighbor.
        /// Choose the neighbor which results in coarsened vertex of minimum degree.
        /// </summary>
        public static int maxmatch9(vtx_data** graph, /* array of vtx data for graph */
            int nvtxs, /* number of vertices in graph */
            int* mflag, /* flag indicating vtx selected or not */
            bool useEdgeWeights /* are edge weights being used? */
        )
        {
            int* order; /* random ordering of vertices */
            int* neighbors; /* scatter array for neighbor list */
            int* iptr; /* loops through integer arrays */
            int* jptr; /* loops through integer arrays */
            float ewgt; /* edge weight */
            int save; /* neighbor vertex if only one active */
            int vtx; /* vertex to process next */
            int neighbor; /* neighbor of a vertex */
            int best; /* best match found so far */
            float same = 0, best_same; /* maximum # neighbors in common so far */
            float best_ewgt; /* edge weight of possible matching edge */
            int nmerged; /* number of edges in matching */
            int i, j, k; /* loop counters */

            /* First, randomly permute the vertices. */
            neighbors = (int*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));
            iptr = order = (int*) Marshal.AllocHGlobal((nvtxs + 1) * sizeof(int));
            jptr = mflag;
            for (i = 1; i <= nvtxs; i++)
            {
                *(++iptr) = i;
                *(++jptr) = 0;
                neighbors[i] = 0;
            }

            randomize(order, nvtxs);

            /*    if (!useEdgeWeights || !HEAVY_MATCH) { */

            nmerged = 0;
            ewgt = 0;
            for (i = 1; i <= nvtxs; i++)
            {
                vtx = order[i];
                if (mflag[vtx] == 0)
                {
                    /* Not already matched. */
                    /* Add up sum of edge weights of neighbors. */
                    save = -1;
                    for (j = 1; j < graph[vtx]->nedges; j++)
                    {
                        neighbor = graph[vtx]->edges[j];
                        neighbors[neighbor] = i;
                        if (mflag[neighbor] == 0)
                        {
                            /* Set flag for single possible neighbor. */
                            if (save == -1)
                            {
                                save = neighbor;
                            }
                            else
                            {
                                save = 0;
                            }
                        }
                        else
                        {
                            neighbors[mflag[neighbor]] = i;
                        }
                    }

                    if (save != -1)
                    {
                        /* Does vertex have contractible edges? */
                        nmerged++;
                        if (save > 0)
                        {
                            /* Only one neighbor, easy special case. */
                            mflag[vtx] = save;
                            mflag[save] = vtx;
                        }
                        else
                        {
                            /* Merge with best neighbor */
                            best = 0;
                            best_same = -1;
                            best_ewgt = -1;
                            for (j = 1; j < graph[vtx]->nedges; j++)
                            {
                                neighbor = graph[vtx]->edges[j];
                                if (mflag[neighbor] == 0)
                                {
                                    if (useEdgeWeights && HEAVY_MATCH)
                                    {
                                        ewgt = graph[vtx]->ewgts[j];
                                    }

                                    if (ewgt > best_ewgt)
                                    {
                                        best = neighbor;
                                        best_same = same;
                                        best_ewgt = ewgt;
                                    }
                                    else if (ewgt == best_ewgt)
                                    {
                                        /* break ties by larger same value */
                                        if (best_same == -1)
                                        {
                                            /* Compute same value for current best vtx. */
                                            best_same = 0;
                                            for (k = 1; k < graph[best]->nedges; k++)
                                            {
                                                if (neighbors[graph[best]->edges[k]] == i)
                                                {
                                                    if (useEdgeWeights)
                                                    {
                                                        best_same += graph[best]->ewgts[k];
                                                    }
                                                    else
                                                    {
                                                        best_same += 1;
                                                    }
                                                }
                                            }
                                        }

                                        /* Now compute same value for candidate vtx. */
                                        same = 0;
                                        for (k = 1; k < graph[neighbor]->nedges; k++)
                                        {
                                            if (neighbors[graph[neighbor]->edges[k]] == i)
                                            {
                                                if (useEdgeWeights)
                                                {
                                                    same += graph[neighbor]->ewgts[k];
                                                }
                                                else
                                                {
                                                    same += 1;
                                                }
                                            }
                                        }

                                        if (same > best_same)
                                        {
                                            best = neighbor;
                                            best_same = same;
                                            best_ewgt = ewgt;
                                        }
                                    }
                                }
                            }

                            mflag[vtx] = best;
                            mflag[best] = vtx;
                        }
                    }
                }
            }

            Marshal.FreeHGlobal((IntPtr) order);
            Marshal.FreeHGlobal((IntPtr) neighbors);
            return (nmerged);
        }
    }
}
