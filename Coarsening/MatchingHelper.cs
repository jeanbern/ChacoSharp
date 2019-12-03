#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
#pragma warning disable HAA0101 // Array allocation for params parameter
using System;
using static ChacoSharp.StaticConstants;
// ReSharper disable InvertIf

namespace ChacoSharp.Coarsening
{
    public static unsafe class MatchingHelper
    {
        /* Despite the formal disclaimer above, this routine is modified from
   code provided by Ed Rothberg at SGI. */

/*
The following code takes a bipartite graph as input (with
n_left+n_right nodes, where node 'i' is adjacent to nodes
'indices[pointers[i]]' through 'indices[pointers[i+1]-1]', all
0-based) and returns a minimum size edge cover (in sep_nodes[]).
*/

        public static void bpcover(int leftVertexCount, /* number of vertices on left side */
            int rightVertexCount, /* number of vertices on right side */
            int* pointers, /* start/stop of adjacency lists */
            int* indices, /* adjacency list for each vertex */
            int* sep_size, /* returned size of separator */
            int* sep_nodes /* list of separator nodes */
        )
        {
            int i; /* loop counter */

            if (DEBUG_COVER != 0)
            {
                Console.WriteLine("-> Entering bpcover, nleft = {0:d}, nright = {1:d}, 2*nedges = {2:d}", leftVertexCount, rightVertexCount, pointers[leftVertexCount + rightVertexCount] - pointers[0]);
            }

            var matching = new int[leftVertexCount + rightVertexCount];
            var touched = new bool[leftVertexCount + rightVertexCount];

            bpmatching(leftVertexCount, rightVertexCount, pointers, indices, matching, touched);

            reachability(leftVertexCount, rightVertexCount, pointers, indices, matching, touched);

            /* Separator includes untouched nodes on left, touched on right. */
            /* Left separator nodes if unconnected to unmatched left node via */
            /* augmenting path, right separator nodes otherwise. */

            *sep_size = 0;
            for (i = 0; i < leftVertexCount; i++)
            {
                if (!touched[i])
                {
                    sep_nodes[(*sep_size)++] = i;
                }
            }

            for (i = leftVertexCount; i < leftVertexCount + rightVertexCount; i++)
            {
                if (touched[i])
                {
                    sep_nodes[(*sep_size)++] = i;
                }
            }

            sep_nodes[*sep_size] = 0;

            if (DEBUG_COVER != 0)
            {
                confirm_match(leftVertexCount, rightVertexCount, pointers, indices, matching, *sep_size, sep_nodes);
            }
        }

        private static void bpmatching(int leftVertexCount, /* number of vertices on left side */
            int rightVertexCount, /* number of vertices on right side */
            int* pointers, /* start/stop of adjacency lists */
            int* indices, /* adjacency list for each vertex */
            int[] matching, /* array to encode matching */
            bool[] touched /* flags for each vertex */
        )
        {
            int i, j; /* loop counters */

            /* First mark all the vertices as unmatched & untouched. */
            for (i = 0; i < leftVertexCount + rightVertexCount; i++)
            {
                matching[i] = -1;
                touched[i] = false;
            }

            /* Now generate a fast, greedy matching to start. */
            for (i = leftVertexCount; i < leftVertexCount + rightVertexCount; i++)
            {
                for (j = pointers[i]; j < pointers[i + 1]; j++)
                {
                    if (matching[indices[j]] == -1)
                    {
                        /* Node not already matched. */
                        matching[i] = indices[j];
                        matching[indices[j]] = i;
                        break;
                    }
                }
            }

            /* Now try to enlarge it via augmenting paths. */

            var seen = new int[leftVertexCount + rightVertexCount];

            /* Look for an augmenting path. */
            for (i = 0; i < leftVertexCount; i++)
            {
                if (matching[i] == -1)
                {
                    augment(i, pointers, indices, matching, touched, seen);
                }
            }
        }

        private static void augment(int node, /* start node in augmenting path */
            int* pointers, /* start/stop of adjacency lists */
            int* indices, /* adjacency list for each vertex */
            int[] matching, /* array to encode matching */
            bool[] touched, /* flags for each vertex */
            int[] seen /* keeps list of vertices encountered */
        )
        {
            /* Look for augmenting path in graph. */
            var verticesSeen = 0;
            var wasEnlarged = touch(node, pointers, indices, matching, touched, seen, &verticesSeen);

            if (wasEnlarged)
            {
                /* Found an augmenting path! */
                /* Free all the vertices encountered in search. */
                /* Otherwise, they can't be involved in augmentation, */
                /* so leave them touched. */
                for (var i = 0; i < verticesSeen; i++)
                {
                    touched[seen[i]] = false;
                }
            }
        }

/* Mark everybody in my alternating path tree, and recursively update */
/* matching if augmenting path found. */
        private static bool touch(int node, int* pointers, /* start/stop of adjacency lists */
            int* indices, /* adjacency list for each vertex */
            int[] matching, /* array to encode matching */
            bool[] touched, /* flags for each vertex */
            int[] seen, /* list of vertices encountered */
            int* nseen /* number of vertices encountered */
        )
        {
            touched[node] = true;
            seen[(*nseen)++] = node;

            for (var j = pointers[node]; j < pointers[node + 1]; j++)
            {
                var neighbor = indices[j]; /* neighbor of a vertex */
                if (!touched[neighbor])
                {
                    touched[neighbor] = true;
                    seen[(*nseen)++] = neighbor;
                    if (matching[neighbor] == -1)
                    {
                        /* Found augmenting path! */
                        matching[neighbor] = node;
                        matching[node] = neighbor;
                        return true;
                    }

                    var result = touch(matching[neighbor], pointers, indices, matching, touched, seen, nseen); /* return node number (or -1) */
                    if (result)
                    {
                        matching[neighbor] = node;
                        matching[node] = neighbor;
                        return true;
                    }
                }
            }

            return false;
        }

        private static void reachability(int n_left, /* number of vertices on left side */
            int n_right, /* number of vertices on right side */
            int* pointers, /* start/stop of adjacency lists */
            int* indices, /* adjacency list for each vertex */
            int[] matching, /* array to encode matching */
            bool[] touched /* flags for each vertex */
        )
        {
            /* Initialize all the vertices to be untouched */
            for (var i = 0; i < n_left + n_right; i++)
            {
                touched[i] = false;
            }

            for (var i = 0; i < n_left; i++)
            {
                if (!touched[i] && matching[i] == -1)
                {
                    touch2(i, pointers, indices, matching, touched);
                }
            }
        }

/* Mark everybody in my alternating path tree, and return vertex at */
/* end of augmenting path if found. */
        private static bool touch2(int node, int* pointers, /* start/stop of adjacency lists */
            int* indices, /* adjacency list for each vertex */
            int[] matching, /* array to encode matching */
            bool[] touched /* flags for each vertex */
        )
        {
            touched[node] = true;
            for (var j = pointers[node]; j < pointers[node + 1]; j++)
            {
                var neighbor = indices[j]; /* neighbor of a vertex */
                if (!touched[neighbor])
                {
                    touched[neighbor] = true;
                    if (matching[neighbor] == -1)
                    {
                        return true;
                    }

                    var result = touch2(matching[neighbor], pointers, indices, matching, touched); /* return node number (or -1) */
                    if (result)
                    {
                        return true;
                    }
                }
            }

            return false;
        }

        private static void confirm_match(int n_left, /* number of vertices on left side */
            int n_right, /* number of vertices on right side */
            int* pointers, /* start/stop of adjacency lists */
            int* indices, /* adjacency list for each vertex */
            int[] matching, /* array to encode matching */
            int sep_size, /* returned size of separator */
            int* sep_nodes /* list of separator nodes */
        )
        {
            var marked = new bool[n_left + n_right];

            for (var i = 0; i < n_left + n_right; i++)
            {
                marked[i] = false;
            }

            for (var i = 0; i < sep_size; i++)
            {
                marked[sep_nodes[i]] = true;
            }

            for (var i = 0; i < n_left; i++)
            {
                if (!marked[i])
                {
                    for (var j = pointers[i]; j < pointers[i + 1]; j++)
                    {
                        var neighbor = indices[j];
                        if (!marked[neighbor])
                        {
                            Console.WriteLine("Edge ({0:d}, {1:d}) not covered", i, neighbor);
                        }
                    }
                }
            }

            var matchSize = match_size(matching, n_left);
            if (sep_size != matchSize)
            {
                Console.WriteLine("ERROR: sep_size = {0:d}, but match_size = {1:d}", sep_size, matchSize);
            }

            for (var i = 0; i < n_left + n_right; i++)
            {
                if (matching[i] != -1 && matching[matching[i]] != i)
                {
                    Console.WriteLine("ERROR: matching[{0:d}] = {1:d}, but matching[{2:d}] = {3:d}", i, matching[i], matching[i], matching[matching[i]]);
                }
            }
        }

        private static int match_size(int[] matching, int nleft)
        {
            var nmatch = 0;
            for (var i = 0; i < nleft; i++)
            {
                if (matching[i] != -1)
                {
                    ++nmatch;
                }
            }

            return nmatch;
        }
    }
}
