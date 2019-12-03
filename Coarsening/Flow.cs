#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
#pragma warning disable HAA0101 // Array allocation for params parameter
using System;
using System.Runtime.InteropServices;
using static ChacoSharp.StaticConstants;
// ReSharper disable PossibleNullReferenceException
// ReSharper disable InvertIf

namespace ChacoSharp.Coarsening
{
    public static unsafe class Flow
    {
/*
    STATUS:
        film.graph random+KL.  before size = 24, should stay 24.
        However, I'm not finding this much flow.

    Start w/ a simple greedy weighted matching.
    For each node on left side w/ some remaining unmatched weight:
        Look for augmenting path via recursive call.
        If my neighbor has uncommitted weight, take it.
        Otherwise, if some committed elsewhere, see who it is committed to.
            If this second order neighbor can get flow elsewhere, patch up augmenting path.

    When no further improvement possible:
        For all unsatisfied nodes on left, include right nodes in their search tree, and exclude left nodes.
        Include all other left nodes.

    This is similar to matching except:
        Each flow edge has a value associated with it.
        Vertex can have several flow neighbors.

    These differences require a different data structure.
        I need to modify the flow associated w/ an edge and see it immediately from either vertex.
        Use a single representation of flow on an edge, and have each edge representation point to it.
*/

/*
    The following code takes as input a weighted bipartite graph
        with n_left+n_right nodes
            where node 'i' is adjacent to nodes 'indices[pointers[i]]' through 'indices[pointers[i+1]-1]' (all 0-based)
    and returns a minimum weight edge cover (in sep_nodes[]).
*/

        public static void wbpcover(int n_left, /* number of vertices on left side */
            int n_right, /* number of vertices on right side */
            int* pointers, /* start/stop of adjacency lists */
            int* indices, /* adjacency list for each vertex */
            int* vweight, /* vertex weights */
            int* psep_size, /* returned size of separator */
            int* psep_weight, /* returned weight of separator */
            int** psep_nodes /* list of separator nodes */
        )
        {
            if (DEBUG_COVER != 0)
            {
                Console.WriteLine("-> Entering wbpcover, nleft = {0:d}, nright = {1:d}, 2*nedges = {2:d}", n_left, n_right, pointers[n_left + n_right] - pointers[0]);

                var wright = 0;
                var wleft = 0;
                var wedges = 0;
                for (var i = 0; i < n_left; i++)
                {
                    wleft += vweight[i];
                    for (var j = pointers[i]; j < pointers[i + 1]; j++)
                    {
                        wedges += vweight[i] * vweight[indices[j]];
                    }
                }

                for (var i = n_left; i < n_left + n_right; i++)
                {
                    wright += vweight[i];
                    for (var j = pointers[i]; j < pointers[i + 1]; j++)
                    {
                        wedges += vweight[i] * vweight[indices[j]];
                    }
                }

                Console.WriteLine("    Corresponds to unweighted, nleft = {0:d}, nright = {1:d}, 2*nedges = {2:d}", wleft, wright, wedges);
            }

            /* number of edges in bipartite graph */
            var nedges = pointers[n_left + n_right] - pointers[0];

            /* remaining, unmatched vertex weight */
            var resid = new int[n_left + n_right];

            /* flags for each vertex */
            var touched = new bool[n_left + n_right];

            /* flow on each right->left edge */
            var flow = new int[nedges + 1];

            /* Not a matching.  I can be connected to multiple nodes. */
            bpflow(n_left, n_right, pointers, indices, vweight, resid, flow, touched);

            reachability(n_left, n_right, pointers, indices, resid, flow, touched);

            /* Separator includes untouched nodes on left, touched on right. */
            /* Left separator nodes if unconnected to unmatched right node via */
            /* augmenting path, right separator nodes otherwise. */

            /* First count the separator size for malloc. */
            var separatorSize = 0; /* returned size of separator */
            for (var i = 0; i < n_left; i++)
            {
                if (!touched[i])
                {
                    separatorSize++;
                }
            }

            for (var i = n_left; i < n_left + n_right; i++)
            {
                if (touched[i])
                {
                    separatorSize++;
                }
            }

            /* list of separator nodes */
            int* separatorNodes = (int*) Marshal.AllocHGlobal((separatorSize + 1) * sizeof(int));

            separatorSize = 0;
            /* returned weight of separator */
            var sepWeight = 0;
            for (var i = 0; i < n_left; i++)
            {
                // ReSharper disable once InvertIf
                if (!touched[i])
                {
                    separatorNodes[separatorSize++] = i;
                    sepWeight += vweight[i];
                }
            }

            for (var i = n_left; i < n_left + n_right; i++)
            {
                // ReSharper disable once InvertIf
                if (touched[i])
                {
                    separatorNodes[separatorSize++] = i;
                    sepWeight += vweight[i];
                }
            }

            separatorNodes[separatorSize] = 0;

            *psep_size = separatorSize;
            *psep_weight = sepWeight;
            *psep_nodes = separatorNodes;

            /* Check the answer. */
            if (DEBUG_COVER != 0)
            {
                confirm_cover(n_left, n_right, pointers, indices, flow, vweight, resid, separatorSize, separatorNodes);
            }
        }

        private static void bpflow(int leftVertexCount, /* number of vertices on left side */
            int rightVertexCount, /* number of vertices on right side */
            int* pointers, /* start/stop of adjacency lists */
            int* indices, /* adjacency list for each vertex */
            int* vweight, /* vertex weights */
            int[] resid, /* residual weight at each vertex */
            int[] flow, /* flow on right->left edges */
            bool[] touched /* flags for each vertex */
        )
        {
            /* First mark all the vertices as unmatched & untouched. */
            for (var i = 0; i < leftVertexCount + rightVertexCount; i++)
            {
                resid[i] = vweight[i];
                touched[i] = false;
            }

            /* Note that I only keep flow in edges from right to left. */
            for (var i = pointers[leftVertexCount]; i < pointers[leftVertexCount + rightVertexCount]; i++)
            {
                flow[i] = 0;
            }

            /* Now generate a fast, greedy flow to start. */
            for (var i = leftVertexCount; i < leftVertexCount + rightVertexCount; i++)
            {
                for (var j = pointers[i]; j < pointers[i + 1] && resid[i] != 0; j++)
                {
                    var neighbor = indices[j]; /* neighboring vertex */
                    if (resid[neighbor] != 0)
                    {
                        var flow1 = Math.Min(resid[i], resid[neighbor]); /* flow through a particular edge */
                        resid[neighbor] -= flow1;
                        resid[i] -= flow1;
                        flow[j] = flow1;
                    }
                }
            }

            /* Now try to enlarge flow via augmenting paths. */

            var seen = new int[leftVertexCount + rightVertexCount];

            /* Look for an augmenting path. */
            for (var i = 0; i < leftVertexCount; i++)
            {
                var modified = true; /* was flow enlarged? */
                while (resid[i] != 0 && modified)
                {
                    modified = augment(i, pointers, indices, resid, flow, touched, seen);
                }
            }
        }

        private static bool augment(int node, /* start node in augmenting path */
            int* pointers, /* start/stop of adjacency lists */
            int* indices, /* adjacency list for each vertex */
            int[] resid, /* residual weight at each vertex */
            int[] flow, /* flow on right->left edges */
            bool[] touched, /* flags for each vertex */
            int[] seen /* keeps list of vertices encountered */
        )
        {
            /* Look for augmenting path in graph. */
            var nseen = 0;
            var flow1 = resid[node];
            touch(node, pointers, indices, resid, flow, touched, &flow1, seen, &nseen);

            if (flow1 != 0)
            {
                /* Found an augmenting path! */
                /* Free all the vertices encountered in search. */
                /* Otherwise they can't be involved in augmentation, */
                /* so leave them touched. */
                for (var i = 0; i < nseen; i++)
                {
                    touched[seen[i]] = false;
                }

                return true;
            }

            return false;
        }

/* Mark everybody in my alternating path tree, and return vertex at */
/* end of augmenting path if found. */
        private static void touch(int node, int* pointers, /* start/stop of adjacency lists */
            int* indices, /* adjacency list for each vertex */
            int[] resid, /* residual weight at each vertex */
            int[] flow, /* flow on right->left edges */
            bool[] touched, /* flags for each vertex */
            int* flow1, /* max flow we are looking for */
            int[] seen, /* list of vertices encountered */
            int* nseen /* number of vertices encountered */
        )
        {
            touched[node] = true;
            seen[(*nseen)++] = node;

            for (var i = pointers[node]; i < pointers[node + 1]; i++)
            {
                var neighbor = indices[i];
                if (!touched[neighbor])
                {
                    /* Not yet considered. */
                    touched[neighbor] = true;
                    seen[(*nseen)++] = neighbor;
                    if (resid[neighbor] != 0)
                    {
                        /* Has flow to spare. */
                        var flow2 = Math.Min(*flow1, resid[neighbor]);

                        /* Adjust flow & resid values. */
                        resid[neighbor] -= flow2;
                        resid[node] -= flow2;
                        /* Perhaps better to compute these upfront once? */
                        int k;
                        for (k = pointers[neighbor]; k < pointers[neighbor + 1] && indices[k] != node; k++) { }

                        flow[k] += flow2;
                        *flow1 = flow2;
                        return;
                    }

                    /* Has no flow to spare. */
                    for (var j = pointers[neighbor]; j < pointers[neighbor + 1]; j++)
                    {
                        if (flow[j] != 0 && !touched[indices[j]])
                        {
                            var flow2 = Math.Min(*flow1, flow[j]);
                            touch(indices[j], pointers, indices, resid, flow, touched, &flow2, seen, nseen);
                            if (flow2 != 0)
                            {
                                /* Found some flow to spare! */
                                /* Adjust flow & resid values. */
                                resid[indices[j]] += flow2;
                                resid[node] -= flow2;
                                flow[j] -= flow2;
                                int k;
                                for (k = pointers[neighbor]; k < pointers[neighbor + 1] && indices[k] != node; k++) { }

                                flow[k] += flow2;
                                *flow1 = flow2;
                                return;
                            }
                        }
                    }
                }
            }

            *flow1 = 0;
        }

        static void reachability(int leftVertexCount, /* number of vertices on left side */
            int rightVertexCount, /* number of vertices on right side */
            int* pointers, /* start/stop of adjacency lists */
            int* indices, /* adjacency list for each vertex */
            int[] resid, /* residual weight at each vertex */
            int[] flow, /* flow on right->left edges */
            bool[] touched /* flags for each vertex */
        )
        {
            int i; /* loop counter */

            /* Initialize all the vertices to be untouched */
            for (i = 0; i < leftVertexCount + rightVertexCount; i++)
            {
                touched[i] = false;
            }

            for (i = 0; i < leftVertexCount; i++)
            {
                if (!touched[i] && resid[i] != 0)
                {
                    touch2(i, pointers, indices, flow, touched);
                }
            }
        }

/* Mark everybody in my alternating path tree, and return vertex at */
/* end of augmenting path if found. */
        private static bool touch2(int node, int* pointers, /* start/stop of adjacency lists */
            int* indices, /* adjacency list for each vertex */
            int[] flow, /* flow on right->left edges */
            bool[] touched /* flags for each vertex */
        )
        {
            touched[node] = true;

            for (var i = pointers[node]; i < pointers[node + 1]; i++)
            {
                var neighbor = indices[i]; /* neighbor of a vertex */
                if (!touched[neighbor])
                {
                    /* Not yet considered. */
                    touched[neighbor] = true;
                    for (var j = pointers[neighbor]; j < pointers[neighbor + 1]; j++)
                    {
                        if (flow[j] != 0 && !touched[indices[j]])
                        {
                            var result = touch2(indices[j], pointers, indices, flow, touched);
                            if (result)
                            {
                                return true;
                            }
                        }
                    }
                }
            }

            return false;
        }

        private static void confirm_cover(int leftVertexCount, int rightVertexCount, int* pointers, int* indices,
            int[] flow, int* vweight,
            int[] resid, /* residual weight at each vertex */
            int sepSize, int* sepNodes)
        {
            var marked = new bool[leftVertexCount + rightVertexCount];

            for (var i = 0; i < leftVertexCount + rightVertexCount; i++)
            {
                marked[i] = false;
            }

            var sepWeight = 0;
            for (var i = 0; i < sepSize; i++)
            {
                marked[sepNodes[i]] = true;
                sepWeight += vweight[sepNodes[i]];
            }

            for (var i = 0; i < leftVertexCount; i++)
            {
                // ReSharper disable once InvertIf
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

            var countFlow = count_flow(leftVertexCount, rightVertexCount, pointers, flow);
            // ReSharper disable once ConvertIfStatementToConditionalTernaryExpression
            if (countFlow != sepWeight)
            {
                Console.WriteLine("ERROR: total_flow = {0:d}, sep_weight = {1:d}, sep_size = {2:d}", countFlow, sepWeight, sepSize);
            }
            else
            {
                Console.WriteLine("total_flow = {0:d}, sep_weight = {1:d}, sep_size = {2:d}", countFlow, sepWeight, sepSize);
            }

            /* Now check if any separator nodes have remaining flow. */
            count_resid(leftVertexCount, rightVertexCount, resid, vweight, marked);

            check_resid(leftVertexCount, rightVertexCount, vweight, resid, pointers, indices, flow);
        }

        static int count_flow(int leftVertexCount, int rightVertexCount, int* pointers, int[] flow)
        {
            var totalFlow = 0;
            for (var i = pointers[leftVertexCount]; i < pointers[leftVertexCount + rightVertexCount]; i++)
            {
                totalFlow += flow[i];
            }

            return totalFlow;
        }

        private static void count_resid(int leftVertexCount, int rightVertexCount, int[] resid, int* vweight, bool[] marked)
        {
            for (var i = 0; i < leftVertexCount + rightVertexCount; i++)
            {
                if (resid[i] < 0 || resid[i] > vweight[i])
                {
                    Console.WriteLine("BAD resid[{0:d}] = {1:d}, vweight = {2:d}", i, resid[i], vweight[i]);
                }
            }

            var leftUsed = 0;
            var rightUsed = 0;
            for (var i = 0; i < leftVertexCount; i++)
            {
                leftUsed += vweight[i] - resid[i];
                if (marked[i] && resid[i] != 0)
                {
                    Console.WriteLine("Vertex {0:d} in separator, but resid = {1:d} (vweight = {2:d})", i, resid[i], vweight[i]);
                }
            }

            for (var i = leftVertexCount; i < leftVertexCount + rightVertexCount; i++)
            {
                rightUsed += vweight[i] - resid[i];
                if (marked[i] && resid[i] != 0)
                {
                    Console.WriteLine("Vertex {0:d} in separator, but resid = {1:d} (vweight = {2:d})", i, resid[i], vweight[i]);
                }
            }

            if (leftUsed != rightUsed)
            {
                Console.WriteLine("left_used = {0:d}, NOT EQUAL TO right_used = {1:d}", leftUsed, rightUsed);
            }
        }

        private static void check_resid(int leftVertexCount, int rightVertexCount, int* vweight, /* vertex weights */
            int[] resid, int* pointers, int* indices, int[] flow)
        {
            for (var i = 0; i < leftVertexCount + rightVertexCount; i++)
            {
                if (resid[i] < 0 || resid[i] > vweight[i])
                {
                    Console.WriteLine("BAD resid[{0:d}] = {1:d}, vweight = {2:d}", i, resid[i], vweight[i]);
                }
            }

            int leftUsed = 0, rightUsed = 0;
            for (var i = 0; i < leftVertexCount; i++)
            {
                leftUsed += vweight[i] - resid[i];
            }

            for (var i = leftVertexCount; i < leftVertexCount + rightVertexCount; i++)
            {
                rightUsed += vweight[i] - resid[i];
            }

            if (leftUsed != rightUsed)
            {
                Console.WriteLine("left_used = {0:d}, NOT EQUAL TO right_resid = {1:d}", leftUsed, rightUsed);
            }

            var diff = new int[leftVertexCount + rightVertexCount];
            for (var i = 0; i < leftVertexCount + rightVertexCount; i++)
            {
                diff[i] = 0;
            }

            for (var i = leftVertexCount; i < leftVertexCount + rightVertexCount; i++)
            {
                int j;
                for (j = pointers[i]; j < pointers[i + 1]; j++)
                {
                    if (flow[j] < 0)
                    {
                        Console.WriteLine("Negative flow ({0:d},{1:d}) = {2:d}", i, indices[j], flow[j]);
                    }

                    diff[i] += flow[j];
                    diff[indices[j]] += flow[j];
                }
            }

            for (var i = 0; i < leftVertexCount + rightVertexCount; i++)
            {
                if (diff[i] != vweight[i] - resid[i])
                {
                    Console.WriteLine("ERROR: diff[{0:d}] = {1:d}, but vweight = {2:d} and resid = {3:d}", i, diff[i], vweight[i], resid[i]);
                }
            }
        }

    }
}
