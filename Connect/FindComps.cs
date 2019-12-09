using System.Runtime.InteropServices;
using static ChacoSharp.Utilities.Randomize;
using System;
using System.Diagnostics;

namespace ChacoSharp.Connect
{
    public static unsafe class FindComps
    {
        /// <summary>
        /// Breadth first search algorithm to find & mark connected components.
        /// </summary>
        /// <param name="graph">graph data structure</param>
        /// <param name="nvtxs">number of vertices in graph</param>
        /// <param name="mark">space for nvtxs+1 ints</param>
        /// <param name="vtxlist">space for nvtxs ints</param>
        /// <returns>The number of connected components.</returns>
        public static int find_comps(vtx_data** graph, int nvtxs, int* mark, int* vtxlist)
        {
            Trace.WriteLine($"<Entering {nameof(find_comps)}>");

            for (var i = 1; i <= nvtxs; i++)
            {
                mark[i] = -1;
            }

            var visitedVertexCount = 0;
            var componentCount = 0;
            // vertex to start the dfs
            var root = (int) (nvtxs * drandom()) + 1;

            bfsearch(graph, root, &visitedVertexCount, mark, vtxlist, componentCount);

            // Are there any remaining vertices?
            while (visitedVertexCount != nvtxs)
            {
                // Find starting vtx for next BFS.
                root = (int) (nvtxs * drandom()) + 1;
                while (mark[root] >= 0)
                {
                    root++;
                    if (root > nvtxs)
                    {
                        root = 1;
                    }
                }

                // Add new edge to list needed for connectivity.
                componentCount++;
                bfsearch(graph, root, &visitedVertexCount, mark, vtxlist, componentCount);
                Trace.WriteLine($"Visited {visitedVertexCount} out of {nvtxs} vertices");
            }

            Trace.WriteLine($"<Exiting {nameof(find_comps)}>");

            return componentCount + 1;
        }

        /// <summary>
        /// Breadth first search algorithm to find & mark connected components.
        /// Returns list of edges to connect them together.
        /// </summary>
        /// <param name="graph">graph data structure</param>
        /// <param name="nvtxs">number of vertices in graph</param>
        /// <param name="mark">space for nvtxs+1 ints</param>
        /// <param name="vtxlist">space for nvtxs ints</param>
        /// <param name="edges">list of edges connecting graph</param>
        /// <returns></returns>
        public static int find_edges(vtx_data** graph, int nvtxs, int* mark, int* vtxlist, edgeslist** edges)
        {
            Trace.WriteLine($"<Entering {nameof(find_edges)}>");
            for (var i = 1; i <= nvtxs; i++)
            {
                mark[i] = -1;
            }

            var visitedVertexCount = 0;
            var edgesAdded = 0;
            *edges = null;
            // vertex to start the dfs
            var root = (int) (nvtxs * drandom()) + 1;

            // last vertex seen in BFS
            var last = bfsearch(graph, root, &visitedVertexCount, mark, vtxlist, edgesAdded);

            while (visitedVertexCount != nvtxs)
            {
                // Are there any remaining vertices?
                // Find starting vtx for next BFS.
                root = (int) (nvtxs * drandom()) + 1;
                //Trace.WriteLine($"looking for next root from: {root}");
                while (mark[root] >= 0)
                {
                    root++;
                    if (root > nvtxs)
                    {
                        root = 1;
                    }
                }



                // Add new edge to list needed for connectivity.
                var newedge = (edgeslist*) Marshal.AllocHGlobal(sizeof(edgeslist));
                newedge->next = *edges;
                newedge->vtx1 = last;
                newedge->vtx2 = root;
                *edges = newedge;
                edgesAdded++;
                last = bfsearch(graph, root, &visitedVertexCount, mark, vtxlist, edgesAdded);
            }

            return edgesAdded;
        }

        /// <summary>
        /// BFS to find connected component.
        /// </summary>
        /// <param name="graph">graph data structure</param>
        /// <param name="root">start vertex for DFS</param>
        /// <param name="count">number of vertices in component</param>
        /// <param name="mark">has vtx been seen?</param>
        /// <param name="vtxlist">space for storing vtxs to search</param>
        /// <param name="currentComponentNumber">current component number</param>
        /// <returns></returns>
        public static int bfsearch(vtx_data** graph, int root, int* count, int* mark, int* vtxlist, int currentComponentNumber)
        {
            var firstVertex = 1;

            var lastVertexInList = 1;
            mark[root] = currentComponentNumber;
            vtxlist[0] = root;

            // Copy root's neighbors to vtxlist, incrementing count
            var edgePointer = graph[root]->edges;

            if (graph[root]->nedges - 1 < 0)
            {
                throw new InvalidOperationException("A node should have at least one edge");
            }
            for (var i = graph[root]->nedges - 1; i != 0; i--)
            {
                // neighbor of vertex
                var neighbor = *(++edgePointer);
                vtxlist[lastVertexInList++] = neighbor;
                mark[neighbor] = currentComponentNumber;
            }

            while (firstVertex < lastVertexInList)
            {
                // vertex being processed
                var vtx = vtxlist[firstVertex++];
                // Loop through neighbors, copying to vtxlist if unmarked.
                var iptr = graph[vtx]->edges;
                if (graph[vtx]->nedges - 1 < 0)
                {
                    throw new InvalidOperationException("A node should have at least one edge");
                }
                for (var i = graph[vtx]->nedges - 1; i != 0; i--)
                {
                    // neighbor of vertex
                    var neighbor = *(++iptr);
                    if (mark[neighbor] != currentComponentNumber)
                    {
                        mark[neighbor] = currentComponentNumber;
                        vtxlist[lastVertexInList++] = neighbor;
                    }
                }
            }

            *count += lastVertexInList;
            return vtxlist[lastVertexInList - 1];
        }
    }
}
