using System.Diagnostics;

namespace ChacoSharp.Graph
{
    public static unsafe class Subgraph
    {
        /// <summary>
        ///
        /// </summary>
        /// <param name="graph">graph data structure</param>
        /// <param name="subgraph">subgraph data structure</param>
        /// <param name="subnvtxs">number of vtxs in subgraph</param>
        /// <param name="psubnedges">ptr to number of edges in subgraph</param>
        /// <param name="assignment">values designating subgraph inclusion</param>
        /// <param name="set">assignment value indicating inclusion</param>
        /// <param name="glob2loc">mapping from graph to subgraph numbering</param>
        /// <param name="loc2glob">mapping from subgraph to graph numbering</param>
        /// <param name="degree">degrees of vertices in graph</param>
        /// <param name="useEdgeWeights">are edge weights being used?</param>
        public static void make_subgraph(vtx_data** graph, vtx_data** subgraph, int subnvtxs, int* psubnedges, int* assignment, int set, int* glob2loc, int* loc2glob, int* degree, bool useEdgeWeights)
        {
            Trace.WriteLine($"<Entering {nameof(make_subgraph)}>");

            var subgraphEdgeCount = 0;
            for (var i = 1; i <= subnvtxs; i++)
            {
                /* First get the vertices organized properly. */
                var subgptr = subgraph[i] = graph[loc2glob[i]]; /* loops through subgraph */
                var newnedges = degree[i] = subgptr->nedges; /* vertex degree in subgraph */

                /* Now work on the edges. */
                subgptr->edges[0] = i;
                double subgraphEdgeWeightTotal = 0; /* sum of weights of subgraph edges */

                /* Move all deactivated edges to the end of the list. */
                int* iptr = subgptr->edges + 1; /* loops through edge list */
                float* fptr = null; /* loops through edge weights */
                if (useEdgeWeights)
                {
                    fptr = subgptr->ewgts + 1;
                }

                for (var j = 1; j < newnedges;)
                {
                    var neighbor = *iptr; /* neighbor vertex in graph */
                    if (assignment[neighbor] == set)
                    {
                        /* Keep vertex in edge list. */
                        subgptr->edges[j] = glob2loc[neighbor];
                        if (useEdgeWeights)
                        {
                            subgraphEdgeWeightTotal += *fptr++;
                        }

                        j++;
                        iptr++;
                    }
                    else
                    {
                        /* Move vertex to back of edge list. */
                        --newnedges;
                        var tempvtx = subgptr->edges[newnedges];
                        subgptr->edges[newnedges] = neighbor;
                        *iptr = tempvtx;
                        if (useEdgeWeights)
                        {
                            var tempwgt = subgptr->ewgts[newnedges];
                            subgptr->ewgts[newnedges] = *fptr;
                            *fptr = tempwgt;
                        }
                    }
                }

                subgptr->nedges = newnedges;
                subgraphEdgeCount += newnedges;
                if (useEdgeWeights)
                {
                    subgptr->ewgts[0] = (float) -subgraphEdgeWeightTotal;
                }
            }

            *psubnedges = (subgraphEdgeCount - subnvtxs) / 2;
        }

        /// <summary>
        /// Undo the construction of the subgraph.
        /// </summary>
        /// <param name="subgraph">subgraph data structure</param>
        /// <param name="subnvtxs">number of vtxs in subgraph</param>
        /// <param name="loc2glob">mapping from subgraph to graph numbering</param>
        /// <param name="degree">degrees of vertices in graph</param>
        /// <param name="useEdgeWeights">are edge weights being used?</param>
        public static void remake_graph(vtx_data** subgraph, int subnvtxs, int* loc2glob, int* degree, bool useEdgeWeights)
        {
            Trace.WriteLine($"<Entering {nameof(remake_graph)}>");

            vtx_data* subgptr; /* loops through subgraph */
            double ewgtsum; /* sum of weights of subgraph edges */
            int nedges; /* vertex degree in subgraph */

            /* For each vertex in subgraph, expand the edge set back out. */
            for (var i = 1; i <= subnvtxs; i++)
            {
                subgptr = subgraph[i];
                subgptr->edges[0] = loc2glob[i];
                nedges = subgptr->nedges;
                /* Change edges back to global numbering system. */
                var iptr = subgptr->edges; /* loops through adjacency list */
                for (var j = nedges - 1; j != 0; j--)
                {
                    iptr++;
                    *iptr = loc2glob[*iptr];
                }

                subgptr->nedges = degree[i];

                /* Now get the diagonal value right. */
                if (useEdgeWeights)
                {
                    ewgtsum = 0;
                    var fptr = subgptr->ewgts; /* loops through edge weights */
                    for (var j = degree[i] - 1; j != 0; j--)
                    {
                        ewgtsum += *(++fptr);
                    }

                    subgptr->ewgts[0] = (float) -ewgtsum;
                }
            }
        }
    }
}
