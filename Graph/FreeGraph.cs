using System;
using System.Runtime.InteropServices;

namespace ChacoSharp.Graph
{
    public static unsafe class FreeGraph
    {
        /// <summary>
        /// Free a graph data structure.
        /// </summary>
        public static void free_graph(vtx_data** graph)
        {
            if (graph == null)
            {
                return;
            }

            if (graph[1] != null)
            {
                if (graph[1]->ewgts != null)
                {
                    Marshal.FreeHGlobal((IntPtr) graph[1]->ewgts);
                }

                if (graph[1]->edges != null)
                {
                    Marshal.FreeHGlobal((IntPtr) graph[1]->edges);
                }

                Marshal.FreeHGlobal((IntPtr) graph[1]);
            }

            Marshal.FreeHGlobal((IntPtr) graph);
        }
    }
}
