// ReSharper disable InconsistentNaming
// ReSharper disable UnusedMember.Global
namespace ChacoSharp
{
    public unsafe struct vtx_data
    {
        public int vwgt; /**< weight of vertex */

        public int nedges; /**< number of neighbors of vertex in subgraph */

        /**< Note: above always includes self-edge first */
        public int* edges; /**< neighbor list in subgraph numbering scheme */

        public float* ewgts; /**< weights of all the edges */
        /**< Note: above 2 fields have self-edge first */
    }

    /**< An Array of lists made of these stores scheduler's message table. */
    public unsafe struct msg_data
    {
        public int dest; /**< destination of the message */
        public double dur; /**< duration of message */
        public double beg; /**< time at which message begins */
        public double end; /**< time at which message end */
        public list* route; /**< linked list of ints stores message route */
        public msg_data* pntr; /**< pointer to next outgoing message from this set */
    }

    /**< A linked list of these stores the selective orthogonalization set */
    public unsafe struct orthlink
    {
        public int depth; /**< bottom of list is 0, previous is 1 etc */
        public int index; /**< position in list of ritz vals (i index) */
        public double ritzval; /**< good ritz value */
        public double betaji; /**< residual bound on good ritz pair */
        public double tau; /**< from orthogonality recursion */
        public double prevtau; /**< from orthogonality recursion */
        public double* vec; /**< vector to orthogonalize against */
        public orthlink* pntr; /**< pointer to next link */
    }

/**< A linked list of these stores the selective orthogonalization set */
    public unsafe struct orthlink_float
    {
        public int depth; /**< bottom of list is 0, previous is 1 etc */
        public int index; /**< position in list of ritz vals (i index) */
        public double ritzval; /**< good ritz value */
        public double betaji; /**< residual bound on good ritz pair */
        public double tau; /**< from orthogonality recursion */
        public double prevtau; /**< from orthogonality recursion */
        public float* vec; /**< vector to orthogonalize against */
        public orthlink_float* pntr; /**< pointer to next link */
    }

/**< Array data structure for heap information */
    public struct heap
    {
        public double val; /**< value being kept in a heap */
        public int tag; /**< info associated with value */
    }

/**< A linked list of these stores the minimum elements of a vector */
    public unsafe struct scanlink
    {
        public double val; /**< value of vector entry */
        public int indx; /**< index of corresponding entry */
        public scanlink* pntr; /**< pointer to next link */
    }

/**< These store the phantom edges needed to keep a subgraph connected */
    public unsafe struct edgeslist
    {
        public int vtx1; /**< first vertex in edge */
        public int vtx2; /**< second vertex in edge */
        public edgeslist* next; /**< pointer to next element in list */
    }

/**< These store all the data needed to modify edges for connectivity. */
    public unsafe struct connect_data
    {
        public ilists* old_edges; /**< overwritten old edges */
        public flists* old_ewgts; /**< overwritten old weights */
        public edgeslist* new_edges; /**< list of new edges */
        public int old_nedges; /**< original number of edges in graph */
    }

/**< Information about subsets of processors is needed in recurse. */
    public unsafe struct set_info
    {
        public int setnum; /**< assignment value for this set */
        public int ndims; /**< log of # processors if hypercube */
        public fixed int low[3]; /**< low limit for grid dimensions if mesh */
        public fixed int span[3]; /**< size of grid dimensions if mesh */
        public set_info* next; /**< pointer to next element in linked list */
    }

/**< Linked list stuff for various uses */
    public unsafe struct list
    {
        /**< linked list of integers */
        public int num; /**< element number */
        public list* next; /**< ptr to next element in list */
    }

    public unsafe struct lists
    {
        /**< linked list of lists */
        public list* begin; /**< pointer to list */
        public lists* nextlist; /**< next list header */
    }

    public unsafe struct bilist
    {
        /**< bidirectional list */
        public bilist* prev; /**< pointer to previous element */
        public bilist* next; /**< ptr to next element in list */
    }

    struct ipairs
    {
        /**< stores pairs of integers */
        public int val1;
        public int val2;
    }

    public unsafe struct ilists
    {
        /**< linked list of integer lists */
        public int* list;
        public ilists* next;
    }

    public unsafe struct flists
    {
        /**< linked list of floating lists */
        public float* list;
        public flists* next;
    }
}