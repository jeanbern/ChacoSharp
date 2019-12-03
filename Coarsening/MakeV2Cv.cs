namespace ChacoSharp.Coarsening
{
    public static unsafe class MakeV2Cv
    {
        public static void makev2cv(
            /* Construct mapping from original graph vtxs to coarsened graph vtxs. */
            int* mflag, /* flag indicating vtx selected or not */
            int nvtxs, /* number of vtxs in original graph */
            int* v2cv /* mapping from vtxs to coarsened vtxs */
        )
        {
            var j = 1;
            for (var i = 1; i <= nvtxs; i++)
            {
                if (mflag[i] == 0 || mflag[i] > i)
                {
                    v2cv[i] = j++;
                }
                else
                {
                    v2cv[i] = v2cv[mflag[i]];
                }
            }
        }
    }
}
