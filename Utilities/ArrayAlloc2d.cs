using System.Runtime.InteropServices;

namespace ChacoSharp.Utilities
{
    public static unsafe class ArrayAlloc2d
    {
        /* Dynamically allocate a 2 dimensional array. */
/* Return instead of dying if out of space. */

        public static T* array_alloc_2D_ret<T>(int dim1, int dim2, int size)
        where T : unmanaged
/* size of first dimension */
/* size of second dimension */
/* size of array elements */
        {
            int total;       /* Total size of the array */
            int aligned_dim; /* dim1 or dim1+1 to ensure data alignment */
            int offset;      /* offset of array elements */
            T * field;       /* The multi-dimensional array */
            T **ptr;         /* Pointer offset */
            T * data;        /* Data offset */
            int j;           /* loop counter */

            aligned_dim = (dim1 % 2) != 0 ? dim1 + 1 : dim1;
            offset      = aligned_dim * sizeof(void *);
            total       = offset + dim1 * dim2 * size;
            field       = (T*)Marshal.AllocHGlobal(total);

            if (field != null) {
                ptr  = (T **)field;
                data = field;
                data += offset;
                for (j = 0; j < dim1; j++) {
                    ptr[j] = data + j * size * dim2;
                }
            }

            return (field);
        }
    }
}
