using System;
using static ChacoSharp.StaticConstants;

namespace ChacoSharp.Utilities
{
    public static unsafe class DivideProcs
    {
        public static bool divide_procs(int              architecture, /* 0 => hypercube, d => d-dimensional mesh */
                 int              ndims,        /* normal dimension of each cut */
                 int              ndims_tot,    /* total number of hypercube dimensions */
                 set_info *info_set,     /* data for all sets */
                 set_info *divide_set,   /* data for set being divided */
                 int[]             subsets,      /* subsets to be created */
                 bool              inert,        /* using inertial method? */
                 int *            pndims_real,  /* actual ndims for this cut */
                 int *            pnsets_real,  /* # sets created by this cut */
                 int *            pstriping,    /* cut in single direction? */
                 int []            cut_dirs,     /* direction of each cut if mesh */
                 int []            mesh_dims,    /* size of full mesh */
                 int[][]              hops_special/*[][MAXSETS]*/ /* hop matrix for nonstandard cases */
)
{
  int nsets_real = -1; /* number of sets to divide into */
  int ndims_real = -1; /* number of eigenvectors to use */
  int striping   = -1; /* cut in single direction? */
  bool flag       = true; /* unusual partition => use special hops */
  int ndim_poss;       /* largest dimensionality possible */
  int idims;           /* true dimensionality of subgrid */
  int i;               /* loop counter */

  if (architecture > 0) { /* Mesh, complicated case. */
    nsets_real = divide_set->span[0] * divide_set->span[1] * divide_set->span[2];
    nsets_real = Math.Min(1 << ndims, nsets_real);
    ndims_real = ndims;
    while (1 << ndims_real > nsets_real) {
      --ndims_real;
    }

    ndim_poss = 0;
    idims     = 0;
    for (i = 0; i < 3; i++) {
      if (divide_set->span[i] >= 2) {
        ndim_poss++;
        idims++;
      }
      if (divide_set->span[i] >= 4) {
        ndim_poss++;
      }
      if (divide_set->span[i] >= 8) {
        ndim_poss++;
      }
    }
    ndims_real = Math.Min(ndim_poss, ndims_real);

    if (idims > 1) {
      nsets_real = 1 << ndims_real;
    }

    flag = define_submeshes(nsets_real, architecture, mesh_dims, divide_set, info_set, subsets,
                            inert, &striping, cut_dirs, hops_special);
    if (striping != 0) {
      ndims_real = 1;
    }
  }

  else if (architecture == 0) { /* Hypercube, easy case. */
    ndims_real = Math.Min(ndims, divide_set->ndims);
    nsets_real = 1 << ndims_real;

    flag = define_subcubes(nsets_real, ndims_tot, ndims_real, divide_set, info_set, subsets, inert,
                           &striping, hops_special);

    if (striping != 0) {
      ndims_real = 1;
    }
  }

  *pndims_real = ndims_real;
  *pnsets_real = nsets_real;
  *pstriping   = striping;

  return (flag);
}

        private static bool define_subcubes(int              nsets_real, /* actual number of sets being created */
                    int              ndims_tot,  /* total hypercube dimensions */
                    int              ndims,      /* # dimension in this cut */
                    set_info *set,        /* data for set being divided */
                    set_info *set_info,   /* data for all sets */
                    int []            subsets,    /* subsets to be created */
                    bool              inert,      /* using inertial method? */
                    int *            pstriping,  /* cut in single direction? */
                    int[][]              hop_mtx_special/*[MAXSETS][MAXSETS]*/ /* nonstandard hop values */
)
{
  bool        hop_flag;  /* use special hop matrix? */
  int        nsets;     /* number of sets being created */
  int        setnum;    /* global number of subset */
  int        bits;      /* number of bits in which two sets differ */
  int        i, j, k;   /* loop counters */

  nsets    = 1 << ndims;
  hop_flag = false;

  for (k = nsets - 1; k >= 0; k--) { /* Backwards to not overwrite current set. */

    setnum                 = set->setnum | (k << (ndims_tot - set->ndims));
    set_info[setnum].ndims = set->ndims - ndims;
    subsets[k]             = setnum;
  }

  *pstriping = (inert && nsets_real > 2) ? 1 : 0;

  if (*pstriping != 0) { /* Gray code for better mapping. */
    for (k = 0; k < nsets; k++) {
      subsets[k] = gray(subsets[k]);
    }

    if (KL_METRIC == KernighanLinMetric.Hops) {
      hop_flag = true;
      for (i = 0; i < nsets; i++) {
        hop_mtx_special[i][i] = 0;
        for (j = 0; j < i; j++) {
          hop_mtx_special[i][j] = 0;
          bits                  = (subsets[i]) ^ (subsets[j]);
          while (bits != 0) {
            if ((bits & 1) != 0) {
              ++hop_mtx_special[i][j];
            }
            bits >>= 1;
          }
          hop_mtx_special[j][i] = hop_mtx_special[i][j];
        }
      }
    }
  }

  return (hop_flag);
}

        /* Figure out how to divide mesh into pieces.  Return true if nonstandard. */
private static bool define_submeshes(int              nsets,        /* number of subsets in this partition */
                     int              cube_or_mesh, /* 0=> hypercube, d=> d-dimensional mesh */
                     int []            mesh_dims,    /* shape of mesh */
                     set_info *set,          /* set data for set I'm partitioning */
                     set_info *set_info,     /* set data for all sets */
                     int []            subsets,      /* subsets being created by partition */
                     bool              inert,        /* using inertial method? */
                     int *            striping,     /* should I partition with parallel cuts? */
                     int []            dir,          /* directions of each cut */
                     int[][]              hop_mtx_special/*[MAXSETS][MAXSETS]*/ /* hops values if unusual */
)
{
  int        ndims;     /* dimension of cut */
  int[]        dims = new int[3];   /* local copy of mesh_dims to modify */
  int        maxdim;    /* longest dimension of the mesh */
  int        mindim;    /* intest dimension of mesh */
  int[]        start = new int[3];  /* start in each index of submesh */
  int[]        width = new int[3];  /* length in each index of submesh */
  int[]        nbits = new int[3];  /* values for computing hops */
  int[]        coords = new int[3]; /* location of set in logical grid */
  int[]        mask = new int[3];   /* values for computing hops */
  int        setnum;    /* number of created set */
  bool        flag;      /* return condition */
  bool        snaking;   /* is single stripe snaking through grid? */
  bool        reverse;   /* should I reverse direction for embedding? */
  int        i, j, k;   /* loop counters */

  dims[0] = set->span[0];
  dims[1] = set->span[1];
  dims[2] = set->span[2];

  ndims = 1;
  while ((2 << ndims) <= nsets) {
    ndims++;
  }

  /* Find the intest and longest directions in mesh. */
  maxdim = -1;
  mindim = dims[0];
  dir[1] = dir[2] = 0;
  for (i = 0; i < cube_or_mesh; i++) {
    if (dims[i] > maxdim) {
      maxdim = dims[i];
      dir[0] = i;
    }
    if (dims[i] < mindim) {
      mindim = dims[i];
    }
  }

  /* Decide whether or not to force striping. */
  i = 0;
  for (j = 0; j < cube_or_mesh; j++) {
    if (set->span[j] > 1) {
      i++;
    }
  }

  *striping = (i <= 1 || nsets == 3 ||
               (maxdim > nsets && (maxdim > .6 * nsets * mindim || (inert && nsets > 2)))) ? 1 : 0;

  snaking = !(*striping != 0) && inert && nsets > 2;

  if (!(*striping != 0)) {
    if (nsets >= 4) { /* Find direction of second & third cuts. */
      dims[dir[0]] /= 2;
      maxdim = -1;
      for (i = 0; i < cube_or_mesh; i++) {
        if (dims[i] > maxdim) {
          maxdim = dims[i];
          dir[1] = i;
        }
      }
    }

    if (nsets == 8) { /* Find a third direction. */
      dims[dir[1]] /= 2;
      maxdim = -1;
      for (i = 0; i < cube_or_mesh; i++) {
        if (dims[i] > maxdim) {
          maxdim = dims[i];
          dir[2] = i;
        }
      }
    }
    nbits[0] = nbits[1] = nbits[2] = 0;
    for (i = 0; i < ndims; i++) {
      ++nbits[dir[i]];
    }
    for (i = 0; i < 3; i++) {
      mask[i] = (1 << nbits[i]) - 1;
    }
    mask[1] <<= nbits[0];
    mask[2] <<= nbits[0] + nbits[1];
  }

  for (k = nsets - 1; k >= 0; k--) { /* Backwards to not overwrite current set. */

    for (i = 0; i < 3; i++) {
      start[i] = 0;
      width[i] = dims[i] = set->span[i];
    }

    if ((*striping) != 0) { /* Use longest direction for all cuts. */
      start[dir[0]] = (k * dims[dir[0]] + nsets - 1) / nsets;
      width[dir[0]] = ((k + 1) * dims[dir[0]] + nsets - 1) / nsets - start[dir[0]];
    }

    else { /* Figure out partition based on cut directions. */
      coords[0] = k & mask[0];
      coords[1] = (k & mask[1]) >> nbits[0];
      coords[2] = (k & mask[2]) >> (nbits[0] + nbits[1]);
      if (snaking) {
        reverse = (coords[1] & 1) != 0;
        if (reverse) {
          coords[0] = mask[0] - coords[0];
        }
        reverse = (coords[2] & 1) != 0;
        if (reverse) {
          coords[1] = (mask[1] >> nbits[0]) - coords[1];
        }
      }

      for (j = 0; j < ndims; j++) {
        --nbits[dir[j]];
        if ((coords[dir[j]] & (1 << nbits[dir[j]])) != 0) {
          /* Right side of partition. */
          start[dir[j]] += (width[dir[j]] + 1) / 2;
          width[dir[j]] /= 2;
        }
        else { /* Left side of partition */
          width[dir[j]] = (width[dir[j]] + 1) / 2;
        }
      }

      /* Now restore nbits values. */
      nbits[0] = nbits[1] = nbits[2] = 0;
      for (i = 0; i < ndims; i++) {
        ++nbits[dir[i]];
      }
    }

    for (i = 0; i < 3; i++) {
      start[i] += set->low[i];
    }

    setnum = (start[2] * mesh_dims[1] + start[1]) * mesh_dims[0] + start[0];

    for (i = 0; i < 3; i++) {
      set_info[setnum].low[i]  = start[i];
      set_info[setnum].span[i] = width[i];
    }
    subsets[k] = setnum;
  }

  /* Check to see if hop_mtx is nonstandard. */
  flag = false;
  if (KL_METRIC == KernighanLinMetric.Hops) {
    if ((*striping) != 0) {
      flag = true;
      for (i = 0; i < nsets; i++) {
        for (j = 0; j < nsets; j++) {
          hop_mtx_special[i][j] = Math.Abs(i - j);
        }
      }
    }

    else if (nsets == 4) {
      if (dir[0] == dir[1] || snaking) {
        flag = true;
        for (i = 0; i < nsets; i++) {
          start[0] = i & mask[0];
          start[1] = (i & mask[1]) >> nbits[0];
          if (snaking) {
            reverse = (start[1] & 1) != 0;
            if (reverse) {
              start[0] = mask[0] - start[0];
            }
          }
          for (j = i; j < nsets; j++) {
            coords[0] = j & mask[0];
            coords[1] = (j & mask[1]) >> nbits[0];
            if (snaking) {
              reverse = (coords[1] & 1) != 0;
              if (reverse) {
                coords[0] = mask[0] - coords[0];
              }
            }

            hop_mtx_special[i][j] = hop_mtx_special[j][i] =
                Math.Abs(start[0] - coords[0]) + Math.Abs(start[1] - coords[1]);
          }
        }
      }
    }
    else if (nsets == 8) {
      if (dir[0] == dir[1] || dir[0] == dir[2] || dir[1] == dir[2] || snaking) {
        flag = true;
        for (i = 0; i < nsets; i++) {
          start[0] = i & mask[0];
          start[1] = (i & mask[1]) >> nbits[0];
          start[2] = (i & mask[2]) >> (nbits[0] + nbits[1]);
          if (snaking) {
            reverse = (start[1] & 1) != 0;
            if (reverse) {
              start[0] = mask[0] - start[0];
            }
            reverse = (start[2] & 1) != 0;
            if (reverse) {
              start[1] = (mask[1] >> nbits[0]) - start[1];
            }
          }
          for (j = i; j < nsets; j++) {
            coords[0] = j & mask[0];
            coords[1] = (j & mask[1]) >> nbits[0];
            coords[2] = (j & mask[2]) >> (nbits[0] + nbits[1]);
            if (snaking) {
              reverse = (coords[1] & 1) != 0;
              if (reverse) {
                coords[0] = mask[0] - coords[0];
              }
              reverse = (coords[2] & 1) != 0;
              if (reverse) {
                coords[1] = (mask[1] >> nbits[0]) - coords[1];
              }
            }

            hop_mtx_special[i][j] = hop_mtx_special[j][i] =
                Math.Abs(start[0] - coords[0]) + Math.Abs(start[1] - coords[1]) + Math.Abs(start[2] - coords[2]);
          }
        }
      }
    }
  }

  *striping = ((*striping != 0) && snaking) ? 1 : 0;

  return (flag);
}

/* Compute the binary reflected Gray code of a value. */
private static int gray(int i) { return ((i >> 1) ^ i); }

/* Compute the inverse of the binary reflected Gray code of a value. */
/*
int       invgray(i)
int       i;
{
    int       k;
    k = i;
    while (k) {
        k >>= 1;
        i ^= k;
    }
    return (i);
}
*/
    }
}
