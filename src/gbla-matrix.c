/* gbla: Gröbner Basis Linear Algebra
 * Copyright (C) 2018 Christian Eder <ederc@mathematik.uni-kl.de>
 * This file is part of gbla.
 * gbla is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * gbla is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with gbla . If not, see <http://www.gnu.org/licenses/>.
 */

#include "gbla-matrix.h"

#define NOT_DENSE_COPYING 0

sm_t *sort_schreyer_matrix(sm_t *M)
{
  ri_t i, j, mi;
  re_t *temp_rows;
  ci_t *temp_pos;
  ci_t temp_rwidth;

  // slow selection sort, needs merge/quick sort
  for (i=0; i<M->nrows-1; ++i) {
    mi  = i;
    for (j=i+1; j<M->nrows; ++j) {
      if (M->pos[j][0]<=M->pos[mi][0])
        mi  = j;
    }
    temp_pos    = M->pos[mi];
    temp_rows   = M->rows[mi];
    temp_rwidth = M->rwidth[mi];

    for (j=mi; j>i; --j) {
      M->pos[j]     = M->pos[j-1];
      M->rows[j]    = M->rows[j-1];
      M->rwidth[j]  = M->rwidth[j-1];
    }
    M->pos[i]     = temp_pos;
    M->rows[i]    = temp_rows;
    M->rwidth[i]  = temp_rwidth;
  }

  return M;
}

#if GBLA_WITH_FFLAS
void copy_block_ml_matrix_to_dns_matrix(sbm_fl_t **source, DNS **destination) {
  sbm_fl_t *src = *source;
  DNS *dst      = *destination;

  const ri_t src_row_idx  = (ri_t) ceil((float) src->nrows / (float)src->bheight);
  const ci_t src_col_idx  = (ci_t) ceil((float) src->ncols / (float)src->bwidth);
  const bi_t ml_bheight   = src->bheight / __GBLA_NROWS_MULTILINE;

  const ri_t max_rows  = dst->row;

  ri_t i;
  ci_t j;
  bi_t k, l;

  bi_t pos;

  ri_t dst_row_idx, dst_col_idx;
  for (i=0; i<src_row_idx; ++i) {
    for (j=0; j<src_col_idx; ++j) {
      if (src->blocks[i][j] == NULL)
        continue;
      dst_row_idx = i * src->bheight;
      dst_col_idx = j * src->bwidth;
      for (k=0; k<ml_bheight; ++k) {
        if (src->blocks[i][j][k].sz == 0) {
          dst_row_idx +=  2;
          continue;
        }
        if (src->blocks[i][j][k].dense == 0) {
          for (l=0; l<src->blocks[i][j][k].sz; ++l) {
            pos = src->blocks[i][j][k].idx[l];
            dst->ptr[dst_row_idx * dst->ld + dst_col_idx + pos] =
              (elemt_t)src->blocks[i][j][k].val[2*l];
          }
          dst_row_idx++;
          if (dst_row_idx<max_rows) {
            for (l=0; l<src->blocks[i][j][k].sz; ++l) {
              pos = src->blocks[i][j][k].idx[l];
              dst->ptr[dst_row_idx * dst->ld + dst_col_idx + pos] =
                (elemt_t)src->blocks[i][j][k].val[2*l+1];
            }
          }
          dst_row_idx++;
        } else {
          for (l=0; l<src->blocks[i][j][k].sz; ++l) {
            dst->ptr[dst_row_idx * dst->ld + dst_col_idx + l] =
              (elemt_t)src->blocks[i][j][k].val[2*l];
          }
          dst_row_idx++;
          if (dst_row_idx<max_rows) {
            for (l=0; l<src->blocks[i][j][k].sz; ++l) {
              /* printf("position %ld\n",dst_row_idx * dst->ld + dst_col_idx + l/2-1); */
              dst->ptr[dst_row_idx * dst->ld + dst_col_idx + l] =
                (elemt_t)src->blocks[i][j][k].val[2*l+1];
            }
          }
          dst_row_idx++;
        }
      }
    }
  }
  /*  free memory in block multiline matrix */
  for (i=0; i<src_row_idx; ++i) {
    for (j=0; j<src_col_idx; ++j) {
      if (src->blocks[i][j] == NULL)
        continue;
      for (k=0; k<ml_bheight; ++k) {
        free(src->blocks[i][j][k].idx);
        free(src->blocks[i][j][k].val);
      }
      free(src->blocks[i][j]);
      src->blocks[i][j] = NULL;
    }
    free(src->blocks[i]);
    src->blocks[i]  = NULL;
  }
  free(src->blocks);
  src->blocks = NULL;
  free(src);
  src = NULL;

  *source       = src;
  *destination  = dst;
}
#endif

void copy_block_ml_matrices_to_sparse_matrix(sbm_fl_t **input_bl,
    sm_fl_ml_t **input_ml, ri_t rank_input_ml, sm_t **output,
    int deleteIn, int nthrds) {

  sbm_fl_t *in_bl   = *input_bl;
  sm_fl_ml_t *in_ml = *input_ml;
  sm_t *out         = *output;

  ri_t i;
  ci_t j;
  bi_t k;
  bi_t l;

	/* ri_t rank = rank_input_ml / __GBLA_NROWS_MULTILINE + rank_input_ml % __GBLA_NROWS_MULTILINE; */
  /*  initialize meta data for multiline matrix out */
  out->nrows    = in_bl->nrows + rank_input_ml; /*  row dimension */
  out->ncols    = in_bl->ncols;                 /*  col dimension */

  /*  allocate memory for blocks */
  ri_t rl = out->nrows;

  out->rows   = (re_t **)malloc(rl * sizeof(re_t *));
  out->pos    = (ci_t **)malloc(rl * sizeof(ci_t *));
  out->rwidth = (ci_t *) malloc(rl * sizeof(ci_t));
  memset(out->rwidth, 0, rl * sizeof(ci_t));
  const ri_t rlin_bl = (ri_t) ceil((float) in_bl->nrows / (float)in_bl->bheight);
  const ci_t clin_bl = (ci_t) ceil((float) in_bl->ncols / (float)in_bl->bwidth);
  /*  we need buffers for all multiline entries since the copying process */
  /*  revisits already filled up multilines.if we reset the buffer to zero, we */
  /*  might realloc only init_buffer memory and lose what we have already in the */
  /*  multiline stored. */

  ci_t block_idx;
  ri_t block_vert;
  ri_t sparse_idx;
  re_t v1, v2;
	/* mli_t crb = 0; */
  const bi_t ml_bheight = in_bl->bheight / __GBLA_NROWS_MULTILINE;
  /*  write D first into output matrix */

  const ri_t rlin_ml = (ri_t) ceil((float) in_ml->nrows / (float)__GBLA_NROWS_MULTILINE);
  const ci_t clin_ml = in_ml->ncols;

  /*  we need buffers for all multiline entries since the copying process */
  /*  revisits already filled up multilines.if we reset the buffer to zero, we */
  /*  might realloc only init_buffer memory and lose what we have already in the */
  /*  multiline stored. */
  sparse_idx  = 0;
  for (i=0; i<rlin_ml; ++i) { /*  loop over multilines */
    if (out->rwidth[sparse_idx] > 0)
      sparse_idx++;
    if (out->rwidth[sparse_idx] > 0)
      sparse_idx++;
    if (in_ml->ml[i].val == NULL || in_ml->ml[i].sz == 0) {
      continue;
    }
    out->rows[sparse_idx]   = (re_t *)malloc(out->ncols * sizeof(re_t));
    out->pos[sparse_idx]    = (ci_t *)malloc(out->ncols * sizeof(ci_t));
    out->rows[sparse_idx+1] = (re_t *)malloc(out->ncols * sizeof(re_t));
    out->pos[sparse_idx+1]  = (ci_t *)malloc(out->ncols * sizeof(ci_t));
    if (in_ml->ml[i].dense == 0) {
      for (j=0; j<clin_ml; ++j) {
        v1  = in_ml->ml[i].val[2*j];
        v2  = in_ml->ml[i].val[2*j+1];
        if (v1 != 0) {
          out->rows[sparse_idx][out->rwidth[sparse_idx]]  = v1;
          out->pos[sparse_idx][out->rwidth[sparse_idx]]   = in_ml->ml[i].idx[j];
          out->rwidth[sparse_idx]++;
        }
        if (v2 != 0) {
          out->rows[sparse_idx+1][out->rwidth[sparse_idx+1]]  = v2;
          out->pos[sparse_idx+1][out->rwidth[sparse_idx+1]]   = in_ml->ml[i].idx[j];
          out->rwidth[sparse_idx+1]++;
        }
      }
    } else {
      for (j=0; j<clin_ml; ++j) {
        v1  = in_ml->ml[i].val[2*j];
        v2  = in_ml->ml[i].val[2*j+1];
        if (v1 != 0) {
          out->rows[sparse_idx][out->rwidth[sparse_idx]]  = v1;
          out->pos[sparse_idx][out->rwidth[sparse_idx]]   = j;
          out->rwidth[sparse_idx]++;
        }
        if (v2 != 0) {
          out->rows[sparse_idx+1][out->rwidth[sparse_idx+1]]  = v2;
          out->pos[sparse_idx+1][out->rwidth[sparse_idx+1]]   = j;
          out->rwidth[sparse_idx+1]++;
        }
      }
    }
    out->rows[sparse_idx] = (re_t*) realloc(out->rows[sparse_idx],
        out->rwidth[sparse_idx]*sizeof(re_t));
    out->pos[sparse_idx]  = (ci_t*) realloc(out->pos[sparse_idx],
        out->rwidth[sparse_idx]*sizeof(ci_t));
    out->rows[sparse_idx+1] = (re_t*) realloc(out->rows[sparse_idx+1],
        out->rwidth[sparse_idx+1]*sizeof(re_t));
    out->pos[sparse_idx+1]  = (ci_t*) realloc(out->pos[sparse_idx+1],
        out->rwidth[sparse_idx+1]*sizeof(ci_t));
  }
  if (deleteIn == 1) {
    /*  free memory for input matrix */
#pragma omp parallel for num_threads(nthrds)
    for (i=0; i<rlin_ml; ++i) {
      free(in_ml->ml[i].idx);
      in_ml->ml[i].idx = NULL;
      free(in_ml->ml[i].val);
      in_ml->ml[i].val = NULL;
    }
    free(in_ml->ml);
    in_ml->ml  = NULL;
    free(in_ml);
    in_ml  = NULL;
  }
  if (out->rwidth[sparse_idx] > 0)
    sparse_idx++;
  if (out->rwidth[sparse_idx] > 0)
    sparse_idx++;

  const ri_t offset = sparse_idx;
  /*  write B second into output matrix */

  for (i=0; i<rlin_bl; ++i) { /*  loop over multilines */
    if (in_bl->blocks[i] == NULL) {
      continue;
    }
    block_vert  = i * in_bl->bheight;
    sparse_idx  = block_vert + offset;
    k = 0;
    while (sparse_idx+k < out->nrows && k<in_bl->bheight) {
      out->rows[sparse_idx+k] = (re_t *)malloc(out->ncols * sizeof(re_t));
      out->pos[sparse_idx+k]  = (ci_t *)malloc(out->ncols * sizeof(ci_t));
      ++k;
    }
    for (j=0; j<clin_bl; ++j) {
      block_idx = j * in_bl->bwidth;
      if (in_bl->blocks[i][j] == NULL) {
        continue;
      }
      for (k=0;k<ml_bheight; ++k) {
        if (in_bl->blocks[i][j][k].sz == 0) {
          continue;
        }
        if (in_bl->blocks[i][j][k].dense == 0) {
          for (l=0; l<in_bl->blocks[i][j][k].sz; ++l) {
            /*  fill in data */
            v1  = in_bl->blocks[i][j][k].val[2*l];
            v2  = in_bl->blocks[i][j][k].val[2*l+1];
            if (v1 != 0) {
              out->rows[sparse_idx+(2*k)][out->rwidth[sparse_idx+(2*k)]]  = v1;
              out->pos[sparse_idx+(2*k)][out->rwidth[sparse_idx+(2*k)]]   = block_idx + in_bl->blocks[i][j][k].idx[l];
              out->rwidth[sparse_idx+(2*k)]++;
            }
            if (v2 != 0) {
              out->rows[sparse_idx+(2*k+1)][out->rwidth[sparse_idx+(2*k+1)]]  = v2;
              out->pos[sparse_idx+(2*k+1)][out->rwidth[sparse_idx+(2*k+1)]]   = block_idx + in_bl->blocks[i][j][k].idx[l];
              out->rwidth[sparse_idx+(2*k+1)]++;
            }
          }
        } else {
          for (l=0; l<in_bl->blocks[i][j][k].sz; ++l) {
            /*  fill in data */
            v1  = in_bl->blocks[i][j][k].val[2*l];
            v2  = in_bl->blocks[i][j][k].val[2*l+1];
            if (v1 != 0) {
              out->rows[sparse_idx+(2*k)][out->rwidth[sparse_idx+(2*k)]]  = v1;
              out->pos[sparse_idx+(2*k)][out->rwidth[sparse_idx+(2*k)]]   = block_idx + l;
              out->rwidth[sparse_idx+(2*k)]++;
            }
            if (v2 != 0) {
              out->rows[sparse_idx+(2*k+1)][out->rwidth[sparse_idx+(2*k+1)]]  = v2;
              out->pos[sparse_idx+(2*k+1)][out->rwidth[sparse_idx+(2*k+1)]]   = block_idx + l;
              out->rwidth[sparse_idx+(2*k+1)]++;
            }
          }
        }
      }
    }
    k = 0;
    while (sparse_idx+k < out->nrows && k<in_bl->bheight) {
      out->rows[sparse_idx+k] = (re_t*) realloc(out->rows[sparse_idx+k],
          out->rwidth[sparse_idx+k]*sizeof(re_t));
      out->pos[sparse_idx+k]  = (ci_t*) realloc(out->pos[sparse_idx+k],
          out->rwidth[sparse_idx+k]*sizeof(ci_t));
      ++k;
    }
  }

  if (deleteIn == 1) {
    /*  free memory for input matrix */
#pragma omp parallel for private(i,j,k) num_threads(nthrds)
    for (i=0; i<rlin_bl; ++i) {
      for (j=0; j<clin_bl; ++j) {
        if (in_bl->blocks[i][j] != NULL) {
          for (k=0; k<ml_bheight; ++k) {
            free(in_bl->blocks[i][j][k].idx);
            in_bl->blocks[i][j][k].idx = NULL;

            free(in_bl->blocks[i][j][k].val);
            in_bl->blocks[i][j][k].val = NULL;
          }
          free(in_bl->blocks[i][j]);
          in_bl->blocks[i][j] = NULL;
        }
      }
      free(in_bl->blocks[i]);
      in_bl->blocks[i] = NULL;
    }
    free(in_bl->blocks);
    in_bl->blocks  = NULL;
    free(in_bl);
    in_bl  = NULL;
    /*  free memory for input matrix */
  }
  *input_bl = in_bl;
  *input_ml = in_ml;

	out->nnz = 0 ;
	for (i=0 ; i < out->nrows ; ++i)
		out->nnz += out->rwidth[i] ;
	out->density = (float)compute_density(out->nnz, out->nrows, out->ncols);

  *output   = out;
}

sm_fl_ml_t *copy_block_matrix_to_multiline_matrix(sbm_fl_t **input,
    sm_fl_ml_t *out, int deleteIn, int nthrds) {

  sbm_fl_t *in  = *input;

  ri_t i, ii;
  ci_t j;
  bi_t k;
  bi_t l;

  /*  initialize meta data for multiline matrix out */
  out->nrows    = in->nrows;  /*  row dimension */
  out->ncols    = in->ncols;  /*  col dimension */
  out->ba       = dtrl;       /*  block alignment */
  out->fe       = 0;          /*  fill empty blocks? */
  out->hr       = 0;          /*  allow hybrid rows? */
  out->nnz      = 0;          /*  number nonzero elements */
  out->density  = (double)0;

  /*  allocate memory for blocks */
  ri_t rl = out->nrows / 2;
  if (out->nrows % 2)
    rl++;

  out->ml = (ml_t *)malloc(rl * sizeof(ml_t));
  for (i=0; i<rl; ++i) {
#if NOT_DENSE_COPYING
    out->ml[i].val    = NULL;
    out->ml[i].idx    = NULL;
    out->ml[i].sz     = 0;
    out->ml[i].dense  = 0;
#else
    out->ml[i].val    = (re_t *)malloc(2 * out->ncols * sizeof(re_t));
    out->ml[i].idx    = NULL;
    out->ml[i].sz     = out->ncols;
    out->ml[i].dense  = 1;
    /* out->ml[i].val  = realloc(out->ml[i].val,2 * out->ncols * sizeof(re_t)); */
#endif
  }

  const ri_t rlin = (ri_t) ceil((float) in->nrows / (float)in->bheight);
  const ci_t clin = (ci_t) ceil((float) in->ncols / (float)in->bwidth);
	/* mli_t init_buffer = 2 * in->bwidth; */
  /*  we need buffers for all multiline entries since the copying process */
  /*  revisits already filled up multilines.if we reset the buffer to zero, we */
  /*  might realloc only init_buffer memory and lose what we have already in the */
  /*  multiline stored. */
  mli_t *buffer = (mli_t *)malloc(rl * sizeof(mli_t));
  memset(buffer, 0, rl * sizeof(mli_t));

  ci_t block_idx;
  re_t v1, v2;
	/* mli_t crb = 0; */
  const bi_t ml_bheight = in->bheight / __GBLA_NROWS_MULTILINE;
  #pragma omp parallel shared(buffer) num_threads(nthrds)
  {
    #pragma omp for private(i,ii,j,k,l,block_idx,v1,v2) nowait
    for (i=0; i<rl; ++i) { /*  loop over multilines */
      memset(out->ml[i].val, 0, 2 * out->ncols * sizeof(re_t));
      ii  = i / ml_bheight;
      k   = i % ml_bheight;
      for (j=0; j<clin; ++j) {
        block_idx = j * in->bwidth;
        if (in->blocks[ii][j] == NULL) {
          continue;
        }
        if (in->blocks[ii][j][k].sz == 0) {
          continue;
        }
        if (in->blocks[ii][j][k].dense == 0) {
          for (l=0; l<in->blocks[ii][j][k].sz; ++l) {
            /*  fill in data */
            v1  = in->blocks[ii][j][k].val[2*l];
            v2  = in->blocks[ii][j][k].val[2*l+1];
            if (v1 != 0 || v2 != 0) {
              out->ml[i].val[2*(block_idx + in->blocks[ii][j][k].idx[l])]   = v1;
              out->ml[i].val[2*(block_idx + in->blocks[ii][j][k].idx[l])+1] = v2;
            }
            buffer[i]++;
          }
        } else {
          for (l=0; l<in->blocks[ii][j][k].sz; ++l) {
            /*  fill in data */
            v1  = in->blocks[ii][j][k].val[2*l];
            v2  = in->blocks[ii][j][k].val[2*l+1];
            if (v1 != 0 || v2 != 0) {
              out->ml[i].val[2*(block_idx+l)]   = v1;
              out->ml[i].val[2*(block_idx+l)+1] = v2;
            }
            buffer[i]++;
          }
        }
      }
      if (buffer[i] == 0) {
        free (out->ml[i].val);
        out->ml[i].val  = NULL;
        out->ml[i].sz   = 0;
      }
    }
  }
  if (deleteIn) {
  /*  free memory for input matrix */
  #pragma omp parallel num_threads(nthrds)
    {
#pragma omp for private(i,j,k) nowait
      for (i=0; i<rlin; ++i) {
        for (j=0; j<clin; ++j) {
          if (in->blocks[i][j] != NULL) {
            for (k=0; k<ml_bheight; ++k) {
              free(in->blocks[i][j][k].idx);
              in->blocks[i][j][k].idx = NULL;

              free(in->blocks[i][j][k].val);
              in->blocks[i][j][k].val = NULL;
            }
            free(in->blocks[i][j]);
            in->blocks[i][j] = NULL;
          }
        }
        free(in->blocks[i]);
        in->blocks[i] = NULL;
      }
    }
    free(in->blocks);
    in->blocks  = NULL;
    free(in);
    in  = NULL;
  }
/*
 *  NOTE:
 *  old and deprecated code: keep for possible tests with sparse-hybrid-dense
 *  representations for D in multiline structures later on
 */
#if 0
  #pragma omp parallel shared(buffer) num_threads(nthrds)
  {
    mli_t crb = 0;
    #pragma omp for private(i,j,k,l,block_idx,crb) nowait
    for (i=0; i<rlin; ++i) {
      /*  curr_row_base in LELA */
      crb  = i * in->bheight / __GBLA_NROWS_MULTILINE;

      for (j=0; j<clin; ++j) {
        if (in->blocks[i][j] == NULL) {
          continue;
        }

        block_idx  = in->bwidth * j;
        for (k=0; k<in->bheight/__GBLA_NROWS_MULTILINE; ++k) {
#if NOT_DENSE_COPYING
          if (in->blocks[i][j][k].sz == 0) {
            continue;
          }
          if (in->blocks[i][j][k].dense == 0) {
            for (l=0; l<in->blocks[i][j][k].sz; ++l) {
              /*  possibly realloc more memory */
              if (out->ml[crb+k].sz == buffer[crb+k]) {
                buffer[crb+k] +=  init_buffer;
                out->ml[crb+k].idx = realloc(out->ml[crb+k].idx,
                    buffer[crb+k] * sizeof(mli_t));
                out->ml[crb+k].val = realloc(out->ml[crb+k].val,
                    2 * buffer[crb+k] * sizeof(re_t));
              }
              /*  fill in data */
              out->ml[crb+k].idx[out->ml[crb+k].sz] =
                block_idx + in->blocks[i][j][k].idx[l];
              out->ml[crb+k].val[2*out->ml[crb+k].sz] =
                in->blocks[i][j][k].val[2*l];
              out->ml[crb+k].val[2*out->ml[crb+k].sz+1] =
                in->blocks[i][j][k].val[2*l+1];
              out->ml[crb+k].sz++;
            }
          } else { /*  block submatrix in is dense */
            for (l=0; l<in->blocks[i][j][k].sz; ++l) {
              /*  possibly realloc more memory */
              if (out->ml[crb+k].sz == buffer[crb+k]) {
                buffer[crb+k]  +=  init_buffer;
                out->ml[crb+k].idx = realloc(out->ml[crb+k].idx,
                    buffer[crb+k] * sizeof(mli_t));
                out->ml[crb+k].val = realloc(out->ml[crb+k].val,
                    2 * buffer[crb+k] * sizeof(re_t));
              }
              /*  fill in data */
              v1  = in->blocks[i][j][k].val[2*l];
              v2  = in->blocks[i][j][k].val[2*l+1];
              if (v1 != 0 || v2 != 0) {
                out->ml[crb+k].idx[out->ml[crb+k].sz] =
                  block_idx + l;
                out->ml[crb+k].val[2*out->ml[crb+k].sz]   = v1;
                out->ml[crb+k].val[2*out->ml[crb+k].sz+1] = v2;
                out->ml[crb+k].sz++;
              }
            }
          }
#else
          if (in->blocks[i][j][k].sz == 0) {
            continue;
          }
          if (in->blocks[i][j][k].dense == 0) {
            for (l=0; l<in->blocks[i][j][k].sz; ++l) {
              printf("%d - %d - %d - %d - %d - %d\n",i,j,k,l,crb,block_idx);
              /*  fill in data */
              v1  = in->blocks[i][j][k].val[2*l];
              v2  = in->blocks[i][j][k].val[2*l+1];
              if (v1 != 0 || v2 != 0) {
                out->ml[crb+k].val[2*(block_idx + in->blocks[i][j][k].idx[l])]   = v1;
                out->ml[crb+k].val[2*(block_idx + in->blocks[i][j][k].idx[l])+1] = v2;
              }
              buffer[crb+k]++;
            }
          } else {
            for (l=0; l<in->blocks[i][j][k].sz; ++l) {
              printf("%d - %d - %d - %d - %d - %d\n",i,j,k,l,crb,block_idx);
              /*  fill in data */
              v1  = in->blocks[i][j][k].val[2*l];
              v2  = in->blocks[i][j][k].val[2*l+1];
              if (v1 != 0 || v2 != 0) {
                out->ml[crb+k].val[2*(block_idx+l)]   = v1;
                out->ml[crb+k].val[2*(block_idx+l)+1] = v2;
              }
              buffer[crb+k]++;
            }
          }
#endif
          /*  destruct input matrix? */
          if (deleteIn) {
            free(in->blocks[i][j][k].idx);
            in->blocks[i][j][k].idx = NULL;

            free(in->blocks[i][j][k].val);
            in->blocks[i][j][k].val = NULL;
          }
        }
        if (deleteIn) {
          free(in->blocks[i][j]);
          in->blocks[i][j] = NULL;
        }
      }
      if (deleteIn) {
        free(in->blocks[i]);
        in->blocks[i] = NULL;
      }
    }
    for (i=0; i<rlin; ++i) {
#if NOT_DENSE_COPYING
      /*  realloc memory, only needed if the multilines are not dense copied */
      if (out->ml[i].sz > 0) {
        if (out->ml[i].sz < out->ncols) {
          out->ml[i].idx  = realloc(out->ml[i].idx,
              out->ml[i].sz * sizeof(mli_t));
        } else {
          free (out->ml[i].idx);
          out->ml[i].idx = NULL;
        }
        out->ml[i].val  = realloc(out->ml[i].val,
            2 * out->ml[i].sz * sizeof(re_t));
      }
#else
      if (buffer[i] == 0) {
        free (out->ml[i].val);
        out->ml[i].val  = NULL;
        out->ml[i].sz   = 0;
      }
#endif
    }
  }
  if (deleteIn) {
    free(in->blocks);
    in->blocks  = NULL;
    printf("addr in %p\n",in);
    free(in);
    in  = NULL;
    printf("addr in %p\n",in);
  }
#endif
  free(buffer);
  *input  = in;
  return out;
}

/**
 * \brief Reallocates memory for the rows of the blocks during the splicing of
 * the input matrix. The buffer size buffer_A is doubled during this process
 *
 * \param block matrix A
 *
 * \param row block index rbi in A
 *
 * \param block index in row bir
 *
 * \param block index in row bir
 *
 * \param line index in block lib
 *
 * \param buffer size buffer_A
 *
 */
void realloc_block_rows_B(sbm_fl_t *A, const ri_t rbi, const ci_t bir,
    const bi_t lib, const bi_t init_buffer_A, bi_t *buffer_A) {
  *buffer_A +=  init_buffer_A;
  A->blocks[rbi][bir][lib].idx = (bi_t*) realloc(
      A->blocks[rbi][bir][lib].idx,
      (*buffer_A) * sizeof(bi_t));
  A->blocks[rbi][bir][lib].val = (re_t*) realloc(
      A->blocks[rbi][bir][lib].val,
      2 * (*buffer_A) * sizeof(re_t));
}

sbm_fl_t *copy_multiline_to_block_matrix_rl(sm_fl_ml_t **A_in,
    ri_t bheight, ci_t bwidth, int free_memory, int nthrds) {

  sm_fl_ml_t *A = *A_in;
  sbm_fl_t *B = (sbm_fl_t *)malloc(sizeof(sbm_fl_t));

  ci_t i;
  ri_t j;
  ri_t k;

  const ri_t rlB            = (ri_t) ceil((float) A->nrows / (float)bheight);
  const ci_t clB            = (ci_t) ceil((float) A->ncols / (float)bwidth);
  const ri_t ml_rlB         = (A->nrows % __GBLA_NROWS_MULTILINE == 0) ?
    A->nrows / __GBLA_NROWS_MULTILINE :
    A->nrows / __GBLA_NROWS_MULTILINE + 1;
  const bi_t ml_bheight     = bheight / __GBLA_NROWS_MULTILINE;

  /*  initialize meta data for block submatrices */
  B->nrows    = A->nrows; /*  row dimension */
  B->ncols    = A->ncols; /*  col dimension */
  B->bheight  = bheight;  /*  block height */
  B->bwidth   = bwidth;   /*  block width */
  B->ba       = dtlr;     /*  block alignment */
  B->fe       = 1;        /*  fill empty blocks? */
  B->hr       = 1;        /*  allow hybrid rows? */
  B->nnz      = 0;        /*  number nonzero elements */
  /*  initialize B */
  B->blocks = (mbl_t ***)malloc(rlB * sizeof(mbl_t **));
  for (i=0; i<rlB; ++i) {
    B->blocks[i]  = (mbl_t **)malloc(clB * sizeof(mbl_t *));
    for (j=0; j<clB; ++j) {
      B->blocks[i][j] = (mbl_t *)malloc(
          ml_bheight * sizeof(mbl_t));
      for (k=0; k<ml_bheight; ++k) {
        B->blocks[i][j][k].val  = NULL;
        B->blocks[i][j][k].idx  = NULL;
        /* B->blocks[i][j][k].val  = (re_t *)malloc(2*bwidth*sizeof(re_t)); */
        /* B->blocks[i][j][k].idx  = (bi_t *)malloc(bwidth*sizeof(bi_t)); */
        B->blocks[i][j][k].sz   = B->blocks[i][j][k].dense  = 0;
      }
    }
  }
  const bi_t init_buffer_B  = (bi_t)(B->bwidth/2);
#pragma omp parallel private(j,k) num_threads(nthrds)
  {
    ri_t rbi;           /*  row block index */
    mli_t idx;          /*  multiline index of values */
    bi_t eil;           /*  element index in line */
    bi_t bir;           /*  block index in row */
    bi_t lib;           /*  line index in block */
    bi_t buffer_B[clB]; /*  buffer status for blocks in B */
    memset(buffer_B, 0, clB * sizeof(bi_t));

#pragma omp for nowait ordered
    for (i=0; i<ml_rlB; ++i) {
      memset(buffer_B, 0, clB * sizeof(bi_t));
      if (A->ml[i].sz == 0)
        continue;

      rbi = i / ml_bheight;
      lib = i % ml_bheight;

      for (j=A->ml[i].sz; j>0; --j) {
        idx = A->ml[i].idx[j-1];
        bir = (A->ncols - 1 - idx) / bwidth;
        eil = (A->ncols - 1 - idx) % bwidth;
        /*  realloc memory if needed */
        if (B->blocks[rbi][bir][lib].sz == buffer_B[bir]) {
          realloc_block_rows_B(B, rbi, bir, lib, init_buffer_B, &buffer_B[bir]);
        }
        /*  set values */
        B->blocks[rbi][bir][lib].idx[B->blocks[rbi][bir][lib].sz]   = eil;
        B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz]   =
          A->ml[i].val[2*(j-1)];
        B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz+1] =
          A->ml[i].val[2*(j-1)+1];
        B->blocks[rbi][bir][lib].sz++;
      }
      if (free_memory == 1) {
        free(A->ml[i].idx);
        /* A->ml[i].idx  = NULL; */
        free(A->ml[i].val);
        /* A->ml[i].val  = NULL; */
      }
    }
  }
  if (free_memory == 1) {
    free(A->ml);
    A->ml = NULL;
    free(A);
    A = NULL;
  }
  A_in  = &A;
  /*  hybrid multirows for the righthand side block matrices? */
  if (B->hr) {
#pragma omp parallel num_threads(nthrds)
    {
      ri_t idx;
      bi_t l;
#pragma omp for private(i,j,k,l) schedule(dynamic) nowait ordered
      /*  TODO: Implement hybrid stuff */
      for (i=0; i<rlB; ++i) {
        for (j=0; j<clB; ++j) {
          for (k=0; k<ml_bheight; ++k) {
            if (B->blocks[i][j][k].sz == 0) {
              free(B->blocks[i][j][k].idx);
              B->blocks[i][j][k].idx  = NULL;
              free(B->blocks[i][j][k].val);
              B->blocks[i][j][k].val  = NULL;
            }
            if ((float)B->blocks[i][j][k].sz / (float)B->bwidth
                < __GBLA_HYBRID_THRESHOLD) {
              B->blocks[i][j][k].idx =  (bi_t*) realloc(
                  B->blocks[i][j][k].idx,
                  B->blocks[i][j][k].sz * sizeof(bi_t));
              B->blocks[i][j][k].val =  (re_t*) realloc(
                  B->blocks[i][j][k].val,
                  2 * B->blocks[i][j][k].sz * sizeof(re_t));
              continue;
            }
            re_t *tmp_val_ptr = (re_t *)malloc(2 * B->bwidth * sizeof(re_t));
            idx  = 0;
            for (l=0; l<B->bwidth; ++l) {
              if (idx < B->blocks[i][j][k].sz && B->blocks[i][j][k].idx[idx] == l) {
                tmp_val_ptr[2*l]    = B->blocks[i][j][k].val[2*idx];
                tmp_val_ptr[2*l+1]  = B->blocks[i][j][k].val[2*idx+1];
                idx++;
              } else {
                tmp_val_ptr[2*l]    = 0;
                tmp_val_ptr[2*l+1]  = 0;
              }
            }
            free(B->blocks[i][j][k].idx);
            B->blocks[i][j][k].idx    = NULL;
            free(B->blocks[i][j][k].val);
            B->blocks[i][j][k].val    = tmp_val_ptr;
            B->blocks[i][j][k].sz     = B->bwidth;
            B->blocks[i][j][k].dense  = 1;
          }
        }
      }
    }
  } else { /*  cut down memory usage */
    /*  Realloc memory usage: */
    /*  Note that A is reallocated during the swapping of the data, so we */
    /*  reallocate only B here. */
#pragma omp parallel num_threads(nthrds)
    {
#pragma omp for schedule(dynamic) nowait ordered
      for (i=0; i<rlB; ++i) {
        int ctr;
        for (j=0; j<clB; ++j) {
          ctr = 0;
          for (k=0; k<ml_bheight; ++k) {
            if (B->blocks[i][j][k].sz > 0) {
              ctr = 1;
              B->blocks[i][j][k].idx = (bi_t*) realloc(
                  B->blocks[i][j][k].idx,
                  B->blocks[i][j][k].sz * sizeof(bi_t));
              B->blocks[i][j][k].val = (re_t*) realloc(
                  B->blocks[i][j][k].val,
                  2 * B->blocks[i][j][k].sz  * sizeof(re_t));
            }
          }
          /*  if full block is empty, remove it! */
          if (ctr == 0) {
            free(B->blocks[i][j]);
            B->blocks[i][j] = NULL;
          }
        }
      }
    }
  }
  return B;
}

double compute_density(nnz_t nnz, ri_t nrows, ri_t ncols) {
	return (double) (nnz * 100) / (double)(nrows * (nnz_t)ncols);
}

/* vim:sts=2:sw=2:ts=2:
 */
