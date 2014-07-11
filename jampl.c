// ==========================================
// Generic Ampl interface for use with Julia.
// Dominique Orban
// Vancouver, April 2014.
// ==========================================

#include "jampl.h"

// Module functions.

void *jampl_init(char *stub) {
    ASL_pfgh *asl = (ASL_pfgh*)ASL_alloc(ASL_read_pfgh);
    if (!asl) return NULL;

    FILE *ampl_file = jac0dim(stub, (fint)strlen(stub));

    // Allocate room to store problem data
    if (! (asl->i.X0_    = (real *)M1alloc(asl->i.n_var_ * sizeof(real)))) return NULL;
    if (! (asl->i.LUv_   = (real *)M1alloc(asl->i.n_var_ * sizeof(real)))) return NULL;
    if (! (asl->i.Uvx_   = (real *)M1alloc(asl->i.n_var_ * sizeof(real)))) return NULL;
    if (! (asl->i.pi0_   = (real *)M1alloc(asl->i.n_con_ * sizeof(real)))) return NULL;
    if (! (asl->i.LUrhs_ = (real *)M1alloc(asl->i.n_con_ * sizeof(real)))) return NULL;
    if (! (asl->i.Urhsx_ = (real* )M1alloc(asl->i.n_con_ * sizeof(real)))) return NULL;

    // Read in ASL structure
    asl->i.want_xpi0_ = 3;        // Read primal and dual estimates
    pfgh_read(ampl_file , 0);     // pfgh_read closes the file.

    return (void *)asl;
}

void jampl_finalize(void *asl) {
    ASL *this_asl = (ASL *)asl;
    ASL_free((ASL **)(&this_asl));
    return;
}

// Problem setup.

int jampl_objtype(void *asl) {
  return ((ASL *)asl)->i.objtype_[0];  // 0 means minimization problem.
}

int jampl_nvar(void *asl) {
  return ((ASL *)asl)->i.n_var_;
}

int jampl_ncon(void *asl) {
  return ((ASL *)asl)->i.n_con_;
}

int jampl_nlc(void *asl) {
  return ((ASL *)asl)->i.nlc_;
}

int jampl_nlnc(void *asl) {
  return ((ASL *)asl)->i.nlnc_;
}

int jampl_nnzj(void *asl) {
  return ((ASL *)asl)->i.nzc_;
}

int jampl_nnzh(void *asl) {
  ASL *this_asl = (ASL *)asl;
  return (int)((*(this_asl->p.Sphset))(this_asl, 0, -1, 1, 1, 1));
}

int jampl_islp(void *asl) {
  ASL *this_asl = (ASL *)asl;
  return ((this_asl->i.nlo_ + this_asl->i.nlc_ + this_asl->i.nlnc_) > 0 ? 0 : 1);
}

double *jampl_x0(void *asl) {
  return ((ASL *)asl)->i.X0_;
}

double *jampl_y0(void *asl) {
  return ((ASL *)asl)->i.pi0_;
}

double *jampl_lvar(void *asl) {
  return ((ASL *)asl)->i.LUv_;
}

double *jampl_uvar(void *asl) {
  return ((ASL *)asl)->i.Uvx_;
}

double *jampl_lcon(void *asl) {
  return ((ASL *)asl)->i.LUrhs_;
}

double *jampl_ucon(void *asl) {
  return ((ASL *)asl)->i.Urhsx_;
}

// Objective.

void jampl_varscale(void *asl, double *s) {
  fint ne;
  int this_nvar = ((ASL *)asl)->i.n_var_;

  for (int i = 0; i < this_nvar; i++)
    varscale_ASL((ASL *)asl, i, s[i], &ne);
  return;
}

double jampl_obj(void *asl, double *x) {
  fint ne;
  return (*((ASL *)asl)->p.Objval)((ASL *)asl, 0, x, &ne);
}

double *jampl_grad(void *asl, double *x) {
  fint ne;
  int this_nvar = ((ASL *)asl)->i.n_var_;
  double *g = (double *)Malloc(this_nvar * sizeof(real));

  (*((ASL *)asl)->p.Objgrd)((ASL *)asl, 0, x, g, &ne);
  return g;
}

// Lagrangian.

void jampl_lagscale(void *asl, double s) {
  fint ne;
  lagscale_ASL((ASL *)asl, s, &ne);
  return;
}

// Constraints and Jacobian.

void jampl_conscale(void *asl, double *s) {
  fint ne;
  int this_ncon = ((ASL *)asl)->i.n_con_;

  for (int j = 0; j < this_ncon; j++)
    conscale_ASL((ASL *)asl, j, s[j], &ne);
  return;
}

double *jampl_cons(void *asl, double *x) {
  fint ne;
  int this_ncon = ((ASL *)asl)->i.n_con_;
  double *c = (double *)Malloc(this_ncon * sizeof(real));

  (*((ASL *)asl)->p.Conval)((ASL *)asl, x, c, &ne);
  return c;
}

double jampl_jcon(void *asl, double *x, int j) {
  fint ne;
  return (*((ASL *)asl)->p.Conival)((ASL *)asl, j, x, &ne);
}

double *jampl_jcongrad(void *asl, double *x, int j) {
  ASL *this_asl = (ASL *)asl;
  int this_nvar = this_asl->i.n_var_;
  double *g = (double *)Malloc(this_nvar * sizeof(real));

  fint ne;
  (*(this_asl->p.Congrd))(this_asl, j, x, g, &ne);
  return g;
}

jl_tuple_t *jampl_sparse_congrad(void *asl, double *x, int j) {
  ASL *this_asl = (ASL *)asl;
  size_t nzgj = 0;
  cgrad *cg;
  for (cg = this_asl->i.Cgrad_[j]; cg; cg = cg->next) nzgj++;

  // Declare double and long Julia array types.
  jl_value_t *float64_array_type = jl_apply_array_type(jl_float64_type, 1);
  jl_value_t *int64_array_type   = jl_apply_array_type(jl_int64_type,   1);

  jl_array_t *jvals = NULL, *jinds = NULL;
  JL_GC_PUSH2(&jvals, &jinds);  // Let Julia worry about these chunks.

  jvals = jl_alloc_array_1d(float64_array_type, nzgj);
  jinds = jl_alloc_array_1d(int64_array_type,   nzgj);

  long   *inds = (long   *)jl_array_data(jinds);
  double *vals = (double *)jl_array_data(jvals);

  int congrd_mode_bkup = this_asl->i.congrd_mode;
  this_asl->i.congrd_mode = 1;  // Sparse gradient mode.

  fint ne;
  (*(this_asl->p.Congrd))(this_asl, j, x, vals, &ne);

  int k = 0;
  for (cg = this_asl->i.Cgrad_[j]; cg; cg = cg->next)
      inds[k++] = (long)(cg->varno) + 1;

  this_asl->i.congrd_mode = congrd_mode_bkup;  // Restore gradient mode.

  jl_tuple_t *tuple = jl_alloc_tuple(2);
  jl_tupleset(tuple, 0, jinds);
  jl_tupleset(tuple, 1, jvals);

  JL_GC_POP();
  return tuple;
}

// Return Jacobian at x in triplet form (rows, vals, cols).
jl_tuple_t *jampl_jac(void *asl, double *x) {
  ASL *this_asl = (ASL *)asl;
  int this_nzc = this_asl->i.nzc_, this_ncon = this_asl->i.n_con_;

  // Declare double and long Julia array types.
  jl_value_t *float64_array_type = jl_apply_array_type(jl_float64_type, 1);
  jl_value_t *int64_array_type   = jl_apply_array_type(jl_int64_type,   1);

  jl_array_t *jvals = NULL, *jrows = NULL, *jcols = NULL;
  JL_GC_PUSH3(&jvals, &jrows, &jcols);  // Let Julia worry about these chunks.

  jvals = jl_alloc_array_1d(float64_array_type, (size_t)(this_nzc));
  jrows = jl_alloc_array_1d(int64_array_type,   (size_t)(this_nzc));
  jcols = jl_alloc_array_1d(int64_array_type,   (size_t)(this_nzc));

  long   *rows = (long   *)jl_array_data(jrows);
  long   *cols = (long   *)jl_array_data(jcols);
  double *vals = (double *)jl_array_data(jvals);

  fint ne;
  (*(this_asl->p.Jacval))(this_asl, x, vals, &ne);

  // Fill in sparsity pattern. Account for 1-based indexing.
  for (int j = 0; j < this_ncon; j++)
    for (cgrad *cg = this_asl->i.Cgrad_[j]; cg; cg = cg->next) {
        rows[cg->goff] = (long)j + 1;
        cols[cg->goff] = (long)(cg->varno) + 1;
    }

#ifdef DEBUG_AMPL_JL
  printf("AMPL says nzc = %d\n", this_nzc);
  printf("array lengths are %d\n", jl_array_len(jrows));
  printf("Error code = %d\n", ne);
  size_t l;
  int nelts = this_nzc/2 < 4 ? this_nzc/2 : 4;
  printf("J: vals = [ ");
  for (l = 0; l < nelts; l++) printf("%8.1e ", vals[l]);
  printf("... ");
  for (l = 0; l < nelts; l++) printf("%8.1e ", vals[this_nzc - nelts + l]);
  printf("]\n");
  printf("J: rows = [ ");
  for (l = 0; l < nelts; l++) printf("%d ", rows[l]);
  printf("... ");
  for (l = 0; l < nelts; l++) printf("%d ", rows[this_nzc - nelts + l]);
  printf("]\n");
  printf("J: cols = [ ");
  for (l = 0; l < nelts; l++) printf("%d ", cols[l]);
  printf("... ");
  for (l = 0; l < nelts; l++) printf("%d ", cols[this_nzc - nelts + l]);
  printf("]\n");
#endif

  jl_tuple_t *tuple = jl_alloc_tuple(3);
  jl_tupleset(tuple, 0, jrows);
  jl_tupleset(tuple, 1, jcols);
  jl_tupleset(tuple, 2, jvals);

  JL_GC_POP();
  return tuple;
}

// Hessian.

double *jampl_hprod(void *asl, double *y, double *v, double w) {
  ASL *this_asl = (ASL *)asl;
  int this_nvar = this_asl->i.n_var_;
  double *hv = (double *)Malloc(this_nvar * sizeof(real));
  double ow[1];  // Objective weight.

  ow[0]  = this_asl->i.objtype_[0] ? -w : w;
  hvpinit_ASL(this_asl, this_asl->p.ihd_limit_, 0, NULL, y);
  (*(this_asl->p.Hvcomp))(this_asl, hv, v, -1, ow, y); // nobj=-1 so ow takes precendence.
  return hv;
}

double *jampl_hvcompd(void *asl, double *v, int nobj) {
  ASL *this_asl = (ASL *)asl;
  int this_nvar = this_asl->i.n_var_;
  double *hv = (double *)Malloc(this_nvar * sizeof(real));
  (*(this_asl->p.Hvcompd))(this_asl, hv, v, nobj);
  return hv;
}

double *jampl_ghjvprod(void *asl, double *g, double *v) {
  ASL *this_asl = (ASL *)asl;
  int this_ncon = this_asl->i.n_con_;
  int this_nvar = this_asl->i.n_var_;
  int this_nlc  = this_asl->i.nlc_;
  double *hv    = (double *)Malloc(this_nvar * sizeof(real));
  double *ghjv  = (double *)Malloc(this_ncon * sizeof(real));

  fint ne;
  double prod;
  int i, j;

  // Process nonlinear constraints.
  for (j = 0 ; j < this_nlc ; j++) {
    (*(this_asl->p.Hvcompd))(this_asl, hv, v, j);

    // Compute dot product g'Hi*v. Should use BLAS.
    for (i = 0, prod = 0 ; i < this_nvar ; i++)
      prod += (hv[i] * g[i]);
    ghjv[j] = prod;
  }
  free(hv);

  // All terms corresponding to linear constraints are zero.
  for (j = this_nlc ; j < this_ncon ; j++) ghjv[j] = 0.;

  return ghjv;
}

// Return Hessian at (x,y) in triplet form (rows, vals, cols).
jl_tuple_t *jampl_hess(void *asl, double *y, double w) {
  ASL *this_asl = (ASL *)asl;
  double ow[1];  // Objective weight.
  size_t nnzh = (size_t)((*(this_asl->p.Sphset))(this_asl, 0, -1, 1, 1, 1)); // nobj=-1 so ow takes precendence.

  // Declare double and long Julia array types.
  jl_value_t *float64_array_type = jl_apply_array_type(jl_float64_type, 1);
  jl_value_t *int64_array_type   = jl_apply_array_type(jl_int64_type,   1);

  jl_array_t *jvals = NULL, *jrows = NULL, *jcols = NULL;
  JL_GC_PUSH3(&jvals, &jrows, &jcols);  // Let Julia worry about these chunks.

  jvals = jl_alloc_array_1d(float64_array_type, nnzh);
  jrows = jl_alloc_array_1d(int64_array_type,   nnzh);
  jcols = jl_alloc_array_1d(int64_array_type,   nnzh);

  long   *rows = (long   *)jl_array_data(jrows);
  long   *cols = (long   *)jl_array_data(jcols);
  double *vals = (double *)jl_array_data(jvals);

  ow[0]  = this_asl->i.objtype_[0] ? -w : w;
  (*(this_asl->p.Sphes))(this_asl, 0, vals, -1, ow, y);

  // Fill in sparsity pattern. Account for 1-based indexing.
  int k = 0;
  for (int i = 0; i < this_asl->i.n_var_; i++)
    for (int j = this_asl->i.sputinfo_->hcolstarts[i]; j < this_asl->i.sputinfo_->hcolstarts[i+1]; j++) {
      rows[k] = this_asl->i.sputinfo_->hrownos[j] + 1;
      cols[k] = i + 1;
      k++;
    }

#ifdef DEBUG_AMPL_JL
  size_t l;
  printf("vals = [ ");
  for (l = 0; l < jl_array_len(jvals); l++) printf("%8.1e ", vals[l]);
  printf("]\n");
  printf("rows = [ ");
  for (l = 0; l < jl_array_len(jrows); l++) printf("%d ", rows[l]);
  printf("]\n");
  printf("cols = [ ");
  for (l = 0; l < jl_array_len(jcols); l++) printf("%d ", cols[l]);
  printf("]\n");
#endif

  jl_tuple_t *tuple = jl_alloc_tuple(3);
  jl_tupleset(tuple, 0, jrows);
  jl_tupleset(tuple, 1, jcols);
  jl_tupleset(tuple, 2, jvals);

  JL_GC_POP();
  return tuple;
}
