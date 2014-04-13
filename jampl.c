/* ==========================================
 * Generic Ampl interface for use with Julia.
 * Dominique Orban
 * Vancouver, April 2014.
 * ==========================================
 */

#include <julia/julia.h>
#include "asl/asl_pfgh.h"              /* Ampl library headers    */
#include "asl/jacpdim.h"               /* For partially-separable structure */

/* ========================================================================== */

/*
 *        P r o t o t y p e s   f o r   m o d u l e   f u n c t i o n s
 */

/* ========================================================================== */

/* Ampl driver specific declarations */

static ASL_pfgh *asl;  /* Main ASL structure */

static FILE *ampl_file  =  NULL;   /* Connection with Ampl nl file */
static int   ampl_file_open = 0;   /* Number of open files counter */
static int   ampl_written_sol = 0; /* Indicates whether solution was written */

/* ========================================================================== */

/*
 *                    M o d u l e   f u n c t i o n s
 */

/* ========================================================================== */

int jampl_init(char *stub) {
    /* Fast return if a file is already open */
    if (ampl_file_open) return 0;

    /* Initialize main ASL structure */
    if (! (asl  = (ASL_pfgh*)ASL_alloc(ASL_read_pfgh))) return -1;

    ampl_file = jac0dim(stub, (fint)strlen(stub));

    /* Allocate room to store problem data */
    if (! (X0    = (real *)M1alloc(n_var*sizeof(real)))) return -2;
    if (! (LUv   = (real *)M1alloc(n_var*sizeof(real)))) return -2;
    if (! (Uvx   = (real *)M1alloc(n_var*sizeof(real)))) return -2;
    if (! (pi0   = (real *)M1alloc(n_con*sizeof(real)))) return -2;
    if (! (LUrhs = (real *)M1alloc(n_con*sizeof(real)))) return -2;
    if (! (Urhsx = (real *)M1alloc(n_con*sizeof(real)))) return -2;

    /* Read in ASL structure */
    want_xpi0 = 3;           /* Read primal and dual estimates */
    pfgh_read(ampl_file , 0);

    ampl_file_open = 1;
    return 0;
}

void jampl_finalize(void) {
    ASL_free((ASL **)(&asl));
    ampl_file_open  = 0;
    return;
}

// Problem setup.

int jampl_objtype(void) {
  return objtype[0];  /* 0 means minimization problem. */
}

int jampl_nvar(void) {
  return asl->i.n_var_;
}

int jampl_ncon(void) {
  return asl->i.n_con_;
}

int jampl_nlc(void) {
  return asl->i.nlc_;
}

int jampl_nlnc(void) {
  return asl->i.nlnc_;
}

int jampl_nnzj(void) {
  return asl->i.nzc_;
}

int jampl_nnzh(void) {
  return (int)sphsetup(-1, 1, 1, 1);
}

int jampl_islp(void) {
  return ((nlo + nlc + nlnc) > 0 ? 0 : 1);
}

double *jampl_x0(void) {
  return asl->i.X0_;
}

double *jampl_y0(void) {
  return asl->i.pi0_;
}

double *jampl_lvar(void) {
  return asl->i.LUv_;
}

double *jampl_uvar(void) {
  return asl->i.Uvx_;
}

double *jampl_lcon(void) {
  return asl->i.LUrhs_;
}

double *jampl_ucon(void) {
  return asl->i.Urhsx_;
}

// Objective.

void jampl_varscale(double *s) {
  fint ne;
  for (int i = 0; i < n_var; i++)
    conscale_ASL((ASL*)asl, i, s[i], &ne);
  return;
}

double jampl_obj(double *x) {
  fint ne;  // Error code. Currently ignored.
  return (*((ASL*)asl)->p.Objval)((ASL*)asl, 0, x, &ne);
}

double *jampl_grad(double *x) {
  fint ne;
  double *g;

  g = (double *)Malloc(n_var * sizeof(real));
  (*((ASL*)asl)->p.Objgrd)((ASL*)asl, 0, x, g, &ne);
  return g;
}

// Lagrangian.

void jampl_lagscale(double s) {
  fint ne;
  lagscale_ASL((ASL*)asl, s, &ne);
  return;
}

// Constraints and Jacobian.

void jampl_conscale(double *s) {
  fint ne;
  for (int i = 0; i < n_con; i++)
    conscale_ASL((ASL*)asl, i, s[i], &ne);
  return;
}

double *jampl_cons(double *x) {
  fint ne;
  double *c;

  c = (double *)Malloc(n_con * sizeof(real));
  (*((ASL*)asl)->p.Conval)((ASL*)asl, x, c, &ne);
  return c;
}

double jampl_jcon(double *x, int j) {
  fint ne;
  return (*((ASL*)asl)->p.Conival)((ASL*)asl, j, x, &ne);
}

double *jampl_jcongrad(double *x, int j) {
  fint ne;
  double *g;

  g = (double *)Malloc(n_var * sizeof(real));
  (*((ASL*)asl)->p.Congrd)((ASL*)asl, j, x, g, &ne);
  return g;
}

jl_tuple_t *jampl_sparse_congrad(double *x, int j) {
  size_t nzgj = 0;
  cgrad *cg;
  for (cg = asl->i.Cgrad_[j]; cg; cg = cg->next) nzgj++;

  // Declare double and long Julia array types.
  jl_value_t *float64_array_type = jl_apply_array_type(jl_float64_type, 1);
  jl_value_t *int64_array_type   = jl_apply_array_type(jl_int64_type,   1);

  jl_array_t *jvals = jl_alloc_array_1d(float64_array_type, nzgj);
  jl_array_t *jinds = jl_alloc_array_1d(int64_array_type,   nzgj);
  JL_GC_PUSH2(&jvals, &jinds);  // Let Julia worry about these chunks.

  long   *inds = (long   *)jl_array_data(jinds);
  double *vals = (double *)jl_array_data(jvals);

  int congrd_mode_bkup = asl->i.congrd_mode;
  asl->i.congrd_mode = 1;  // Sparse gradient mode.

  fint ne;
  congrd(j, x, vals, &ne);

  int k = 0;
  for (cg = Cgrad[j]; cg; cg = cg->next)
      inds[k++] = (long)(cg->varno) + 1;

  asl->i.congrd_mode = congrd_mode_bkup;  // Restore gradient mode.

  jl_tuple_t *tuple = jl_alloc_tuple(2);
  jl_tupleset(tuple, 0, jinds);
  jl_tupleset(tuple, 1, jvals);

  return tuple;
}

// Return Jacobian at x in triplet form (rows, vals, cols).
jl_tuple_t *jampl_jac(double *x) {
  // Declare double and long Julia array types.
  jl_value_t *float64_array_type = jl_apply_array_type(jl_float64_type, 1);
  jl_value_t *int64_array_type   = jl_apply_array_type(jl_int64_type,   1);

  jl_array_t *jvals = jl_alloc_array_1d(float64_array_type, (size_t)(asl->i.nzc_));
  jl_array_t *jrows = jl_alloc_array_1d(int64_array_type,   (size_t)(asl->i.nzc_));
  jl_array_t *jcols = jl_alloc_array_1d(int64_array_type,   (size_t)(asl->i.nzc_));
  JL_GC_PUSH3(&jvals, &jrows, &jcols);  // Let Julia worry about these chunks.

  long   *rows = (long   *)jl_array_data(jrows);
  long   *cols = (long   *)jl_array_data(jcols);
  double *vals = (double *)jl_array_data(jvals);

  fint ne;
  jacval(x, vals, &ne);

  // Fill in sparsity pattern. Account for 1-based indexing.
  for (int i = 0; i < n_con; i++)
    for (cgrad *cg = Cgrad[i]; cg; cg = cg->next) {
        rows[cg->goff] = (long)i + 1;
        cols[cg->goff] = (long)(cg->varno) + 1;
    }

#ifdef DEBUG_AMPL_JL
  size_t l;
  printf("J: vals = [ ");
  for (l = 0; l < jl_array_len(jvals); l++) printf("%8.1e ", vals[l]);
  printf("]\n");
  printf("J: rows = [ ");
  for (l = 0; l < jl_array_len(jrows); l++) printf("%d ", rows[l]);
  printf("]\n");
  printf("J: cols = [ ");
  for (l = 0; l < jl_array_len(jcols); l++) printf("%d ", cols[l]);
  printf("]\n");
#endif

  jl_tuple_t *tuple = jl_alloc_tuple(3);
  jl_tupleset(tuple, 0, jrows);
  jl_tupleset(tuple, 1, jcols);
  jl_tupleset(tuple, 2, jvals);

  return tuple;
}

// Hessian.

double *jampl_hprod(double *x, double *y, double *v, double w) {
  double *hv;
  double ow[1];  // Objective weight.

  ow[0]  = objtype[0] ? -w : w;
  hv = (double *)Malloc(n_var * sizeof(real));
  hvpinit_ASL((ASL*)asl, ihd_limit, 0, NULL, y);
  (*((ASL*)asl)->p.Hvcomp)((ASL*)asl, hv, v, -1, ow, y); // nobj=-1 so ow takes precendence.
  return hv;
}

double *jampl_ghjvprod(double *x, double *g, double *v) {
  fint ne;
  double *y, *hv, *ghjv, prod;
  int i, j;

  y = (double *)Malloc(n_con * sizeof(real));
  hv = (double *)Malloc(n_var * sizeof(real));
  ghjv = (double *)Malloc(n_con * sizeof(real));

  for (j = 0 ; j < n_con ; j++) y[j] = 0.;

  // Process nonlinear constraints.
  for (j = 0 ; j < nlc ; j++) {
    // Set vector of multipliers to (0, 0, ..., -1, ..., 0).
    y[j] = -1.;

    // Compute hv = Hj * v by setting OW to NULL.
    hvpinit_ASL((ASL*)asl, ihd_limit, 0, NULL, y);
    hvcomp(hv, v, 0, NULL, y);

    // Compute dot product g'Hi*v. Should use BLAS.
    for (i = 0, prod = 0 ; i < n_var ; i++)
      prod += (hv[i] * g[i]);
    ghjv[j] = prod;

    // Reset j-th multiplier.
    y[j] = 0.;
  }
  free(y);
  free(hv);

  // All terms corresponding to linear constraints are zero.
  for (j = nlc ; j < n_con ; j++) ghjv[j] = 0.;

  return ghjv;
}

// Return Hessian at (x,y) in triplet form (rows, vals, cols).
jl_tuple_t *jampl_hess(double *x, double *y, double w) {
  double ow[1];  // Objective weight.
  ow[0]  = objtype[0] ? -w : w;
  size_t nnzh = (size_t)sphsetup(-1, 1, 1, 1); // nobj=-1 so ow takes precendence.

  // Declare double and long Julia array types.
  jl_value_t *float64_array_type = jl_apply_array_type(jl_float64_type, 1);
  jl_value_t *int64_array_type   = jl_apply_array_type(jl_int64_type,   1);

  jl_array_t *jvals = jl_alloc_array_1d(float64_array_type, nnzh);
  jl_array_t *jrows = jl_alloc_array_1d(int64_array_type,   nnzh);
  jl_array_t *jcols = jl_alloc_array_1d(int64_array_type,   nnzh);
  JL_GC_PUSH3(&jvals, &jrows, &jcols);  // Let Julia worry about these chunks.

  long   *rows = (long   *)jl_array_data(jrows);
  long   *cols = (long   *)jl_array_data(jcols);
  double *vals = (double *)jl_array_data(jvals);

  sphes(vals, -1, ow, y);

  // Fill in sparsity pattern. Account for 1-based indexing.
  int k = 0;
  for (int i = 0; i < n_var; i++)
    for (int j = sputinfo->hcolstarts[i]; j < sputinfo->hcolstarts[i+1]; j++) {
      rows[k] = sputinfo->hrownos[j] + 1;
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

  return tuple;
}
