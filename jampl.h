#include <julia/julia.h>
#include "asl/asl_pfgh.h"              // Ampl library headers
#include "asl/jacpdim.h"               // For partially-separable structure

// ==========================================================================

//
//        P r o t o t y p e s   f o r   m o d u l e   f u n c t i o n s

// ==========================================================================

void *jampl_init(char *stub);
void jampl_finalize(void *asl);

int jampl_objtype(void *asl);
int jampl_nvar(   void *asl);
int jampl_ncon(   void *asl);
int jampl_nlc(    void *asl);
int jampl_nlnc(   void *asl);
int jampl_nnzj(   void *asl);
int jampl_nnzh(   void *asl);
int jampl_islp(   void *asl);

double *jampl_x0(  void *asl);
double *jampl_y0(  void *asl);
double *jampl_lvar(void *asl);
double *jampl_uvar(void *asl);
double *jampl_lcon(void *asl);
double *jampl_ucon(void *asl);

void jampl_varscale(void *asl, double *s);
void jampl_lagscale(void *asl, double  s);
void jampl_conscale(void *asl, double *s);

double  jampl_obj(     void *asl, double *x);
double *jampl_grad(    void *asl, double *x);
double *jampl_cons(    void *asl, double *x);
double  jampl_jcon(    void *asl, double *x, int j);
double *jampl_jcongrad(void *asl, double *x, int j);
double *jampl_hprod(   void *asl, double *y, double *v, double w);
double *jampl_hvcompd( void *asl, double *v, int nobj);
double *jampl_ghjvprod(void *asl, double *g, double *v);

// Functions that return a Julia-specific value.
jl_tuple_t *jampl_sparse_congrad(void *asl, double *x, int j);
jl_tuple_t *jampl_jac( void *asl, double *x);
jl_tuple_t *jampl_hess(void *asl, double *y, double w);
