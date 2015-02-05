#include <sys/types.h>                 // For ssize_t
#include "asl/asl_pfgh.h"              // Ampl library headers

// ==========================================================================

//
//        P r o t o t y p e s   f o r   m o d u l e   f u n c t i o n s

// ==========================================================================

void *jampl_init(char *stub);
void jampl_finalize(void *asl);
void jampl_write_sol(void *asl, const char *msg, double *x, double *y);

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
void    jampl_grad(    void *asl, double *x, double *g);
void    jampl_cons(    void *asl, double *x, double *c);
double  jampl_jcon(    void *asl, double *x, int j);
void    jampl_jcongrad(void *asl, double *x, double *g, int j);
void    jampl_hprod(   void *asl, double *y, double *v, double *hv, double w);
void    jampl_hvcompd( void *asl, double *v, double *hv, int nobj);
void    jampl_ghjvprod(void *asl, double *g, double *v, double *ghjv);

size_t jampl_sparse_congrad_nnz(void *asl, int j);
void jampl_sparse_congrad(void *asl, double *x, int j, int64_t *inds, double *vals);
void jampl_jac( void *asl, double *x, int64_t *rows, int64_t *cols, double *vals);
void jampl_hess(void *asl, double *y, double w, int64_t *rows, int64_t *cols, double *vals);
