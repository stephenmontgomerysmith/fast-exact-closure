#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <ctype.h>

/* The position of the components of A and B in X. */
#define A11 0
#define A12 1
#define A13 2
#define A22 3
#define A23 4
#define B11 5
#define B12 6
#define B13 7
#define B22 8
#define B23 9
#define length 10

typedef struct {
  int verbose;
  int verbose_print;
  const char *outfilename;
  double h;
  double tol;
  int print_every;
  double tstart;
  double tend;
  double w[3];
  double gamm[9];
  double lambda;
  double CI;
  int ode_rk_4;
  int do_ard;
  double b1, b2, b3, b4, b5;
  int do_rsc;
  double kappa;
  int do_reset;
  int debug;
} param_list_t;

void derivs(double t, double *X, double *Xdot, param_list_t *param);
void reset(double *X);
void ode_rk_4_solve(double *t, double *x, double h,
                    void derivs(double t, double *x, double *diffx, param_list_t *param),
                    param_list_t *param);
void ode_rkf_45_solve(double *t, double *x, double *h_use,
                      void derivs(double t, double *x, double *diffx, param_list_t *param),
                      param_list_t *param);
void diagonalize_sym(int n, double *A, double *evec);
void diagonalize_sym_3(double *A, double *evec);
void rotate(double *A, double *evec);
void rotate_diag(double *Adiag, double *A, double *evec);
void reverse_rotate(double *A, double *evec);
void inverse_sym(double *B, double *inv);
void get_parameters(int argc, const char **argv, param_list_t *param);
int param_bool(const char *p);
int param_int(const char *p);
double param_double(const char *p);
int param_choice(const char *p, ...);
char *param_string(const char *p);
double evaluate_string(char *s, char var, double x);
void set_param_filename(const char *f);
void set_param_verbose_level(int v);
void done_with_param();
void param_ignore(const char *p);
int check_param(const char *p);
