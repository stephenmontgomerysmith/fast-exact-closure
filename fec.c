#include "fec.h"

int no_close, two_close[3], three_close;

int main(int argc, const char **argv) {
  FILE *sout;
  double *X;
  double t;
  int iteration;
  param_list_t param;

  get_parameters(argc, argv, &param);

  no_close = 0;
  memset(two_close,0,sizeof(two_close));
  three_close = 0;

  X = (double*)malloc(sizeof(double)*length);
/* Start with isotropic distribution. */
  memset(X,0,sizeof(double)*length);
  X[A11] = X[A22] = 1./3;
  X[B11] = X[B22] = 1;

  t = param.tstart;
  iteration=0;

  sout = fopen(param.outfilename,"w");
  if (sout==NULL) {
    perror("unable to open output file");
    exit(1);
  }

  while (1) {
    if (iteration%param.print_every==0) {
      fprintf(sout,"%g",t);
      fprintf(sout," %g",X[A11]);
      fprintf(sout," %g",X[A12]);
      fprintf(sout," %g",X[A13]);
      fprintf(sout," %g",X[A22]);
      fprintf(sout," %g",X[A23]);
      fprintf(sout," %g",1-X[A11]-X[A22]);
      fprintf(sout,"\n");
      fflush(sout);
      if (param.verbose_print!=0 && (iteration/param.print_every)%param.verbose_print==0) {
        printf("%g",t);
        printf(" %g",X[A11]);
        printf(" %g",X[A12]);
        printf(" %g",X[A13]);
        printf(" %g",X[A22]);
        printf(" %g",X[A23]);
        printf(" %g",1-X[A11]-X[A22]);
        printf("\n");
      }
    }

    if (t>=param.tend) break;
    if (param.do_reset)
      reset(X);
    if (param.ode_rk_4)
      ode_rk_4_solve(&t,X,param.h,derivs,&param);
    else
      ode_rkf_45_solve(&t,X,&(param.h),derivs,&param);

    if (t+param.h>param.tend) param.h = (1+1e-7)*(param.tend-t);
    iteration++;
  }

  if (param.verbose)
    printf("Total number of iterations = %d\n",iteration);
  if (param.debug) {
    printf("Three close %d\n",three_close);
    printf("Two close 0 1 %d\n",two_close[2]);
    printf("Two close 0 2 %d\n",two_close[1]);
    printf("Two close 1 2 %d\n",two_close[0]);
    printf("No close %d\n",no_close);
  }

  exit(0);
}
