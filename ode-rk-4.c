#include "fec.h"

static double *k1, *k2, *k3, *k4;
static double *tempx;
static int first = 1;

void ode_rk_4_solve(double *t, double *x, double h,
                    void derivs(double t, double *x, double *diffx, param_list_t *param),
                    param_list_t *param) {
  int i;

  if (first) {
    first = 0;
    k1 = (double*)malloc(sizeof(double)*length);
    k2 = (double*)malloc(sizeof(double)*length);
    k3 = (double*)malloc(sizeof(double)*length);
    k4 = (double*)malloc(sizeof(double)*length);
    tempx = (double*)malloc(sizeof(double)*length);
  }
  derivs(*t,x,k1,param);
  for (i=0;i<length;i++) tempx[i] = x[i] + h*k1[i]/2;
  derivs(*t+h/2,tempx,k2,param);
  for (i=0;i<length;i++) tempx[i] = x[i] + h*k2[i]/2;
  derivs(*t+h/2,tempx,k3,param);
  for (i=0;i<length;i++) tempx[i] = x[i] + h*k3[i];
  derivs(*t+h,tempx,k4,param);
  for (i=0;i<length;i++) x[i] += h/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
  *t += h;
}
