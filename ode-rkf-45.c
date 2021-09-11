#include "fec.h"

/* Taken from "Numerical Recipes in C" by Press, Teukolsky, Vetterling, and Flannery. */

static double **k;
static double *tempx;
static int first = 1;

static double
a[6] = {0,1./5,3./10,3./5,1,7./8},
b[6][6] = {{0, 0, 0, 0, 0, 0},
{1./5, 0, 0, 0, 0},
{3./40, 9./40, 0, 0, 0},
{3./10, -9./10, 6./5, 0, 0},
{-11./54, 5./2, -70./27, 35./27, 0},
{1631./55296, 175./512, 575./13824, 44275./110592, 253./4096}},
c5[6] = {37./378, 0, 250./621, 125./594, 0, 512./1771},
c4[6] = {2825./27648, 0, 18575./48384, 13525./55296, 277./14336, 1./4};

void ode_rkf_45_solve(double *t, double *x, double *h_use,
                      void derivs(double t, double *x, double *diffx, param_list_t *param),
                      param_list_t *param) {
  int i,j,m;
  double error, h, delta, temp1;

  if (first) {
    first = 0;
    k = (double**)malloc(sizeof(double*)*6);
    for (i=0;i<6;i++)
      k[i] = (double*)malloc(sizeof(double)*length);
    tempx = (double*)malloc(sizeof(double)*length);
  }

  h = *h_use;
  while (1) {
    for (i=0;i<6;i++) {
      memcpy(tempx,x,sizeof(double)*length);
      for (j=0;j<i;j++)
        for (m=0;m<length;m++)
          tempx[m] += h*b[i][j]*k[j][m];
      derivs(*t+a[i]*h,tempx,k[i],param);
    }
/* l_infinity error between RK4 and RK5 */
    error = 0;
    for (m=0;m<length;m++) {
      temp1 = 0;
      for (i=0;i<6;i++)
        temp1 += (c5[i]-c4[i])*k[i][m];
      if (fabs(temp1)>error) error = fabs(temp1);
    }
    error = h*error/param->tol;
    if (error<=1) break;
    delta = 0.9*pow(error,-0.25);
    if (delta>0.1)
      h = delta*h;
    else
      h = 0.1*h;
  }

/* RK order 5 */
  for (i=0;i<6;i++)
    for (m=0;m<length;m++)
      x[m] += h*c5[i]*k[i][m];
  *t += h;

  if (error>1.89e-4)
    *h_use = 0.9*h*pow(error,-0.2);
  else
    *h_use = 5*h;
}
