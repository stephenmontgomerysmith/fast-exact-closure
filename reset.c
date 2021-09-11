#include "fec.h"

void reset(double *X) {
  double B[9];
  double evec[9];
#define Beval(i) B[(i)*3+(i)]
  double A[9];
  double c1, c2, c3, b0, e, b3, b03, temp, I[6];
  int i;

  B[0*3+0] = X[B11];
  B[0*3+1] = B[1*3+0] = X[B12];
  B[0*3+2] = B[2*3+0] = X[B13];
  B[1*3+1] = X[B22];
  B[1*3+2] = B[2*3+1] = X[B23];
  B[2*3+2] = (1+pow(X[B13],2)*X[B22]-2*X[B12]*X[B13]*X[B23]+X[B11]*pow(X[B23],2))
              /(X[B11]*X[B22]-pow(X[B12],2));

  diagonalize_sym_3(B, evec);

  if (fabs(Beval(0)-Beval(1)) + fabs(Beval(0)-Beval(2)) + fabs(Beval(1)-Beval(2)) < 3e-2) {
    c1 = Beval(0)-1;
    c2 = Beval(1)-1;
    c3 = Beval(2)-1;
#define calc_A_three_close(c1,c2,c3) 0.3333333333333333 \
   - (3*c1)/10. - c2/10. - c3/10. \
   +(15*pow(c1,2))/56. + (3*c1*c2)/28. + (3*pow(c2,2))/56. + \
    (3*c1*c3)/28. + (c2*c3)/28. + (3*pow(c3,2))/56. \
   -(35*pow(c1,3))/144. - (5*pow(c1,2)*c2)/48. - \
    (c1*pow(c2,2))/16. - (5*pow(c2,3))/144. - \
    (5*pow(c1,2)*c3)/48. - (c1*c2*c3)/24. - (pow(c2,2)*c3)/48. - \
    (c1*pow(c3,2))/16. - (c2*pow(c3,2))/48. - (5*pow(c3,3))/144. \
   +(315*pow(c1,4))/1408. + (35*pow(c1,3)*c2)/352. + \
    (45*pow(c1,2)*pow(c2,2))/704. + (15*c1*pow(c2,3))/352. + \
    (35*pow(c2,4))/1408. + (35*pow(c1,3)*c3)/352. + \
    (15*pow(c1,2)*c2*c3)/352. + (9*c1*pow(c2,2)*c3)/352. + \
    (5*pow(c2,3)*c3)/352. + (45*pow(c1,2)*pow(c3,2))/704. + \
    (9*c1*c2*pow(c3,2))/352. + (9*pow(c2,2)*pow(c3,2))/704. + \
    (15*c1*pow(c3,3))/352. + (5*c2*pow(c3,3))/352. + \
    (35*pow(c3,4))/1408.
    memset(A,0,sizeof(A));
    A[0*3+0] = calc_A_three_close(c1,c2,c3);
    A[1*3+1] = calc_A_three_close(c2,c1,c3);
    A[2*3+2] = calc_A_three_close(c3,c1,c2);
  }
#define else_if_two_close(arg0,arg1,arg2) \
  else if (fabs(Beval(arg0)-Beval(arg1)) < 1e-2) { \
    b0 = 0.5*(Beval(arg0)+Beval(arg1)); \
    e = 0.5*(Beval(arg0)-Beval(arg1)); \
    b3 = Beval(arg2); \
    b03 = b0-b3; \
    temp = sqrt(b3)/b03; \
    I[0] = (b03>0) ? \
           2./sqrt(b03)*acos(sqrt(b3/b0)) : \
           2./sqrt(-b03)*acosh(sqrt(b3/b0)); \
    for (i=1;i<=5;i++) { \
      temp /= b0; \
      I[i] = (2*i-1)/2./i/b03*I[i-1] - temp/i; \
    } \
    memset(A,0,sizeof(A)); \
    A[arg0*3+arg0] = 1./2*I[1] - 1./2*I[2]*e + 3./4*I[3]*e*e \
                      - 3./4*I[4]*e*e*e + 15./16*I[5]*e*e*e*e; \
    A[arg1*3+arg1] = 1./2*I[1] + 1./2*I[2]*e + 3./4*I[3]*e*e \
                      + 3./4*I[4]*e*e*e + 15./16*I[5]*e*e*e*e; \
    A[arg2*3+arg2] = 1 - A[arg0*3+arg0] - A[arg1*3+arg1]; \
  }
  else_if_two_close(0,1,2)
  else_if_two_close(0,2,1)
  else_if_two_close(1,2,0)
  else {
    return;
  }

  reverse_rotate(A,evec);

  X[A11] = A[0*3+0];
  X[A12] = A[0*3+1];
  X[A13] = A[0*3+2];
  X[A22] = A[1*3+1];
  X[A23] = A[1*3+2];
}
