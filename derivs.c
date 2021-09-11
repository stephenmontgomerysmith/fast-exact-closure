#include "fec.h"

static void check_for_ip_rights();

extern int no_close, two_close[3], three_close;

void derivs(double t, double *X, double *Xdot, param_list_t *param) {
  double A[9];
  double B[9];
  double Aeval[3];
#define Beval(i) B[(i)*3+(i)]
  double evec[9];
  double dA[9];
  double dB[9];
  double g; // sqrt(0.5 Gamma:Gamma)
  double H[9]; // Omega + lambda Gamma

/* Used for ARD */
  double Dr[9];
  double gamm[9];
  double gamm_2[9]; // gamm^2
  double trace_Dr;
  double F1[9]; // F1 = 5(B.Dr + Dr.B)
  double F2[9]; // F2 = 2 Dr + 3 tr(Dr) A

/* Used for Folgar-Tucker */
  double F[3]; // F = CI g (2I - 6A), diagonal entries

/* Used to calculate ciijj */
  double ciijj[9];
  double diijj[9];
  double c1, c2, c3, b0, e, b3, b03, temp, I[5];

  int i,j,k;

/* Extract the tensors A and B */

  B[0*3+0] = X[B11];
  B[0*3+1] = B[1*3+0] = X[B12];
  B[0*3+2] = B[2*3+0] = X[B13];
  B[1*3+1] = X[B22];
  B[1*3+2] = B[2*3+1] = X[B23];
// Det(B) = 1
  B[2*3+2] = (1+pow(X[B13],2)*X[B22]-2*X[B12]*X[B13]*X[B23]+X[B11]*pow(X[B23],2))
             /(X[B11]*X[B22]-pow(X[B12],2));

  A[0*3+0] = X[A11];
  A[0*3+1] = A[1*3+0] = X[A12];
  A[0*3+2] = A[2*3+0] = X[A13];
  A[1*3+1] = X[A22];
  A[1*3+2] = A[2*3+1] = X[A23];
// Tr(A) = 1
  A[2*3+2] = 1-X[A11]-X[A22];

/* Diagonalize B and use this to diagonalize A */

  diagonalize_sym_3(B, evec);
  rotate_diag(Aeval,A,evec);

/* Calculate C, but only C_{iijj}, since C_{ijij}=C_{iijj} */

  if (fabs(Beval(0)-Beval(1)) + fabs(Beval(0)-Beval(2)) + fabs(Beval(1)-Beval(2)) < 3e-2) {
    three_close++;
    c1 = Beval(0)-1;
    c2 = Beval(1)-1;
    c3 = Beval(2)-1;
#define calc_c_three_close(c1,c2,c3) 0.1 - (3*c1)/28. - (3*c2)/28. - c3/28. \
  + (5*pow(c1,2))/48. + (c1*c2)/8. + (5*pow(c2,2))/48. \
  + (c1*c3)/24. + (c2*c3)/24. + pow(c3,2)/48. \
  - (35*pow(c1,3))/352. - (45*pow(c1,2)*c2)/352. - (45*c1*pow(c2,2))/352. \
  - (35*pow(c2,3))/352. - (15*pow(c1,2)*c3)/352. \
  - (9*c1*c2*c3)/176. - (15*pow(c2,2)*c3)/352. - (9*c1*pow(c3,2))/352. \
  - (9*c2*pow(c3,2))/352. - (5*pow(c3,3))/352.
    ciijj[0*3+1] = ciijj[1*3+0] = calc_c_three_close(c1,c2,c3);
    ciijj[0*3+2] = ciijj[2*3+0] = calc_c_three_close(c1,c3,c2);
    ciijj[1*3+2] = ciijj[2*3+1] = calc_c_three_close(c2,c3,c1);
  }
#define else_if_two_close(arg0,arg1,arg2) \
  else if (fabs(Beval(arg0)-Beval(arg1)) < 1e-2) { \
    two_close[arg2]++; \
    b0 = 0.5*(Beval(arg0)+Beval(arg1)); \
    e = 0.5*(Beval(arg0)-Beval(arg1)); \
    b3 = Beval(arg2); \
    b03 = b0-b3; \
    temp = sqrt(b3)/b03; \
    I[0] = (b03>0) ? \
           2./sqrt(b03)*acos(sqrt(b3/b0)) : \
           2./sqrt(-b03)*acosh(sqrt(b3/b0)); \
    for (i=1;i<5;i++) { \
      temp /= b0; \
      I[i] = (2*i-1)/2./i/b03*I[i-1] - temp/i; \
    } \
    ciijj[arg0*3+arg1] = ciijj[arg1*3+arg0] = 0.25*I[2] + 3/8.*I[4]*e*e; \
    ciijj[arg0*3+arg2] = ciijj[arg2*3+arg0] \
                       = (Aeval[arg0]-Aeval[arg2])/2/(Beval(arg2)-Beval(arg0)); \
    ciijj[arg1*3+arg2] = ciijj[arg2*3+arg1] \
                       = (Aeval[arg1]-Aeval[arg2])/2/(Beval(arg2)-Beval(arg1)); \
  }
  else_if_two_close(0,1,2)
  else_if_two_close(0,2,1)
  else_if_two_close(1,2,0)
  else {
    no_close++;
#define do_calc_c(arg0,arg1) \
  ciijj[arg0*3+arg1] = ciijj[arg1*3+arg0] \
                     = (Aeval[arg0]-Aeval[arg1])/2/(Beval(arg1)-Beval(arg0))
    do_calc_c(0,1);
    do_calc_c(0,2);
    do_calc_c(1,2);
  }
  ciijj[0*3+0] = 0.5/Beval(0) - ciijj[0*3+1] - ciijj[0*3+2];
  ciijj[1*3+1] = 0.5/Beval(1) - ciijj[1*3+0] - ciijj[1*3+2];
  ciijj[2*3+2] = 0.5/Beval(2) - ciijj[2*3+0] - ciijj[2*3+1];

/* Calculate D, but only D_{iijj} (since D_{ijij}=1/4/C_{ijij}) */

  inverse_sym(ciijj,diijj);

/* g = sqrt(0.5 Gamma:Gamma) */

  g = 0;
  for (i=0;i<3;i++) for (j=0;j<3;j++)
    g += pow(param->gamm[i*3+j],2);
  g = sqrt(g/2);

/* H = Omega + lambda Gamma */

  for (i=0;i<9;i++)
    H[i] = param->lambda*param->gamm[i];

  H[0*3+1] += -param->w[2];
  H[1*3+0] +=  param->w[2];
  H[0*3+2] +=  param->w[1];
  H[2*3+0] += -param->w[1];
  H[1*3+2] += -param->w[0];
  H[2*3+1] +=  param->w[0];

  rotate(H,evec);

/* Start computations of DA/DT and DB/DT */

/* First the Jeffery's part */

// DB/Dt = - (B.H + H^T.B)/2

  for (i=0;i<3;i++) for (j=0;j<3;j++)
    dB[i*3+j] = -(H[i*3+j]*Beval(i) + H[j*3+i]*Beval(j))/2;

// DA/Dt = C:(B.H + H^T.B)/2

  for (i=0;i<3;i++) {
    dA[i*3+i] = 0;
    for (j=0;j<3;j++)
      dA[i*3+i] -= ciijj[i*3+j]*dB[j*3+j];
  }
  dA[0*3+1] = dA[1*3+0] = - 2*ciijj[0*3+1]*dB[0*3+1];
  dA[0*3+2] = dA[2*3+0] = - 2*ciijj[0*3+2]*dB[0*3+2];
  dA[1*3+2] = dA[2*3+1] = - 2*ciijj[1*3+2]*dB[1*3+2];

/* Then the diffusion part, either ARD ... */

  if (param->do_ard) {

// Dr = g * (b1 I + b2 A + b3 A^2 + b4/2/g Gamma + b5/4/g^2 Gamma^2)

    memcpy(gamm,param->gamm,sizeof(gamm));
    rotate(gamm,evec);
    for (i=0;i<3;i++) for (j=0;j<3;j++) {
      gamm_2[i*3+j] = 0;
      for (k=0;k<3;k++)
        gamm_2[i*3+k] += gamm[i*3+j]*gamm[j*3+k];
    }

    for (i=0;i<3;i++) for (j=0;j<3;j++)
      Dr[i*3+j] = param->b4/2*param->gamm[i*3+j]+param->b5/4/g*gamm_2[i*3+j];
    for (i=0;i<3;i++)
      Dr[i*3+i] += g*(param->b1+param->b2*Aeval[i]+param->b3*pow(Aeval[i],2));

    trace_Dr = 0;
    for (i=0;i<3;i++)
      trace_Dr += Dr[i*3+i];

// F1 = 5(B.Dr + Dr.B)
// F2 = 2 Dr + 3 tr(Dr) A

    for (i=0;i<3;i++) for (j=0;j<3;j++)
      F1[i*3+j] = 5*Dr[i*3+j]*(Beval(i)+Beval(j));
    for (i=0;i<3;i++) for (j=0;j<3;j++)
      F2[i*3+j] = 2*Dr[i*3+j];
    for (i=0;i<3;i++)
      F2[i*3+i] += 3*trace_Dr*Aeval[i];

// DA/Dt += F2
// DB/Dt += F1

    for (i=0;i<9;i++) {
      dA[i] += F2[i];
      dB[i] += F1[i];
    }

// DA/Dt -= C:F1
// DB/Dt -= D:F2

    for (i=0;i<3;i++) for (j=0;j<3;j++) {
      dA[i*3+i] -= ciijj[i*3+j]*F1[j*3+j];
      dB[i*3+i] -= diijj[i*3+j]*F2[j*3+j];
    }
    dA[0*3+1] = dA[1*3+0] -= 2*ciijj[0*3+1]*F1[0*3+1];
    dA[0*3+2] = dA[2*3+0] -= 2*ciijj[0*3+2]*F1[0*3+2];
    dA[1*3+2] = dA[2*3+1] -= 2*ciijj[1*3+2]*F1[1*3+2];
    dB[0*3+1] = dB[1*3+0] -= 0.5/ciijj[0*3+1]*F2[0*3+1];
    dB[0*3+2] = dB[2*3+0] -= 0.5/ciijj[0*3+2]*F2[0*3+2];
    dB[1*3+2] = dB[2*3+1] -= 0.5/ciijj[1*3+2]*F2[1*3+2];

  } else {

/* ... or Folgar-Tucker */

// F = CI g (2I - 6A)
// DA/Dt += F
// DB/Dt -= D:F

    for (i=0;i<3;i++) {
      F[i] = param->CI*g*(2-6*Aeval[i]);
      dA[i*3+i] += F[i];
    }
    for (i=0;i<3;i++) for (j=0;j<3;j++)
      dB[i*3+i] -= diijj[i*3+j]*F[j];
  }

/* Perform RSC */

// Multiply diagonal entries of DA/Dt and DB/Dt by kappa

  if (param->do_rsc) {
    check_for_ip_rights();
    for (i=0;i<3;i++) {
      dA[i*3+i] *= param->kappa;
      dB[i*3+i] *= param->kappa;
    }
  }

/* Rotate back into the standard frame */

  reverse_rotate(dA,evec);
  reverse_rotate(dB,evec);

/* and return the values of DA/Dt and DB/Dt */

  Xdot[A11] = dA[0*3+0];
  Xdot[A12] = dA[0*3+1];
  Xdot[A13] = dA[0*3+2];
  Xdot[A22] = dA[1*3+1];
  Xdot[A23] = dA[1*3+2];
  Xdot[B11] = dB[0*3+0];
  Xdot[B12] = dB[0*3+1];
  Xdot[B13] = dB[0*3+2];
  Xdot[B22] = dB[1*3+1];
  Xdot[B23] = dB[1*3+2];
}

/* Miscellaneous functions */

/* Compute the inverse of a 3x3 symmetric matrix */

void inverse_sym(double *A, double *inv) {
  int i;
  double det;

  inv[0*3+0] = A[1*3+1]*A[2*3+2] - pow(A[1*3+2],2);
  inv[0*3+1] = inv[1*3+0] = A[0*3+2]*A[1*3+2] - A[0*3+1]*A[2*3+2];
  inv[0*3+2] = inv[2*3+0] = A[0*3+1]*A[1*3+2] - A[0*3+2]*A[1*3+1];
  inv[1*3+1] = A[0*3+0]*A[2*3+2] - pow(A[0*3+2],2);
  inv[1*3+2] = inv[2*3+1] = A[0*3+1]*A[0*3+2] - A[0*3+0]*A[1*3+2];
  inv[2*3+2] = A[0*3+0]*A[1*3+1] - pow(A[0*3+1],2);
  det = A[0*3+0]*inv[0*3+0] + A[1*3+0]*inv[0*3+1] + A[2*3+0]*inv[0*3+2];
  for (i=0;i<9;i++) inv[i] /= det;
}

/* Rotate rank 2 tensor into the principal frame */

void rotate(double *A, double *evec) {
  double temp[9];
  int i,j,k;

  for (i=0;i<3;i++) for (j=0;j<3;j++) {
    temp[3*i+j] = 0;
    for (k=0;k<3;k++)
      temp[3*i+j] += A[3*k+j]*evec[i*3+k];
  }

  for (i=0;i<3;i++) for (j=0;j<3;j++) {
    A[3*i+j] = 0;
    for (k=0;k<3;k++)
      A[3*i+j] += temp[3*i+k]*evec[j*3+k];
  }
}

/* Rotate rank 2 tensor into the principal frame, assuming result is diagonal */

void rotate_diag(double *Adiag, double *A, double *evec) {
  double temp[9];
  int i,j,k;

  for (i=0;i<3;i++) for (j=0;j<3;j++) {
    temp[3*i+j] = 0;
    for (k=0;k<3;k++)
      temp[3*i+j] += A[3*k+j]*evec[i*3+k];
  }

  for (i=0;i<3;i++) {
    Adiag[i] = 0;
    for (k=0;k<3;k++)
      Adiag[i] += temp[3*i+k]*evec[i*3+k];
  }
}

/* Rotate rank 2 tensor out of the principal frame, assuming result is symmetric */

void reverse_rotate(double *A, double *evec) {
  double temp[9];
  int i,j,k;

  for (i=0;i<3;i++) for (j=0;j<3;j++) {
    temp[3*i+j] = 0;
    for (k=0;k<3;k++)
      temp[3*i+j] += A[3*k+j]*evec[k*3+i];
  }

  for (i=0;i<3;i++) for (j=i;j<3;j++) {
    A[i*3+j] = 0;
    for (k=0;k<3;k++)
      A[3*i+j] += temp[3*i+k]*evec[k*3+j];
  }
  for (i=1;i<3;i++) for (j=0;j<i;j++)
    A[3*i+j] = A[3*j+i];
}

/* Check IP rights for RSC */

static int first_time=1;
static void check_for_ip_rights() {
  char answer[10];

  if (first_time) {
    first_time = 0;
    if (getenv("MAY_USE_PATENT_7266469")==NULL) {
      printf("Do you understand that this program uses intellectual property protected by U.S. Patent No. 7,266,469? ");
      if (fgets(answer,10,stdin) != NULL && strcmp(answer,"yes\n")!=0) {
        printf("This program will not run if given an answer other than \"yes\" (written in full).\n");
        exit(1);
      }
    }
  }
}
