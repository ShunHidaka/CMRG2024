// solve (A + \sigma_{k}B)\vb{x} = \vb{b}
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include "function_util.h"
#include "function_blas.h"


int main(int argc, char *argv[]){
  fprintf(stdout, "# Generalizd shifted MINRES method\n");
  int j, k;
  // prepare Matrix A, B
  int N, DATASIZE;
  int *A_row, *A_col; double complex *A_ele;
  int *B_row, *B_col; double complex *B_ele;
  read_csr(ANAME, &N, &DATASIZE, &A_row, &A_col, &A_ele);
  read_csr(BNAME, &N, &DATASIZE, &B_row, &B_col, &B_ele);
  fprintf(stdout, "# A = %s\n# B = %s\n", ANAME, BNAME);

  switch( atoi(argv[1]) ){
  case 0:
    EPS_INNER = 1e-12;
    break;
  case 1:
    EPS_INNER = 1e-10;
    break;
  case 2:
    EPS_INNER = 1e-8;
    break;
  case 3:
    EPS_INNER = 1e-6;
    break;
  default:
    fprintf(stderr, "./%s {0,1,2,3} // {1e-13,1e-12,1e-10,1e-8}\n", argv[0]);
  }
  fprintf(stdout, "# inner=%e, outet=%e\n", EPS_INNER, EPS_OUTER);

  // prepare rhs vector b
  double complex *b;
  double bnrm;
  set_rhsVector(N, &b, &bnrm);
  // prepare shifts \sigma
  int M=1;
  double complex *sigma;
  set_shifts(&M, &sigma);

  // declare variables
  double complex *w, *u;
  double         *alpha, *beta;
  double complex *T;
  double         **G1;
  double complex **G2;
  double complex **p, **x, *f;
  double         *h, *res, r0nrm;
  // allocate memory
  w     = (double complex  *)calloc(N*3, sizeof(double complex));
  u     = (double complex  *)calloc(N*3, sizeof(double complex));
  alpha = (double          *)calloc(1,   sizeof(double));
  beta  = (double          *)calloc(2,   sizeof(double));
  T     = (double complex  *)calloc(4,   sizeof(double complex));
  G1    = (double         **)calloc(M,   sizeof(double *));
  G2    = (double complex **)calloc(M,   sizeof(double complex *));
  p     = (double complex **)calloc(M,   sizeof(double complex *));
  x     = (double complex **)calloc(M,   sizeof(double complex *));
  f     = (double complex  *)calloc(M,   sizeof(double complex));
  h     = (double          *)calloc(M,   sizeof(double));
  res   = (double          *)calloc(M,   sizeof(double));
  for(k=0; k<M; k++){
    G1[k] = (double         *)calloc(3,   sizeof(double));
    G2[k] = (double complex *)calloc(3,   sizeof(double complex));
    p[k]  = (double complex *)calloc(N*3, sizeof(double complex));
    x[k]  = (double complex *)calloc(N,   sizeof(double complex));
  }
  int conv_num = 0, ret_CG;
  int *is_conv = (int *)calloc(M, sizeof(int));

  double complex *r_tmp1, *r_tmp2;
  r_tmp1 = (double complex *)calloc(N, sizeof(double complex));
  r_tmp2 = (double complex *)calloc(N, sizeof(double complex));

  // Generalized shfited MINRES method
  // preprocess
  ret_CG = CG_method(B_row,B_col,B_ele, &(u[N]), b, N, EPS_INNER);
  if(ret_CG == 0){ fprintf(stdout, "# unconverged in preprocess\n"); exit(1);}
  r0nrm = sqrt( creal(zdotc_(&N, b, &ONE, &(u[N]), &ONE)) );
  dTMP = 1 / r0nrm;
  zcopy_(&N, &(u[N]), &ONE, &(w[N]), &ONE);
  zcopy_(&N,       b, &ONE, &(u[N]), &ONE);
  zdscal_(&N, &dTMP, &(w[N]), &ONE);
  zdscal_(&N, &dTMP, &(u[N]), &ONE);
  beta[0] = 0;
  for(k=0; k<M; k++){
    f[k] = 1.0;
    h[k] = r0nrm;
  }

  // main loop
  time_t start_time, end_time;
  start_time = time(NULL);
  for(j=1; j<=MAX_ITR; j++){
    // generalized lanczos process
    SpMV(A_row,A_col,A_ele, &(w[N]), &(u[2*N]), N);
    alpha[0] = creal( zdotc_(&N, &(w[N]), &ONE, &(u[2*N]), &ONE) );
    cTMP = -alpha[0];
    zaxpy_(&N, &cTMP, &(u[N]), &ONE, &(u[2*N]), &ONE);
    cTMP = -beta[0];
    zaxpy_(&N, &cTMP, &(u[0]), &ONE, &(u[2*N]), &ONE);
    ret_CG = CG_method(B_row,B_col,B_ele, &(w[2*N]), &(u[2*N]), N, EPS_INNER);
    if(ret_CG == 0){ fprintf(stderr, "# unconverged in %d inner solution\n", j); exit(1);}
    beta[1] = sqrt( creal(zdotc_(&N, &(u[2*N]), &ONE, &(w[2*N]), &ONE)) );
    dTMP = 1 / beta[1];
    zdscal_(&N, &dTMP, &(u[2*N]), &ONE);
    zdscal_(&N, &dTMP, &(w[2*N]), &ONE);
    // shift systems
    for(k=0; k<M; k++){
      if(is_conv[k] != 0)
        continue;
      T[0]=0.0;
      T[1]=beta[0]; T[2]=alpha[0]+sigma[k]; T[3]=beta[1];
      if(j >= 3){ zrot_(&ONE, &(T[0]), &ONE, &(T[1]), &ONE, &(G1[k][0]), &(G2[k][0]));}
      if(j >= 2){ zrot_(&ONE, &(T[1]), &ONE, &(T[2]), &ONE, &(G1[k][1]), &(G2[k][1]));}
      zlartg_(&(T[2]), &(T[3]), &(G1[k][2]), &(G2[k][2]), &cTMP);

      zcopy_(&N, &(p[k][N]), &ONE, &(p[k][0]), &ONE);
      zcopy_(&N, &(p[k][2*N]), &ONE, &(p[k][N]), &ONE);
      zcopy_(&N, &(w[N]), &ONE, &(p[k][2*N]), &ONE);
      cTMP = -T[0];
      zaxpy_(&N, &cTMP, &(p[k][0]), &ONE, &(p[k][2*N]), &ONE);
      cTMP = -T[1];
      zaxpy_(&N, &cTMP, &(p[k][N]), &ONE, &(p[k][2*N]), &ONE);
      cTMP = 1.0 / T[2];
      zscal_(&N, &cTMP, &(p[k][2*N]), &ONE);

      cTMP = r0nrm * G1[k][2] * f[k];
      zaxpy_(&N, &cTMP, &(p[k][2*N]), &ONE, &(x[k][0]), &ONE);
      
      f[k] = -conj(G2[k][2])*f[k];
      h[k]    = cabs( -conj(G2[k][2]) ) * h[k];
      G1[k][0]=G1[k][1]; G1[k][1]=G1[k][2];
      G2[k][0]=G2[k][1]; G2[k][1]=G2[k][2];

      res[k] = h[k];
      if(res[k]/r0nrm < EPS_OUTER){
        conv_num++;
        is_conv[k] = j;
        continue;
      }
    }
    beta[0] = beta[1];
    zcopy_(&N, &(w[N]),   &ONE, &(w[0]), &ONE);
    zcopy_(&N, &(w[2*N]), &ONE, &(w[N]), &ONE);
    zcopy_(&N, &(u[N]),   &ONE, &(u[0]), &ONE);
    zcopy_(&N, &(u[2*N]), &ONE, &(u[N]), &ONE);

    // Determine Convergence
    if(conv_num == M){
      break;
    }
  }
  end_time = time(NULL);

  fprintf(stdout, "# status=-%d, time=%ld\n", j, end_time-start_time);
  fprintf(stdout, "# converged %d / %d\n", conv_num, M);
  for(k=0; k<M; k++){
    zcopy_(&N, b, &ONE, r_tmp1, &ONE);
    SpMV(A_row,A_col,A_ele, &(x[k][0]), r_tmp2, N);
    cTMP = -1.0;
    zaxpy_(&N, &cTMP, r_tmp2, &ONE, r_tmp1, &ONE);
    SpMV(B_row,B_col,B_ele, &(x[k][0]), r_tmp2, N);
    cTMP = -sigma[k];
    zaxpy_(&N, &cTMP, r_tmp2, &ONE, r_tmp1, &ONE);
    fprintf(stdout, "%d %lf %lf %d %e %e\n",
	    k, creal(sigma[k]), cimag(sigma[k]),
	    is_conv[k], res[k]/r0nrm, dznrm2_(&N,r_tmp1,&ONE)/bnrm);
  }
  return 0;
}
