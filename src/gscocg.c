#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include "function_util.h"
#include "function_blas.h"


int main(int argc, char *argv[])
{
  fprintf(stdout, "# Generalizd shifted COCG method\n");
  int j, k;
  // prepare Matrix A
  int N, DATASIZE;
  int *A_row, *A_col; double complex *A_ele;
  int *B_row, *B_col; double complex *B_ele;
  read_csr(ANAME, &N, &DATASIZE, &A_row, &A_col, &A_ele);
  read_csr(BNAME, &N, &DATASIZE, &B_row, &B_col, &B_ele);
  fprintf(stdout, "# A = %s\n# B = %s\n", ANAME, BNAME);
  fprintf(stdout, "# inner=%e, outer=%e\n", EPS_INNER, EPS_OUTER);
  // prepare right-hand side vector b
  double complex *b;
  double bnrm;
  set_rhsVector(N, &b, &bnrm);
  // prepare shifts sigma
  int M;
  double complex *sigma;
  set_shifts(&M, &sigma);

  // declare variables
  int s, t;
  double complex **x, **r, **p, **pi;
  double complex *Ap, *Bp, *Binv_r;
  double complex *alpha, *alpha_old, *beta;
  double complex rr, rr_old;
  double         *res, r0nrm;
  // allocate memory
  x         = (double complex **)calloc(M, sizeof(double complex *));
  r         = (double complex **)calloc(M, sizeof(double complex *));
  p         = (double complex **)calloc(M, sizeof(double complex *));
  pi        = (double complex **)calloc(M, sizeof(double complex *));
  Ap        = (double complex  *)calloc(N, sizeof(double complex));
  Bp        = (double complex  *)calloc(N, sizeof(double complex));
  Binv_r    = (double complex  *)calloc(N, sizeof(double complex));
  alpha     = (double complex  *)calloc(M, sizeof(double complex));
  alpha_old = (double complex  *)calloc(M, sizeof(double complex));
  beta      = (double complex  *)calloc(M, sizeof(double complex));
  res       = (double          *)calloc(M, sizeof(double));
  for(k=0; k<M; k++){
    x[k]  = (double complex *)calloc(N, sizeof(double complex));
    r[k]  = (double complex *)calloc(N, sizeof(double complex));
    p[k]  = (double complex *)calloc(N, sizeof(double complex));
    pi[k] = (double complex *)calloc(3, sizeof(double complex));
  }
  int conv_num = 0, ret_CG;
  int *is_conv = (int *)calloc(M, sizeof(int));

  double tmp1, tmp2;
  double complex *r_tmp1, *r_tmp2;
  r_tmp1 = (double complex *)calloc(N, sizeof(double complex));
  r_tmp2 = (double complex *)calloc(N, sizeof(double complex));

  // Generalized shifted COCG method
  // preprocess
  for(k=0; k<M; k++){
    zcopy_(&N, &(b[0]), &ONE, &(r[k][0]), &ONE);
    pi[k][0] = pi[k][1] = 1;
  }
  r0nrm = bnrm;
  s = 0;
  ret_CG = CG_method(B_row,B_col,B_ele, Binv_r, r[s], N, EPS_INNER);
  if(ret_CG == 0){ fprintf(stdout, "# unconverged in : preprocess\n"); exit(1);}
  rr = zdotu_(&N, &(r[s][0]), &ONE, Binv_r, &ONE);
  alpha[s] = 1; beta[s] = 0;

  {
    fprintf(stderr, "0");
    for(k=0; k<M; k+=OUTPUT_K){
      zcopy_(&N, b, &ONE, r_tmp1, &ONE);
      SpMV(A_row,A_col,A_ele, &(x[k][0]), r_tmp2, N);
      cTMP = -1.0;
      zaxpy_(&N, &cTMP, r_tmp2, &ONE, r_tmp1, &ONE);
      SpMV(B_row,B_col,B_ele, &(x[k][0]), r_tmp2, N);
      cTMP = -sigma[k];
      zaxpy_(&N, &cTMP, r_tmp2, &ONE, r_tmp1, &ONE);
      tmp1=r0nrm; tmp2=dznrm2_(&N, r_tmp1, &ONE);
      fprintf(stderr, " %e %e", tmp1/r0nrm, tmp2/r0nrm);
    }
    fprintf(stderr, "\n");
  }

  // main loop
  time_t start_time, end_time;
  start_time = time(NULL);
  for(j=1; j<=MAX_ITR; j++){
    // seed system
    cTMP = beta[s];
    zscal_(&N, &cTMP, &(p[s][0]), &ONE);
    cTMP = 1.0 + 0.0I;
    zaxpy_(&N, &cTMP, Binv_r, &ONE, &(p[s][0]), &ONE);
    alpha_old[s] = alpha[s];
    SpMV(A_row,A_col,A_ele, p[s], Ap, N);
    SpMV(B_row,B_col,B_ele, p[s], Bp, N);
    cTMP = sigma[s];
    zaxpy_(&N, &cTMP, Bp, &ONE, Ap, &ONE);
    alpha[s] = rr / zdotu_(&N, &(p[s][0]), &ONE, &(Ap[0]), &ONE);
    cTMP = alpha[s];
    zaxpy_(&N, &cTMP, &(p[s][0]),&ONE, &(x[s][0]), &ONE);
    cTMP = -alpha[s];
    zaxpy_(&N, &cTMP, &(Ap[0]), &ONE, &(r[s][0]), &ONE);
    res[s] = dznrm2_(&N, &(r[s][0]), &ONE);
    if(res[s]/r0nrm < EPS_OUTER){
      conv_num++;
      is_conv[s] = j+1;
    }
    // shift systems
    for(k=0; k<M; k++){
      if(is_conv[k] != 0 || k == s)
        continue;
      pi[k][2] = (1 + (beta[s]/alpha_old[s])*alpha[s] + alpha[s]*(sigma[k]-sigma[s]))*pi[k][1] - (beta[s]/alpha_old[s])*alpha[s]*pi[k][0];
      alpha[k] = (pi[k][1]/pi[k][2])*alpha[s];
      beta[k] = (pi[k][0]/pi[k][1])*(pi[k][0]/pi[k][1])*beta[s];

      cTMP = beta[k];
      zscal_(&N, &cTMP, &(p[k][0]), &ONE);
      cTMP = 1.0 / pi[k][1];
      zaxpy_(&N, &cTMP, Binv_r, &ONE, &(p[k][0]), &ONE);
      cTMP = alpha[k];
      zaxpy_(&N, &cTMP, &(p[k][0]), &ONE, &(x[k][0]), &ONE);
      zcopy_(&N, &(r[s][0]), &ONE, &(r[k][0]), &ONE);

      cTMP = 1 / pi[k][2];
      zscal_(&N, &cTMP, &(r[k][0]), &ONE);

      res[k] = dznrm2_(&N, &(r[k][0]), &ONE);
      if(res[k]/r0nrm < EPS_OUTER){
        conv_num++;
        is_conv[k] = j+1;
      }
      pi[k][0] = pi[k][1];
      pi[k][1] = pi[k][2];
    }
    rr_old = rr;
    ret_CG = CG_method(B_row,B_col,B_ele, Binv_r, r[s], N, EPS_INNER);
    if(ret_CG == 0){ fprintf(stdout, "# unconverged in %d : inner solution\n", j); exit(1);}
    rr = zdotu_(&N, &(r[s][0]), &ONE, Binv_r, &ONE);
    beta[s] = rr / rr_old;
    // seed switching
    if(is_conv[s] != 0){
      t = s;
      for(k=0; k<M; k++) if(res[k] > res[t] && k != s) t = k;
      beta[t]  = (pi[t][0]/pi[t][1])*(pi[t][0]/pi[t][1])*beta[s];
      for(k=0; k<M; k++){
        if(k == t) continue;
        pi[k][0] = pi[k][0] / pi[t][0];
        pi[k][1] = pi[k][1] / pi[t][1];
      }
      fprintf(stdout, "# SWITCH [%d] to [%d] in %d : %lf %lf\n", s, t, j, log10(res[s]), log10(res[t]));
      ret_CG = CG_method(B_row,B_col,B_ele, Binv_r, r[t], N, EPS_INNER);
      if(ret_CG == 0){ fprintf(stdout, "# unconverged in %d : seed switching\n", j); exit(1);}
      rr = zdotu_(&N, &(r[t][0]), &ONE, Binv_r, &ONE);
      s = t;
    }

    {
      if(j % OUTPUT_J == 0){
        fprintf(stderr, "%d", j);
        for(k=0; k<M; k+=OUTPUT_K){
          zcopy_(&N, b, &ONE, r_tmp1, &ONE);
          SpMV(A_row,A_col,A_ele, &(x[k][0]), r_tmp2, N);
          cTMP = -1.0;
          zaxpy_(&N, &cTMP, r_tmp2, &ONE, r_tmp1, &ONE);
          SpMV(B_row,B_col,B_ele, &(x[k][0]), r_tmp2, N);
          cTMP = -sigma[k];
          zaxpy_(&N, &cTMP, r_tmp2, &ONE, r_tmp1, &ONE);
          double tmp1=res[k], tmp2=dznrm2_(&N, r_tmp1, &ONE);
          fprintf(stderr, " %e %e", tmp1/r0nrm, tmp2/r0nrm);
        }
        fprintf(stderr, "\n");
      }
    }

    // Determine Convergence
    if(conv_num == M){
      break;
    }
  }
  end_time = time(NULL);

  fprintf(stdout, "# status=-%d, time=%ld\n", j+1, end_time-start_time);
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
	    is_conv[k], res[k]/r0nrm, dznrm2_(&N,r_tmp1,&ONE)/r0nrm);
  }
  return 0;
}
