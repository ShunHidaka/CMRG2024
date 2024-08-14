#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "function_util.h"
#include "function_blas.h"

//char *ANAME="/home/u11674/DATA/VCNT90000_A.csr",   *BNAME="/home/u11674/DATA/VCNT90000_B.csr";
char *ANAME="/home/u11674/DATA/VCNT10800h_A.csr", *BNAME="/home/u11674/DATA/VCNT10800h_B.csr";
//char *ANAME="/home/u11674/DATA/PPE3594_A.csr",   *BNAME="/home/u11674/DATA/PPE3594_B.csr";
double EPS_INNER = 1e-13;
double EPS_OUTER = 1e-13;
int OUTPUT_J = 10;
int OUTPUT_K = 1;
int MAX_ITR = 10000;
int ZERO = 0;
int ONE  = 1;
double dTMP = 0;
double complex cTMP = 0.0 + 0.0I;

void set_rhsVector(int N, double complex **b, double *norm)
{
  *b = (double complex *)calloc(N, sizeof(double complex));
  for(int i=0; i<N; i++)
    (*b)[i] = 1.0;
  *norm = dznrm2_(&N, *b, &ONE);
}

void set_shifts(int *M, double complex **sigma)
{
  /*
  *M = 2;
  *sigma = (double complex *)calloc(*M, sizeof(double complex));
  for(int j=1; j<=(*M); j++)
    (*sigma)[j-1] = 0.001 * cexp( 2 * M_PI * I * (j - 0.5) / (*M) );
  */

  *M = 50;
  *sigma = (double complex *)calloc(*M, sizeof(double complex));
  for(int j=1; j<=(*M); j++)
    (*sigma)[j-1] = 0.01 * cexp( 2 * M_PI * I * (j - 0.5) / (*M) );

  /*
  *M = 1001;
  *sigma = (double complex *)calloc(*M, sizeof(double complex));
  for(int j=0; j<(*M); j++)
    (*sigma)[j] = (0.4 + 0.001*j) + 0.01I;
  */
}

void SpMV(const int *A_row, const int *A_col, const double complex *A_ele,
          const double complex *x, double complex *b, int N)
{
  int i, j;
  double complex tmp;
  //#pragma omp prallel for private(j, tmp)
  for(i=0; i<N; i++){
    tmp = 0.0+0.0I;
    for(j=A_row[i]; j<A_row[i+1]; j++)
      tmp += A_ele[j]*x[A_col[j]];
    b[i] = tmp;
  }
}

int CG_method(const int *B_row, const int *B_col, const double complex *B_ele,
              double complex *x, const double complex *b, int N, const double tol)
{
  int status = 0;
  double complex r[N], p[N], Bp[N];
  double complex alpha, beta, rr, rr_old;
  for(int i=0; i<N; i++){
    x[i] = 0.0 + 0.0I;
    r[i] = b[i];
    p[i] = r[i];
  }
  rr = zdotc_(&N, r, &ONE, r, &ONE);
  double r0_nrm = dznrm2_(&N, r, &ONE);
  for(int j=0; j<MAX_ITR; j++){
    SpMV(B_row,B_col,B_ele, p, Bp, N);
    alpha = rr / zdotc_(&N, p, &ONE, Bp, &ONE);
    zaxpy_(&N, &alpha, p, &ONE, x, &ONE);
    alpha = -alpha;
    zaxpy_(&N, &alpha, Bp, &ONE, r, &ONE);
    if(dznrm2_(&N, r, &ONE)/r0_nrm <= tol){
      //fprintf(stderr, "[CG] %d %lf %lf\n", j, dznrm2_(&N, r, &ONE), r0_nrm);
      status = 1;
      break;
    }
    alpha = -alpha;
    rr_old = rr;
    rr = zdotc_(&N, r, &ONE, r, &ONE);
    beta = rr / rr_old;
    zscal_(&N, &beta, p, &ONE);
    cTMP = 1.0 + 0.0I;
    zaxpy_(&N, &cTMP, r, &ONE, p, &ONE);
  }
  return status;
}

FILE *fopen_mtx(const char *fname, const char *mode, int *row_size, int *col_size, int *ele_size)
{
  FILE *fp = NULL;
  fp = fopen(fname, mode);
  if(fp == NULL){
    fprintf(stderr, "Can not open file : %s\n", fname);
    exit(1);
  }
  char chr;
  while( (chr = fgetc(fp)) != EOF && (chr=='%' || chr=='#') )
    while( (chr = fgetc(fp)) != EOF )
      if(chr == '\n') break;
  fseek(fp, -sizeof(char), SEEK_CUR);
  int num, tmp1, tmp2, tmp3;
  num = fscanf(fp, "%d %d %d", &tmp1, &tmp2, &tmp3);
  if(num != 3){ fprintf(stderr, "fopen err\n"); exit(1);}
  if(col_size != NULL) *row_size = tmp1;
  if(row_size != NULL) *col_size = tmp2;
  if(ele_size != NULL) *ele_size = tmp3;
  return fp;
}
void read_csr(const char *fname, int *N, int *DATASIZE,
              int **row_ptr, int **col_ind, double complex **element)
{
  FILE *fp;
  int rsize, csize, esize;
  fp = fopen_mtx(fname, "r", &rsize, &csize, &esize);
  *N = rsize-1; *DATASIZE = esize;

  *row_ptr = (int *)calloc(rsize, sizeof(int));
  *col_ind = (int *)calloc(csize, sizeof(int));
  *element = (double complex *)calloc(esize, sizeof(double complex));

  int row, col, i=0;
  double real, imag=0;
  while( fscanf(fp, "%d %d %lf %lf", &row, &col, &real, &imag) != EOF ){
    if(i < rsize) (*row_ptr)[i] = row;
    if(i < csize) (*col_ind)[i] = col;
    if(i < esize) (*element)[i] = real + imag*I;
    i = i + 1;
  }
  fclose(fp);
}
