#ifndef FUNCTION_BLAS_H
#define FUNCTION_BLAS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>


// 倍精度複素数ベクトルに倍精度(実数/複素数)のスカラ倍を行う
// x <= a*x
void zdscal_(int *n, double *a, double complex *x, int *incx);
void zscal_(int *n, double complex *a, double complex *x, int *incx);

// 倍精度複素数ベクトル同士の加算
// y <= alpha*x + y
void zaxpy_(int *n, double complex *alpha, double complex *x, int *incx, double complex *y, int *incy);

// 倍精度複素数ベクトルのコピー
// y <= x
void zcopy_(int *n, double complex *x, int *incx, double complex *y, int *incy);

// 倍精度複素数ベクトルどうしの内積（共役転置）
// ret <= <x, y> = x^H y
double complex zdotc_(int *n, double complex *x, int *incx, double complex *y, int *incy);

// 倍精度複素数ベクトルどうしの疑似的な内積（転置）
// ret <= <x, y>' = x^T y
double complex zdotu_(int *n, double complex *x, int *incx, double complex *y, int *incy);

// 倍精度複素数ベクトルのノルム
// ret <= ||x||
double dznrm2_(int *n, double complex *x, int *incx);

// 倍精度複素数で与えられた点に対してy座標をゼロにするGivens回転を求める
void zrotg_(double complex *a,double complex *b, double *c, double complex *s);
void zlartg_(double complex *f, double complex *g, double *c, double complex *s, double complex *r);

// 倍精度複素数で生成されたGivens回転を適用する
void zrot_(int *n, double complex *x, int *incx, double complex *y, int *incy, double *c, double complex *s);

// 倍精度複素数エルミート行列とベクトルの積
// y <= alpha*A*x + beta*y
void zgemv_(char *trans, int *m, int *n,
            double complex *alpha, double complex *A, int *ldA, double complex *x, int *incx,
            double complex *beta , double complex *y, int *incy);
void zhpmv_(char *uplo, int *n, double complex *alpha, double complex *A,
            double complex *x, int *incx, double complex *beta, double complex *y, int *incy);


#endif
