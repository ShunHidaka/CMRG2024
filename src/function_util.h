#ifndef FUNCTION_UTIL_H
#define FUNCTION_UTIL_H

// 行列ファイル名
extern char *ANAME, *BNAME;
// ループの収束判定用定数
extern double EPS_INNER;
extern double EPS_OUTER;
// 残差の出力の間隔用定数
extern int OUTPUT_J;
extern int OUTPUT_K;
// メインループの最大値
extern int MAX_ITR;
// BLAS用定数
extern int ZERO;
extern int ONE;
extern double dTMP;
extern double complex cTMP;
// M_PIが未定義の場合定義する
#ifndef M_PI
#define M_PI 3.14159265358979
#endif

// Ax=bのベクトルbを作成する
void set_rhsVector(int N, double complex **b, double *norm);

// M個の純虚数を要素に持つ配列sigmaを作成する
void set_shifts(int *M, double complex **sigma);

// CSR形式による疎行列ベクトル積
void SpMV(const int *A_row, const int *A_col, const double complex *A_ele,
          const double complex *x, double complex *b, int N);

// 返り値が0ならば収束していない
int CG_method(const int *B_row, const int *B_col, const double complex *B_ele,
              double complex *x, const double complex *b,
              int N, const double tol);


// mtxファイルfnameを開き、コメントアウト部分を読み飛ばし、FILEポインタを返す。
// 行列の行数をrow_sizeに、列数をcol_sizeに、非ゼロ要素数をele_sizeに格納する。
FILE* fopen_mtx(const char *fname, const char *mode,
                int *row_size, int *col_size, int *ele_size);

// CSR形式のデータが格納されたファイル"fname"を開き、
// 行数N, 非ゼロ要素数DATASIZE, CRS形式の行列を返す。
void read_csr(const char *fname,
              int *N, int *DATASIZE,
              int **row_ptr, int **col_ind, double complex **element);

#endif
