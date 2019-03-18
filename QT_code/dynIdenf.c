/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: dynIdenf.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 24-Oct-2018 21:35:46
 */

/* Include Files */
#include <math.h>
#include <string.h>
#include "dynIdenf.h"
#include "dynIdenf_emxutil.h"

/* Variable Definitions */
static emxArray_real_T *a;
static emxArray_real_T *alpha;
static emxArray_real_T *d;
static double g[3];
static emxArray_real_T *phi_r0;
static double pfilt;
static double noise_err;
static double cond_max;
static double lambda;
static emxArray_real_T *segErr;
static double count;
static bool count_not_empty;
static emxArray_real_T *phi_pre;
static emxArray_real_T *setResample;
static emxArray_real_T *setEvalTorq;
static emxArray_real_T *phi_r;
static emxArray_real_T *Rpk;
static emxArray_real_T *Rphi;
static emxArray_real_T *num;
static double den;

/* Function Declarations */
static void assertValidIndexArg(const emxArray_real_T *s, emxArray_int32_T *sint);
static void b_abs(const emxArray_real_T *x, emxArray_real_T *y);
static void b_eye(double varargin_1, emxArray_real_T *I);
static void b_heapsort(emxArray_int32_T *x, int xstart, int xend, const
  cell_wrap_4 cmp_tunableEnvironment[2]);
static void b_rdivide(const emxArray_real_T *x, const emxArray_real_T *y,
                      emxArray_real_T *z);
static void b_rms(const emxArray_real_T *x, emxArray_real_T *y);
static void b_sqrt(double *x);
static void b_xaxpy(int n, double b_a, const emxArray_real_T *x, int ix0,
                    emxArray_real_T *y, int iy0);
static double b_xnrm2(int n, const emxArray_real_T *x, int ix0);
static void c_rdivide(const emxArray_real_T *x, const emxArray_real_T *y,
                      emxArray_real_T *z);
static void c_xaxpy(int n, double b_a, const emxArray_real_T *x, int ix0,
                    emxArray_real_T *y, int iy0);
static double combineVectorElements(const emxArray_real_T *x);
static void diag(const emxArray_real_T *v, emxArray_real_T *b_d);
static int div_s32_floor(int numerator, int denominator);
static void do_vectors(const emxArray_real_T *b_a, const emxArray_real_T *b,
  emxArray_real_T *c, emxArray_int32_T *ia, emxArray_int32_T *ib);
static void dynIdenf_free(void);
static void dynIdenf_init(void);
static void eml_signed_integer_colon(int b, emxArray_int32_T *y);
static void evalTorque(const emxArray_real_T *K, const emxArray_real_T *phi,
  const emxArray_real_T *b_phi_pre, const emxArray_real_T *tau_filt, const
  emxArray_real_T *b_segErr, emxArray_real_T *tau_pos, emxArray_real_T *tau_pre,
  emxArray_real_T *distErr, emxArray_real_T *convErr);
static void eye(double varargin_1, double varargin_2, emxArray_real_T *I);
static void filter(emxArray_real_T *b, emxArray_real_T *b_a, const
                   emxArray_real_T *x, const emxArray_real_T *zi,
                   emxArray_real_T *y);
static void fordKinematics(const emxArray_real_T *theta, const emxArray_real_T
  *theta_dot, const emxArray_real_T *theta_ddot, const double b_g[3], const
  emxArray_real_T *b_a, const emxArray_real_T *b_alpha, const emxArray_real_T
  *b_d, const emxArray_real_T *b_Rpk, emxArray_real_T *K);
static void getCoeffs(double fpass, emxArray_real_T *b_num, double *b_den);
static void heapify(emxArray_int32_T *x, int idx, int xstart, int xend, const
                    cell_wrap_4 cmp_tunableEnvironment[2]);
static void insertionsort(emxArray_int32_T *x, int xstart, int xend, const
  cell_wrap_4 cmp_tunableEnvironment[2]);
static void introsort(emxArray_int32_T *x, int xend, const cell_wrap_4
                      cmp_tunableEnvironment[2]);
static void lsqSVD(emxArray_real_T *K, emxArray_real_T *tau, const
                   emxArray_real_T *b_phi_pre, double b_count, const
                   emxArray_real_T *b_phi_r, double b_noise_err, double
                   b_cond_max, double b_lambda, emxArray_real_T *phi, double *rk);
static void merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int np,
                  int nq, emxArray_int32_T *iwork, emxArray_real_T *xwork);
static void merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
  int n, int preSortLevel, emxArray_int32_T *iwork, emxArray_real_T *xwork);
static void mldivide(const emxArray_real_T *A, emxArray_real_T *B);
static void mrdivide(emxArray_real_T *A, const emxArray_real_T *B);
static void myfiltfilt(emxArray_real_T *b, double b_a, const emxArray_real_T *x,
  emxArray_real_T *y_out);
static void permuteVector(const emxArray_int32_T *idx, emxArray_int32_T *y);
static int rankFromQR(const emxArray_real_T *A);
static void rdivide(const emxArray_real_T *x, const emxArray_real_T *y,
                    emxArray_real_T *z);
static void regroup(bool isLevel, double nparJoint, double nparMinSet, const
                    emxArray_real_T *b_a, const emxArray_real_T *b_alpha, const
                    emxArray_real_T *b_d, emxArray_real_T *b_Rpk,
                    emxArray_real_T *b_Rphi);
static void repmat(const emxArray_real_T *b_a, double varargin_1,
                   emxArray_real_T *b);
static double rms(const emxArray_real_T *x);
static double rt_hypotd(double u0, double u1);
static double skip_to_last_equal_value(int *k, const emxArray_real_T *x);
static void sort(emxArray_real_T *x);
static void sortIdx(emxArray_real_T *x, emxArray_int32_T *idx);
static void sparse(const emxArray_real_T *varargin_1, const emxArray_real_T
                   *varargin_2, const emxArray_real_T *varargin_3,
                   coder_internal_sparse *y);
static void sparse_full(const emxArray_real_T *this_d, const emxArray_int32_T
  *this_colidx, const emxArray_int32_T *this_rowidx, int this_m, int this_n,
  emxArray_real_T *y);
static void svd(const emxArray_real_T *A, emxArray_real_T *U, emxArray_real_T *s,
                emxArray_real_T *V);
static void vecnorm(const emxArray_real_T *x, emxArray_real_T *y);
static void xaxpy(int n, double b_a, int ix0, emxArray_real_T *y, int iy0);
static double xdotc(int n, const emxArray_real_T *x, int ix0, const
                    emxArray_real_T *y, int iy0);
static void xgeqp3(emxArray_real_T *A, emxArray_real_T *tau, emxArray_int32_T
                   *jpvt);
static void xgetrf(int m, int n, emxArray_real_T *A, int lda, emxArray_int32_T
                   *ipiv, int *info);
static double xnrm2(int n, const emxArray_real_T *x, int ix0);
static void xrot(int n, emxArray_real_T *x, int ix0, int iy0, double c, double s);
static void xrotg(double *b_a, double *b, double *c, double *s);
static void xscal(int n, double b_a, emxArray_real_T *x, int ix0);
static void xswap(int n, emxArray_real_T *x, int ix0, int iy0);

/* Function Definitions */

/*
 * Arguments    : const emxArray_real_T *s
 *                emxArray_int32_T *sint
 * Return Type  : void
 */
static void assertValidIndexArg(const emxArray_real_T *s, emxArray_int32_T *sint)
{
  int k;
  k = sint->size[0];
  sint->size[0] = s->size[1];
  emxEnsureCapacity_int32_T(sint, k);
  for (k = 0; k < s->size[1]; k++) {
    sint->data[k] = (int)s->data[k];
  }
}

/*
 * Arguments    : const emxArray_real_T *x
 *                emxArray_real_T *y
 * Return Type  : void
 */
static void b_abs(const emxArray_real_T *x, emxArray_real_T *y)
{
  int i17;
  unsigned int uv1[2];
  int k;
  int i18;
  int subs_idx_1;
  int b_k;
  for (i17 = 0; i17 < 2; i17++) {
    uv1[i17] = (unsigned int)x->size[i17];
  }

  i17 = y->size[0] * y->size[1];
  y->size[0] = (int)uv1[0];
  y->size[1] = (int)uv1[1];
  emxEnsureCapacity_real_T(y, i17);
  if (!((x->size[0] == 0) || (x->size[1] == 0))) {
    i17 = y->size[1];
    for (k = 1; k <= i17; k++) {
      i18 = y->size[0];
      if (1 <= i18) {
        subs_idx_1 = k;
      }

      for (b_k = 1; b_k <= i18; b_k++) {
        y->data[(b_k + y->size[0] * (subs_idx_1 - 1)) - 1] = fabs(x->data[(b_k +
          x->size[0] * (subs_idx_1 - 1)) - 1]);
      }
    }
  }
}

/*
 * Arguments    : double varargin_1
 *                emxArray_real_T *I
 * Return Type  : void
 */
static void b_eye(double varargin_1, emxArray_real_T *I)
{
  double t;
  int k;
  int loop_ub;
  int b_loop_ub;
  int i22;
  if (varargin_1 < 0.0) {
    t = 0.0;
  } else {
    t = varargin_1;
  }

  k = I->size[0] * I->size[1];
  I->size[0] = (int)t;
  I->size[1] = (int)t;
  emxEnsureCapacity_real_T(I, k);
  loop_ub = (int)t;
  for (k = 0; k < loop_ub; k++) {
    b_loop_ub = (int)t;
    for (i22 = 0; i22 < b_loop_ub; i22++) {
      I->data[i22 + I->size[0] * k] = 0.0;
    }
  }

  if ((int)t > 0) {
    for (k = 0; k < (int)t; k++) {
      I->data[k + I->size[0] * k] = 1.0;
    }
  }
}

/*
 * Arguments    : emxArray_int32_T *x
 *                int xstart
 *                int xend
 *                const cell_wrap_4 cmp_tunableEnvironment[2]
 * Return Type  : void
 */
static void b_heapsort(emxArray_int32_T *x, int xstart, int xend, const
  cell_wrap_4 cmp_tunableEnvironment[2])
{
  int n;
  int idx;
  int t;
  n = xend - xstart;
  for (idx = n + 1; idx > 0; idx--) {
    heapify(x, idx, xstart, xend, cmp_tunableEnvironment);
  }

  for (idx = 1; idx <= n; idx++) {
    t = x->data[xend - 1];
    x->data[xend - 1] = x->data[xstart - 1];
    x->data[xstart - 1] = t;
    xend--;
    heapify(x, 1, xstart, xend, cmp_tunableEnvironment);
  }
}

/*
 * Arguments    : const emxArray_real_T *x
 *                const emxArray_real_T *y
 *                emxArray_real_T *z
 * Return Type  : void
 */
static void b_rdivide(const emxArray_real_T *x, const emxArray_real_T *y,
                      emxArray_real_T *z)
{
  int i23;
  int loop_ub;
  i23 = z->size[0];
  z->size[0] = x->size[0];
  emxEnsureCapacity_real_T1(z, i23);
  loop_ub = x->size[0];
  for (i23 = 0; i23 < loop_ub; i23++) {
    z->data[i23] = x->data[i23] / y->data[i23];
  }
}

/*
 * Arguments    : const emxArray_real_T *x
 *                emxArray_real_T *y
 * Return Type  : void
 */
static void b_rms(const emxArray_real_T *x, emxArray_real_T *y)
{
  emxArray_real_T *b_x;
  int i26;
  int loop_ub;
  int vlen;
  int b_loop_ub;
  unsigned int sz[2];
  emxInit_real_T(&b_x, 2);
  i26 = b_x->size[0] * b_x->size[1];
  b_x->size[0] = x->size[0];
  b_x->size[1] = x->size[1];
  emxEnsureCapacity_real_T(b_x, i26);
  loop_ub = x->size[1];
  for (i26 = 0; i26 < loop_ub; i26++) {
    b_loop_ub = x->size[0];
    for (vlen = 0; vlen < b_loop_ub; vlen++) {
      b_x->data[vlen + b_x->size[0] * i26] = x->data[vlen + x->size[0] * i26] *
        x->data[vlen + x->size[0] * i26];
    }
  }

  vlen = b_x->size[0];
  if ((b_x->size[0] == 0) || (b_x->size[1] == 0)) {
    for (i26 = 0; i26 < 2; i26++) {
      sz[i26] = (unsigned int)b_x->size[i26];
    }

    i26 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = (int)sz[1];
    emxEnsureCapacity_real_T(y, i26);
    loop_ub = (int)sz[1];
    for (i26 = 0; i26 < loop_ub; i26++) {
      y->data[y->size[0] * i26] = 0.0;
    }
  } else {
    for (i26 = 0; i26 < 2; i26++) {
      sz[i26] = (unsigned int)b_x->size[i26];
    }

    i26 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = (int)sz[1];
    emxEnsureCapacity_real_T(y, i26);
    for (loop_ub = 0; loop_ub < b_x->size[1]; loop_ub++) {
      y->data[y->size[0] * loop_ub] = b_x->data[b_x->size[0] * loop_ub];
      for (b_loop_ub = 1; b_loop_ub < vlen; b_loop_ub++) {
        if (vlen >= 2) {
          y->data[y->size[0] * loop_ub] += b_x->data[b_loop_ub + b_x->size[0] *
            loop_ub];
        }
      }
    }
  }

  b_loop_ub = b_x->size[0];
  i26 = y->size[0] * y->size[1];
  y->size[0] = 1;
  emxEnsureCapacity_real_T(y, i26);
  loop_ub = y->size[1];
  emxFree_real_T(&b_x);
  for (i26 = 0; i26 < loop_ub; i26++) {
    y->data[y->size[0] * i26] /= (double)b_loop_ub;
  }

  i26 = y->size[1];
  for (loop_ub = 1; loop_ub <= i26; loop_ub++) {
    y->data[y->size[0] * (loop_ub - 1)] = sqrt(y->data[y->size[0] * (loop_ub - 1)]);
  }
}

/*
 * Arguments    : double *x
 * Return Type  : void
 */
static void b_sqrt(double *x)
{
  *x = sqrt(*x);
}

/*
 * Arguments    : int n
 *                double b_a
 *                const emxArray_real_T *x
 *                int ix0
 *                emxArray_real_T *y
 *                int iy0
 * Return Type  : void
 */
static void b_xaxpy(int n, double b_a, const emxArray_real_T *x, int ix0,
                    emxArray_real_T *y, int iy0)
{
  int ix;
  int iy;
  int k;
  if ((n < 1) || (b_a == 0.0)) {
  } else {
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k < n; k++) {
      y->data[iy] += b_a * x->data[ix];
      ix++;
      iy++;
    }
  }
}

/*
 * Arguments    : int n
 *                const emxArray_real_T *x
 *                int ix0
 * Return Type  : double
 */
static double b_xnrm2(int n, const emxArray_real_T *x, int ix0)
{
  double y;
  double scale;
  int kend;
  int k;
  double absxk;
  double t;
  y = 0.0;
  if (!(n < 1)) {
    if (n == 1) {
      y = fabs(x->data[ix0 - 1]);
    } else {
      scale = 3.3121686421112381E-170;
      kend = (ix0 + n) - 1;
      for (k = ix0; k <= kend; k++) {
        absxk = fabs(x->data[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          y = 1.0 + y * t * t;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * sqrt(y);
    }
  }

  return y;
}

/*
 * Arguments    : const emxArray_real_T *x
 *                const emxArray_real_T *y
 *                emxArray_real_T *z
 * Return Type  : void
 */
static void c_rdivide(const emxArray_real_T *x, const emxArray_real_T *y,
                      emxArray_real_T *z)
{
  int i27;
  int loop_ub;
  i27 = z->size[0] * z->size[1];
  z->size[0] = 1;
  z->size[1] = x->size[1];
  emxEnsureCapacity_real_T(z, i27);
  loop_ub = x->size[1];
  for (i27 = 0; i27 < loop_ub; i27++) {
    z->data[z->size[0] * i27] = x->data[x->size[0] * i27] / y->data[y->size[0] *
      i27];
  }
}

/*
 * Arguments    : int n
 *                double b_a
 *                const emxArray_real_T *x
 *                int ix0
 *                emxArray_real_T *y
 *                int iy0
 * Return Type  : void
 */
static void c_xaxpy(int n, double b_a, const emxArray_real_T *x, int ix0,
                    emxArray_real_T *y, int iy0)
{
  int ix;
  int iy;
  int k;
  if ((n < 1) || (b_a == 0.0)) {
  } else {
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k < n; k++) {
      y->data[iy] += b_a * x->data[ix];
      ix++;
      iy++;
    }
  }
}

/*
 * Arguments    : const emxArray_real_T *x
 * Return Type  : double
 */
static double combineVectorElements(const emxArray_real_T *x)
{
  double y;
  int vlen;
  int k;
  vlen = x->size[0];
  if (x->size[0] == 0) {
    y = 0.0;
  } else {
    y = x->data[0];
    for (k = 2; k <= vlen; k++) {
      if (vlen >= 2) {
        y += x->data[k - 1];
      }
    }
  }

  return y;
}

/*
 * Arguments    : const emxArray_real_T *v
 *                emxArray_real_T *b_d
 * Return Type  : void
 */
static void diag(const emxArray_real_T *v, emxArray_real_T *b_d)
{
  int u0;
  int u1;
  if ((v->size[0] == 1) && (v->size[1] == 1)) {
    u0 = b_d->size[0];
    b_d->size[0] = 1;
    emxEnsureCapacity_real_T1(b_d, u0);
    b_d->data[0] = v->data[0];
  } else {
    u0 = v->size[0];
    u1 = v->size[1];
    if (u0 < u1) {
      u1 = u0;
    }

    if (0 < v->size[1]) {
    } else {
      u1 = 0;
    }

    u0 = b_d->size[0];
    b_d->size[0] = u1;
    emxEnsureCapacity_real_T1(b_d, u0);
    for (u0 = 0; u0 < u1; u0++) {
      b_d->data[u0] = v->data[u0 + v->size[0] * u0];
    }
  }
}

/*
 * Arguments    : int numerator
 *                int denominator
 * Return Type  : int
 */
static int div_s32_floor(int numerator, int denominator)
{
  int quotient;
  unsigned int absNumerator;
  unsigned int absDenominator;
  bool quotientNeedsNegation;
  unsigned int tempAbsQuotient;
  if (denominator == 0) {
    if (numerator >= 0) {
      quotient = MAX_int32_T;
    } else {
      quotient = MIN_int32_T;
    }
  } else {
    if (numerator < 0) {
      absNumerator = ~(unsigned int)numerator + 1U;
    } else {
      absNumerator = (unsigned int)numerator;
    }

    if (denominator < 0) {
      absDenominator = ~(unsigned int)denominator + 1U;
    } else {
      absDenominator = (unsigned int)denominator;
    }

    quotientNeedsNegation = ((numerator < 0) != (denominator < 0));
    tempAbsQuotient = absNumerator / absDenominator;
    if (quotientNeedsNegation) {
      absNumerator %= absDenominator;
      if (absNumerator > 0U) {
        tempAbsQuotient++;
      }

      quotient = -(int)tempAbsQuotient;
    } else {
      quotient = (int)tempAbsQuotient;
    }
  }

  return quotient;
}

/*
 * Arguments    : const emxArray_real_T *b_a
 *                const emxArray_real_T *b
 *                emxArray_real_T *c
 *                emxArray_int32_T *ia
 *                emxArray_int32_T *ib
 * Return Type  : void
 */
static void do_vectors(const emxArray_real_T *b_a, const emxArray_real_T *b,
  emxArray_real_T *c, emxArray_int32_T *ia, emxArray_int32_T *ib)
{
  int na;
  int iafirst;
  int iblast;
  int nc;
  int nia;
  int ialast;
  int b_ialast;
  double ak;
  double bk;
  double absxk;
  emxArray_int32_T *b_ia;
  int exponent;
  na = b_a->size[1];
  iafirst = b_a->size[1];
  iblast = c->size[0] * c->size[1];
  c->size[0] = 1;
  c->size[1] = iafirst;
  emxEnsureCapacity_real_T(c, iblast);
  iblast = ia->size[0];
  ia->size[0] = b_a->size[1];
  emxEnsureCapacity_int32_T(ia, iblast);
  iblast = ib->size[0];
  ib->size[0] = 0;
  emxEnsureCapacity_int32_T(ib, iblast);
  nc = 0;
  nia = 0;
  iafirst = 0;
  ialast = 1;
  iblast = 1;
  while ((ialast <= na) && (iblast <= b->size[1])) {
    b_ialast = ialast;
    ak = skip_to_last_equal_value(&b_ialast, b_a);
    ialast = b_ialast;
    bk = skip_to_last_equal_value(&iblast, b);
    absxk = fabs(bk / 2.0);
    if (absxk <= 2.2250738585072014E-308) {
      absxk = 4.94065645841247E-324;
    } else {
      frexp(absxk, &exponent);
      absxk = ldexp(1.0, exponent - 53);
    }

    if (fabs(bk - ak) < absxk) {
      ialast = b_ialast + 1;
      iafirst = b_ialast;
      iblast++;
    } else if (ak < bk) {
      nc++;
      nia++;
      c->data[nc - 1] = ak;
      ia->data[nia - 1] = iafirst + 1;
      ialast = b_ialast + 1;
      iafirst = b_ialast;
    } else {
      iblast++;
    }
  }

  while (ialast <= na) {
    iafirst = ialast;
    ak = skip_to_last_equal_value(&iafirst, b_a);
    nc++;
    nia++;
    c->data[nc - 1] = ak;
    ia->data[nia - 1] = ialast;
    ialast = iafirst + 1;
  }

  if (b_a->size[1] > 0) {
    if (1 > nia) {
      iafirst = 0;
    } else {
      iafirst = nia;
    }

    emxInit_int32_T1(&b_ia, 1);
    iblast = b_ia->size[0];
    b_ia->size[0] = iafirst;
    emxEnsureCapacity_int32_T(b_ia, iblast);
    for (iblast = 0; iblast < iafirst; iblast++) {
      b_ia->data[iblast] = ia->data[iblast];
    }

    iblast = ia->size[0];
    ia->size[0] = b_ia->size[0];
    emxEnsureCapacity_int32_T(ia, iblast);
    iafirst = b_ia->size[0];
    for (iblast = 0; iblast < iafirst; iblast++) {
      ia->data[iblast] = b_ia->data[iblast];
    }

    emxFree_int32_T(&b_ia);
    iblast = c->size[0] * c->size[1];
    if (1 > nc) {
      c->size[1] = 0;
    } else {
      c->size[1] = nc;
    }

    emxEnsureCapacity_real_T(c, iblast);
  }
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void dynIdenf_free(void)
{
  emxFree_real_T(&a);
  emxFree_real_T(&alpha);
  emxFree_real_T(&d);
  emxFree_real_T(&phi_r0);
  emxFree_real_T(&segErr);
  emxFree_real_T(&phi_pre);
  emxFree_real_T(&setResample);
  emxFree_real_T(&setEvalTorq);
  emxFree_real_T(&phi_r);
  emxFree_real_T(&Rpk);
  emxFree_real_T(&Rphi);
  emxFree_real_T(&num);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void dynIdenf_init(void)
{
  emxInit_real_T(&num, 2);
  emxInit_real_T(&Rphi, 2);
  emxInit_real_T(&Rpk, 2);
  emxInit_real_T1(&phi_r, 1);
  emxInit_real_T(&setEvalTorq, 2);
  emxInit_real_T(&setResample, 2);
  emxInit_real_T1(&phi_pre, 1);
  emxInit_real_T(&segErr, 2);
  emxInit_real_T(&phi_r0, 2);
  emxInit_real_T(&d, 2);
  emxInit_real_T(&alpha, 2);
  emxInit_real_T(&a, 2);
}

/*
 * Arguments    : int b
 *                emxArray_int32_T *y
 * Return Type  : void
 */
static void eml_signed_integer_colon(int b, emxArray_int32_T *y)
{
  int n;
  int yk;
  int k;
  if (b < 1) {
    n = 0;
  } else {
    n = b;
  }

  yk = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = n;
  emxEnsureCapacity_int32_T1(y, yk);
  if (n > 0) {
    y->data[0] = 1;
    yk = 1;
    for (k = 2; k <= n; k++) {
      yk++;
      y->data[k - 1] = yk;
    }
  }
}

/*
 * function [tau_pos, tau_pre, distErr, convErr] = evalTorque(K, phi, phi_pre, tau_filt, segErr)
 * this function evaluate the torques predicted by the model and evaluate
 *  the error of torque
 * Arguments    : const emxArray_real_T *K
 *                const emxArray_real_T *phi
 *                const emxArray_real_T *b_phi_pre
 *                const emxArray_real_T *tau_filt
 *                const emxArray_real_T *b_segErr
 *                emxArray_real_T *tau_pos
 *                emxArray_real_T *tau_pre
 *                emxArray_real_T *distErr
 *                emxArray_real_T *convErr
 * Return Type  : void
 */
static void evalTorque(const emxArray_real_T *K, const emxArray_real_T *phi,
  const emxArray_real_T *b_phi_pre, const emxArray_real_T *tau_filt, const
  emxArray_real_T *b_segErr, emxArray_real_T *tau_pos, emxArray_real_T *tau_pre,
  emxArray_real_T *distErr, emxArray_real_T *convErr)
{
  int peval;
  emxArray_real_T *b_tau_pos;
  int m;
  int i24;
  int idx;
  int i;
  int k;
  int i25;
  emxArray_real_T *tau_err;
  emxArray_real_T *r9;
  emxArray_boolean_T *tau_flag;
  emxArray_int32_T *ii;
  double b_idx;
  int j;
  emxArray_real_T *tau_err_pos;
  emxArray_real_T *r10;
  emxArray_real_T *r11;
  bool exitg1;

  /*  inputs: */
  /*  K - the K matrix */
  /*  phi - the currently resolved dynamics paramters */
  /*  phi_pre - the previously resolved dynamics paramters */
  /*  tau_filt - the filtered torque */
  /*  segErr - the segmentation of the percentage error */
  /*  parameters: */
  /*  n - the number of joints */
  /*  outputs: */
  /*  tau_pre - the torque calculated by previously solved dynamic parameters */
  /*  tau_pos - the torque calculated by currently solved dynamic parameters */
  /*  distErr - the distribution of error */
  /*  convErr - the convergence error */
  /*  obtain the number of joints and the number of error segmentations */
  /* 'dynIdenf:512' [peval, n] = size(tau_filt); */
  peval = tau_filt->size[0];

  /* 'dynIdenf:513' nSeg = length(segErr); */
  /*  calculate tau_pos */
  /* 'dynIdenf:515' tau_pos = K * phi; */
  emxInit_real_T1(&b_tau_pos, 1);
  if ((K->size[1] == 1) || (phi->size[0] == 1)) {
    i24 = b_tau_pos->size[0];
    b_tau_pos->size[0] = K->size[0];
    emxEnsureCapacity_real_T1(b_tau_pos, i24);
    idx = K->size[0];
    for (i24 = 0; i24 < idx; i24++) {
      b_tau_pos->data[i24] = 0.0;
      m = K->size[1];
      for (i25 = 0; i25 < m; i25++) {
        b_tau_pos->data[i24] += K->data[i24 + K->size[0] * i25] * phi->data[i25];
      }
    }
  } else {
    m = K->size[0];
    i24 = b_tau_pos->size[0];
    b_tau_pos->size[0] = K->size[0];
    emxEnsureCapacity_real_T1(b_tau_pos, i24);
    for (i = 1; i <= m; i++) {
      b_tau_pos->data[i - 1] = 0.0;
    }

    for (k = 0; k < K->size[1]; k++) {
      if (phi->data[k] != 0.0) {
        for (i = 0; i < m; i++) {
          b_tau_pos->data[i] += phi->data[k] * K->data[i + K->size[0] * k];
        }
      }
    }
  }

  /* 'dynIdenf:516' tau_pos = reshape(tau_pos,n,peval)'; */
  m = tau_filt->size[1];
  idx = tau_filt->size[0];
  i24 = tau_pos->size[0] * tau_pos->size[1];
  tau_pos->size[0] = idx;
  tau_pos->size[1] = m;
  emxEnsureCapacity_real_T(tau_pos, i24);
  for (i24 = 0; i24 < m; i24++) {
    for (i25 = 0; i25 < idx; i25++) {
      tau_pos->data[i25 + tau_pos->size[0] * i24] = b_tau_pos->data[i24 + m *
        i25];
    }
  }

  /*  calculate tau_pre */
  /* 'dynIdenf:518' tau_pre = K * phi_pre; */
  if ((K->size[1] == 1) || (b_phi_pre->size[0] == 1)) {
    i24 = b_tau_pos->size[0];
    b_tau_pos->size[0] = K->size[0];
    emxEnsureCapacity_real_T1(b_tau_pos, i24);
    idx = K->size[0];
    for (i24 = 0; i24 < idx; i24++) {
      b_tau_pos->data[i24] = 0.0;
      m = K->size[1];
      for (i25 = 0; i25 < m; i25++) {
        b_tau_pos->data[i24] += K->data[i24 + K->size[0] * i25] *
          b_phi_pre->data[i25];
      }
    }
  } else {
    m = K->size[0];
    i24 = b_tau_pos->size[0];
    b_tau_pos->size[0] = K->size[0];
    emxEnsureCapacity_real_T1(b_tau_pos, i24);
    for (i = 1; i <= m; i++) {
      b_tau_pos->data[i - 1] = 0.0;
    }

    for (k = 0; k < K->size[1]; k++) {
      if (b_phi_pre->data[k] != 0.0) {
        for (i = 0; i < m; i++) {
          b_tau_pos->data[i] += b_phi_pre->data[k] * K->data[i + K->size[0] * k];
        }
      }
    }
  }

  /* 'dynIdenf:519' tau_pre = reshape(tau_pre,n,peval)'; */
  m = tau_filt->size[1];
  idx = tau_filt->size[0];
  i24 = tau_pre->size[0] * tau_pre->size[1];
  tau_pre->size[0] = idx;
  tau_pre->size[1] = m;
  emxEnsureCapacity_real_T(tau_pre, i24);
  for (i24 = 0; i24 < m; i24++) {
    for (i25 = 0; i25 < idx; i25++) {
      tau_pre->data[i25 + tau_pre->size[0] * i24] = b_tau_pos->data[i24 + m *
        i25];
    }
  }

  emxInit_real_T(&tau_err, 2);

  /*  calculate the distribution of prediction error */
  /* 'dynIdenf:521' distErr = zeros(nSeg,n); */
  /* 'dynIdenf:522' tau_err = sort(abs((tau_pos - tau_filt)./tau_filt),'descend'); */
  i24 = tau_err->size[0] * tau_err->size[1];
  tau_err->size[0] = tau_pos->size[0];
  tau_err->size[1] = tau_pos->size[1];
  emxEnsureCapacity_real_T(tau_err, i24);
  idx = tau_pos->size[1];
  for (i24 = 0; i24 < idx; i24++) {
    m = tau_pos->size[0];
    for (i25 = 0; i25 < m; i25++) {
      tau_err->data[i25 + tau_err->size[0] * i24] = tau_pos->data[i25 +
        tau_pos->size[0] * i24] - tau_filt->data[i25 + tau_filt->size[0] * i24];
    }
  }

  emxInit_real_T(&r9, 2);
  rdivide(tau_err, tau_filt, r9);
  b_abs(r9, tau_err);
  sort(tau_err);

  /* 'dynIdenf:523' for i = 1:nSeg */
  i24 = distErr->size[0] * distErr->size[1];
  distErr->size[0] = b_segErr->size[1];
  distErr->size[1] = tau_filt->size[1];
  emxEnsureCapacity_real_T(distErr, i24);
  i = 0;
  emxFree_real_T(&r9);
  emxInit_boolean_T1(&tau_flag, 2);
  emxInit_int32_T1(&ii, 1);
  while (i <= b_segErr->size[1] - 1) {
    /* 'dynIdenf:524' tau_flag = tau_err <= segErr(i); */
    b_idx = b_segErr->data[i];
    i24 = tau_flag->size[0] * tau_flag->size[1];
    tau_flag->size[0] = tau_err->size[0];
    tau_flag->size[1] = tau_err->size[1];
    emxEnsureCapacity_boolean_T(tau_flag, i24);
    idx = tau_err->size[1];
    for (i24 = 0; i24 < idx; i24++) {
      m = tau_err->size[0];
      for (i25 = 0; i25 < m; i25++) {
        tau_flag->data[i25 + tau_flag->size[0] * i24] = (tau_err->data[i25 +
          tau_err->size[0] * i24] <= b_idx);
      }
    }

    /* 'dynIdenf:525' for j = 1:n */
    for (j = 0; j < tau_filt->size[1]; j++) {
      /* 'dynIdenf:526' idx = sum(find(tau_flag(:,j),1,'first')); */
      i24 = tau_flag->size[0];
      i25 = tau_flag->size[0];
      if (1 <= i25) {
        k = 1;
      } else {
        k = tau_flag->size[0];
      }

      idx = 0;
      i25 = ii->size[0];
      ii->size[0] = k;
      emxEnsureCapacity_int32_T(ii, i25);
      m = 1;
      exitg1 = false;
      while ((!exitg1) && (m <= i24)) {
        if (tau_flag->data[(m + tau_flag->size[0] * j) - 1]) {
          idx++;
          ii->data[idx - 1] = m;
          if (idx >= k) {
            exitg1 = true;
          } else {
            m++;
          }
        } else {
          m++;
        }
      }

      if (k == 1) {
        if (idx == 0) {
          i24 = ii->size[0];
          ii->size[0] = 0;
          emxEnsureCapacity_int32_T(ii, i24);
        }
      } else {
        i24 = ii->size[0];
        if (1 > idx) {
          ii->size[0] = 0;
        } else {
          ii->size[0] = idx;
        }

        emxEnsureCapacity_int32_T(ii, i24);
      }

      i24 = b_tau_pos->size[0];
      b_tau_pos->size[0] = ii->size[0];
      emxEnsureCapacity_real_T1(b_tau_pos, i24);
      idx = ii->size[0];
      for (i24 = 0; i24 < idx; i24++) {
        b_tau_pos->data[i24] = ii->data[i24];
      }

      b_idx = combineVectorElements(b_tau_pos);

      /* 'dynIdenf:527' if isempty(idx) == 0 */
      /* 'dynIdenf:528' distErr(i,j) = (peval - idx + (idx~=peval)) / peval; */
      distErr->data[i + distErr->size[0] * j] = (((double)peval - b_idx) +
        (double)(b_idx != peval)) / (double)peval;
    }

    i++;
  }

  emxFree_int32_T(&ii);
  emxFree_real_T(&b_tau_pos);
  emxFree_boolean_T(&tau_flag);

  /*  calculate the convergence error */
  /* 'dynIdenf:533' tau_err_pos = rms(tau_pos - tau_filt) ./ rms(tau_filt); */
  i24 = tau_err->size[0] * tau_err->size[1];
  tau_err->size[0] = tau_pos->size[0];
  tau_err->size[1] = tau_pos->size[1];
  emxEnsureCapacity_real_T(tau_err, i24);
  idx = tau_pos->size[1];
  for (i24 = 0; i24 < idx; i24++) {
    m = tau_pos->size[0];
    for (i25 = 0; i25 < m; i25++) {
      tau_err->data[i25 + tau_err->size[0] * i24] = tau_pos->data[i25 +
        tau_pos->size[0] * i24] - tau_filt->data[i25 + tau_filt->size[0] * i24];
    }
  }

  emxInit_real_T(&tau_err_pos, 2);
  emxInit_real_T(&r10, 2);
  emxInit_real_T(&r11, 2);
  b_rms(tau_err, r10);
  b_rms(tau_filt, r11);
  c_rdivide(r10, r11, tau_err_pos);

  /* 'dynIdenf:534' tau_err_pre = rms(tau_pre - tau_filt) ./ rms(tau_filt); */
  i24 = tau_err->size[0] * tau_err->size[1];
  tau_err->size[0] = tau_pre->size[0];
  tau_err->size[1] = tau_pre->size[1];
  emxEnsureCapacity_real_T(tau_err, i24);
  idx = tau_pre->size[1];
  for (i24 = 0; i24 < idx; i24++) {
    m = tau_pre->size[0];
    for (i25 = 0; i25 < m; i25++) {
      tau_err->data[i25 + tau_err->size[0] * i24] = tau_pre->data[i25 +
        tau_pre->size[0] * i24] - tau_filt->data[i25 + tau_filt->size[0] * i24];
    }
  }

  b_rms(tau_err, r10);
  b_rms(tau_filt, r11);
  c_rdivide(r10, r11, convErr);

  /* 'dynIdenf:535' convErr = tau_err_pre - tau_err_pos; */
  i24 = convErr->size[0] * convErr->size[1];
  convErr->size[0] = 1;
  emxEnsureCapacity_real_T(convErr, i24);
  idx = convErr->size[1];
  emxFree_real_T(&r11);
  emxFree_real_T(&r10);
  emxFree_real_T(&tau_err);
  for (i24 = 0; i24 < idx; i24++) {
    convErr->data[convErr->size[0] * i24] -= tau_err_pos->data[tau_err_pos->
      size[0] * i24];
  }

  emxFree_real_T(&tau_err_pos);
}

/*
 * Arguments    : double varargin_1
 *                double varargin_2
 *                emxArray_real_T *I
 * Return Type  : void
 */
static void eye(double varargin_1, double varargin_2, emxArray_real_T *I)
{
  int m;
  double t;
  int b_d;
  int i4;
  int loop_ub;
  int i5;
  if (varargin_1 < 0.0) {
    m = 0;
  } else {
    m = (int)varargin_1;
  }

  if (varargin_2 < 0.0) {
    t = 0.0;
  } else {
    t = varargin_2;
  }

  if (m <= (int)t) {
    b_d = m;
  } else {
    b_d = (int)t;
  }

  i4 = I->size[0] * I->size[1];
  I->size[0] = m;
  I->size[1] = (int)t;
  emxEnsureCapacity_real_T(I, i4);
  loop_ub = (int)t;
  for (i4 = 0; i4 < loop_ub; i4++) {
    for (i5 = 0; i5 < m; i5++) {
      I->data[i5 + I->size[0] * i4] = 0.0;
    }
  }

  if (b_d > 0) {
    for (m = 0; m < b_d; m++) {
      I->data[m + I->size[0] * m] = 1.0;
    }
  }
}

/*
 * Arguments    : emxArray_real_T *b
 *                emxArray_real_T *b_a
 *                const emxArray_real_T *x
 *                const emxArray_real_T *zi
 *                emxArray_real_T *y
 * Return Type  : void
 */
static void filter(emxArray_real_T *b, emxArray_real_T *b_a, const
                   emxArray_real_T *x, const emxArray_real_T *zi,
                   emxArray_real_T *y)
{
  int na;
  int nb;
  int u0;
  int ndbuffer;
  signed char zicase;
  double a1;
  int k;
  unsigned int b_x[2];
  unsigned int uv0[2];
  int nx;
  na = b_a->size[1];
  nb = b->size[1];
  u0 = b_a->size[1];
  ndbuffer = b->size[1];
  if (u0 > ndbuffer) {
    ndbuffer = u0;
  }

  if (zi->size[0] == ndbuffer - 1) {
    zicase = 1;
  } else {
    zicase = 2;
  }

  a1 = b_a->data[0];
  if ((!(b_a->data[0] == 0.0)) && (b_a->data[0] != 1.0)) {
    for (k = 0; k < nb; k++) {
      b->data[k] /= a1;
    }

    for (k = 1; k < na; k++) {
      b_a->data[k] /= a1;
    }

    b_a->data[0] = 1.0;
  }

  b_x[0] = (unsigned int)x->size[0];
  b_x[1] = 1U;
  for (u0 = 0; u0 < 2; u0++) {
    uv0[u0] = b_x[u0];
  }

  u0 = y->size[0];
  y->size[0] = (int)uv0[0];
  emxEnsureCapacity_real_T1(y, u0);
  nx = x->size[0];
  u0 = x->size[0];
  ndbuffer--;
  if (u0 < ndbuffer) {
    ndbuffer = u0;
  }

  if (zicase == 2) {
    for (k = 0; k < ndbuffer; k++) {
      y->data[k] = zi->data[k];
    }

    for (k = ndbuffer + 1; k <= nx; k++) {
      y->data[k - 1] = 0.0;
    }
  } else {
    for (k = 0; k < ndbuffer; k++) {
      y->data[k] = zi->data[k];
    }

    for (k = ndbuffer + 1; k <= nx; k++) {
      y->data[k - 1] = 0.0;
    }
  }

  if ((na == 1) && (x->size[0] >= (nb << 1))) {
    for (k = 1; k <= nb; k++) {
      for (u0 = k; u0 <= nx; u0++) {
        y->data[u0 - 1] += b->data[k - 1] * x->data[u0 - k];
      }
    }
  } else {
    for (k = 0; k < nx; k++) {
      ndbuffer = nx - k;
      if (!(ndbuffer < nb)) {
        ndbuffer = nb;
      }

      for (u0 = 0; u0 < ndbuffer; u0++) {
        y->data[k + u0] += x->data[k] * b->data[u0];
      }

      u0 = (nx - k) - 1;
      ndbuffer = na - 1;
      if (u0 < ndbuffer) {
        ndbuffer = u0;
      }

      a1 = -y->data[k];
      for (u0 = 1; u0 <= ndbuffer; u0++) {
        y->data[k + u0] += a1 * b_a->data[u0];
      }
    }
  }
}

/*
 * function K = fordKinematics(theta, theta_dot, theta_ddot,...
 *     g, a, alpha, d, Rpk)
 * this function calculates the K kinematics, noting that the
 *  modified D-H model is used
 * Arguments    : const emxArray_real_T *theta
 *                const emxArray_real_T *theta_dot
 *                const emxArray_real_T *theta_ddot
 *                const double b_g[3]
 *                const emxArray_real_T *b_a
 *                const emxArray_real_T *b_alpha
 *                const emxArray_real_T *b_d
 *                const emxArray_real_T *b_Rpk
 *                emxArray_real_T *K
 * Return Type  : void
 */
static void fordKinematics(const emxArray_real_T *theta, const emxArray_real_T
  *theta_dot, const emxArray_real_T *theta_ddot, const double b_g[3], const
  emxArray_real_T *b_a, const emxArray_real_T *b_alpha, const emxArray_real_T
  *b_d, const emxArray_real_T *b_Rpk, emxArray_real_T *K)
{
  int n;
  int i;
  int m;
  int i12;
  int i13;
  int k;
  emxArray_real_T *o_X_i;
  emxArray_int32_T *r5;
  emxArray_real_T *c_a;
  int loop_ub;
  int i14;
  int inner;
  double i_X_o[6][6];
  int b_i;
  double i_V_i[6];
  double i_A_i[6];
  int j;
  double i1_R_i[3][3];
  double i1_p_i[3];
  double dv11[3][3];
  double dv12[3][3];
  double i1_X_i[6][6];
  double i_X_i1[6][6];
  double b_i1_R_i[3][3];
  double i_V_i1_i[6];
  double i_A_i1_i[6];
  double b_o_X_i[6][6];
  double c_o_X_i;
  double i_V_i1[6];
  double b_i_X_i1[6][6];
  double c_i_X_i1[6];
  double dv13[6];
  double b_i_V_i[6][3];
  double b_i_A_i[6][3];
  double dv14[3][3];
  double i_Am_i[10][6];
  double b_i_X_o[3];
  double dv15[6][3];
  double c_i_X_o;
  double d_i_X_o[12];
  double e_i_X_o[10][6];

  /*  inputs: */
  /*  theta - the angular position of joints */
  /*  theta_dot - the angular velocities of joints */
  /*  theta_ddot - the angular acceleration of joints */
  /*  parameters: */
  /*  g - the gravity */
  /*  a, alpha, d - the D-H parameters */
  /*  Rpk - the matrix for regrouping K */
  /*  outputs */
  /*  K - the matrix that makes tau = K*phi. */
  /*  obtain dimension of data */
  /* 'dynIdenf:280' [p, n] = size(theta); */
  n = theta->size[1];

  /* 'dynIdenf:281' m = length(g); */
  /* 'dynIdenf:282' nparJoint = (m^2+3*m+6)/2; */
  /*  initialize K */
  /* 'dynIdenf:284' K = zeros(n*p,n*nparJoint); */
  i = (int)((double)theta->size[1] * (double)theta->size[0]);
  m = (int)((double)theta->size[1] * 12.0);
  i12 = K->size[0] * K->size[1];
  K->size[0] = i;
  K->size[1] = m;
  emxEnsureCapacity_real_T(K, i12);
  for (i12 = 0; i12 < m; i12++) {
    for (i13 = 0; i13 < i; i13++) {
      K->data[i13 + K->size[0] * i12] = 0.0;
    }
  }

  /*  calculate the K */
  /*  the outer-most loop indicates the number of sample points */
  /* 'dynIdenf:287' for k = 1:p */
  k = 0;
  emxInit_real_T2(&o_X_i, 3);
  emxInit_int32_T1(&r5, 1);
  while (k <= theta->size[0] - 1) {
    /*  the transformation matrix from reference i to reference 0, */
    /*  and its reverse */
    /* 'dynIdenf:290' o_X_i = zeros(2*m,2*m,n); */
    i12 = o_X_i->size[0] * o_X_i->size[1] * o_X_i->size[2];
    o_X_i->size[0] = 6;
    o_X_i->size[1] = 6;
    o_X_i->size[2] = n;
    emxEnsureCapacity_real_T2(o_X_i, i12);
    for (i12 = 0; i12 < n; i12++) {
      for (i13 = 0; i13 < 6; i13++) {
        for (i14 = 0; i14 < 6; i14++) {
          o_X_i->data[(i14 + o_X_i->size[0] * i13) + o_X_i->size[0] *
            o_X_i->size[1] * i12] = 0.0;
        }
      }
    }

    /* 'dynIdenf:291' i_X_o = zeros(2*m,2*m); */
    /*  the absolute velocity and acceleration of reference i with components in */
    /*  itself */
    /* 'dynIdenf:294' i_V_i = zeros(2*m,1); */
    /* 'dynIdenf:295' i_A_i = zeros(2*m,1); */
    memset(&i_X_o[0][0], 0, 36U * sizeof(double));
    for (b_i = 0; b_i < 6; b_i++) {
      i_V_i[b_i] = 0.0;
      i_A_i[b_i] = 0.0;
    }

    /*  the secon loop indicates the number of joints */
    /* 'dynIdenf:297' for i = 1:n */
    for (b_i = 0; b_i < n; b_i++) {
      /*  the rotation matrix from reference i to i-1 */
      /* 'dynIdenf:299' i1_R_i = [cos(theta(k,i)), -sin(theta(k,i)), 0;... */
      /* 'dynIdenf:300'             sin(theta(k,i))*cos(alpha(i)), cos(theta(k,i))*cos(alpha(i)), -sin(alpha(i));... */
      /* 'dynIdenf:301'             sin(theta(k,i))*sin(alpha(i)), cos(theta(k,i))*sin(alpha(i)), cos(alpha(i))]; */
      i1_R_i[0][0] = cos(theta->data[k + theta->size[0] * b_i]);
      i1_R_i[1][0] = -sin(theta->data[k + theta->size[0] * b_i]);
      i1_R_i[2][0] = 0.0;
      i1_R_i[0][1] = sin(theta->data[k + theta->size[0] * b_i]) * cos
        (b_alpha->data[b_i]);
      i1_R_i[1][1] = cos(theta->data[k + theta->size[0] * b_i]) * cos
        (b_alpha->data[b_i]);
      i1_R_i[2][1] = -sin(b_alpha->data[b_i]);
      i1_R_i[0][2] = sin(theta->data[k + theta->size[0] * b_i]) * sin
        (b_alpha->data[b_i]);
      i1_R_i[1][2] = cos(theta->data[k + theta->size[0] * b_i]) * sin
        (b_alpha->data[b_i]);
      i1_R_i[2][2] = cos(b_alpha->data[b_i]);

      /*  the offset vector from reference i to i-1 */
      /* 'dynIdenf:303' i1_p_i = [a(i); -d(i+1)*sin(alpha(i)); d(i+1)*cos(alpha(i))]; */
      i1_p_i[0] = b_a->data[b_i];
      i1_p_i[1] = -b_d->data[b_i + 1] * sin(b_alpha->data[b_i]);
      i1_p_i[2] = b_d->data[b_i + 1] * cos(b_alpha->data[b_i]);

      /*  the spatial form of the homogeneous transformation */
      /*  matrix, and its reverse */
      /* 'dynIdenf:306' i1_X_i = [i1_R_i, zeros(m,m); skew(i1_p_i)*i1_R_i, i1_R_i]; */
      /*  transform a R3 vector into a skew matrix */
      /* 'dynIdenf:352' y = [0, -x(3), x(2); x(3), 0, -x(1); -x(2) x(1) 0]; */
      dv11[0][0] = 0.0;
      dv11[1][0] = -i1_p_i[2];
      dv11[2][0] = i1_p_i[1];
      dv11[0][1] = i1_p_i[2];
      dv11[1][1] = 0.0;
      dv11[2][1] = -i1_p_i[0];
      dv11[0][2] = -i1_p_i[1];
      dv11[1][2] = i1_p_i[0];
      dv11[2][2] = 0.0;
      for (i12 = 0; i12 < 3; i12++) {
        for (i13 = 0; i13 < 3; i13++) {
          dv12[i13][i12] = 0.0;
          for (i14 = 0; i14 < 3; i14++) {
            dv12[i13][i12] += dv11[i14][i12] * i1_R_i[i13][i14];
          }

          i1_X_i[i12][i13] = i1_R_i[i12][i13];
          i1_X_i[i12 + 3][i13] = 0.0;
        }
      }

      /* 'dynIdenf:307' i_X_i1 = [i1_R_i', zeros(m,m); i1_R_i'*skew(i1_p_i)', i1_R_i']; */
      /*  transform a R3 vector into a skew matrix */
      /* 'dynIdenf:352' y = [0, -x(3), x(2); x(3), 0, -x(1); -x(2) x(1) 0]; */
      dv11[0][0] = 0.0;
      dv11[0][1] = -i1_p_i[2];
      dv11[0][2] = i1_p_i[1];
      dv11[1][0] = i1_p_i[2];
      dv11[1][1] = 0.0;
      dv11[1][2] = -i1_p_i[0];
      dv11[2][0] = -i1_p_i[1];
      dv11[2][1] = i1_p_i[0];
      dv11[2][2] = 0.0;
      for (i12 = 0; i12 < 3; i12++) {
        for (i13 = 0; i13 < 3; i13++) {
          i1_X_i[i12][i13 + 3] = dv12[i12][i13];
          i1_X_i[i12 + 3][i13 + 3] = i1_R_i[i12][i13];
          b_i1_R_i[i13][i12] = 0.0;
          for (i14 = 0; i14 < 3; i14++) {
            b_i1_R_i[i13][i12] += i1_R_i[i12][i14] * dv11[i13][i14];
          }

          i_X_i1[i12][i13] = i1_R_i[i13][i12];
          i_X_i1[i12 + 3][i13] = 0.0;
        }
      }

      for (i12 = 0; i12 < 3; i12++) {
        for (i13 = 0; i13 < 3; i13++) {
          i_X_i1[i12][i13 + 3] = b_i1_R_i[i12][i13];
          i_X_i1[i12 + 3][i13 + 3] = i1_R_i[i13][i12];
        }
      }

      /*  the relative velocity and acceleration of reference */
      /*  i to i-1 with components in reference i */
      /* 'dynIdenf:310' i_V_i1_i = [zeros(m-1,1); theta_dot(k,i); zeros(m,1)]; */
      for (i12 = 0; i12 < 2; i12++) {
        i_V_i1_i[i12] = 0.0;
      }

      i_V_i1_i[2] = theta_dot->data[k + theta_dot->size[0] * b_i];
      for (i12 = 0; i12 < 3; i12++) {
        i_V_i1_i[i12 + 3] = 0.0;
      }

      /* 'dynIdenf:311' i_A_i1_i = [zeros(m-1,1); theta_ddot(k,i); zeros(m,1)]; */
      for (i12 = 0; i12 < 2; i12++) {
        i_A_i1_i[i12] = 0.0;
      }

      i_A_i1_i[2] = theta_ddot->data[k + theta_ddot->size[0] * b_i];
      for (i12 = 0; i12 < 3; i12++) {
        i_A_i1_i[i12 + 3] = 0.0;
      }

      /*  if it is the frist joint, the abolute transformation, */
      /*  velocity and acceleration are equal to the relative ones */
      /* 'dynIdenf:314' if i == 1 */
      if (1 + b_i == 1) {
        /* 'dynIdenf:315' o_X_i(:,:,i) = i1_X_i; */
        /* 'dynIdenf:316' i_X_o = i_X_i1; */
        /* 'dynIdenf:317' i_V_i = i_V_i1_i; */
        /* 'dynIdenf:318' i_A_i = i_A_i1_i; */
        for (i = 0; i < 6; i++) {
          for (i12 = 0; i12 < 6; i12++) {
            o_X_i->data[(i12 + o_X_i->size[0] * i) + o_X_i->size[0] *
              o_X_i->size[1] * b_i] = i1_X_i[i][i12];
            i_X_o[i][i12] = i_X_i1[i][i12];
          }

          i_V_i[i] = i_V_i1_i[i];
          i_A_i[i] = i_A_i1_i[i];
        }

        /*  otherwise, obtain the absolute transformation, */
        /*  velocity and acceleration by iteration */
      } else {
        /* 'dynIdenf:321' else */
        /* 'dynIdenf:322' o_X_i(:,:,i) = o_X_i(:,:,i-1) * i1_X_i; */
        for (i12 = 0; i12 < 6; i12++) {
          for (i13 = 0; i13 < 6; i13++) {
            b_o_X_i[i13][i12] = 0.0;
            for (i14 = 0; i14 < 6; i14++) {
              c_o_X_i = b_o_X_i[i13][i12] + o_X_i->data[(i12 + o_X_i->size[0] *
                i14) + o_X_i->size[0] * o_X_i->size[1] * (b_i - 1)] * i1_X_i[i13]
                [i14];
              b_o_X_i[i13][i12] = c_o_X_i;
            }
          }
        }

        /* 'dynIdenf:323' i_X_o = i_X_i1 * i_X_o; */
        for (i12 = 0; i12 < 6; i12++) {
          for (i13 = 0; i13 < 6; i13++) {
            o_X_i->data[(i13 + o_X_i->size[0] * i12) + o_X_i->size[0] *
              o_X_i->size[1] * b_i] = b_o_X_i[i12][i13];
            b_i_X_i1[i13][i12] = 0.0;
            for (i14 = 0; i14 < 6; i14++) {
              b_i_X_i1[i13][i12] += i_X_i1[i14][i12] * i_X_o[i13][i14];
            }
          }
        }

        /* 'dynIdenf:324' i_V_i1 = i_X_i1 * i_V_i; */
        for (i12 = 0; i12 < 6; i12++) {
          i_V_i1[i12] = 0.0;
          for (i13 = 0; i13 < 6; i13++) {
            i_X_o[i12][i13] = b_i_X_i1[i12][i13];
            i_V_i1[i12] += i_X_i1[i13][i12] * i_V_i[i13];
          }
        }

        /* 'dynIdenf:325' i_V_i = i_V_i1 + i_V_i1_i; */
        for (i = 0; i < 6; i++) {
          i_V_i[i] = i_V_i1[i] + i_V_i1_i[i];
        }

        /* 'dynIdenf:326' i_A_i = i_X_i1 * i_A_i - skew6(i_V_i1_i) * i_V_i1 + i_A_i1_i; */
        /*  transform a R6 vector into a skew matrix */
        /* 'dynIdenf:356' y = [skew(x(1:3)), zeros(3,3); skew(x(4:6)), skew(x(1:3))]; */
        /*  transform a R3 vector into a skew matrix */
        /* 'dynIdenf:352' y = [0, -x(3), x(2); x(3), 0, -x(1); -x(2) x(1) 0]; */
        /*  transform a R3 vector into a skew matrix */
        /* 'dynIdenf:352' y = [0, -x(3), x(2); x(3), 0, -x(1); -x(2) x(1) 0]; */
        /*  transform a R3 vector into a skew matrix */
        /* 'dynIdenf:352' y = [0, -x(3), x(2); x(3), 0, -x(1); -x(2) x(1) 0]; */
        b_o_X_i[0][0] = 0.0;
        b_o_X_i[1][0] = -i_V_i1_i[2];
        b_o_X_i[2][0] = 0.0;
        b_o_X_i[0][1] = i_V_i1_i[2];
        b_o_X_i[1][1] = 0.0;
        b_o_X_i[2][1] = -0.0;
        b_o_X_i[0][2] = -0.0;
        b_o_X_i[1][2] = 0.0;
        b_o_X_i[2][2] = 0.0;
        for (i12 = 0; i12 < 3; i12++) {
          for (i13 = 0; i13 < 3; i13++) {
            b_o_X_i[i12 + 3][i13] = 0.0;
          }
        }

        b_o_X_i[0][3] = 0.0;
        b_o_X_i[1][3] = -0.0;
        b_o_X_i[2][3] = 0.0;
        b_o_X_i[0][4] = 0.0;
        b_o_X_i[1][4] = 0.0;
        b_o_X_i[2][4] = -0.0;
        b_o_X_i[0][5] = -0.0;
        b_o_X_i[1][5] = 0.0;
        b_o_X_i[2][5] = 0.0;
        b_o_X_i[3][3] = 0.0;
        b_o_X_i[4][3] = -i_V_i1_i[2];
        b_o_X_i[5][3] = 0.0;
        b_o_X_i[3][4] = i_V_i1_i[2];
        b_o_X_i[4][4] = 0.0;
        b_o_X_i[5][4] = -0.0;
        b_o_X_i[3][5] = -0.0;
        b_o_X_i[4][5] = 0.0;
        b_o_X_i[5][5] = 0.0;
        for (i12 = 0; i12 < 6; i12++) {
          c_i_X_i1[i12] = 0.0;
          dv13[i12] = 0.0;
          for (i13 = 0; i13 < 6; i13++) {
            c_i_X_i1[i12] += i_X_i1[i13][i12] * i_A_i[i13];
            dv13[i12] += b_o_X_i[i13][i12] * i_V_i1[i13];
          }
        }

        for (i12 = 0; i12 < 6; i12++) {
          i_A_i[i12] = (c_i_X_i1[i12] - dv13[i12]) + i_A_i1_i[i12];
        }
      }

      /*  the absolute angular velocity and acceleration (different from classical one) of reference i with */
      /*  components in itself */
      /* 'dynIdenf:330' i_omega_i = i_V_i(1:m); */
      /* 'dynIdenf:331' i_omega_dot_i = i_A_i(1:m); */
      /*  the classical absolute acceleration of reference i with */
      /*  components in itself */
      /* 'dynIdenf:334' i_d_ddot_i = i_A_i(m+1:2*m) + skew(i_V_i(1:m)) * i_V_i(m+1:2*m) - i_X_o(1:m,1:m) * g; */
      /*  transform a R3 vector into a skew matrix */
      /* 'dynIdenf:352' y = [0, -x(3), x(2); x(3), 0, -x(1); -x(2) x(1) 0]; */
      dv11[0][0] = 0.0;
      dv11[1][0] = -i_V_i[2];
      dv11[2][0] = i_V_i[1];
      dv11[0][1] = i_V_i[2];
      dv11[1][1] = 0.0;
      dv11[2][1] = -i_V_i[0];
      dv11[0][2] = -i_V_i[1];
      dv11[1][2] = i_V_i[0];
      dv11[2][2] = 0.0;

      /*  i_Am_i is a intermediate matrix, please refer to the book */
      /*  Springer Hand of Robotics, */
      /* 'dynIdenf:337' i_Am_i = [zeros(m,1),-skew(i_d_ddot_i),lin(i_omega_dot_i)+skew(i_omega_i)*lin(i_omega_i);... */
      /* 'dynIdenf:338'             i_d_ddot_i,skew(i_omega_dot_i)+skew(i_omega_i)^2,zeros(m,2*m)]; */
      /*  transform a R3 vector into a skew matrix */
      /* 'dynIdenf:352' y = [0, -x(3), x(2); x(3), 0, -x(1); -x(2) x(1) 0]; */
      /*  transform the angular velocity and acceleration vector into a */
      /*  matrix for constructing the linear form of dynamics */
      /* 'dynIdenf:361' y = [x(1),x(2),x(3),0,0,0;... */
      /* 'dynIdenf:362'             0,x(1),0,x(2),x(3),0;... */
      /* 'dynIdenf:363'             0,0,x(1),0,x(2),x(3)]; */
      /*  transform a R3 vector into a skew matrix */
      /* 'dynIdenf:352' y = [0, -x(3), x(2); x(3), 0, -x(1); -x(2) x(1) 0]; */
      /*  transform the angular velocity and acceleration vector into a */
      /*  matrix for constructing the linear form of dynamics */
      /* 'dynIdenf:361' y = [x(1),x(2),x(3),0,0,0;... */
      /* 'dynIdenf:362'             0,x(1),0,x(2),x(3),0;... */
      /* 'dynIdenf:363'             0,0,x(1),0,x(2),x(3)]; */
      /*  transform a R3 vector into a skew matrix */
      /* 'dynIdenf:352' y = [0, -x(3), x(2); x(3), 0, -x(1); -x(2) x(1) 0]; */
      /*  transform a R3 vector into a skew matrix */
      /* 'dynIdenf:352' y = [0, -x(3), x(2); x(3), 0, -x(1); -x(2) x(1) 0]; */
      i1_R_i[0][0] = 0.0;
      i1_R_i[1][0] = -i_V_i[2];
      i1_R_i[2][0] = i_V_i[1];
      i1_R_i[0][1] = i_V_i[2];
      i1_R_i[1][1] = 0.0;
      i1_R_i[2][1] = -i_V_i[0];
      i1_R_i[0][2] = -i_V_i[1];
      i1_R_i[1][2] = i_V_i[0];
      i1_R_i[2][2] = 0.0;
      dv12[0][0] = 0.0;
      dv12[1][0] = -i_V_i[2];
      dv12[2][0] = i_V_i[1];
      dv12[0][1] = i_V_i[2];
      dv12[1][1] = 0.0;
      dv12[2][1] = -i_V_i[0];
      dv12[0][2] = -i_V_i[1];
      dv12[1][2] = i_V_i[0];
      dv12[2][2] = 0.0;
      b_i_V_i[0][0] = i_V_i[0];
      b_i_V_i[1][0] = i_V_i[1];
      b_i_V_i[2][0] = i_V_i[2];
      b_i_V_i[3][0] = 0.0;
      b_i_V_i[4][0] = 0.0;
      b_i_V_i[5][0] = 0.0;
      b_i_V_i[0][1] = 0.0;
      b_i_V_i[1][1] = i_V_i[0];
      b_i_V_i[2][1] = 0.0;
      b_i_V_i[3][1] = i_V_i[1];
      b_i_V_i[4][1] = i_V_i[2];
      b_i_V_i[5][1] = 0.0;
      b_i_V_i[0][2] = 0.0;
      b_i_V_i[1][2] = 0.0;
      b_i_V_i[2][2] = i_V_i[0];
      b_i_V_i[3][2] = 0.0;
      b_i_V_i[4][2] = i_V_i[1];
      b_i_V_i[5][2] = i_V_i[2];
      b_i_A_i[0][0] = i_A_i[0];
      b_i_A_i[1][0] = i_A_i[1];
      b_i_A_i[2][0] = i_A_i[2];
      b_i_A_i[3][0] = 0.0;
      b_i_A_i[4][0] = 0.0;
      b_i_A_i[5][0] = 0.0;
      b_i_A_i[0][1] = 0.0;
      b_i_A_i[1][1] = i_A_i[0];
      b_i_A_i[2][1] = 0.0;
      b_i_A_i[3][1] = i_A_i[1];
      b_i_A_i[4][1] = i_A_i[2];
      b_i_A_i[5][1] = 0.0;
      b_i_A_i[0][2] = 0.0;
      b_i_A_i[1][2] = 0.0;
      b_i_A_i[2][2] = i_A_i[0];
      b_i_A_i[3][2] = 0.0;
      b_i_A_i[4][2] = i_A_i[1];
      b_i_A_i[5][2] = i_A_i[2];
      dv14[0][0] = 0.0;
      dv14[1][0] = -i_A_i[2];
      dv14[2][0] = i_A_i[1];
      dv14[0][1] = i_A_i[2];
      dv14[1][1] = 0.0;
      dv14[2][1] = -i_A_i[0];
      dv14[0][2] = -i_A_i[1];
      dv14[1][2] = i_A_i[0];
      dv14[2][2] = 0.0;
      for (i12 = 0; i12 < 3; i12++) {
        c_o_X_i = 0.0;
        for (i13 = 0; i13 < 3; i13++) {
          c_o_X_i += dv11[i13][i12] * i_V_i[3 + i13];
        }

        b_i_X_o[i12] = 0.0;
        for (i13 = 0; i13 < 3; i13++) {
          b_i_X_o[i12] += i_X_o[i13][i12] * b_g[i13];
        }

        i1_p_i[i12] = (i_A_i[3 + i12] + c_o_X_i) - b_i_X_o[i12];
        for (i13 = 0; i13 < 6; i13++) {
          dv15[i13][i12] = 0.0;
          for (i14 = 0; i14 < 3; i14++) {
            dv15[i13][i12] += dv12[i14][i12] * b_i_V_i[i13][i14];
          }
        }

        for (i13 = 0; i13 < 3; i13++) {
          b_i1_R_i[i13][i12] = 0.0;
          for (i14 = 0; i14 < 3; i14++) {
            b_i1_R_i[i13][i12] += i1_R_i[i14][i12] * i1_R_i[i13][i14];
          }
        }

        i_Am_i[0][i12] = 0.0;
      }

      i_Am_i[1][0] = -0.0;
      i_Am_i[2][0] = -(-i1_p_i[2]);
      i_Am_i[3][0] = -i1_p_i[1];
      i_Am_i[1][1] = -i1_p_i[2];
      i_Am_i[2][1] = -0.0;
      i_Am_i[3][1] = -(-i1_p_i[0]);
      i_Am_i[1][2] = -(-i1_p_i[1]);
      i_Am_i[2][2] = -i1_p_i[0];
      i_Am_i[3][2] = -0.0;
      for (i12 = 0; i12 < 6; i12++) {
        for (i13 = 0; i13 < 3; i13++) {
          i_Am_i[i12 + 4][i13] = b_i_A_i[i12][i13] + dv15[i12][i13];
        }
      }

      for (i12 = 0; i12 < 3; i12++) {
        i_Am_i[0][i12 + 3] = i1_p_i[i12];
        for (i13 = 0; i13 < 3; i13++) {
          i_Am_i[i12 + 1][i13 + 3] = dv14[i12][i13] + b_i1_R_i[i12][i13];
        }
      }

      for (i12 = 0; i12 < 6; i12++) {
        for (i13 = 0; i13 < 3; i13++) {
          i_Am_i[i12 + 4][i13 + 3] = 0.0;
        }
      }

      /*  calculate the matrix K by iterations */
      /* 'dynIdenf:340' for j = 1:i */
      c_o_X_i = 12.0 * ((1.0 + (double)b_i) - 1.0) + 1.0;
      c_i_X_o = 12.0 * (1.0 + (double)b_i);
      if (c_o_X_i > c_i_X_o) {
        i12 = 0;
        i13 = 0;
      } else {
        i12 = (int)c_o_X_i - 1;
        i13 = (int)c_i_X_o;
      }

      loop_ub = i13 - i12;
      for (j = 0; j <= b_i; j++) {
        /* 'dynIdenf:341' temp = (i_X_o * o_X_i(:,:,j))' * i_Am_i; */
        /* 'dynIdenf:342' K(n*(k-1)+j,nparJoint*(i-1)+1:nparJoint*i) = [temp(m,:), (i==j)*sign(theta_dot(k,i)), (i==j)*theta_dot(k,i)]; */
        i14 = (int)((double)n * ((1.0 + (double)k) - 1.0) + (1.0 + (double)j)) -
          1;
        inner = r5->size[0];
        r5->size[0] = i13 - i12;
        emxEnsureCapacity_int32_T(r5, inner);
        for (inner = 0; inner < loop_ub; inner++) {
          r5->data[inner] = i12 + inner;
        }

        c_o_X_i = theta_dot->data[k + theta_dot->size[0] * b_i];
        if (theta_dot->data[k + theta_dot->size[0] * b_i] < 0.0) {
          c_o_X_i = -1.0;
        } else {
          if (theta_dot->data[k + theta_dot->size[0] * b_i] > 0.0) {
            c_o_X_i = 1.0;
          }
        }

        for (inner = 0; inner < 6; inner++) {
          for (i = 0; i < 6; i++) {
            b_o_X_i[i][inner] = 0.0;
            for (m = 0; m < 6; m++) {
              c_i_X_o = b_o_X_i[i][inner] + i_X_o[m][i] * o_X_i->data[(m +
                o_X_i->size[0] * inner) + o_X_i->size[0] * o_X_i->size[1] * j];
              b_o_X_i[i][inner] = c_i_X_o;
            }
          }

          for (i = 0; i < 10; i++) {
            e_i_X_o[i][inner] = 0.0;
            for (m = 0; m < 6; m++) {
              e_i_X_o[i][inner] += b_o_X_i[m][inner] * i_Am_i[i][m];
            }
          }
        }

        for (inner = 0; inner < 10; inner++) {
          d_i_X_o[inner] = e_i_X_o[inner][2];
        }

        d_i_X_o[10] = (double)(1 + b_i == 1 + j) * c_o_X_i;
        d_i_X_o[11] = (double)(1 + b_i == 1 + j) * theta_dot->data[k +
          theta_dot->size[0] * b_i];
        i = r5->size[0];
        for (inner = 0; inner < i; inner++) {
          K->data[i14 + K->size[0] * r5->data[inner]] = d_i_X_o[inner];
        }
      }
    }

    k++;
  }

  emxFree_int32_T(&r5);
  emxFree_real_T(&o_X_i);
  emxInit_real_T(&c_a, 2);

  /*  regroup the K corresponding to the minimum set of dynamic */
  /*  parameters */
  /* 'dynIdenf:348' K = K * Rpk; */
  i12 = c_a->size[0] * c_a->size[1];
  c_a->size[0] = K->size[0];
  c_a->size[1] = K->size[1];
  emxEnsureCapacity_real_T(c_a, i12);
  loop_ub = K->size[1];
  for (i12 = 0; i12 < loop_ub; i12++) {
    m = K->size[0];
    for (i13 = 0; i13 < m; i13++) {
      c_a->data[i13 + c_a->size[0] * i12] = K->data[i13 + K->size[0] * i12];
    }
  }

  if ((K->size[1] == 1) || (b_Rpk->size[0] == 1)) {
    i12 = c_a->size[0] * c_a->size[1];
    c_a->size[0] = K->size[0];
    c_a->size[1] = b_Rpk->size[1];
    emxEnsureCapacity_real_T(c_a, i12);
    loop_ub = K->size[0];
    for (i12 = 0; i12 < loop_ub; i12++) {
      m = b_Rpk->size[1];
      for (i13 = 0; i13 < m; i13++) {
        c_a->data[i12 + c_a->size[0] * i13] = 0.0;
        i = K->size[1];
        for (i14 = 0; i14 < i; i14++) {
          c_a->data[i12 + c_a->size[0] * i13] += K->data[i12 + K->size[0] * i14]
            * b_Rpk->data[i14 + b_Rpk->size[0] * i13];
        }
      }
    }

    i12 = K->size[0] * K->size[1];
    K->size[0] = c_a->size[0];
    K->size[1] = c_a->size[1];
    emxEnsureCapacity_real_T(K, i12);
    loop_ub = c_a->size[1];
    for (i12 = 0; i12 < loop_ub; i12++) {
      m = c_a->size[0];
      for (i13 = 0; i13 < m; i13++) {
        K->data[i13 + K->size[0] * i12] = c_a->data[i13 + c_a->size[0] * i12];
      }
    }
  } else {
    m = K->size[0];
    inner = K->size[1];
    i = K->size[0];
    i12 = K->size[0] * K->size[1];
    K->size[0] = i;
    K->size[1] = b_Rpk->size[1];
    emxEnsureCapacity_real_T(K, i12);
    for (j = 0; j < b_Rpk->size[1]; j++) {
      for (b_i = 1; b_i <= m; b_i++) {
        K->data[(b_i + K->size[0] * j) - 1] = 0.0;
      }

      for (k = 0; k < inner; k++) {
        if (b_Rpk->data[k + b_Rpk->size[0] * j] != 0.0) {
          for (b_i = 0; b_i < m; b_i++) {
            K->data[b_i + K->size[0] * j] += b_Rpk->data[k + b_Rpk->size[0] * j]
              * c_a->data[b_i + c_a->size[0] * k];
          }
        }
      }
    }
  }

  emxFree_real_T(&c_a);
}

/*
 * function [num, den] = getCoeffs(fpass)
 * the function generates the numerator and denominator of the filter
 *  according to pass frequency
 * Arguments    : double fpass
 *                emxArray_real_T *b_num
 *                double *b_den
 * Return Type  : void
 */
static void getCoeffs(double fpass, emxArray_real_T *b_num, double *b_den)
{
  int i6;
  static const double dv1[684] = { -0.000520978002302694, -5.20169990936757E-5,
    -5.45841105354002E-5, -5.72195149660382E-5, -5.99216529332766E-5,
    -6.26902993572831E-5, -6.5527757231746E-5, -6.84327924698304E-5,
    -7.14024783352774E-5, -7.44382315108702E-5, -7.75419955514318E-5,
    -8.07125864750608E-5, -8.3949220270288E-5, -8.72527281450113E-5,
    -9.06217587629433E-5, -9.40534210975644E-5, -9.75483334878834E-5,
    -0.000101109902849751, -0.000104738576954986, -0.000108431175299353,
    -0.000112185059058074, -0.000115999352356773, -0.000119874199376914,
    -0.000123811981761444, -0.000127815442979604, -0.000131881513112167,
    -0.000136002977132762, -0.0001401798588683, -0.000144420119293729,
    -0.000148723321214873, -0.000153076275831861, -0.000157476132364005,
    -0.000161941418859757, -0.000166463101250286, -0.000171015020968334,
    -0.000175636719599759, -0.000180296194896843, -0.000185006087018483,
    -0.000189757036251856, -0.000194554159433506, -0.000199389160913883,
    -0.000204258823938341, -0.000209166018405914, -0.000214108569106062,
    -0.000219080375140543, -0.000224079825588825, -0.000229108242018683,
    -0.000234163373811686, -0.000239241113430869, -0.000244340284774501,
    -0.000249460464775636, -0.000254598436236085, -0.000259750812540093,
    -0.000264916528117114, -0.000270093766899914, -0.000275278552826386,
    -0.000280468430692387, -0.000285662688010565, -0.00029085760600134,
    -0.000296047191551294, -0.00030122915946475, -0.000306403210484201,
    -0.000311563497335765, -0.000316701436270235, -0.000321814827623077,
    -0.000326902443782036, -0.000331953435387113, -0.000336972143859848,
    -0.000341937352236693, -0.000346885828652329, -0.000351751509861389,
    -0.000356584499208047, -0.00036137547994058, -0.000366096777583531,
    -0.000370757272781832, -0.000375368791800677, -0.000379922791039298,
    -0.000384404674196856, -0.000388811574995772, -0.00039314670724098,
    -0.000397406523568648, -0.000401580477775939, -0.000405657064973403,
    -0.000409625775268277, -0.000413478116422792, -0.000417209639330773,
    -0.00042081780637044, -0.000424298109597335, -0.000427646677871602,
    -0.000430866002181443, -0.000433963703733845, -0.000436946687593806,
    -0.00043981922138183, -0.000442581391775045, -0.000445221558359002,
    -0.000447710843745925, -0.000450008300141517, -0.000452070267800783,
    -0.000453863926632936, -0.000455404190153653, -0.000456810882109912,
    -0.000458350037646638, -0.000460310821393207, -0.000460348968389753,
    -0.000461459608392154, -0.000462016902456549, -0.000462378433406761,
    -0.000462524399925029, -0.000462451346355499, -0.000462151708972203,
    -0.000461618145097995, -0.000460847998429803, -0.000459837145402885,
    -0.000458578890725335, -0.000457068727987613, -0.000455302981268372,
    -0.000453274535611701, -0.000450975934296376, -0.000448404400488288,
    -0.000445558845249149, -0.00044243445974399, -0.000439023959439499,
    -0.000435321003830649, -0.000431319729545988, -0.000427014842399289,
    -0.000422403822988828, -0.000417484095055867, -0.000412247617334492,
    -0.000406685229515458, -0.00040079601020735, -0.000394582525310948,
    -0.000388037089025544, -0.000381146478738261, -0.000373911374121919,
    -0.000366340474934905, -0.000358414692150382, -0.000350122619270888,
    -0.000341486217704825, -0.000332478179656225, -0.000323106326908167,
    -0.000313362143318375, -0.00030324909832825, -0.000292754542035043,
    -0.00028187712436533, -0.000270618581875681, -0.000258972238338168,
    -0.000246930199298563, -0.000234491181193772, -0.000221654608059194,
    -0.000208415138473624, -0.000194767852636238, -0.000180712190136185,
    -0.000166247400760134, -0.000151369883516121, -0.000136077105983035,
    -0.000120368829513933, -0.00010424343758073, -8.76979772704483E-5,
    -7.07315042938449E-5, -5.33432778794504E-5, -3.55290954195747E-5,
    -1.72848723039368E-5, 1.38872191841378E-6, 2.04911018494417E-5,
    4.00271837304823E-5, 6.00012002257589E-5, 8.04110281968471E-5,
    0.000101256579820492, 0.000122543061367604, 0.000144261922674221,
    0.000166430525578065, 0.000188999271198846, 0.000212060291713194,
    0.000235509463682494, 0.000259391814847139, 0.00028373168109523,
    0.000308500057039596, 0.000333684825244486, 0.000359302512586066,
    0.000385367215620059, 0.000411873717316901, 0.000438810480214042,
    0.000466172926860413, 0.000493961235466057, 0.000522174877221049,
    0.000550811482289791, 0.00057986556184007, 0.000609328804955685,
    0.000639196293728137, 0.000669471252782896, 0.000700160521508816,
    0.000731267776747155, 0.000762792692118473, 0.000794730701659777,
    0.000827069207154016, 0.000859788963085936, 0.000892874456342517,
    0.000926322448581773, 0.000960143220924291, 0.000994359339319614,
    0.00102899356216643, 0.00106403359253597, 0.00109940448183309,
    0.0011350127622632, 0.00117092633398118, 0.00120756603038993,
    0.00124414286807791, 0.00128124465925545, 0.00131863665906098,
    0.00135635578729621, 0.00139439095544912, 0.00143273921326114,
    0.00147139468118457, 0.00151034727681911, 0.0015495900608145,
    0.00158911709434638, 0.00162891923891444, 0.00166898846896255,
    0.00170932107190954, 0.00174991207709709, 0.00179075141906086,
    0.00183182857807123, 0.00187313619701259, 0.00191466727083169,
    0.0019564134765558, 0.00199836703727341, 0.00204051967345753,
    0.00208286050291843, 0.00212538001142598, 0.00216807371366821,
    0.00221093576143788, 0.00225395281482843, 0.00229711200656131,
    0.00234040938008887, 0.00238383968702757, 0.00242738493017718,
    0.00247102931641753, 0.00251477773243179, 0.00255861614240247,
    0.00260252203198765, 0.00264650059615073, 0.00269053584052438,
    0.00273461819792636, 0.00277873398090667, 0.00282288392309461,
    0.00286705106750894, 0.00291122100593392, 0.00295538889026983,
    0.00299954644597872, 0.00304367821272944, 0.0030877714892921,
    0.00313181876280088, 0.00317580991187943, 0.0032197307991366,
    0.00326356957138265, 0.00330731742211808, 0.00335096383756065,
    0.00339449670702713, 0.00343790560431774, 0.00348118036145555,
    0.00352430894741731, 0.00356728063999266, 0.00361008749181655,
    0.00365271900135355, 0.0036951601983629, 0.0037373983268239,
    0.00377942516996933, 0.00382122980465069, 0.00386279739394695,
    0.00390411765143349, 0.00394518232614754, 0.00398597405354187,
    0.00402649854545029, 0.00406669350872107, 0.00410667702228297,
    0.00414625801694625, 0.00418555412451145, 0.0042245533013057,
    0.00426318977685051, 0.0043014485984287, 0.00433934592051402,
    0.00437688313874072, 0.00441403882587412, 0.0044507906652811,
    0.00448712803027301, 0.00452304763244816, 0.0045585475568763,
    0.00459362319520842, 0.00462826233752556, 0.00466244580837228,
    0.00469615596427826, 0.00472938226584746, 0.00476211864665913,
    0.00479436106160584, 0.0048261082246571, 0.00485735774695412,
    0.00488809853788567, 0.00491831081285634, 0.00494797438450371,
    0.00497707538941737, 0.00500561055815462, 0.00503359021864047,
    0.00506102608993381, 0.00508790129864639, 0.00511416372355433,
    0.00513978126324414, 0.00516482484920594, 0.0051893798127028,
    0.00521313023739376, 0.0052363822494749, 0.00525895958423336,
    0.00528090765629976, 0.00530220995844425, 0.00532286619204487,
    0.00534287022373767, 0.00536221380126212, 0.00538089270680401,
    0.00539890107887976, 0.00541622949372172, 0.00543287183790025,
    0.00544882648850096, 0.00546408964030538, 0.00547865411988842,
    0.00549251489038215, 0.00550566997871806, 0.00551811665232018,
    0.00552985145136765, 0.00554087212638495, 0.00555117526124427,
    0.0055607549582844, 0.00556960768651819, 0.00557773394998625,
    0.00558513146790002, 0.00559179314635169, 0.00559771613805072,
    0.005602904459173, 0.00560735723206231, 0.00561106560533546,
    0.00561402913964094, 0.0056162588894292, 0.00561774501983524,
    0.00561848476850506, 0.00561848476850506, 0.00561774501983524,
    0.0056162588894292, 0.00561402913964094, 0.00561106560533546,
    0.00560735723206231, 0.005602904459173, 0.00559771613805072,
    0.00559179314635169, 0.00558513146790002, 0.00557773394998625,
    0.00556960768651819, 0.0055607549582844, 0.00555117526124427,
    0.00554087212638495, 0.00552985145136765, 0.00551811665232018,
    0.00550566997871806, 0.00549251489038215, 0.00547865411988842,
    0.00546408964030538, 0.00544882648850096, 0.00543287183790025,
    0.00541622949372172, 0.00539890107887976, 0.00538089270680401,
    0.00536221380126212, 0.00534287022373767, 0.00532286619204487,
    0.00530220995844425, 0.00528090765629976, 0.00525895958423336,
    0.0052363822494749, 0.00521313023739376, 0.0051893798127028,
    0.00516482484920594, 0.00513978126324414, 0.00511416372355433,
    0.00508790129864639, 0.00506102608993381, 0.00503359021864047,
    0.00500561055815462, 0.00497707538941737, 0.00494797438450371,
    0.00491831081285634, 0.00488809853788567, 0.00485735774695412,
    0.0048261082246571, 0.00479436106160584, 0.00476211864665913,
    0.00472938226584746, 0.00469615596427826, 0.00466244580837228,
    0.00462826233752556, 0.00459362319520842, 0.0045585475568763,
    0.00452304763244816, 0.00448712803027301, 0.0044507906652811,
    0.00441403882587412, 0.00437688313874072, 0.00433934592051402,
    0.0043014485984287, 0.00426318977685051, 0.0042245533013057,
    0.00418555412451145, 0.00414625801694625, 0.00410667702228297,
    0.00406669350872107, 0.00402649854545029, 0.00398597405354187,
    0.00394518232614754, 0.00390411765143349, 0.00386279739394695,
    0.00382122980465069, 0.00377942516996933, 0.0037373983268239,
    0.0036951601983629, 0.00365271900135355, 0.00361008749181655,
    0.00356728063999266, 0.00352430894741731, 0.00348118036145555,
    0.00343790560431774, 0.00339449670702713, 0.00335096383756065,
    0.00330731742211808, 0.00326356957138265, 0.0032197307991366,
    0.00317580991187943, 0.00313181876280088, 0.0030877714892921,
    0.00304367821272944, 0.00299954644597872, 0.00295538889026983,
    0.00291122100593392, 0.00286705106750894, 0.00282288392309461,
    0.00277873398090667, 0.00273461819792636, 0.00269053584052438,
    0.00264650059615073, 0.00260252203198765, 0.00255861614240247,
    0.00251477773243179, 0.00247102931641753, 0.00242738493017718,
    0.00238383968702757, 0.00234040938008887, 0.00229711200656131,
    0.00225395281482843, 0.00221093576143788, 0.00216807371366821,
    0.00212538001142598, 0.00208286050291843, 0.00204051967345753,
    0.00199836703727341, 0.0019564134765558, 0.00191466727083169,
    0.00187313619701259, 0.00183182857807123, 0.00179075141906086,
    0.00174991207709709, 0.00170932107190954, 0.00166898846896255,
    0.00162891923891444, 0.00158911709434638, 0.0015495900608145,
    0.00151034727681911, 0.00147139468118457, 0.00143273921326114,
    0.00139439095544912, 0.00135635578729621, 0.00131863665906098,
    0.00128124465925545, 0.00124414286807791, 0.00120756603038993,
    0.00117092633398118, 0.0011350127622632, 0.00109940448183309,
    0.00106403359253597, 0.00102899356216643, 0.000994359339319614,
    0.000960143220924291, 0.000926322448581773, 0.000892874456342517,
    0.000859788963085936, 0.000827069207154016, 0.000794730701659777,
    0.000762792692118473, 0.000731267776747155, 0.000700160521508816,
    0.000669471252782896, 0.000639196293728137, 0.000609328804955685,
    0.00057986556184007, 0.000550811482289791, 0.000522174877221049,
    0.000493961235466057, 0.000466172926860413, 0.000438810480214042,
    0.000411873717316901, 0.000385367215620059, 0.000359302512586066,
    0.000333684825244486, 0.000308500057039596, 0.00028373168109523,
    0.000259391814847139, 0.000235509463682494, 0.000212060291713194,
    0.000188999271198846, 0.000166430525578065, 0.000144261922674221,
    0.000122543061367604, 0.000101256579820492, 8.04110281968471E-5,
    6.00012002257589E-5, 4.00271837304823E-5, 2.04911018494417E-5,
    1.38872191841378E-6, -1.72848723039368E-5, -3.55290954195747E-5,
    -5.33432778794504E-5, -7.07315042938449E-5, -8.76979772704483E-5,
    -0.00010424343758073, -0.000120368829513933, -0.000136077105983035,
    -0.000151369883516121, -0.000166247400760134, -0.000180712190136185,
    -0.000194767852636238, -0.000208415138473624, -0.000221654608059194,
    -0.000234491181193772, -0.000246930199298563, -0.000258972238338168,
    -0.000270618581875681, -0.00028187712436533, -0.000292754542035043,
    -0.00030324909832825, -0.000313362143318375, -0.000323106326908167,
    -0.000332478179656225, -0.000341486217704825, -0.000350122619270888,
    -0.000358414692150382, -0.000366340474934905, -0.000373911374121919,
    -0.000381146478738261, -0.000388037089025544, -0.000394582525310948,
    -0.00040079601020735, -0.000406685229515458, -0.000412247617334492,
    -0.000417484095055867, -0.000422403822988828, -0.000427014842399289,
    -0.000431319729545988, -0.000435321003830649, -0.000439023959439499,
    -0.00044243445974399, -0.000445558845249149, -0.000448404400488288,
    -0.000450975934296376, -0.000453274535611701, -0.000455302981268372,
    -0.000457068727987613, -0.000458578890725335, -0.000459837145402885,
    -0.000460847998429803, -0.000461618145097995, -0.000462151708972203,
    -0.000462451346355499, -0.000462524399925029, -0.000462378433406761,
    -0.000462016902456549, -0.000461459608392154, -0.000460348968389753,
    -0.000460310821393207, -0.000458350037646638, -0.000456810882109912,
    -0.000455404190153653, -0.000453863926632936, -0.000452070267800783,
    -0.000450008300141517, -0.000447710843745925, -0.000445221558359002,
    -0.000442581391775045, -0.00043981922138183, -0.000436946687593806,
    -0.000433963703733845, -0.000430866002181443, -0.000427646677871602,
    -0.000424298109597335, -0.00042081780637044, -0.000417209639330773,
    -0.000413478116422792, -0.000409625775268277, -0.000405657064973403,
    -0.000401580477775939, -0.000397406523568648, -0.00039314670724098,
    -0.000388811574995772, -0.000384404674196856, -0.000379922791039298,
    -0.000375368791800677, -0.000370757272781832, -0.000366096777583531,
    -0.00036137547994058, -0.000356584499208047, -0.000351751509861389,
    -0.000346885828652329, -0.000341937352236693, -0.000336972143859848,
    -0.000331953435387113, -0.000326902443782036, -0.000321814827623077,
    -0.000316701436270235, -0.000311563497335765, -0.000306403210484201,
    -0.00030122915946475, -0.000296047191551294, -0.00029085760600134,
    -0.000285662688010565, -0.000280468430692387, -0.000275278552826386,
    -0.000270093766899914, -0.000264916528117114, -0.000259750812540093,
    -0.000254598436236085, -0.000249460464775636, -0.000244340284774501,
    -0.000239241113430869, -0.000234163373811686, -0.000229108242018683,
    -0.000224079825588825, -0.000219080375140543, -0.000214108569106062,
    -0.000209166018405914, -0.000204258823938341, -0.000199389160913883,
    -0.000194554159433506, -0.000189757036251856, -0.000185006087018483,
    -0.000180296194896843, -0.000175636719599759, -0.000171015020968334,
    -0.000166463101250286, -0.000161941418859757, -0.000157476132364005,
    -0.000153076275831861, -0.000148723321214873, -0.000144420119293729,
    -0.0001401798588683, -0.000136002977132762, -0.000131881513112167,
    -0.000127815442979604, -0.000123811981761444, -0.000119874199376914,
    -0.000115999352356773, -0.000112185059058074, -0.000108431175299353,
    -0.000104738576954986, -0.000101109902849751, -9.75483334878834E-5,
    -9.40534210975644E-5, -9.06217587629433E-5, -8.72527281450113E-5,
    -8.3949220270288E-5, -8.07125864750608E-5, -7.75419955514318E-5,
    -7.44382315108702E-5, -7.14024783352774E-5, -6.84327924698304E-5,
    -6.5527757231746E-5, -6.26902993572831E-5, -5.99216529332766E-5,
    -5.72195149660382E-5, -5.45841105354002E-5, -5.20169990936757E-5,
    -0.000520978002302694 };

  static const double dv2[723] = { 0.000524469736335486, 5.00586518032772E-5,
    5.23815850259164E-5, 5.47674639762573E-5, 5.7167569336772E-5,
    5.96287108106926E-5, 6.21002944473329E-5, 6.46252421778573E-5,
    6.71538849242156E-5, 6.97343785522134E-5, 7.23168342051082E-5,
    7.49510402638714E-5, 7.75832509210173E-5, 8.02551921214269E-5,
    8.29159682127781E-5, 8.56184309621784E-5, 8.83159101463416E-5,
    9.10537173189225E-5, 9.37690208083129E-5, 9.6501764344623E-5,
    9.92159696337733E-5, 0.000101967709425591, 0.000104680883398715,
    0.000107392581745282, 0.000110108789624883, 0.000112800194579151,
    0.000115468005236082, 0.000118132310510059, 0.000120750142158453,
    0.000123362723034365, 0.000125918478995377, 0.000128458153510528,
    0.000130933154839445, 0.000133391093408533, 0.000135776345148193,
    0.000138133925206034, 0.00014040937231192, 0.000142648282511158,
    0.000144799075700251, 0.000146905389952095, 0.000148909574313128,
    0.000150860785731569, 0.000152708088942437, 0.000154490393439512,
    0.000156153810147281, 0.000157751952845765, 0.000159227356866869,
    0.000160617329536024, 0.000161881818578263, 0.000163049900688652,
    0.000164080891749098, 0.000165012261291061, 0.000165793466249417,
    0.000166463912498038, 0.000166977599482095, 0.000167371723214612,
    0.000167596209337449, 0.00016769546361647, 0.000167612685064898,
    0.000167395706802922, 0.000166989668474572, 0.000166435386587554,
    0.000165681305169779, 0.000164776375709852, 0.000163662215298123,
    0.000162382142477805, 0.000160885927205586, 0.000159218674345988,
    0.000157324536310168, 0.000155251833704162, 0.000152947120543503,
    0.00015045433430142, 0.000147727128052281, 0.000144804508496339,
    0.00014163973860547, 0.000138276287207714, 0.000134668974387507,
    0.000130850035109437, 0.000126785208378037, 0.000122502135281101,
    0.000117960899953838, 0.000113210405380578, 0.00010818882774188,
    0.000102945185568431, 9.74338509893395E-5, 9.16972189395753E-5,
    8.5679575283788E-5, 7.94321200124609E-5, 7.2904053299632E-5,
    6.61402958740868E-5, 5.9088983210048E-5, 5.17992381686803E-5,
    4.42231681996236E-5, 3.64139273452416E-5, 2.83209585120617E-5,
    1.99914518409441E-5, 1.1382115264166E-5, 2.54292500196256E-6,
    -6.56931307847807E-6, -1.59122093433135E-5, -2.55059336623784E-5,
    -3.53685936727762E-5, -4.54400540476134E-5, -5.57696549960879E-5,
    -6.63771677159197E-5, -7.7184238904831E-5, -8.82595020734168E-5,
    -9.95753499456616E-5, -0.000111172538004551, -0.000122985841755254,
    -0.000135054315559792, -0.000147323340028758, -0.000159835502885448,
    -0.000172526661615128, -0.000185423796157919, -0.000198461680990827,
    -0.000211686103776358, -0.000225058353634012, -0.000238655847061092,
    -0.000252470845359277, -0.000266585641496571, -0.000280950954741698,
    -0.000295511902559795, -0.000309919414694471, -0.000323659926009783,
    -0.000339293525390544, -0.0003539411416973, -0.000368940000057705,
    -0.000384067514638573, -0.000399266350183937, -0.000414577645869965,
    -0.000429943114160157, -0.000445402498835314, -0.000460893486056534,
    -0.000476456042465615, -0.000492024291725858, -0.000507642405972521,
    -0.000523252423509115, -0.000538894259811931, -0.000554498478923738,
    -0.000570099611785535, -0.000585633180585017, -0.000601147785670214,
    -0.000616582149430934, -0.000631961341646842, -0.000647217512299709,
    -0.000662409644078374, -0.0006774691855177, -0.00069239830826689,
    -0.000707187922078236, -0.000721827244014612, -0.000736275102735415,
    -0.000750555856424403, -0.00076460715899914, -0.000778464178550299,
    -0.000792065298356889, -0.000805441953162336, -0.000818525790246944,
    -0.000831355249152882, -0.000843861392426555, -0.000856081097708123,
    -0.000867944860467514, -0.000879485936148014, -0.000890637435221489,
    -0.000901438137545911, -0.000911814674348623, -0.000921802153901256,
    -0.000931338658585429, -0.000940458537303633, -0.000949086060018919,
    -0.000957261986294197, -0.000964921530360863, -0.000972089771447915,
    -0.000978711973604225, -0.000984810266900323, -0.000990324037377398,
    -0.000995284735363989, -0.000999628829097016, -0.00100338311657393,
    -0.00100648840624336, -0.00100897383676334, -0.00101077273655536,
    -0.00101192337158328, -0.00101235351098358, -0.00101210187328732,
    -0.00101110505554398, -0.00100939628756502, -0.00100690493368262,
    -0.00100367346039378, -0.000999635478971382, -0.000994825353485726,
    -0.000989177746920166, -0.000982732992102768, -0.000975422972123854,
    -0.00096728806560114, -0.000958264039244553, -0.000948385818031927,
    -0.000937596065739797, -0.000925928377171722, -0.000913320924087787,
    -0.000899806758135586, -0.000885336927490748, -0.000869930983136534,
    -0.00085355194217766, -0.000836221973888521, -0.000817881480076714,
    -0.000798585263593298, -0.00077826718196216, -0.000756965390748311,
    -0.000734626928834411, -0.000711299284505064, -0.000686917624821226,
    -0.000661526837259463, -0.000635071485390282, -0.000607599524831393,
    -0.000579052925211444, -0.000549475774246958, -0.000518807513898672,
    -0.000487098981962924, -0.000454297513830191, -0.000420448886202128,
    -0.00038550104935414, -0.000349502134594296, -0.000312408368692861,
    -0.000274245563333016, -0.000235052758243438, -0.000194716827169965,
    -0.000153368898219367, -0.000110976175697161, -6.74507400684875E-5,
    -2.28794438080612E-5, 2.27625335533632E-5, 6.94323009513472E-5,
    0.000117191616338152, 0.000165983944139365, 0.000215841768100563,
    0.00026668962419746, 0.000318565611876379, 0.000371418819152364,
    0.000425307568836569, 0.000480188155160603, 0.000536114684986784,
    0.000593018858295846, 0.000650913251078771, 0.000709697566461656,
    0.000769388541859999, 0.000829970522307472, 0.000891610805735191,
    0.000954285968006834, 0.00101744016505917, 0.00108179886173261,
    0.00114690729131758, 0.00121285144909704, 0.00127964433977491,
    0.00134721656050543, 0.00141559926007632, 0.00148472334593358,
    0.00155462231617105, 0.00162522635776969, 0.00169657224734207,
    0.00176858246215739, 0.00184128292924902, 0.00191459755274841,
    0.00198855971401751, 0.00206309700201124, 0.00213823884606977,
    0.00221390072824991, 0.00229011116259731, 0.00236680486462171,
    0.00244400829749066, 0.00252162755198125, 0.00259969797970943,
    0.00267815865187951, 0.00275699238707736, 0.00283615837128198,
    0.00291565612145279, 0.00299541841485039, 0.00307546448692997,
    0.00315571496135086, 0.00323618817147553, 0.00331680810872238,
    0.00339759728367135, 0.00347847095160714, 0.00355945302176549,
    0.00364045801355713, 0.00372150918167257, 0.00380252443146705,
    0.00388352332653804, 0.00396441939294586, 0.00404523950814305,
    0.00412589779852082, 0.00420640941352418, 0.00428669196063306,
    0.00436677183715796, 0.00444655832219488, 0.00452607142285969,
    0.00460523299812716, 0.00468405241174943, 0.00476245700439543,
    0.00484045896527067, 0.00491797738602819, 0.00499502806288758,
    0.00507153449172557, 0.00514750647263878, 0.00522286825169297,
    0.00529763742455926, 0.00537172869163968, 0.00544516922139953,
    0.00551787118353883, 0.0055898529352443, 0.00566103527103497,
    0.00573144265123906, 0.00580098564991597, 0.00586968765119188,
    0.00593746838994351, 0.0060043498380394, 0.00607024842906037,
    0.0061351908702529, 0.0061990943126447, 0.00626198447081014,
    0.00632378326795535, 0.00638450943149585, 0.00644408895040196,
    0.00650254590547329, 0.00655980744711141, 0.00661588627447452,
    0.00667072799918385, 0.00672432819037122, 0.00677663947138492,
    0.00682768236374154, 0.00687736905905883, 0.00692574000472932,
    0.00697272478704275, 0.00701834618360901, 0.00706252998982802,
    0.00710531700968379, 0.00714663590495778, 0.00718651793431919,
    0.00722488904769741, 0.00726178404355198, 0.00729713536227356,
    0.00733098161601232, 0.00736324911867787, 0.00739397283653819,
    0.00742308861567995, 0.00745063588114617, 0.00747655471033473,
    0.0075008835193354, 0.00752358594716239, 0.00754458981333914,
    0.00756409696501951, 0.00758179825183088, 0.00759784259849564,
    0.00761229933332533, 0.00762504302649761, 0.00763608902775804,
    0.00764540147060972, 0.0076530561615182, 0.00765901092419252,
    0.00766330901046047, 0.00766587753821448, 0.00766675273805942,
    0.00766587753821448, 0.00766330901046047, 0.00765901092419252,
    0.0076530561615182, 0.00764540147060972, 0.00763608902775804,
    0.00762504302649761, 0.00761229933332533, 0.00759784259849564,
    0.00758179825183088, 0.00756409696501951, 0.00754458981333914,
    0.00752358594716239, 0.0075008835193354, 0.00747655471033473,
    0.00745063588114617, 0.00742308861567995, 0.00739397283653819,
    0.00736324911867787, 0.00733098161601232, 0.00729713536227356,
    0.00726178404355198, 0.00722488904769741, 0.00718651793431919,
    0.00714663590495778, 0.00710531700968379, 0.00706252998982802,
    0.00701834618360901, 0.00697272478704275, 0.00692574000472932,
    0.00687736905905883, 0.00682768236374154, 0.00677663947138492,
    0.00672432819037122, 0.00667072799918385, 0.00661588627447452,
    0.00655980744711141, 0.00650254590547329, 0.00644408895040196,
    0.00638450943149585, 0.00632378326795535, 0.00626198447081014,
    0.0061990943126447, 0.0061351908702529, 0.00607024842906037,
    0.0060043498380394, 0.00593746838994351, 0.00586968765119188,
    0.00580098564991597, 0.00573144265123906, 0.00566103527103497,
    0.0055898529352443, 0.00551787118353883, 0.00544516922139953,
    0.00537172869163968, 0.00529763742455926, 0.00522286825169297,
    0.00514750647263878, 0.00507153449172557, 0.00499502806288758,
    0.00491797738602819, 0.00484045896527067, 0.00476245700439543,
    0.00468405241174943, 0.00460523299812716, 0.00452607142285969,
    0.00444655832219488, 0.00436677183715796, 0.00428669196063306,
    0.00420640941352418, 0.00412589779852082, 0.00404523950814305,
    0.00396441939294586, 0.00388352332653804, 0.00380252443146705,
    0.00372150918167257, 0.00364045801355713, 0.00355945302176549,
    0.00347847095160714, 0.00339759728367135, 0.00331680810872238,
    0.00323618817147553, 0.00315571496135086, 0.00307546448692997,
    0.00299541841485039, 0.00291565612145279, 0.00283615837128198,
    0.00275699238707736, 0.00267815865187951, 0.00259969797970943,
    0.00252162755198125, 0.00244400829749066, 0.00236680486462171,
    0.00229011116259731, 0.00221390072824991, 0.00213823884606977,
    0.00206309700201124, 0.00198855971401751, 0.00191459755274841,
    0.00184128292924902, 0.00176858246215739, 0.00169657224734207,
    0.00162522635776969, 0.00155462231617105, 0.00148472334593358,
    0.00141559926007632, 0.00134721656050543, 0.00127964433977491,
    0.00121285144909704, 0.00114690729131758, 0.00108179886173261,
    0.00101744016505917, 0.000954285968006834, 0.000891610805735191,
    0.000829970522307472, 0.000769388541859999, 0.000709697566461656,
    0.000650913251078771, 0.000593018858295846, 0.000536114684986784,
    0.000480188155160603, 0.000425307568836569, 0.000371418819152364,
    0.000318565611876379, 0.00026668962419746, 0.000215841768100563,
    0.000165983944139365, 0.000117191616338152, 6.94323009513472E-5,
    2.27625335533632E-5, -2.28794438080612E-5, -6.74507400684875E-5,
    -0.000110976175697161, -0.000153368898219367, -0.000194716827169965,
    -0.000235052758243438, -0.000274245563333016, -0.000312408368692861,
    -0.000349502134594296, -0.00038550104935414, -0.000420448886202128,
    -0.000454297513830191, -0.000487098981962924, -0.000518807513898672,
    -0.000549475774246958, -0.000579052925211444, -0.000607599524831393,
    -0.000635071485390282, -0.000661526837259463, -0.000686917624821226,
    -0.000711299284505064, -0.000734626928834411, -0.000756965390748311,
    -0.00077826718196216, -0.000798585263593298, -0.000817881480076714,
    -0.000836221973888521, -0.00085355194217766, -0.000869930983136534,
    -0.000885336927490748, -0.000899806758135586, -0.000913320924087787,
    -0.000925928377171722, -0.000937596065739797, -0.000948385818031927,
    -0.000958264039244553, -0.00096728806560114, -0.000975422972123854,
    -0.000982732992102768, -0.000989177746920166, -0.000994825353485726,
    -0.000999635478971382, -0.00100367346039378, -0.00100690493368262,
    -0.00100939628756502, -0.00101110505554398, -0.00101210187328732,
    -0.00101235351098358, -0.00101192337158328, -0.00101077273655536,
    -0.00100897383676334, -0.00100648840624336, -0.00100338311657393,
    -0.000999628829097016, -0.000995284735363989, -0.000990324037377398,
    -0.000984810266900323, -0.000978711973604225, -0.000972089771447915,
    -0.000964921530360863, -0.000957261986294197, -0.000949086060018919,
    -0.000940458537303633, -0.000931338658585429, -0.000921802153901256,
    -0.000911814674348623, -0.000901438137545911, -0.000890637435221489,
    -0.000879485936148014, -0.000867944860467514, -0.000856081097708123,
    -0.000843861392426555, -0.000831355249152882, -0.000818525790246944,
    -0.000805441953162336, -0.000792065298356889, -0.000778464178550299,
    -0.00076460715899914, -0.000750555856424403, -0.000736275102735415,
    -0.000721827244014612, -0.000707187922078236, -0.00069239830826689,
    -0.0006774691855177, -0.000662409644078374, -0.000647217512299709,
    -0.000631961341646842, -0.000616582149430934, -0.000601147785670214,
    -0.000585633180585017, -0.000570099611785535, -0.000554498478923738,
    -0.000538894259811931, -0.000523252423509115, -0.000507642405972521,
    -0.000492024291725858, -0.000476456042465615, -0.000460893486056534,
    -0.000445402498835314, -0.000429943114160157, -0.000414577645869965,
    -0.000399266350183937, -0.000384067514638573, -0.000368940000057705,
    -0.0003539411416973, -0.000339293525390544, -0.000323659926009783,
    -0.000309919414694471, -0.000295511902559795, -0.000280950954741698,
    -0.000266585641496571, -0.000252470845359277, -0.000238655847061092,
    -0.000225058353634012, -0.000211686103776358, -0.000198461680990827,
    -0.000185423796157919, -0.000172526661615128, -0.000159835502885448,
    -0.000147323340028758, -0.000135054315559792, -0.000122985841755254,
    -0.000111172538004551, -9.95753499456616E-5, -8.82595020734168E-5,
    -7.7184238904831E-5, -6.63771677159197E-5, -5.57696549960879E-5,
    -4.54400540476134E-5, -3.53685936727762E-5, -2.55059336623784E-5,
    -1.59122093433135E-5, -6.56931307847807E-6, 2.54292500196256E-6,
    1.1382115264166E-5, 1.99914518409441E-5, 2.83209585120617E-5,
    3.64139273452416E-5, 4.42231681996236E-5, 5.17992381686803E-5,
    5.9088983210048E-5, 6.61402958740868E-5, 7.2904053299632E-5,
    7.94321200124609E-5, 8.5679575283788E-5, 9.16972189395753E-5,
    9.74338509893395E-5, 0.000102945185568431, 0.00010818882774188,
    0.000113210405380578, 0.000117960899953838, 0.000122502135281101,
    0.000126785208378037, 0.000130850035109437, 0.000134668974387507,
    0.000138276287207714, 0.00014163973860547, 0.000144804508496339,
    0.000147727128052281, 0.00015045433430142, 0.000152947120543503,
    0.000155251833704162, 0.000157324536310168, 0.000159218674345988,
    0.000160885927205586, 0.000162382142477805, 0.000163662215298123,
    0.000164776375709852, 0.000165681305169779, 0.000166435386587554,
    0.000166989668474572, 0.000167395706802922, 0.000167612685064898,
    0.00016769546361647, 0.000167596209337449, 0.000167371723214612,
    0.000166977599482095, 0.000166463912498038, 0.000165793466249417,
    0.000165012261291061, 0.000164080891749098, 0.000163049900688652,
    0.000161881818578263, 0.000160617329536024, 0.000159227356866869,
    0.000157751952845765, 0.000156153810147281, 0.000154490393439512,
    0.000152708088942437, 0.000150860785731569, 0.000148909574313128,
    0.000146905389952095, 0.000144799075700251, 0.000142648282511158,
    0.00014040937231192, 0.000138133925206034, 0.000135776345148193,
    0.000133391093408533, 0.000130933154839445, 0.000128458153510528,
    0.000125918478995377, 0.000123362723034365, 0.000120750142158453,
    0.000118132310510059, 0.000115468005236082, 0.000112800194579151,
    0.000110108789624883, 0.000107392581745282, 0.000104680883398715,
    0.000101967709425591, 9.92159696337733E-5, 9.6501764344623E-5,
    9.37690208083129E-5, 9.10537173189225E-5, 8.83159101463416E-5,
    8.56184309621784E-5, 8.29159682127781E-5, 8.02551921214269E-5,
    7.75832509210173E-5, 7.49510402638714E-5, 7.23168342051082E-5,
    6.97343785522134E-5, 6.71538849242156E-5, 6.46252421778573E-5,
    6.21002944473329E-5, 5.96287108106926E-5, 5.7167569336772E-5,
    5.47674639762573E-5, 5.23815850259164E-5, 5.00586518032772E-5,
    0.000524469736335486 };

  static const double dv3[703] = { 0.00044218856248694, -0.000110657522425279,
    -9.98455797241236E-5, -9.14897478869632E-5, -8.52361666827068E-5,
    -8.06959517311162E-5, -7.76411712637086E-5, -7.57899038105713E-5,
    -7.49994662714565E-5, -7.50603085922735E-5, -7.58901010430049E-5,
    -7.73273137334636E-5, -7.9324611916843E-5, -8.17496582065717E-5,
    -8.45883226012555E-5, -8.77382876114558E-5, -9.11922203705935E-5,
    -9.48664849214488E-5, -9.87761589957704E-5, -0.00010284417598808,
    -0.000107090784794005, -0.000111440909688352, -0.000115928813200131,
    -0.000120477160220656, -0.000125117715117357, -0.000129794392580865,
    -0.000134536155367815, -0.000139275765199167, -0.000144062787810614,
    -0.000148824518200912, -0.000153595780710076, -0.000158325434209622,
    -0.000163036019066201, -0.000167679678011588, -0.000172283561123639,
    -0.000176796250544662, -0.000181250958338245, -0.00018558383578915,
    -0.000189833189429288, -0.000193936404147956, -0.000197929142837483,
    -0.000201752638548829, -0.000205448168113863, -0.000208952400045264,
    -0.000212298526996463, -0.000215431951730521, -0.000218370006534548,
    -0.000221087779149871, -0.000223568899355575, -0.000225816461603436,
    -0.00022782334706094, -0.000229527953926879, -0.000230995790176611,
    -0.000232140202409209, -0.00023300079340031, -0.000233539688331506,
    -0.00023378633527504, -0.000233659708473604, -0.000233206053130778,
    -0.000232378805458336, -0.000231219358034977, -0.000229656785942769,
    -0.000227727315820203, -0.000225380645567152, -0.00022265859544934,
    -0.00021950906281289, -0.000215967665766597, -0.000211968448558799,
    -0.000207561615470521, -0.000202704140518313, -0.000197438201725542,
    -0.00019168680001509, -0.000185498343418373, -0.000178850781994576,
    -0.000171780962719214, -0.00016419797671807, -0.000156190026970194,
    -0.000147717289209701, -0.00013876270798387, -0.000129375123024395,
    -0.000119521441536058, -0.000109202912399607, -9.84484818706821E-5,
    -8.72242658070373E-5, -7.55832297779589E-5, -6.3480682068005E-5,
    -5.0974193619924E-5, -3.80211853750399E-5, -2.46764348809454E-5,
    -1.09007929801398E-5, 3.24223471168255E-6, 1.77916103388651E-5,
    3.26898421180593E-5, 4.79656255396375E-5, 6.3566002270918E-5,
    7.95144815924272E-5, 9.57553359401E-5, 0.000112311171021065,
    0.000129115507243385, 0.000146208523407364, 0.000163510711781697,
    0.000181055990914089, 0.000198779793426531, 0.000216700488791655,
    0.000234745799216152, 0.00025294821508132, 0.00027122522848054,
    0.00028960230791491, 0.00030800663658308, 0.000326458055947502,
    0.000344882743063236, 0.000363296018691483, 0.000381627256894602,
    0.000399890129212652, 0.000418002629144219, 0.000435987929285981,
    0.000453767065420231, 0.000471355470721825, 0.000488669514310265,
    0.000505726649129089, 0.000522444131355595, 0.000538838194008873,
    0.000554824556931667, 0.000570381891315664, 0.00058574363951865,
    0.000599712098309406, 0.000614748659898585, 0.000628041748771862,
    0.000640655606610034, 0.000653112783043952, 0.000665075491543359,
    0.000676258910887898, 0.000686503321973215, 0.000695927732178446,
    0.000704573751514583, 0.000712526700146433, 0.000719690675126741,
    0.000726032690090996, 0.000731413762698865, 0.000735829824257724,
    0.000739201036584629, 0.000741561040118516, 0.000742864552642016,
    0.00074316335541823, 0.000742393457652286, 0.000740580422933247,
    0.000737638852671666, 0.000733592725282046, 0.00072834588241453,
    0.000721918005128956, 0.000714242597601836, 0.000705348970537505,
    0.000695168304813553, 0.000683764815943311, 0.000671073326074807,
    0.000657138705447006, 0.000641920386361616, 0.000625443291086012,
    0.000607663334750444, 0.000588610658995366, 0.000568235949272277,
    0.000546583747008114, 0.000523600383315921, 0.000499344029128501,
    0.000473767964040842, 0.000446928977698078, 0.000418784086169086,
    0.000389399293243541, 0.000358735991847034, 0.000326863787781629,
    0.000293764207557915, 0.000259496360454215, 0.000224062935335397,
    0.00018749540502519, 0.00014981249517778, 0.000111067243459504,
    7.12362989078123E-5, 3.04355904575503E-5, -1.13692218873472E-5,
    -5.40934650222404E-5, -9.77379655364106E-5, -0.000142218358284181,
    -0.000187553808130795, -0.000233629127780396, -0.000280437394148562,
    -0.000327882387234898, -0.000375979311793093, -0.000424620534910132,
    -0.00047378396698472, -0.000523354431241435, -0.000573322387778651,
    -0.000623584537328208, -0.000674136790271356, -0.000724845720752713,
    -0.000775679872534808, -0.000826525753613857, -0.00087739182426909,
    -0.00092814237722969, -0.000978713698804689, -0.00102899532653899,
    -0.00107902169607718, -0.00112860084065069, -0.0011776693081902,
    -0.00122622838046955, -0.00127408587839297, -0.00132123592440593,
    -0.00136757653727756, -0.00141299760787155, -0.00145745599673916,
    -0.0015008180684154, -0.00154305312456678, -0.00158402096786208,
    -0.00162368591796622, -0.00166191133409511, -0.00169866257719637,
    -0.00173379937250004, -0.00176729056400442, -0.00179900352340403,
    -0.00182889587715962, -0.0018568395050691, -0.00188279319885088,
    -0.00190663066128216, -0.00192831551866106, -0.00194771439987576,
    -0.00196480771427863, -0.00197945287606819, -0.00199162656239649,
    -0.00200120800579529, -0.00200816746672251, -0.00201238080748018,
    -0.00201383794301468, -0.00201240896567734, -0.00200808023856785,
    -0.00200073871634229, -0.00199036976433981, -0.00197686642391883,
    -0.00196021582547589, -0.0019403192293783, -0.00191717081627876,
    -0.00189066708440931, -0.00186081656382531, -0.00182752251490315,
    -0.00179079178272426, -0.00175053388265721, -0.00170676987426391,
    -0.00165941766709796, -0.00160850114906796, -0.00155394758862404,
    -0.00149574927976897, -0.00143399143098721, -0.00136831299612247,
    -0.00129932561355133, -0.00122639037844494, -0.00114972371991373,
    -0.00106955335669744, -0.000985669737727666, -0.000898040252148599,
    -0.000806691221734975, -0.000711774350863481, -0.000613254881865684,
    -0.000511158423926788, -0.00040542401537393, -0.000296148158273841,
    -0.000183347113779057, -6.71406343458254E-5, 5.24827105073237E-5,
    0.000175453268935852, 0.00030178436157789, 0.000431367973228352,
    0.000564188884829856, 0.000700117036511363, 0.000839124467772486,
    0.000981073903589326, 0.00112595807659575, 0.00127365802791523,
    0.00142414162003024, 0.00157727829340564, 0.0017330385315845,
    0.00189124867543978, 0.00205186680827992, 0.0022147422178615,
    0.00237979773402877, 0.00254690057435, 0.00271597700409001,
    0.00288688153935119, 0.00305954623686013, 0.00323380501568201,
    0.00340958986361004, 0.00358671627162477, 0.00376510452014456,
    0.00394456924668566, 0.00412502877604286, 0.00430629431747686,
    0.00448828951085397, 0.00467083359079199, 0.00485383705373487,
    0.00503712657977312, 0.00522058688840571, 0.00540405069335647,
    0.00558738694575325, 0.00577040966215496, 0.00595303526933757,
    0.00613504619542937, 0.00631636903698379, 0.00649681040951683,
    0.00667625932295169, 0.00685451323038891, 0.00703148439372521,
    0.00720696612105708, 0.00738087270312673, 0.00755301874920811,
    0.0077233165064257, 0.00789155222299408, 0.00805761456478647,
    0.00822131875322062, 0.00838259356390575, 0.00854125907927656,
    0.0086972234431769, 0.00885027046921507, 0.00900031283405329,
    0.00914719655546697, 0.00929088009437544, 0.00943115585333465,
    0.00956789639043661, 0.00970096259291177, 0.0098303833523068,
    0.00995589630104736, 0.0100773480018851, 0.0101948546792788,
    0.0103080638144166, 0.0104170283034853, 0.0105216018540367,
    0.0106216692956179, 0.0107171859212136, 0.010808016098072,
    0.0108941318936983, 0.0109753996703859, 0.0110517981558641,
    0.0111232094631204, 0.0111896227087287, 0.0112509221691577,
    0.0113071058431389, 0.0113580739117312, 0.0114038233874121,
    0.0114442687045384, 0.0114794193958816, 0.0115091971892624,
    0.0115336216751115, 0.0115526155187787, 0.0115662222486294,
    0.0115743673395176, 0.0115770980653248, 0.0115743673395176,
    0.0115662222486294, 0.0115526155187787, 0.0115336216751115,
    0.0115091971892624, 0.0114794193958816, 0.0114442687045384,
    0.0114038233874121, 0.0113580739117312, 0.0113071058431389,
    0.0112509221691577, 0.0111896227087287, 0.0111232094631204,
    0.0110517981558641, 0.0109753996703859, 0.0108941318936983,
    0.010808016098072, 0.0107171859212136, 0.0106216692956179,
    0.0105216018540367, 0.0104170283034853, 0.0103080638144166,
    0.0101948546792788, 0.0100773480018851, 0.00995589630104736,
    0.0098303833523068, 0.00970096259291177, 0.00956789639043661,
    0.00943115585333465, 0.00929088009437544, 0.00914719655546697,
    0.00900031283405329, 0.00885027046921507, 0.0086972234431769,
    0.00854125907927656, 0.00838259356390575, 0.00822131875322062,
    0.00805761456478647, 0.00789155222299408, 0.0077233165064257,
    0.00755301874920811, 0.00738087270312673, 0.00720696612105708,
    0.00703148439372521, 0.00685451323038891, 0.00667625932295169,
    0.00649681040951683, 0.00631636903698379, 0.00613504619542937,
    0.00595303526933757, 0.00577040966215496, 0.00558738694575325,
    0.00540405069335647, 0.00522058688840571, 0.00503712657977312,
    0.00485383705373487, 0.00467083359079199, 0.00448828951085397,
    0.00430629431747686, 0.00412502877604286, 0.00394456924668566,
    0.00376510452014456, 0.00358671627162477, 0.00340958986361004,
    0.00323380501568201, 0.00305954623686013, 0.00288688153935119,
    0.00271597700409001, 0.00254690057435, 0.00237979773402877,
    0.0022147422178615, 0.00205186680827992, 0.00189124867543978,
    0.0017330385315845, 0.00157727829340564, 0.00142414162003024,
    0.00127365802791523, 0.00112595807659575, 0.000981073903589326,
    0.000839124467772486, 0.000700117036511363, 0.000564188884829856,
    0.000431367973228352, 0.00030178436157789, 0.000175453268935852,
    5.24827105073237E-5, -6.71406343458254E-5, -0.000183347113779057,
    -0.000296148158273841, -0.00040542401537393, -0.000511158423926788,
    -0.000613254881865684, -0.000711774350863481, -0.000806691221734975,
    -0.000898040252148599, -0.000985669737727666, -0.00106955335669744,
    -0.00114972371991373, -0.00122639037844494, -0.00129932561355133,
    -0.00136831299612247, -0.00143399143098721, -0.00149574927976897,
    -0.00155394758862404, -0.00160850114906796, -0.00165941766709796,
    -0.00170676987426391, -0.00175053388265721, -0.00179079178272426,
    -0.00182752251490315, -0.00186081656382531, -0.00189066708440931,
    -0.00191717081627876, -0.0019403192293783, -0.00196021582547589,
    -0.00197686642391883, -0.00199036976433981, -0.00200073871634229,
    -0.00200808023856785, -0.00201240896567734, -0.00201383794301468,
    -0.00201238080748018, -0.00200816746672251, -0.00200120800579529,
    -0.00199162656239649, -0.00197945287606819, -0.00196480771427863,
    -0.00194771439987576, -0.00192831551866106, -0.00190663066128216,
    -0.00188279319885088, -0.0018568395050691, -0.00182889587715962,
    -0.00179900352340403, -0.00176729056400442, -0.00173379937250004,
    -0.00169866257719637, -0.00166191133409511, -0.00162368591796622,
    -0.00158402096786208, -0.00154305312456678, -0.0015008180684154,
    -0.00145745599673916, -0.00141299760787155, -0.00136757653727756,
    -0.00132123592440593, -0.00127408587839297, -0.00122622838046955,
    -0.0011776693081902, -0.00112860084065069, -0.00107902169607718,
    -0.00102899532653899, -0.000978713698804689, -0.00092814237722969,
    -0.00087739182426909, -0.000826525753613857, -0.000775679872534808,
    -0.000724845720752713, -0.000674136790271356, -0.000623584537328208,
    -0.000573322387778651, -0.000523354431241435, -0.00047378396698472,
    -0.000424620534910132, -0.000375979311793093, -0.000327882387234898,
    -0.000280437394148562, -0.000233629127780396, -0.000187553808130795,
    -0.000142218358284181, -9.77379655364106E-5, -5.40934650222404E-5,
    -1.13692218873472E-5, 3.04355904575503E-5, 7.12362989078123E-5,
    0.000111067243459504, 0.00014981249517778, 0.00018749540502519,
    0.000224062935335397, 0.000259496360454215, 0.000293764207557915,
    0.000326863787781629, 0.000358735991847034, 0.000389399293243541,
    0.000418784086169086, 0.000446928977698078, 0.000473767964040842,
    0.000499344029128501, 0.000523600383315921, 0.000546583747008114,
    0.000568235949272277, 0.000588610658995366, 0.000607663334750444,
    0.000625443291086012, 0.000641920386361616, 0.000657138705447006,
    0.000671073326074807, 0.000683764815943311, 0.000695168304813553,
    0.000705348970537505, 0.000714242597601836, 0.000721918005128956,
    0.00072834588241453, 0.000733592725282046, 0.000737638852671666,
    0.000740580422933247, 0.000742393457652286, 0.00074316335541823,
    0.000742864552642016, 0.000741561040118516, 0.000739201036584629,
    0.000735829824257724, 0.000731413762698865, 0.000726032690090996,
    0.000719690675126741, 0.000712526700146433, 0.000704573751514583,
    0.000695927732178446, 0.000686503321973215, 0.000676258910887898,
    0.000665075491543359, 0.000653112783043952, 0.000640655606610034,
    0.000628041748771862, 0.000614748659898585, 0.000599712098309406,
    0.00058574363951865, 0.000570381891315664, 0.000554824556931667,
    0.000538838194008873, 0.000522444131355595, 0.000505726649129089,
    0.000488669514310265, 0.000471355470721825, 0.000453767065420231,
    0.000435987929285981, 0.000418002629144219, 0.000399890129212652,
    0.000381627256894602, 0.000363296018691483, 0.000344882743063236,
    0.000326458055947502, 0.00030800663658308, 0.00028960230791491,
    0.00027122522848054, 0.00025294821508132, 0.000234745799216152,
    0.000216700488791655, 0.000198779793426531, 0.000181055990914089,
    0.000163510711781697, 0.000146208523407364, 0.000129115507243385,
    0.000112311171021065, 9.57553359401E-5, 7.95144815924272E-5,
    6.3566002270918E-5, 4.79656255396375E-5, 3.26898421180593E-5,
    1.77916103388651E-5, 3.24223471168255E-6, -1.09007929801398E-5,
    -2.46764348809454E-5, -3.80211853750399E-5, -5.0974193619924E-5,
    -6.3480682068005E-5, -7.55832297779589E-5, -8.72242658070373E-5,
    -9.84484818706821E-5, -0.000109202912399607, -0.000119521441536058,
    -0.000129375123024395, -0.00013876270798387, -0.000147717289209701,
    -0.000156190026970194, -0.00016419797671807, -0.000171780962719214,
    -0.000178850781994576, -0.000185498343418373, -0.00019168680001509,
    -0.000197438201725542, -0.000202704140518313, -0.000207561615470521,
    -0.000211968448558799, -0.000215967665766597, -0.00021950906281289,
    -0.00022265859544934, -0.000225380645567152, -0.000227727315820203,
    -0.000229656785942769, -0.000231219358034977, -0.000232378805458336,
    -0.000233206053130778, -0.000233659708473604, -0.00023378633527504,
    -0.000233539688331506, -0.00023300079340031, -0.000232140202409209,
    -0.000230995790176611, -0.000229527953926879, -0.00022782334706094,
    -0.000225816461603436, -0.000223568899355575, -0.000221087779149871,
    -0.000218370006534548, -0.000215431951730521, -0.000212298526996463,
    -0.000208952400045264, -0.000205448168113863, -0.000201752638548829,
    -0.000197929142837483, -0.000193936404147956, -0.000189833189429288,
    -0.00018558383578915, -0.000181250958338245, -0.000176796250544662,
    -0.000172283561123639, -0.000167679678011588, -0.000163036019066201,
    -0.000158325434209622, -0.000153595780710076, -0.000148824518200912,
    -0.000144062787810614, -0.000139275765199167, -0.000134536155367815,
    -0.000129794392580865, -0.000125117715117357, -0.000120477160220656,
    -0.000115928813200131, -0.000111440909688352, -0.000107090784794005,
    -0.00010284417598808, -9.87761589957704E-5, -9.48664849214488E-5,
    -9.11922203705935E-5, -8.77382876114558E-5, -8.45883226012555E-5,
    -8.17496582065717E-5, -7.9324611916843E-5, -7.73273137334636E-5,
    -7.58901010430049E-5, -7.50603085922735E-5, -7.49994662714565E-5,
    -7.57899038105713E-5, -7.76411712637086E-5, -8.06959517311162E-5,
    -8.52361666827068E-5, -9.14897478869632E-5, -9.98455797241236E-5,
    -0.000110657522425279, 0.00044218856248694 };

  static const double dv4[699] = { 0.000532723527805931, 6.90845510945624E-5,
    7.32098567809705E-5, 7.7282881735165E-5, 8.1227965161615E-5,
    8.50730286220366E-5, 8.87432006630032E-5, 9.22632831417983E-5,
    9.55562989757991E-5, 9.86445413311909E-5, 0.000101452080284138,
    0.000104006770866074, 0.000106234517995692, 0.000108156839530409,
    0.00010969767408967, 0.000110875128450003, 0.000111621547646505,
    0.000111970687838987, 0.00011184988429155, 0.000111282242197444,
    0.000110192908129537, 0.000108623081333228, 0.000106515908723708,
    0.000103892652709569, 0.00010067735487843, 9.69389345979157E-5,
    9.26226107885545E-5, 8.77215523829069E-5, 8.22521094557147E-5,
    7.62066042336095E-5, 6.95679622187701E-5, 6.23782262138511E-5,
    5.4595564470071E-5, 4.62793390590956E-5, 3.73910309205438E-5,
    2.79965806474123E-5, 1.80639059740427E-5, 7.66317268571757E-6,
    -3.2342105979732E-6, -1.45519362022193E-5, -2.63058607512367E-5,
    -3.84220186694711E-5, -5.09117413208545E-5, -6.36905034503717E-5,
    -7.67661277548565E-5, -9.00437851570947E-5, -0.000103535649615391,
    -0.000117140359067039, -0.000130853687142025, -0.000144579619052287,
    -0.000158316330340026, -0.000171948811459607, -0.000185476850131979,
    -0.000198794093871137, -0.000211878668856251, -0.000224633417558677,
    -0.000237034781155883, -0.000248978615566231, -0.000260448698693465,
    -0.000271332401327684, -0.000281620445417402, -0.000291202506939356,
    -0.000300068275069311, -0.000308108733299146, -0.000315319515874973,
    -0.000321597411331245, -0.000326934577519601, -0.000331241670082228,
    -0.000334505681801338, -0.00033665013086722, -0.000337671660988197,
    -0.000337481977392415, -0.000336106485795569, -0.000333455311220054,
    -0.000329550165379097, -0.000324325703452802, -0.00031781195342509,
    -0.000309939057334787, -0.000300754241979969, -0.00029020554346551,
    -0.000278344626220354, -0.000265132115339844, -0.000250628169585483,
    -0.000234812182294999, -0.000217759861217467, -0.000199454100181876,
    -0.000179986263437099, -0.000159353716962562, -0.000137659970668377,
    -0.000114912452474471, -9.12241339784792E-5, -6.66119812977303E-5,
    -4.11936497765334E-5, -1.50167581043864E-5, 1.18423295007193E-5,
    3.92565201720357E-5, 6.71914198631196E-5, 9.55422954605711E-5,
    0.000124142970867157, 0.0001529738906533, 0.0001818706532063,
    0.000210742122882789, 0.000239413971555361, 0.000267824245622387,
    0.000295802383646803, 0.000323268902375579, 0.000350044687543034,
    0.000376064507731952, 0.000401179043921861, 0.000425332695313922,
    0.000448370380215483, 0.000470221483334989, 0.000490702014082546,
    0.000509715635440249, 0.000527055147364876, 0.000542619714948337,
    0.000556249759298917, 0.000567969650086771, 0.000577807639663091,
    0.000585912834301797, 0.000591969798291343, 0.000594305253361897,
    0.000595860696550544, 0.000594399332766709, 0.000590571980802667,
    0.000584290296204713, 0.00057546378154099, 0.000564118894623614,
    0.000550184412602726, 0.000533706380462012, 0.000514638074566508,
    0.000493050047060811, 0.000468913400443583, 0.000442313253832272,
    0.00041324132878444, 0.000381803274818962, 0.000348023547052433,
    0.000312034702218758, 0.000273870106753618, 0.000233681198056428,
    0.000191528676530613, 0.000147596411187783, 0.00010196028738839,
    5.48016986829699E-5, 6.22557969660877E-6, -4.3536285407423E-5,
    -9.4388231957352E-5, -0.000146110116296165, -0.000198526522947237,
    -0.00025144485083223, -0.000304678802873628, -0.000357993619188933,
    -0.000411225582184539, -0.000464109785038842, -0.000516482314140015,
    -0.000568068863470274, -0.000618700274194817, -0.000668100604982624,
    -0.000716098863454397, -0.000762415967824261, -0.000806880267618045,
    -0.000849226128650888, -0.000889282065659916, -0.000926788086023849,
    -0.000961585620209408, -0.000993422610739849, -0.00102216328155141,
    -0.00104756127527116, -0.00106949414293076, -0.00108774785205145,
    -0.00110221455827652, -0.00111269268098981, -0.00111911216097625,
    -0.00112129581941446, -0.00111918690073059, -0.00111265273925094,
    -0.00110165907175774, -0.00108610324359315, -0.0010659873359159,
    -0.00104123445629234, -0.00101189035018557, -0.00097790972708705,
    -0.000939371894413249, -0.000896270612403707, -0.000848725929165545,
    -0.000796768602540623, -0.000740552848045237, -0.000680156633333393,
    -0.000615760383245651, -0.000547494664223412, -0.000475569076402148,
    -0.000400142017380367, -0.000321485002025148, -0.000239764850790327,
    -0.000155292447891851, -6.82857751072828E-5, 2.09289179787751E-5,
    0.00011210929762225, 0.000204883420977705, 0.000298981574344912,
    0.000394014750820599, 0.0004896859683268, 0.000585583218912135,
    0.000681382407825788, 0.000776660024773496, 0.000871085500740466,
    0.000964216924345655, 0.0010557159369007, 0.00114513675352477,
    0.0012321401912473, 0.00131628101314346, 0.00139722080626986,
    0.00147451648209021, 0.00154782628278128, 0.0016167564868366,
    0.00168089004027324, 0.00173997814549448, 0.00179355111402787,
    0.00184129196834712, 0.00188297511181281, 0.00191820056831365,
    0.00194671777139435, 0.00196823189611401, 0.0019825664814632,
    0.00198943348273961, 0.0019886606874511, 0.00198000124960408,
    0.00196335543997715, 0.00193855339959633, 0.00190555323215255,
    0.00186424362293629, 0.00181463879312254, 0.00175667524340481,
    0.00169041467832458, 0.00161583185530486, 0.00153304327567335,
    0.00144212179286118, 0.00134331731841221, 0.00123681071188208,
    0.00112283496099542, 0.00100142115160551, 0.000872866227640752,
    0.000738010318540991, 0.000596478847215025, 0.000449101772719016,
    0.000296158967640106, 0.000138068317705179, -2.46513410261154E-5,
    -0.000191559202568725, -0.000362075082993075, -0.000535698168265624,
    -0.000711795700978104, -0.000889822901721817, -0.00106910560605327,
    -0.00124905648482067, -0.00142896209429875, -0.00160819122016243,
    -0.00178600345365117, -0.00196174870426848, -0.00213465415341023,
    -0.00230404356324101, -0.00246912606313768, -0.00262922616543297,
    -0.00278355295102357, -0.00293140840277156, -0.00307200120028402,
    -0.00320467215722319, -0.00332862998939508, -0.0034432003952259,
    -0.00354766165289471, -0.00364133959867959, -0.00372353791259458,
    -0.00379364301321045, -0.0038509756885853, -0.00389498262148901,
    -0.0039250291459149, -0.00394061946397722, -0.00394117996376844,
    -0.00392627630440225, -0.00389539996021668, -0.00384818977002379,
    -0.00378421825818613, -0.00370319322404993, -0.00360477304508528,
    -0.00348875245919541, -0.00335487603222865, -0.00320303723850247,
    -0.00303306169439496, -0.00284494561585071, -0.00263862046510879,
    -0.00241416824095201, -0.0021716226242081, -0.00191117876624268,
    -0.00163296066456587, -0.00133726271460791, -0.00102432543856643,
    -0.000694533535606485, -0.000348233195879792, 1.40895738088142E-5,
    0.000391991534902159, 0.000784879339972465, 0.00119222116184219,
    0.00161332794293357, 0.00204757197540031, 0.00249417437950656,
    0.00295242425888869, 0.00342146101503524, 0.0039004895427804,
    0.00438858491491148, 0.00488486434194883, 0.00538835584449848,
    0.00589810685900003, 0.00641307176552732, 0.0069322805889733,
    0.00745460951227381, 0.00797904473983032, 0.00850445114532105,
    0.00902977164682324, 0.00955383211944169, 0.0100755711107632,
    0.010593809222233, 0.0111074713522684, 0.0116153757129336,
    0.0121164515317209, 0.0126095392988345, 0.0130935850412763,
    0.0135674431488346, 0.0140300933926742, 0.0144804281527954,
    0.0149174721511036, 0.0153401670861071, 0.0157475906337756,
    0.0161387441993727, 0.0165127735255988, 0.0168687716270319,
    0.0172058854778318, 0.0175234683086044, 0.0178204989017597,
    0.0180965132230385, 0.0183508338529866, 0.0185827142905703,
    0.0187916926310988, 0.0189772718951223, 0.0191390675550726,
    0.0192766147002736, 0.0193896183087324, 0.0194777543833272,
    0.0195408769319307, 0.0195787891742262, 0.0195914531025123,
    0.0195787891742262, 0.0195408769319307, 0.0194777543833272,
    0.0193896183087324, 0.0192766147002736, 0.0191390675550726,
    0.0189772718951223, 0.0187916926310988, 0.0185827142905703,
    0.0183508338529866, 0.0180965132230385, 0.0178204989017597,
    0.0175234683086044, 0.0172058854778318, 0.0168687716270319,
    0.0165127735255988, 0.0161387441993727, 0.0157475906337756,
    0.0153401670861071, 0.0149174721511036, 0.0144804281527954,
    0.0140300933926742, 0.0135674431488346, 0.0130935850412763,
    0.0126095392988345, 0.0121164515317209, 0.0116153757129336,
    0.0111074713522684, 0.010593809222233, 0.0100755711107632,
    0.00955383211944169, 0.00902977164682324, 0.00850445114532105,
    0.00797904473983032, 0.00745460951227381, 0.0069322805889733,
    0.00641307176552732, 0.00589810685900003, 0.00538835584449848,
    0.00488486434194883, 0.00438858491491148, 0.0039004895427804,
    0.00342146101503524, 0.00295242425888869, 0.00249417437950656,
    0.00204757197540031, 0.00161332794293357, 0.00119222116184219,
    0.000784879339972465, 0.000391991534902159, 1.40895738088142E-5,
    -0.000348233195879792, -0.000694533535606485, -0.00102432543856643,
    -0.00133726271460791, -0.00163296066456587, -0.00191117876624268,
    -0.0021716226242081, -0.00241416824095201, -0.00263862046510879,
    -0.00284494561585071, -0.00303306169439496, -0.00320303723850247,
    -0.00335487603222865, -0.00348875245919541, -0.00360477304508528,
    -0.00370319322404993, -0.00378421825818613, -0.00384818977002379,
    -0.00389539996021668, -0.00392627630440225, -0.00394117996376844,
    -0.00394061946397722, -0.0039250291459149, -0.00389498262148901,
    -0.0038509756885853, -0.00379364301321045, -0.00372353791259458,
    -0.00364133959867959, -0.00354766165289471, -0.0034432003952259,
    -0.00332862998939508, -0.00320467215722319, -0.00307200120028402,
    -0.00293140840277156, -0.00278355295102357, -0.00262922616543297,
    -0.00246912606313768, -0.00230404356324101, -0.00213465415341023,
    -0.00196174870426848, -0.00178600345365117, -0.00160819122016243,
    -0.00142896209429875, -0.00124905648482067, -0.00106910560605327,
    -0.000889822901721817, -0.000711795700978104, -0.000535698168265624,
    -0.000362075082993075, -0.000191559202568725, -2.46513410261154E-5,
    0.000138068317705179, 0.000296158967640106, 0.000449101772719016,
    0.000596478847215025, 0.000738010318540991, 0.000872866227640752,
    0.00100142115160551, 0.00112283496099542, 0.00123681071188208,
    0.00134331731841221, 0.00144212179286118, 0.00153304327567335,
    0.00161583185530486, 0.00169041467832458, 0.00175667524340481,
    0.00181463879312254, 0.00186424362293629, 0.00190555323215255,
    0.00193855339959633, 0.00196335543997715, 0.00198000124960408,
    0.0019886606874511, 0.00198943348273961, 0.0019825664814632,
    0.00196823189611401, 0.00194671777139435, 0.00191820056831365,
    0.00188297511181281, 0.00184129196834712, 0.00179355111402787,
    0.00173997814549448, 0.00168089004027324, 0.0016167564868366,
    0.00154782628278128, 0.00147451648209021, 0.00139722080626986,
    0.00131628101314346, 0.0012321401912473, 0.00114513675352477,
    0.0010557159369007, 0.000964216924345655, 0.000871085500740466,
    0.000776660024773496, 0.000681382407825788, 0.000585583218912135,
    0.0004896859683268, 0.000394014750820599, 0.000298981574344912,
    0.000204883420977705, 0.00011210929762225, 2.09289179787751E-5,
    -6.82857751072828E-5, -0.000155292447891851, -0.000239764850790327,
    -0.000321485002025148, -0.000400142017380367, -0.000475569076402148,
    -0.000547494664223412, -0.000615760383245651, -0.000680156633333393,
    -0.000740552848045237, -0.000796768602540623, -0.000848725929165545,
    -0.000896270612403707, -0.000939371894413249, -0.00097790972708705,
    -0.00101189035018557, -0.00104123445629234, -0.0010659873359159,
    -0.00108610324359315, -0.00110165907175774, -0.00111265273925094,
    -0.00111918690073059, -0.00112129581941446, -0.00111911216097625,
    -0.00111269268098981, -0.00110221455827652, -0.00108774785205145,
    -0.00106949414293076, -0.00104756127527116, -0.00102216328155141,
    -0.000993422610739849, -0.000961585620209408, -0.000926788086023849,
    -0.000889282065659916, -0.000849226128650888, -0.000806880267618045,
    -0.000762415967824261, -0.000716098863454397, -0.000668100604982624,
    -0.000618700274194817, -0.000568068863470274, -0.000516482314140015,
    -0.000464109785038842, -0.000411225582184539, -0.000357993619188933,
    -0.000304678802873628, -0.00025144485083223, -0.000198526522947237,
    -0.000146110116296165, -9.4388231957352E-5, -4.3536285407423E-5,
    6.22557969660877E-6, 5.48016986829699E-5, 0.00010196028738839,
    0.000147596411187783, 0.000191528676530613, 0.000233681198056428,
    0.000273870106753618, 0.000312034702218758, 0.000348023547052433,
    0.000381803274818962, 0.00041324132878444, 0.000442313253832272,
    0.000468913400443583, 0.000493050047060811, 0.000514638074566508,
    0.000533706380462012, 0.000550184412602726, 0.000564118894623614,
    0.00057546378154099, 0.000584290296204713, 0.000590571980802667,
    0.000594399332766709, 0.000595860696550544, 0.000594305253361897,
    0.000591969798291343, 0.000585912834301797, 0.000577807639663091,
    0.000567969650086771, 0.000556249759298917, 0.000542619714948337,
    0.000527055147364876, 0.000509715635440249, 0.000490702014082546,
    0.000470221483334989, 0.000448370380215483, 0.000425332695313922,
    0.000401179043921861, 0.000376064507731952, 0.000350044687543034,
    0.000323268902375579, 0.000295802383646803, 0.000267824245622387,
    0.000239413971555361, 0.000210742122882789, 0.0001818706532063,
    0.0001529738906533, 0.000124142970867157, 9.55422954605711E-5,
    6.71914198631196E-5, 3.92565201720357E-5, 1.18423295007193E-5,
    -1.50167581043864E-5, -4.11936497765334E-5, -6.66119812977303E-5,
    -9.12241339784792E-5, -0.000114912452474471, -0.000137659970668377,
    -0.000159353716962562, -0.000179986263437099, -0.000199454100181876,
    -0.000217759861217467, -0.000234812182294999, -0.000250628169585483,
    -0.000265132115339844, -0.000278344626220354, -0.00029020554346551,
    -0.000300754241979969, -0.000309939057334787, -0.00031781195342509,
    -0.000324325703452802, -0.000329550165379097, -0.000333455311220054,
    -0.000336106485795569, -0.000337481977392415, -0.000337671660988197,
    -0.00033665013086722, -0.000334505681801338, -0.000331241670082228,
    -0.000326934577519601, -0.000321597411331245, -0.000315319515874973,
    -0.000308108733299146, -0.000300068275069311, -0.000291202506939356,
    -0.000281620445417402, -0.000271332401327684, -0.000260448698693465,
    -0.000248978615566231, -0.000237034781155883, -0.000224633417558677,
    -0.000211878668856251, -0.000198794093871137, -0.000185476850131979,
    -0.000171948811459607, -0.000158316330340026, -0.000144579619052287,
    -0.000130853687142025, -0.000117140359067039, -0.000103535649615391,
    -9.00437851570947E-5, -7.67661277548565E-5, -6.36905034503717E-5,
    -5.09117413208545E-5, -3.84220186694711E-5, -2.63058607512367E-5,
    -1.45519362022193E-5, -3.2342105979732E-6, 7.66317268571757E-6,
    1.80639059740427E-5, 2.79965806474123E-5, 3.73910309205438E-5,
    4.62793390590956E-5, 5.4595564470071E-5, 6.23782262138511E-5,
    6.95679622187701E-5, 7.62066042336095E-5, 8.22521094557147E-5,
    8.77215523829069E-5, 9.26226107885545E-5, 9.69389345979157E-5,
    0.00010067735487843, 0.000103892652709569, 0.000106515908723708,
    0.000108623081333228, 0.000110192908129537, 0.000111282242197444,
    0.00011184988429155, 0.000111970687838987, 0.000111621547646505,
    0.000110875128450003, 0.00010969767408967, 0.000108156839530409,
    0.000106234517995692, 0.000104006770866074, 0.000101452080284138,
    9.86445413311909E-5, 9.55562989757991E-5, 9.22632831417983E-5,
    8.87432006630032E-5, 8.50730286220366E-5, 8.1227965161615E-5,
    7.7282881735165E-5, 7.32098567809705E-5, 6.90845510945624E-5,
    0.000532723527805931 };

  static const double dv5[861] = { 0.000541983718811927, 9.39155210670181E-5,
    0.000100958414377355, 0.000107595463015365, 0.000113688228537891,
    0.0001191770935068, 0.000123929134733962, 0.000127894687581549,
    0.000130945415051525, 0.000133044030326802, 0.000134078026482925,
    0.000134025874269021, 0.0001327911298552, 0.000130373124118671,
    0.000126701511593129, 0.000121799681582674, 0.000115621965373439,
    0.000108220396625604, 9.95750422745516E-5, 8.97652851564876E-5,
    7.88062325959248E-5, 6.68078245004558E-5, 5.38098574168716E-5,
    3.99522311415092E-5, 2.53040869264146E-5, 1.0030902985957E-5,
    -5.77254316079002E-6, -2.19231126541713E-5, -3.83000510956379E-5,
    -5.4711660316172E-5, -7.10064533753939E-5, -8.70064897013418E-5,
    -0.000102529979431761, -0.000117407979761088, -0.000131481469000169,
    -0.000144534880486384, -0.000156470017010974, -0.000167084529934335,
    -0.000176261149179698, -0.00018383717508416, -0.000189747215029418,
    -0.000193841772251799, -0.000196073541936795, -0.000196346119718796,
    -0.000194665923081674, -0.000190965804054552, -0.000185281259751758,
    -0.000177600347515248, -0.000168020586519177, -0.000156570952883738,
    -0.000143383749251984, -0.000128534875844589, -0.000112203089934051,
    -9.45098568138851E-5, -7.56833010219895E-5, -5.58768211555122E-5,
    -3.53418834782659E-5, -1.42714666447916E-5, 7.06086733384884E-6,
    2.84329633486768E-5, 4.95590136185014E-5, 7.02071138914894E-5,
    9.00914671901327E-5, 0.000108987167847066, 0.000126612317432528,
    0.000142760685517128, 0.000157159733488689, 0.000169645517325974,
    0.000179985578337734, 0.000188051801111294, 0.000193659121415496,
    0.000196723579657171, 0.000197109088634444, 0.000194802675707469,
    0.000189736391522837, 0.000181945744112124, 0.000171419148763144,
    0.000158270235612336, 0.000142560624355696, 0.000124455756549172,
    0.00010407408097609, 8.16637642732029E-5, 5.74156605995404E-5,
    3.16134719413095E-5, 4.49760535647356E-6, -2.35688788073929E-5,
    -5.23127671566118E-5, -8.13640577031261E-5, -0.000110381212543967,
    -0.000139001076920169, -0.00016687275536317, -0.000193622797261251,
    -0.000218930579462653, -0.000242428006148879, -0.000263829521034822,
    -0.000282792065964986, -0.00029907288616221, -0.000312368664535981,
    -0.00032248522846086, -0.000329171873533593, -0.000332299375483446,
    -0.000331700761792967, -0.000327336417669472, -0.000319136839763799,
    -0.000307173721793938, -0.00029148787236693, -0.000272250469485773,
    -0.000249588877743172, -0.000223736461708317, -0.000194860168732485,
    -0.000163234645663197, -0.000129096464368348, -9.28427841098224E-5,
    -5.48945053903947E-5, -1.58536847988518E-5, 2.37618534986138E-5,
    6.35611116325544E-5, 0.000103589939605899, 0.000143704940542049,
    0.000181025976607809, 0.000217931456555117, 0.000252357113836307,
    0.000284257828711435, 0.000313182758504491, 0.000338680260437589,
    0.000360413706777267, 0.000377998077818869, 0.000391181527305812,
    0.000399671906379765, 0.000403315939413201, 0.000401925594888096,
    0.000395460692994347, 0.000383854155585259, 0.00036718765573372,
    0.00034551939954402, 0.000319057917826926, 0.000287981876933402,
    0.000252620897756928, 0.000213277183864515, 0.000170392162781004,
    0.000124370213198492, 7.57533418379797E-5, 2.50361530035493E-5,
    -2.7161145425949E-5, -8.02768371814556E-5, -0.000133640898306514,
    -0.000186647951479871, -0.000238612044359073, -0.000288902070780961,
    -0.00033686907410968, -0.000381865684983965, -0.000423315587749875,
    -0.000460629099728361, -0.000493239085464758, -0.000520728588172605,
    -0.000542597316576951, -0.000558500829399108, -0.000568103825078294,
    -0.000571215822744036, -0.000567612524697722, -0.000557245516783168,
    -0.000540076800221215, -0.000516232780295362, -0.000485819662611941,
    -0.000449113416335854, -0.000406397565443762, -0.000358124919507632,
    -0.000304728818471515, -0.000246799554346799, -0.000184906233867251,
    -0.000119764102369112, -5.20614307239963E-5, 1.73789449575496E-5,
    8.7795776565017E-5, 0.000158312370733294, 0.000228118320448514,
    0.000296318447014287, 0.00036208896972677, 0.000424554418003875,
    0.000482925129093528, 0.000536391099262198, 0.000584245077633135,
    0.000625772979312749, 0.000660385014627503, 0.000687485682784789,
    0.000706653899194518, 0.000717462680602415, 0.000719672891575588,
    0.000713058339562528, 0.000697573224095843, 0.00067319585265931,
    0.000640110647994494, 0.000598517083474895, 0.000548801415044521,
    0.000491374900037892, 0.000426844595346898, 0.00035581891628597,
    0.000279073495233852, 0.0001973870440229, 0.000111712026596431,
    2.29692786579564E-5, -6.77874598097355E-5, -0.000159529303918018,
    -0.000251102810739723, -0.000341447274286012, -0.000429399162350645,
    -0.000513872111887188, -0.000593763701806401, -0.000668037007475755,
    -0.000735672733369139, -0.000795757970971713, -0.000847393382183875,
    -0.000889830045262367, -0.000922339566245509, -0.000944383237710939,
    -0.00095545499522514, -0.000955264249053625, -0.000943564158373563,
    -0.000920339162179491, -0.000885626875170628, -0.000839683905906734,
    -0.00078281274902258, -0.000715536243565728, -0.000638420749418111,
    -0.000552260439891737, -0.000457903260615247, -0.000356409314647111,
    -0.000248851915276288, -0.000136464578717182, -2.04444307204847E-5,
    9.78643173315946E-5, 0.000217112766350226, 0.00033576160747592,
    0.0004523378855645, 0.000565409273805852, 0.000673758101284352,
    0.000775880149805064, 0.000870008572946229, 0.000955579555048026,
    0.00103084522155815, 0.00109493959888634, 0.00114684740229226,
    0.00118568928711056, 0.00121079692202438, 0.00122156944464255,
    0.00121766042770094, 0.00119880438152662, 0.00116500239871139,
    0.00111634393692054, 0.00105319209384983, 0.000976000891026337,
    0.000885492683406536, 0.000782475751035737, 0.000668013805650651,
    0.000543234859783321, 0.000409507246153594, 0.000268243354692785,
    0.000121059821395169, -3.04152432307783E-5, -0.000184379347411724,
    -0.000339059179704988, -0.000492548025635855, -0.000643010218370379,
    -0.000788522783484884, -0.000927273197598017, -0.00105741642025622,
    -0.00117724186294024, -0.00128509963403905, -0.00137943108525202,
    -0.00145889954833536, -0.00152219868611472, -0.0015682621132935,
    -0.00159627209149843, -0.00160546004472087, -0.00159540022296035,
    -0.00156583223311964, -0.00151678445534883, -0.00144841035655361,
    -0.00136121045709383, -0.00125585035904654, -0.00113330657899723,
    -0.000994673305073961, -0.000841350754803934, -0.000674867067745653,
    -0.00049703143195459, -0.000309729372712078, -0.000115084270895235,
    8.47444798042395E-5, 0.000287402789594611, 0.000490536127404496,
    0.000691646548084861, 0.000888309185380098, 0.00107801310643884,
    0.0012583731149921, 0.0014269852599288, 0.00158161297276232,
    0.00172006918497051, 0.00184038233686247, 0.00194070288808355,
    0.00201944312775646, 0.00207518637052296, 0.00210681645356649,
    0.00211341152631996, 0.00209441885249221, 0.00204948439910114,
    0.00197866066579206, 0.0018822202455726, 0.0017608276313106,
    0.00161537857650242, 0.0014471716527866, 0.0012577031505066,
    0.00104883948423773, 0.000822638973055532, 0.000581505818928723,
    0.000327975991355641, 6.48740710262127E-5, -0.000204878269979111,
    -0.000478120672898087, -0.00075168575015615, -0.00102224662003798,
    -0.00128651260628978, -0.00154110992486032, -0.00178281151318399,
    -0.00200836474362844, -0.00221472015192204, -0.00239892427229273,
    -0.00255826248660913, -0.00269021370473722, -0.00279257189999333,
    -0.00286337500792373, -0.00290105581462172, -0.0029043318032904,
    -0.00287237823233321, -0.00280467899129824, -0.00270121258812355,
    -0.00256228356681253, -0.00238872096189341, -0.00218169844396906,
    -0.00194290595192604, -0.00167436411638681, -0.00137857899476873,
    -0.00105833699746803, -0.000716870711910641, -0.000357662081547121,
    1.54069234904593E-5, 0.000398269941917141, 0.000786554082636692,
    0.0011758230289691, 0.00156143613107137, 0.00193873592622316,
    0.00230289725078119, 0.00264921375954365, 0.00297305133436034,
    0.00327001605698777, 0.00353565573027526, 0.00376581900389305,
    0.00395697189829865, 0.00410539360304297, 0.00420817378241546,
    0.00426264813419011, 0.00426664918670783, 0.00421857702165025,
    0.00411729676153594, 0.0039623370132401, 0.00375373166757278,
    0.00349221769263537, 0.00317906140829903, 0.00281624784751191,
    0.00240629391400286, 0.0019524231096476, 0.00145836996294566,
    0.00092853544645223, 0.000367776971615253, -0.000218437652524706,
    -0.000824257773321898, -0.00144331686525717, -0.00206896141715662,
    -0.00269412620596608, -0.00331158028138154, -0.00391380697172656,
    -0.00449325707533294, -0.00504223586171052, -0.00555315588129342,
    -0.00601843603020415, -0.00643073558805818, -0.00678291423619245,
    -0.00706812188302244, -0.00728002984015002, -0.00741252860038501,
    -0.00746025640048564, -0.00741838763488326, -0.00728259643511458,
    -0.00704943210024094, -0.00671611785605616, -0.00628069922505195,
    -0.00574193764251256, -0.00509956676718777, -0.0043541195370051,
    -0.00350708574003026, -0.00256073033453729, -0.00151829096644007,
    -0.00038379958925631, 0.000837765270517561, 0.00214069365337124,
    0.00351840339682411, 0.0049636681477122, 0.0064684801757253,
    0.00802427702998046, 0.00962183372004456, 0.0112515176906817,
    0.0129031792166922, 0.0145664200115191, 0.0162304954604694,
    0.0178845643813347, 0.0195176178574509, 0.021118715628301,
    0.0226769516122276, 0.0241816643667645, 0.0256224202753077,
    0.0269891976615814, 0.0282723656290463, 0.0294628998851187,
    0.0305522889367809, 0.0315328063885665, 0.0323973751703697,
    0.033139808318974, 0.0337547063219542, 0.0342376726977091, 0.034585156603869,
    0.034794675320555, 0.0348646677825939, 0.034794675320555, 0.034585156603869,
    0.0342376726977091, 0.0337547063219542, 0.033139808318974,
    0.0323973751703697, 0.0315328063885665, 0.0305522889367809,
    0.0294628998851187, 0.0282723656290463, 0.0269891976615814,
    0.0256224202753077, 0.0241816643667645, 0.0226769516122276,
    0.021118715628301, 0.0195176178574509, 0.0178845643813347,
    0.0162304954604694, 0.0145664200115191, 0.0129031792166922,
    0.0112515176906817, 0.00962183372004456, 0.00802427702998046,
    0.0064684801757253, 0.0049636681477122, 0.00351840339682411,
    0.00214069365337124, 0.000837765270517561, -0.00038379958925631,
    -0.00151829096644007, -0.00256073033453729, -0.00350708574003026,
    -0.0043541195370051, -0.00509956676718777, -0.00574193764251256,
    -0.00628069922505195, -0.00671611785605616, -0.00704943210024094,
    -0.00728259643511458, -0.00741838763488326, -0.00746025640048564,
    -0.00741252860038501, -0.00728002984015002, -0.00706812188302244,
    -0.00678291423619245, -0.00643073558805818, -0.00601843603020415,
    -0.00555315588129342, -0.00504223586171052, -0.00449325707533294,
    -0.00391380697172656, -0.00331158028138154, -0.00269412620596608,
    -0.00206896141715662, -0.00144331686525717, -0.000824257773321898,
    -0.000218437652524706, 0.000367776971615253, 0.00092853544645223,
    0.00145836996294566, 0.0019524231096476, 0.00240629391400286,
    0.00281624784751191, 0.00317906140829903, 0.00349221769263537,
    0.00375373166757278, 0.0039623370132401, 0.00411729676153594,
    0.00421857702165025, 0.00426664918670783, 0.00426264813419011,
    0.00420817378241546, 0.00410539360304297, 0.00395697189829865,
    0.00376581900389305, 0.00353565573027526, 0.00327001605698777,
    0.00297305133436034, 0.00264921375954365, 0.00230289725078119,
    0.00193873592622316, 0.00156143613107137, 0.0011758230289691,
    0.000786554082636692, 0.000398269941917141, 1.54069234904593E-5,
    -0.000357662081547121, -0.000716870711910641, -0.00105833699746803,
    -0.00137857899476873, -0.00167436411638681, -0.00194290595192604,
    -0.00218169844396906, -0.00238872096189341, -0.00256228356681253,
    -0.00270121258812355, -0.00280467899129824, -0.00287237823233321,
    -0.0029043318032904, -0.00290105581462172, -0.00286337500792373,
    -0.00279257189999333, -0.00269021370473722, -0.00255826248660913,
    -0.00239892427229273, -0.00221472015192204, -0.00200836474362844,
    -0.00178281151318399, -0.00154110992486032, -0.00128651260628978,
    -0.00102224662003798, -0.00075168575015615, -0.000478120672898087,
    -0.000204878269979111, 6.48740710262127E-5, 0.000327975991355641,
    0.000581505818928723, 0.000822638973055532, 0.00104883948423773,
    0.0012577031505066, 0.0014471716527866, 0.00161537857650242,
    0.0017608276313106, 0.0018822202455726, 0.00197866066579206,
    0.00204948439910114, 0.00209441885249221, 0.00211341152631996,
    0.00210681645356649, 0.00207518637052296, 0.00201944312775646,
    0.00194070288808355, 0.00184038233686247, 0.00172006918497051,
    0.00158161297276232, 0.0014269852599288, 0.0012583731149921,
    0.00107801310643884, 0.000888309185380098, 0.000691646548084861,
    0.000490536127404496, 0.000287402789594611, 8.47444798042395E-5,
    -0.000115084270895235, -0.000309729372712078, -0.00049703143195459,
    -0.000674867067745653, -0.000841350754803934, -0.000994673305073961,
    -0.00113330657899723, -0.00125585035904654, -0.00136121045709383,
    -0.00144841035655361, -0.00151678445534883, -0.00156583223311964,
    -0.00159540022296035, -0.00160546004472087, -0.00159627209149843,
    -0.0015682621132935, -0.00152219868611472, -0.00145889954833536,
    -0.00137943108525202, -0.00128509963403905, -0.00117724186294024,
    -0.00105741642025622, -0.000927273197598017, -0.000788522783484884,
    -0.000643010218370379, -0.000492548025635855, -0.000339059179704988,
    -0.000184379347411724, -3.04152432307783E-5, 0.000121059821395169,
    0.000268243354692785, 0.000409507246153594, 0.000543234859783321,
    0.000668013805650651, 0.000782475751035737, 0.000885492683406536,
    0.000976000891026337, 0.00105319209384983, 0.00111634393692054,
    0.00116500239871139, 0.00119880438152662, 0.00121766042770094,
    0.00122156944464255, 0.00121079692202438, 0.00118568928711056,
    0.00114684740229226, 0.00109493959888634, 0.00103084522155815,
    0.000955579555048026, 0.000870008572946229, 0.000775880149805064,
    0.000673758101284352, 0.000565409273805852, 0.0004523378855645,
    0.00033576160747592, 0.000217112766350226, 9.78643173315946E-5,
    -2.04444307204847E-5, -0.000136464578717182, -0.000248851915276288,
    -0.000356409314647111, -0.000457903260615247, -0.000552260439891737,
    -0.000638420749418111, -0.000715536243565728, -0.00078281274902258,
    -0.000839683905906734, -0.000885626875170628, -0.000920339162179491,
    -0.000943564158373563, -0.000955264249053625, -0.00095545499522514,
    -0.000944383237710939, -0.000922339566245509, -0.000889830045262367,
    -0.000847393382183875, -0.000795757970971713, -0.000735672733369139,
    -0.000668037007475755, -0.000593763701806401, -0.000513872111887188,
    -0.000429399162350645, -0.000341447274286012, -0.000251102810739723,
    -0.000159529303918018, -6.77874598097355E-5, 2.29692786579564E-5,
    0.000111712026596431, 0.0001973870440229, 0.000279073495233852,
    0.00035581891628597, 0.000426844595346898, 0.000491374900037892,
    0.000548801415044521, 0.000598517083474895, 0.000640110647994494,
    0.00067319585265931, 0.000697573224095843, 0.000713058339562528,
    0.000719672891575588, 0.000717462680602415, 0.000706653899194518,
    0.000687485682784789, 0.000660385014627503, 0.000625772979312749,
    0.000584245077633135, 0.000536391099262198, 0.000482925129093528,
    0.000424554418003875, 0.00036208896972677, 0.000296318447014287,
    0.000228118320448514, 0.000158312370733294, 8.7795776565017E-5,
    1.73789449575496E-5, -5.20614307239963E-5, -0.000119764102369112,
    -0.000184906233867251, -0.000246799554346799, -0.000304728818471515,
    -0.000358124919507632, -0.000406397565443762, -0.000449113416335854,
    -0.000485819662611941, -0.000516232780295362, -0.000540076800221215,
    -0.000557245516783168, -0.000567612524697722, -0.000571215822744036,
    -0.000568103825078294, -0.000558500829399108, -0.000542597316576951,
    -0.000520728588172605, -0.000493239085464758, -0.000460629099728361,
    -0.000423315587749875, -0.000381865684983965, -0.00033686907410968,
    -0.000288902070780961, -0.000238612044359073, -0.000186647951479871,
    -0.000133640898306514, -8.02768371814556E-5, -2.7161145425949E-5,
    2.50361530035493E-5, 7.57533418379797E-5, 0.000124370213198492,
    0.000170392162781004, 0.000213277183864515, 0.000252620897756928,
    0.000287981876933402, 0.000319057917826926, 0.00034551939954402,
    0.00036718765573372, 0.000383854155585259, 0.000395460692994347,
    0.000401925594888096, 0.000403315939413201, 0.000399671906379765,
    0.000391181527305812, 0.000377998077818869, 0.000360413706777267,
    0.000338680260437589, 0.000313182758504491, 0.000284257828711435,
    0.000252357113836307, 0.000217931456555117, 0.000181025976607809,
    0.000143704940542049, 0.000103589939605899, 6.35611116325544E-5,
    2.37618534986138E-5, -1.58536847988518E-5, -5.48945053903947E-5,
    -9.28427841098224E-5, -0.000129096464368348, -0.000163234645663197,
    -0.000194860168732485, -0.000223736461708317, -0.000249588877743172,
    -0.000272250469485773, -0.00029148787236693, -0.000307173721793938,
    -0.000319136839763799, -0.000327336417669472, -0.000331700761792967,
    -0.000332299375483446, -0.000329171873533593, -0.00032248522846086,
    -0.000312368664535981, -0.00029907288616221, -0.000282792065964986,
    -0.000263829521034822, -0.000242428006148879, -0.000218930579462653,
    -0.000193622797261251, -0.00016687275536317, -0.000139001076920169,
    -0.000110381212543967, -8.13640577031261E-5, -5.23127671566118E-5,
    -2.35688788073929E-5, 4.49760535647356E-6, 3.16134719413095E-5,
    5.74156605995404E-5, 8.16637642732029E-5, 0.00010407408097609,
    0.000124455756549172, 0.000142560624355696, 0.000158270235612336,
    0.000171419148763144, 0.000181945744112124, 0.000189736391522837,
    0.000194802675707469, 0.000197109088634444, 0.000196723579657171,
    0.000193659121415496, 0.000188051801111294, 0.000179985578337734,
    0.000169645517325974, 0.000157159733488689, 0.000142760685517128,
    0.000126612317432528, 0.000108987167847066, 9.00914671901327E-5,
    7.02071138914894E-5, 4.95590136185014E-5, 2.84329633486768E-5,
    7.06086733384884E-6, -1.42714666447916E-5, -3.53418834782659E-5,
    -5.58768211555122E-5, -7.56833010219895E-5, -9.45098568138851E-5,
    -0.000112203089934051, -0.000128534875844589, -0.000143383749251984,
    -0.000156570952883738, -0.000168020586519177, -0.000177600347515248,
    -0.000185281259751758, -0.000190965804054552, -0.000194665923081674,
    -0.000196346119718796, -0.000196073541936795, -0.000193841772251799,
    -0.000189747215029418, -0.00018383717508416, -0.000176261149179698,
    -0.000167084529934335, -0.000156470017010974, -0.000144534880486384,
    -0.000131481469000169, -0.000117407979761088, -0.000102529979431761,
    -8.70064897013418E-5, -7.10064533753939E-5, -5.4711660316172E-5,
    -3.83000510956379E-5, -2.19231126541713E-5, -5.77254316079002E-6,
    1.0030902985957E-5, 2.53040869264146E-5, 3.99522311415092E-5,
    5.38098574168716E-5, 6.68078245004558E-5, 7.88062325959248E-5,
    8.97652851564876E-5, 9.95750422745516E-5, 0.000108220396625604,
    0.000115621965373439, 0.000121799681582674, 0.000126701511593129,
    0.000130373124118671, 0.0001327911298552, 0.000134025874269021,
    0.000134078026482925, 0.000133044030326802, 0.000130945415051525,
    0.000127894687581549, 0.000123929134733962, 0.0001191770935068,
    0.000113688228537891, 0.000107595463015365, 0.000100958414377355,
    9.39155210670181E-5, 0.000541983718811927 };

  static const double dv6[431] = { 0.00059108930639301, 0.000203135102293136,
    0.000228702342098083, 0.000249343617535587, 0.0002635069762735,
    0.000269934205063887, 0.000267468276220013, 0.000255402479651796,
    0.000233273621757195, 0.000201219812933285, 0.000159682886147079,
    0.000109724050959972, 5.27085192408347E-5, -9.40957449163658E-6,
    -7.44675389209232E-5, -0.000139841547166302, -0.000202885801604255,
    -0.000260742843701583, -0.000310787524703063, -0.000350412679782949,
    -0.000377500244119596, -0.000390236144409474, -0.000387546998560333,
    -0.000368850380778751, -0.000334438591051295, -0.00028518190011347,
    -0.00022286939488848, -0.000149852306955432, -6.92610818611087E-5,
    1.54578194554596E-5, 0.000100328930107147, 0.000181312682430688,
    0.000254228552368588, 0.000315289227791654, 0.000360922754044401,
    0.000388245755760958, 0.000395072723587488, 0.000380257274688715,
    0.000343456936944303, 0.000285548262117333, 0.000208452985533168,
    0.00011504912925973, 9.22773619128056E-6, -0.000104352698013944,
    -0.000220482509693878, -0.000333510603247651, -0.000437774350204889,
    -0.000527652824120159, -0.000598146628842731, -0.000644809129631383,
    -0.000664324261454271, -0.000654418809838929, -0.000614326456707724,
    -0.000544556021427374, -0.000447273028712333, -0.000325932862966167,
    -0.000185480230425671, -3.18423145588313E-5, 0.000127976430894124,
    0.000286536638081843, 0.000435981691500021, 0.00056873011946875,
    0.000677529452718088, 0.000756169921240225, 0.00079949138541286,
    0.00080400917405841, 0.000767835944931364, 0.00069116925972533,
    0.000576024658008797, 0.000426553487272439, 0.000248616157020845,
    4.98649714080354E-5, -0.000160884175849547, -0.000373698214029173,
    -0.000578284175225528, -0.000764238148896663, -0.000921730736103272,
    -0.00104195550586188, -0.00111743946208126, -0.00114289576019665,
    -0.00111495567630999, -0.00103287630883038, -0.000898679009764979,
    -0.00071671588526, -0.000494043479284734, -0.000239947494484822,
    3.43194002594898E-5, 0.000316273875969929, 0.00059235552217791,
    0.000848893105460828, 0.00107254336564001, 0.00125127031009951,
    0.00137460646711804, 0.0014345399251153, 0.00142569620683721,
    0.0013460326251024, 0.00119672309564059, 0.000982528974077999,
    0.00071137182731132, 0.000394431615214206, 4.54340164485031E-5,
    -0.00031961189899161, -0.000683490853806374, -0.00102828762426528,
    -0.00133657646482561, -0.00159194769783486, -0.00178011022049491,
    -0.0018892781771403, -0.00191116859766644, -0.00184130095806268,
    -0.00167978390347217, -0.00143121554544908, -0.00110479063604924,
    -0.000713421120918264, -0.000273458733759465, 0.000195673531830978,
    0.000671366376562233, 0.0011302146092135, 0.00155149140592462,
    0.00191065109022734, 0.00218952949178247, 0.00237100500122057,
    0.00244271393319439, 0.00239715681481954, 0.00223215536250292,
    0.00195142503243853, 0.00156428896227672, 0.00108576770776838,
    0.000535720184246959, -6.16152130609063E-5, -0.000678958280353827,
    -0.00128686237729246, -0.00185542519534156, -0.00235536537880196,
    -0.00275979008980575, -0.00304525170506758, -0.00319337838493626,
    -0.00319160458969291, -0.00303434112247834, -0.00272315976994321,
    -0.00226733239506472, -0.00168333400111216, -0.000994682377620397,
    -0.000230714554582418, 0.000574271356711954, 0.00138280136505085,
    0.00215554933784529, 0.00285353542133235, 0.00343968761494919,
    0.00388099477089373, 0.00414992270721462, 0.00422634082660383,
    0.00409845181704085, 0.00376392649007728, 0.00323021244990824,
    0.00251482837201619, 0.00164466215588651, 0.000655298804813174,
    -0.000410399793453483, -0.00150401349518631, -0.00257368112952168,
    -0.00356627803023731, -0.00443007902205741, -0.0051170846663053,
    -0.00558567080162908, -0.00580256924799627, -0.0057451743710917,
    -0.0054027330838167, -0.00477772903587716, -0.0038860874195769,
    -0.00275741920859663, -0.00143393039765822, 3.05684537051032E-5,
    0.00157291354971703, 0.00312270604093914, 0.00460555918146983,
    0.00594595608074391, 0.0070710055803154, 0.00791357328678103,
    0.00841599967138872, 0.0085328783812285, 0.00823414920527687,
    0.00750696640123341, 0.00635761364944972, 0.00481202042811583,
    0.00291616685206864, 0.000734945152006006, -0.00164911359246092,
    -0.00413854912587785, -0.00662376584037583, -0.00898713211466898,
    -0.0111068908850509, -0.0128620309563478, -0.0141368054464089,
    -0.0148254937011961, -0.0148372599010804, -0.0140992200430354,
    -0.012561790818122, -0.0101995897012186, -0.0070144801308704,
    -0.0030368526015946, 0.00167519934115319, 0.00703647678686433,
    0.0129367099188529, 0.0192433657583256, 0.0258060198556134,
    0.0324606280440578, 0.0390349134111635, 0.0453535974623025,
    0.0512445768528114, 0.0565444150513016, 0.0611042135615666,
    0.0647943156695586, 0.0675090341357033, 0.0691699572095659,
    0.0697290376973575, 0.0691699572095659, 0.0675090341357033,
    0.0647943156695586, 0.0611042135615666, 0.0565444150513016,
    0.0512445768528114, 0.0453535974623025, 0.0390349134111635,
    0.0324606280440578, 0.0258060198556134, 0.0192433657583256,
    0.0129367099188529, 0.00703647678686433, 0.00167519934115319,
    -0.0030368526015946, -0.0070144801308704, -0.0101995897012186,
    -0.012561790818122, -0.0140992200430354, -0.0148372599010804,
    -0.0148254937011961, -0.0141368054464089, -0.0128620309563478,
    -0.0111068908850509, -0.00898713211466898, -0.00662376584037583,
    -0.00413854912587785, -0.00164911359246092, 0.000734945152006006,
    0.00291616685206864, 0.00481202042811583, 0.00635761364944972,
    0.00750696640123341, 0.00823414920527687, 0.0085328783812285,
    0.00841599967138872, 0.00791357328678103, 0.0070710055803154,
    0.00594595608074391, 0.00460555918146983, 0.00312270604093914,
    0.00157291354971703, 3.05684537051032E-5, -0.00143393039765822,
    -0.00275741920859663, -0.0038860874195769, -0.00477772903587716,
    -0.0054027330838167, -0.0057451743710917, -0.00580256924799627,
    -0.00558567080162908, -0.0051170846663053, -0.00443007902205741,
    -0.00356627803023731, -0.00257368112952168, -0.00150401349518631,
    -0.000410399793453483, 0.000655298804813174, 0.00164466215588651,
    0.00251482837201619, 0.00323021244990824, 0.00376392649007728,
    0.00409845181704085, 0.00422634082660383, 0.00414992270721462,
    0.00388099477089373, 0.00343968761494919, 0.00285353542133235,
    0.00215554933784529, 0.00138280136505085, 0.000574271356711954,
    -0.000230714554582418, -0.000994682377620397, -0.00168333400111216,
    -0.00226733239506472, -0.00272315976994321, -0.00303434112247834,
    -0.00319160458969291, -0.00319337838493626, -0.00304525170506758,
    -0.00275979008980575, -0.00235536537880196, -0.00185542519534156,
    -0.00128686237729246, -0.000678958280353827, -6.16152130609063E-5,
    0.000535720184246959, 0.00108576770776838, 0.00156428896227672,
    0.00195142503243853, 0.00223215536250292, 0.00239715681481954,
    0.00244271393319439, 0.00237100500122057, 0.00218952949178247,
    0.00191065109022734, 0.00155149140592462, 0.0011302146092135,
    0.000671366376562233, 0.000195673531830978, -0.000273458733759465,
    -0.000713421120918264, -0.00110479063604924, -0.00143121554544908,
    -0.00167978390347217, -0.00184130095806268, -0.00191116859766644,
    -0.0018892781771403, -0.00178011022049491, -0.00159194769783486,
    -0.00133657646482561, -0.00102828762426528, -0.000683490853806374,
    -0.00031961189899161, 4.54340164485031E-5, 0.000394431615214206,
    0.00071137182731132, 0.000982528974077999, 0.00119672309564059,
    0.0013460326251024, 0.00142569620683721, 0.0014345399251153,
    0.00137460646711804, 0.00125127031009951, 0.00107254336564001,
    0.000848893105460828, 0.00059235552217791, 0.000316273875969929,
    3.43194002594898E-5, -0.000239947494484822, -0.000494043479284734,
    -0.00071671588526, -0.000898679009764979, -0.00103287630883038,
    -0.00111495567630999, -0.00114289576019665, -0.00111743946208126,
    -0.00104195550586188, -0.000921730736103272, -0.000764238148896663,
    -0.000578284175225528, -0.000373698214029173, -0.000160884175849547,
    4.98649714080354E-5, 0.000248616157020845, 0.000426553487272439,
    0.000576024658008797, 0.00069116925972533, 0.000767835944931364,
    0.00080400917405841, 0.00079949138541286, 0.000756169921240225,
    0.000677529452718088, 0.00056873011946875, 0.000435981691500021,
    0.000286536638081843, 0.000127976430894124, -3.18423145588313E-5,
    -0.000185480230425671, -0.000325932862966167, -0.000447273028712333,
    -0.000544556021427374, -0.000614326456707724, -0.000654418809838929,
    -0.000664324261454271, -0.000644809129631383, -0.000598146628842731,
    -0.000527652824120159, -0.000437774350204889, -0.000333510603247651,
    -0.000220482509693878, -0.000104352698013944, 9.22773619128056E-6,
    0.00011504912925973, 0.000208452985533168, 0.000285548262117333,
    0.000343456936944303, 0.000380257274688715, 0.000395072723587488,
    0.000388245755760958, 0.000360922754044401, 0.000315289227791654,
    0.000254228552368588, 0.000181312682430688, 0.000100328930107147,
    1.54578194554596E-5, -6.92610818611087E-5, -0.000149852306955432,
    -0.00022286939488848, -0.00028518190011347, -0.000334438591051295,
    -0.000368850380778751, -0.000387546998560333, -0.000390236144409474,
    -0.000377500244119596, -0.000350412679782949, -0.000310787524703063,
    -0.000260742843701583, -0.000202885801604255, -0.000139841547166302,
    -7.44675389209232E-5, -9.40957449163658E-6, 5.27085192408347E-5,
    0.000109724050959972, 0.000159682886147079, 0.000201219812933285,
    0.000233273621757195, 0.000255402479651796, 0.000267468276220013,
    0.000269934205063887, 0.0002635069762735, 0.000249343617535587,
    0.000228702342098083, 0.000203135102293136, 0.00059108930639301 };

  static const double dv7[216] = { 0.000703874537800322, 0.000467770792417846,
    0.000539242180785141, 0.000548714935447943, 0.000481522011841131,
    0.000335028055626118, 0.000121223276456724, -0.000133359726248962,
    -0.000390715672681121, -0.000607205143402996, -0.000741333556204524,
    -0.000762141773286633, -0.000656572560307136, -0.000434098146732689,
    -0.000127505132900021, 0.000211213864529279, 0.000518476440093037,
    0.000731495160441613, 0.000799153119402637, 0.000695246620796837,
    0.00042431661526213, 2.47090795265547E-5, -0.000435793659924257,
    -0.000871359045494241, -0.00119303737495923, -0.0013260633211882,
    -0.00122651508944986, -0.000892783882200666, -0.000369395944141702,
    0.000257559799410267, 0.000873723084736679, 0.00135682421064802,
    0.00160057389483701, 0.00153681996160046, 0.00115260727769925,
    0.000497054237143859, -0.000322657933790586, -0.00115812613514566,
    -0.00184572049236201, -0.00223771577754856, -0.00223288551008581,
    -0.00180018865702504, -0.000990574744246535, 6.65230524788864E-5,
    0.00118273010625971, 0.00214338417588515, 0.00274778393008912,
    0.00284990682126778, 0.00239141435582874, 0.00141991192268376,
    8.7532957294945E-5, -0.00137083732615, -0.00267793941285247,
    -0.00356533814117648, -0.00382699219718207, -0.00336439837467214,
    -0.00221403160456867, -0.00055020650992204, 0.00133894443554213,
    0.00309954977047986, 0.00437648940085857, 0.00488292289402572,
    0.00446165445896134, 0.00312562937150274, 0.00106802814472366,
    -0.00136197880230862, -0.00371559282089712, -0.00552480925577379,
    -0.00639220079314186, -0.00607420576815403, -0.00454007720411235,
    -0.0019945287417268, 0.0011437349696741, 0.00430664656510399,
    0.00687535360020624, 0.00829614461523189, 0.00819326036702013,
    0.00645661505509234, 0.00328519219444469, -0.000825238131840536,
    -0.0051519826191329, -0.00886492643018391, -0.0111763865229752,
    -0.0114955598573719, -0.00956066216535001, -0.00551970592973219,
    5.66013232584069E-5, 0.00624100851967568, 0.0118874959665648,
    0.0158228554581239, 0.0170617328212118, 0.0150099407000975,
    0.0096196946673462, 0.00146523736031698, -0.00828158313067491,
    -0.0179786968695017, -0.025728978155432, -0.029655972719434,
    -0.0282029907072548, -0.0204040677918449, -0.00607822692183737,
    0.0140682470909551, 0.0384820357992972, 0.0649166378856776,
    0.0907022576471758, 0.113083794771166, 0.129583642532735, 0.138334705982925,
    0.138334705982925, 0.129583642532735, 0.113083794771166, 0.0907022576471758,
    0.0649166378856776, 0.0384820357992972, 0.0140682470909551,
    -0.00607822692183737, -0.0204040677918449, -0.0282029907072548,
    -0.029655972719434, -0.025728978155432, -0.0179786968695017,
    -0.00828158313067491, 0.00146523736031698, 0.0096196946673462,
    0.0150099407000975, 0.0170617328212118, 0.0158228554581239,
    0.0118874959665648, 0.00624100851967568, 5.66013232584069E-5,
    -0.00551970592973219, -0.00956066216535001, -0.0114955598573719,
    -0.0111763865229752, -0.00886492643018391, -0.0051519826191329,
    -0.000825238131840536, 0.00328519219444469, 0.00645661505509234,
    0.00819326036702013, 0.00829614461523189, 0.00687535360020624,
    0.00430664656510399, 0.0011437349696741, -0.0019945287417268,
    -0.00454007720411235, -0.00607420576815403, -0.00639220079314186,
    -0.00552480925577379, -0.00371559282089712, -0.00136197880230862,
    0.00106802814472366, 0.00312562937150274, 0.00446165445896134,
    0.00488292289402572, 0.00437648940085857, 0.00309954977047986,
    0.00133894443554213, -0.00055020650992204, -0.00221403160456867,
    -0.00336439837467214, -0.00382699219718207, -0.00356533814117648,
    -0.00267793941285247, -0.00137083732615, 8.7532957294945E-5,
    0.00141991192268376, 0.00239141435582874, 0.00284990682126778,
    0.00274778393008912, 0.00214338417588515, 0.00118273010625971,
    6.65230524788864E-5, -0.000990574744246535, -0.00180018865702504,
    -0.00223288551008581, -0.00223771577754856, -0.00184572049236201,
    -0.00115812613514566, -0.000322657933790586, 0.000497054237143859,
    0.00115260727769925, 0.00153681996160046, 0.00160057389483701,
    0.00135682421064802, 0.000873723084736679, 0.000257559799410267,
    -0.000369395944141702, -0.000892783882200666, -0.00122651508944986,
    -0.0013260633211882, -0.00119303737495923, -0.000871359045494241,
    -0.000435793659924257, 2.47090795265547E-5, 0.00042431661526213,
    0.000695246620796837, 0.000799153119402637, 0.000731495160441613,
    0.000518476440093037, 0.000211213864529279, -0.000127505132900021,
    -0.000434098146732689, -0.000656572560307136, -0.000762141773286633,
    -0.000741333556204524, -0.000607205143402996, -0.000390715672681121,
    -0.000133359726248962, 0.000121223276456724, 0.000335028055626118,
    0.000481522011841131, 0.000548714935447943, 0.000539242180785141,
    0.000467770792417846, 0.000703874537800322 };

  static const double dv8[109] = { 0.000701399010505143, 0.00140839928277612,
    0.00111499214990723, 0.000816359876385092, -0.000344933864345112,
    -0.00108141865642759, -0.00124645997896249, -0.000328553861170572,
    0.000953439235176911, 0.00179667825296385, 0.00135345312991716,
    -0.000219458433219591, -0.00194482183672888, -0.0024630397849623,
    -0.00115724040777242, 0.00129882662458625, 0.0031785526935836,
    0.00291164222037781, 0.000326935175146764, -0.00295096520726589,
    -0.00448154193949628, -0.00277888085477348, 0.00136069032614137,
    0.00510846062072828, 0.00548563343107216, 0.00166478513201453,
    -0.00404551341994001, -0.00750445686013723, -0.00567597805715282,
    0.000839555412816575, 0.00771068709028747, 0.00965911058655907,
    0.00440494762676344, -0.00510107364052596, -0.0121479368907136,
    -0.0108638160387891, -0.000888239685836671, 0.0114580689600834,
    0.016954062785566, 0.0101044811810528, -0.00598125671977792,
    -0.0204435882733855, -0.0215909420459943, -0.0057131514725697,
    0.0184513279393411, 0.0336982429435121, 0.0254652063742116,
    -0.00656621578009933, -0.0444036515025244, -0.059329538916181,
    -0.0280398177296805, 0.0517667833144886, 0.156163072823834,
    0.244442499000399, 0.278942485275155, 0.244442499000399, 0.156163072823834,
    0.0517667833144886, -0.0280398177296805, -0.059329538916181,
    -0.0444036515025244, -0.00656621578009933, 0.0254652063742116,
    0.0336982429435121, 0.0184513279393411, -0.0057131514725697,
    -0.0215909420459943, -0.0204435882733855, -0.00598125671977792,
    0.0101044811810528, 0.016954062785566, 0.0114580689600834,
    -0.000888239685836671, -0.0108638160387891, -0.0121479368907136,
    -0.00510107364052596, 0.00440494762676344, 0.00965911058655907,
    0.00771068709028747, 0.000839555412816575, -0.00567597805715282,
    -0.00750445686013723, -0.00404551341994001, 0.00166478513201453,
    0.00548563343107216, 0.00510846062072828, 0.00136069032614137,
    -0.00277888085477348, -0.00448154193949628, -0.00295096520726589,
    0.000326935175146764, 0.00291164222037781, 0.0031785526935836,
    0.00129882662458625, -0.00115724040777242, -0.0024630397849623,
    -0.00194482183672888, -0.000219458433219591, 0.00135345312991716,
    0.00179667825296385, 0.000953439235176911, -0.000328553861170572,
    -0.00124645997896249, -0.00108141865642759, -0.000344933864345112,
    0.000816359876385092, 0.00111499214990723, 0.00140839928277612,
    0.000701399010505143 };

  static const double dv9[55] = { 0.00136160164139038, 0.00267369398929317,
    -0.000697117089159221, -0.00234837763528271, 0.00191581346362666,
    0.00269900377250914, -0.00389621556543998, -0.00225279038429998,
    0.00633783534093455, 0.000605005371784838, -0.00893041710195627,
    0.00276678504202593, 0.0109185935758989, -0.00812278495383672,
    -0.0112857210024778, 0.0154360766314228, 0.00873812519907081,
    -0.0242867022296137, -0.0017025900826685, 0.0338768389503649,
    -0.0120304430378272, -0.0431322973914692, 0.0369550170539418,
    0.0508632874364068, -0.0888408228750916, -0.0560010731682915,
    0.312337552871639, 0.557802460840468, 0.312337552871639, -0.0560010731682915,
    -0.0888408228750916, 0.0508632874364068, 0.0369550170539418,
    -0.0431322973914692, -0.0120304430378272, 0.0338768389503649,
    -0.0017025900826685, -0.0242867022296137, 0.00873812519907081,
    0.0154360766314228, -0.0112857210024778, -0.00812278495383672,
    0.0109185935758989, 0.00276678504202593, -0.00893041710195627,
    0.000605005371784838, 0.00633783534093455, -0.00225279038429998,
    -0.00389621556543998, 0.00269900377250914, 0.00191581346362666,
    -0.00234837763528271, -0.000697117089159221, 0.00267369398929317,
    0.00136160164139038 };

  static const double dv10[34] = { 0.00396227519738595, -0.00432580986830259,
    0.00537645334783375, -0.00525174059537846, 0.00329357485741764,
    0.000923114411272612, -0.00738284264798429, 0.0154672131632766,
    -0.0238706513909533, 0.0306418755845219, -0.0333198455184126,
    0.0290689001544471, -0.0146804244699826, -0.0141610568835452,
    0.0662458874783531, -0.172852103294901, 0.623027351924843, 0.623027351924843,
    -0.172852103294901, 0.0662458874783531, -0.0141610568835452,
    -0.0146804244699826, 0.0290689001544471, -0.0333198455184126,
    0.0306418755845219, -0.0238706513909533, 0.0154672131632766,
    -0.00738284264798429, 0.000923114411272612, 0.00329357485741764,
    -0.00525174059537846, 0.00537645334783375, -0.00432580986830259,
    0.00396227519738595 };

  /*  inputs: */
  /*  fpass - the pass frequency */
  /*  outputs: */
  /*  num - the numerator of the filter */
  /*  den - the numerator and denominator of the filter */
  /*  set the discrete frequences */
  /*  freq = [0.25, 0.50, 1.00, 2.00, 4.00, 8.00, 16.00, 32.00, 64.00, 100.00]; */
  /* 'dynIdenf:551' fdisc = [0.25, 0.35, 0.70, 1.41, 2.82, 5.65, 11.31, 22.62, 45.25, 80.00, 100.00]; */
  /*  get the coefficents */
  /* 'dynIdenf:553' if fpass < fdisc(1) */
  if (fpass < 0.25) {
    /* 'dynIdenf:554' num = [0 0]; */
    i6 = b_num->size[0] * b_num->size[1];
    b_num->size[0] = 1;
    b_num->size[1] = 2;
    emxEnsureCapacity_real_T(b_num, i6);
    for (i6 = 0; i6 < 2; i6++) {
      b_num->data[b_num->size[0] * i6] = 0.0;
    }

    /* 'dynIdenf:555' den = 1; */
  } else if ((0.25 <= fpass) && (fpass < 0.35)) {
    /* 'dynIdenf:556' elseif fdisc(1) <= fpass && fpass < fdisc(2) */
    /* 'dynIdenf:557' num = [-0.000520978002302694,-5.20169990936757e-05,-5.45841105354002e-05,-5.72195149660382e-05,-5.99216529332766e-05,-6.26902993572831e-05,-6.55277572317460e-05,-6.84327924698304e-05,-7.14024783352774e-05,-7.44382315108702e-05,-7.75419955514318e-05,-8.07125864750608e-05,-8.39492202702880e-05,-8.72527281450113e-05,-9.06217587629433e-05,-9.40534210975644e-05,-9.75483334878834e-05,-0.000101109902849751,-0.000104738576954986,-0.000108431175299353,-0.000112185059058074,-0.000115999352356773,-0.000119874199376914,-0.000123811981761444,-0.000127815442979604,-0.000131881513112167,-0.000136002977132762,-0.000140179858868300,-0.000144420119293729,-0.000148723321214873,-0.000153076275831861,-0.000157476132364005,-0.000161941418859757,-0.000166463101250286,-0.000171015020968334,-0.000175636719599759,-0.000180296194896843,-0.000185006087018483,-0.000189757036251856,-0.000194554159433506,-0.000199389160913883,-0.000204258823938341,-0.000209166018405914,-0.000214108569106062,-0.000219080375140543,-0.000224079825588825,-0.000229108242018683,-0.000234163373811686,-0.000239241113430869,-0.000244340284774501,-0.000249460464775636,-0.000254598436236085,-0.000259750812540093,-0.000264916528117114,-0.000270093766899914,-0.000275278552826386,-0.000280468430692387,-0.000285662688010565,-0.000290857606001340,-0.000296047191551294,-0.000301229159464750,-0.000306403210484201,-0.000311563497335765,-0.000316701436270235,-0.000321814827623077,-0.000326902443782036,-0.000331953435387113,-0.000336972143859848,-0.000341937352236693,-0.000346885828652329,-0.000351751509861389,-0.000356584499208047,-0.000361375479940580,-0.000366096777583531,-0.000370757272781832,-0.000375368791800677,-0.000379922791039298,-0.000384404674196856,-0.000388811574995772,-0.000393146707240980,-0.000397406523568648,-0.000401580477775939,-0.000405657064973403,-0.000409625775268277,-0.000413478116422792,-0.000417209639330773,-0.000420817806370440,-0.000424298109597335,-0.000427646677871602,-0.000430866002181443,-0.000433963703733845,-0.000436946687593806,-0.000439819221381830,-0.000442581391775045,-0.000445221558359002,-0.000447710843745925,-0.000450008300141517,-0.000452070267800783,-0.000453863926632936,-0.000455404190153653,-0.000456810882109912,-0.000458350037646638,-0.000460310821393207,-0.000460348968389753,-0.000461459608392154,-0.000462016902456549,-0.000462378433406761,-0.000462524399925029,-0.000462451346355499,-0.000462151708972203,-0.000461618145097995,-0.000460847998429803,-0.000459837145402885,-0.000458578890725335,-0.000457068727987613,-0.000455302981268372,-0.000453274535611701,-0.000450975934296376,-0.000448404400488288,-0.000445558845249149,-0.000442434459743990,-0.000439023959439499,-0.000435321003830649,-0.000431319729545988,-0.000427014842399289,-0.000422403822988828,-0.000417484095055867,-0.000412247617334492,-0.000406685229515458,-0.000400796010207350,-0.000394582525310948,-0.000388037089025544,-0.000381146478738261,-0.000373911374121919,-0.000366340474934905,-0.000358414692150382,-0.000350122619270888,-0.000341486217704825,-0.000332478179656225,-0.000323106326908167,-0.000313362143318375,-0.000303249098328250,-0.000292754542035043,-0.000281877124365330,-0.000270618581875681,-0.000258972238338168,-0.000246930199298563,-0.000234491181193772,-0.000221654608059194,-0.000208415138473624,-0.000194767852636238,-0.000180712190136185,-0.000166247400760134,-0.000151369883516121,-0.000136077105983035,-0.000120368829513933,-0.000104243437580730,-8.76979772704483e-05,-7.07315042938449e-05,-5.33432778794504e-05,-3.55290954195747e-05,-1.72848723039368e-05,1.38872191841378e-06,2.04911018494417e-05,4.00271837304823e-05,6.00012002257589e-05,8.04110281968471e-05,0.000101256579820492,0.000122543061367604,0.000144261922674221,0.000166430525578065,0.000188999271198846,0.000212060291713194,0.000235509463682494,0.000259391814847139,0.000283731681095230,0.000308500057039596,0.000333684825244486,0.000359302512586066,0.000385367215620059,0.000411873717316901,0.000438810480214042,0.000466172926860413,0.000493961235466057,0.000522174877221049,0.000550811482289791,0.000579865561840070,0.000609328804955685,0.000639196293728137,0.000669471252782896,0.000700160521508816,0.000731267776747155,0.000762792692118473,0.000794730701659777,0.000827069207154016,0.000859788963085936,0.000892874456342517,0.000926322448581773,0.000960143220924291,0.000994359339319614,0.00102899356216643,0.00106403359253597,0.00109940448183309,0.00113501276226320,0.00117092633398118,0.00120756603038993,0.00124414286807791,0.00128124465925545,0.00131863665906098,0.00135635578729621,0.00139439095544912,0.00143273921326114,0.00147139468118457,0.00151034727681911,0.00154959006081450,0.00158911709434638,0.00162891923891444,0.00166898846896255,0.00170932107190954,0.00174991207709709,0.00179075141906086,0.00183182857807123,0.00187313619701259,0.00191466727083169,0.00195641347655580,0.00199836703727341,0.00204051967345753,0.00208286050291843,0.00212538001142598,0.00216807371366821,0.00221093576143788,0.00225395281482843,0.00229711200656131,0.00234040938008887,0.00238383968702757,0.00242738493017718,0.00247102931641753,0.00251477773243179,0.00255861614240247,0.00260252203198765,0.00264650059615073,0.00269053584052438,0.00273461819792636,0.00277873398090667,0.00282288392309461,0.00286705106750894,0.00291122100593392,0.00295538889026983,0.00299954644597872,0.00304367821272944,0.00308777148929210,0.00313181876280088,0.00317580991187943,0.00321973079913660,0.00326356957138265,0.00330731742211808,0.00335096383756065,0.00339449670702713,0.00343790560431774,0.00348118036145555,0.00352430894741731,0.00356728063999266,0.00361008749181655,0.00365271900135355,0.00369516019836290,0.00373739832682390,0.00377942516996933,0.00382122980465069,0.00386279739394695,0.00390411765143349,0.00394518232614754,0.00398597405354187,0.00402649854545029,0.00406669350872107,0.00410667702228297,0.00414625801694625,0.00418555412451145,0.00422455330130570,0.00426318977685051,0.00430144859842870,0.00433934592051402,0.00437688313874072,0.00441403882587412,0.00445079066528110,0.00448712803027301,0.00452304763244816,0.00455854755687630,0.00459362319520842,0.00462826233752556,0.00466244580837228,0.00469615596427826,0.00472938226584746,0.00476211864665913,0.00479436106160584,0.00482610822465710,0.00485735774695412,0.00488809853788567,0.00491831081285634,0.00494797438450371,0.00497707538941737,0.00500561055815462,0.00503359021864047,0.00506102608993381,0.00508790129864639,0.00511416372355433,0.00513978126324414,0.00516482484920594,0.00518937981270280,0.00521313023739376,0.00523638224947490,0.00525895958423336,0.00528090765629976,0.00530220995844425,0.00532286619204487,0.00534287022373767,0.00536221380126212,0.00538089270680401,0.00539890107887976,0.00541622949372172,0.00543287183790025,0.00544882648850096,0.00546408964030538,0.00547865411988842,0.00549251489038215,0.00550566997871806,0.00551811665232018,0.00552985145136765,0.00554087212638495,0.00555117526124427,0.00556075495828440,0.00556960768651819,0.00557773394998625,0.00558513146790002,0.00559179314635169,0.00559771613805072,0.00560290445917300,0.00560735723206231,0.00561106560533546,0.00561402913964094,0.00561625888942920,0.00561774501983524,0.00561848476850506,0.00561848476850506,0.00561774501983524,0.00561625888942920,0.00561402913964094,0.00561106560533546,0.00560735723206231,0.00560290445917300,0.00559771613805072,0.00559179314635169,0.00558513146790002,0.00557773394998625,0.00556960768651819,0.00556075495828440,0.00555117526124427,0.00554087212638495,0.00552985145136765,0.00551811665232018,0.00550566997871806,0.00549251489038215,0.00547865411988842,0.00546408964030538,0.00544882648850096,0.00543287183790025,0.00541622949372172,0.00539890107887976,0.00538089270680401,0.00536221380126212,0.00534287022373767,0.00532286619204487,0.00530220995844425,0.00528090765629976,0.00525895958423336,0.00523638224947490,0.00521313023739376,0.00518937981270280,0.00516482484920594,0.00513978126324414,0.00511416372355433,0.00508790129864639,0.00506102608993381,0.00503359021864047,0.00500561055815462,0.00497707538941737,0.00494797438450371,0.00491831081285634,0.00488809853788567,0.00485735774695412,0.00482610822465710,0.00479436106160584,0.00476211864665913,0.00472938226584746,0.00469615596427826,0.00466244580837228,0.00462826233752556,0.00459362319520842,0.00455854755687630,0.00452304763244816,0.00448712803027301,0.00445079066528110,0.00441403882587412,0.00437688313874072,0.00433934592051402,0.00430144859842870,0.00426318977685051,0.00422455330130570,0.00418555412451145,0.00414625801694625,0.00410667702228297,0.00406669350872107,0.00402649854545029,0.00398597405354187,0.00394518232614754,0.00390411765143349,0.00386279739394695,0.00382122980465069,0.00377942516996933,0.00373739832682390,0.00369516019836290,0.00365271900135355,0.00361008749181655,0.00356728063999266,0.00352430894741731,0.00348118036145555,0.00343790560431774,0.00339449670702713,0.00335096383756065,0.00330731742211808,0.00326356957138265,0.00321973079913660,0.00317580991187943,0.00313181876280088,0.00308777148929210,0.00304367821272944,0.00299954644597872,0.00295538889026983,0.00291122100593392,0.00286705106750894,0.00282288392309461,0.00277873398090667,0.00273461819792636,0.00269053584052438,0.00264650059615073,0.00260252203198765,0.00255861614240247,0.00251477773243179,0.00247102931641753,0.00242738493017718,0.00238383968702757,0.00234040938008887,0.00229711200656131,0.00225395281482843,0.00221093576143788,0.00216807371366821,0.00212538001142598,0.00208286050291843,0.00204051967345753,0.00199836703727341,0.00195641347655580,0.00191466727083169,0.00187313619701259,0.00183182857807123,0.00179075141906086,0.00174991207709709,0.00170932107190954,0.00166898846896255,0.00162891923891444,0.00158911709434638,0.00154959006081450,0.00151034727681911,0.00147139468118457,0.00143273921326114,0.00139439095544912,0.00135635578729621,0.00131863665906098,0.00128124465925545,0.00124414286807791,0.00120756603038993,0.00117092633398118,0.00113501276226320,0.00109940448183309,0.00106403359253597,0.00102899356216643,0.000994359339319614,0.000960143220924291,0.000926322448581773,0.000892874456342517,0.000859788963085936,0.000827069207154016,0.000794730701659777,0.000762792692118473,0.000731267776747155,0.000700160521508816,0.000669471252782896,0.000639196293728137,0.000609328804955685,0.000579865561840070,0.000550811482289791,0.000522174877221049,0.000493961235466057,0.000466172926860413,0.000438810480214042,0.000411873717316901,0.000385367215620059,0.000359302512586066,0.000333684825244486,0.000308500057039596,0.000283731681095230,0.000259391814847139,0.000235509463682494,0.000212060291713194,0.000188999271198846,0.000166430525578065,0.000144261922674221,0.000122543061367604,0.000101256579820492,8.04110281968471e-05,6.00012002257589e-05,4.00271837304823e-05,2.04911018494417e-05,1.38872191841378e-06,-1.72848723039368e-05,-3.55290954195747e-05,-5.33432778794504e-05,-7.07315042938449e-05,-8.76979772704483e-05,-0.000104243437580730,-0.000120368829513933,-0.000136077105983035,-0.000151369883516121,-0.000166247400760134,-0.000180712190136185,-0.000194767852636238,-0.000208415138473624,-0.000221654608059194,-0.000234491181193772,-0.000246930199298563,-0.000258972238338168,-0.000270618581875681,-0.000281877124365330,-0.000292754542035043,-0.000303249098328250,-0.000313362143318375,-0.000323106326908167,-0.000332478179656225,-0.000341486217704825,-0.000350122619270888,-0.000358414692150382,-0.000366340474934905,-0.000373911374121919,-0.000381146478738261,-0.000388037089025544,-0.000394582525310948,-0.000400796010207350,-0.000406685229515458,-0.000412247617334492,-0.000417484095055867,-0.000422403822988828,-0.000427014842399289,-0.000431319729545988,-0.000435321003830649,-0.000439023959439499,-0.000442434459743990,-0.000445558845249149,-0.000448404400488288,-0.000450975934296376,-0.000453274535611701,-0.000455302981268372,-0.000457068727987613,-0.000458578890725335,-0.000459837145402885,-0.000460847998429803,-0.000461618145097995,-0.000462151708972203,-0.000462451346355499,-0.000462524399925029,-0.000462378433406761,-0.000462016902456549,-0.000461459608392154,-0.000460348968389753,-0.000460310821393207,-0.000458350037646638,-0.000456810882109912,-0.000455404190153653,-0.000453863926632936,-0.000452070267800783,-0.000450008300141517,-0.000447710843745925,-0.000445221558359002,-0.000442581391775045,-0.000439819221381830,-0.000436946687593806,-0.000433963703733845,-0.000430866002181443,-0.000427646677871602,-0.000424298109597335,-0.000420817806370440,-0.000417209639330773,-0.000413478116422792,-0.000409625775268277,-0.000405657064973403,-0.000401580477775939,-0.000397406523568648,-0.000393146707240980,-0.000388811574995772,-0.000384404674196856,-0.000379922791039298,-0.000375368791800677,-0.000370757272781832,-0.000366096777583531,-0.000361375479940580,-0.000356584499208047,-0.000351751509861389,-0.000346885828652329,-0.000341937352236693,-0.000336972143859848,-0.000331953435387113,-0.000326902443782036,-0.000321814827623077,-0.000316701436270235,-0.000311563497335765,-0.000306403210484201,-0.000301229159464750,-0.000296047191551294,-0.000290857606001340,-0.000285662688010565,-0.000280468430692387,-0.000275278552826386,-0.000270093766899914,-0.000264916528117114,-0.000259750812540093,-0.000254598436236085,-0.000249460464775636,-0.000244340284774501,-0.000239241113430869,-0.000234163373811686,-0.000229108242018683,-0.000224079825588825,-0.000219080375140543,-0.000214108569106062,-0.000209166018405914,-0.000204258823938341,-0.000199389160913883,-0.000194554159433506,-0.000189757036251856,-0.000185006087018483,-0.000180296194896843,-0.000175636719599759,-0.000171015020968334,-0.000166463101250286,-0.000161941418859757,-0.000157476132364005,-0.000153076275831861,-0.000148723321214873,-0.000144420119293729,-0.000140179858868300,-0.000136002977132762,-0.000131881513112167,-0.000127815442979604,-0.000123811981761444,-0.000119874199376914,-0.000115999352356773,-0.000112185059058074,-0.000108431175299353,-0.000104738576954986,-0.000101109902849751,-9.75483334878834e-05,-9.40534210975644e-05,-9.06217587629433e-05,-8.72527281450113e-05,-8.39492202702880e-05,-8.07125864750608e-05,-7.75419955514318e-05,-7.44382315108702e-05,-7.14024783352774e-05,-6.84327924698304e-05,-6.55277572317460e-05,-6.26902993572831e-05,-5.99216529332766e-05,-5.72195149660382e-05,-5.45841105354002e-05,-5.20169990936757e-05,-0.000520978002302694]; */
    i6 = b_num->size[0] * b_num->size[1];
    b_num->size[0] = 1;
    b_num->size[1] = 684;
    emxEnsureCapacity_real_T(b_num, i6);
    for (i6 = 0; i6 < 684; i6++) {
      b_num->data[b_num->size[0] * i6] = dv1[i6];
    }

    /* 'dynIdenf:558' den = 1; */
    /*  transient width 5 */
  } else if ((0.35 <= fpass) && (fpass < 0.7)) {
    /* 'dynIdenf:559' elseif fdisc(2) <= fpass && fpass < fdisc(3) */
    /* 'dynIdenf:560' num = [0.000524469736335486,5.00586518032772e-05,5.23815850259164e-05,5.47674639762573e-05,5.71675693367720e-05,5.96287108106926e-05,6.21002944473329e-05,6.46252421778573e-05,6.71538849242156e-05,6.97343785522134e-05,7.23168342051082e-05,7.49510402638714e-05,7.75832509210173e-05,8.02551921214269e-05,8.29159682127781e-05,8.56184309621784e-05,8.83159101463416e-05,9.10537173189225e-05,9.37690208083129e-05,9.65017643446230e-05,9.92159696337733e-05,0.000101967709425591,0.000104680883398715,0.000107392581745282,0.000110108789624883,0.000112800194579151,0.000115468005236082,0.000118132310510059,0.000120750142158453,0.000123362723034365,0.000125918478995377,0.000128458153510528,0.000130933154839445,0.000133391093408533,0.000135776345148193,0.000138133925206034,0.000140409372311920,0.000142648282511158,0.000144799075700251,0.000146905389952095,0.000148909574313128,0.000150860785731569,0.000152708088942437,0.000154490393439512,0.000156153810147281,0.000157751952845765,0.000159227356866869,0.000160617329536024,0.000161881818578263,0.000163049900688652,0.000164080891749098,0.000165012261291061,0.000165793466249417,0.000166463912498038,0.000166977599482095,0.000167371723214612,0.000167596209337449,0.000167695463616470,0.000167612685064898,0.000167395706802922,0.000166989668474572,0.000166435386587554,0.000165681305169779,0.000164776375709852,0.000163662215298123,0.000162382142477805,0.000160885927205586,0.000159218674345988,0.000157324536310168,0.000155251833704162,0.000152947120543503,0.000150454334301420,0.000147727128052281,0.000144804508496339,0.000141639738605470,0.000138276287207714,0.000134668974387507,0.000130850035109437,0.000126785208378037,0.000122502135281101,0.000117960899953838,0.000113210405380578,0.000108188827741880,0.000102945185568431,9.74338509893395e-05,9.16972189395753e-05,8.56795752837880e-05,7.94321200124609e-05,7.29040532996320e-05,6.61402958740868e-05,5.90889832100480e-05,5.17992381686803e-05,4.42231681996236e-05,3.64139273452416e-05,2.83209585120617e-05,1.99914518409441e-05,1.13821152641660e-05,2.54292500196256e-06,-6.56931307847807e-06,-1.59122093433135e-05,-2.55059336623784e-05,-3.53685936727762e-05,-4.54400540476134e-05,-5.57696549960879e-05,-6.63771677159197e-05,-7.71842389048310e-05,-8.82595020734168e-05,-9.95753499456616e-05,-0.000111172538004551,-0.000122985841755254,-0.000135054315559792,-0.000147323340028758,-0.000159835502885448,-0.000172526661615128,-0.000185423796157919,-0.000198461680990827,-0.000211686103776358,-0.000225058353634012,-0.000238655847061092,-0.000252470845359277,-0.000266585641496571,-0.000280950954741698,-0.000295511902559795,-0.000309919414694471,-0.000323659926009783,-0.000339293525390544,-0.000353941141697300,-0.000368940000057705,-0.000384067514638573,-0.000399266350183937,-0.000414577645869965,-0.000429943114160157,-0.000445402498835314,-0.000460893486056534,-0.000476456042465615,-0.000492024291725858,-0.000507642405972521,-0.000523252423509115,-0.000538894259811931,-0.000554498478923738,-0.000570099611785535,-0.000585633180585017,-0.000601147785670214,-0.000616582149430934,-0.000631961341646842,-0.000647217512299709,-0.000662409644078374,-0.000677469185517700,-0.000692398308266890,-0.000707187922078236,-0.000721827244014612,-0.000736275102735415,-0.000750555856424403,-0.000764607158999140,-0.000778464178550299,-0.000792065298356889,-0.000805441953162336,-0.000818525790246944,-0.000831355249152882,-0.000843861392426555,-0.000856081097708123,-0.000867944860467514,-0.000879485936148014,-0.000890637435221489,-0.000901438137545911,-0.000911814674348623,-0.000921802153901256,-0.000931338658585429,-0.000940458537303633,-0.000949086060018919,-0.000957261986294197,-0.000964921530360863,-0.000972089771447915,-0.000978711973604225,-0.000984810266900323,-0.000990324037377398,-0.000995284735363989,-0.000999628829097016,-0.00100338311657393,-0.00100648840624336,-0.00100897383676334,-0.00101077273655536,-0.00101192337158328,-0.00101235351098358,-0.00101210187328732,-0.00101110505554398,-0.00100939628756502,-0.00100690493368262,-0.00100367346039378,-0.000999635478971382,-0.000994825353485726,-0.000989177746920166,-0.000982732992102768,-0.000975422972123854,-0.000967288065601140,-0.000958264039244553,-0.000948385818031927,-0.000937596065739797,-0.000925928377171722,-0.000913320924087787,-0.000899806758135586,-0.000885336927490748,-0.000869930983136534,-0.000853551942177660,-0.000836221973888521,-0.000817881480076714,-0.000798585263593298,-0.000778267181962160,-0.000756965390748311,-0.000734626928834411,-0.000711299284505064,-0.000686917624821226,-0.000661526837259463,-0.000635071485390282,-0.000607599524831393,-0.000579052925211444,-0.000549475774246958,-0.000518807513898672,-0.000487098981962924,-0.000454297513830191,-0.000420448886202128,-0.000385501049354140,-0.000349502134594296,-0.000312408368692861,-0.000274245563333016,-0.000235052758243438,-0.000194716827169965,-0.000153368898219367,-0.000110976175697161,-6.74507400684875e-05,-2.28794438080612e-05,2.27625335533632e-05,6.94323009513472e-05,0.000117191616338152,0.000165983944139365,0.000215841768100563,0.000266689624197460,0.000318565611876379,0.000371418819152364,0.000425307568836569,0.000480188155160603,0.000536114684986784,0.000593018858295846,0.000650913251078771,0.000709697566461656,0.000769388541859999,0.000829970522307472,0.000891610805735191,0.000954285968006834,0.00101744016505917,0.00108179886173261,0.00114690729131758,0.00121285144909704,0.00127964433977491,0.00134721656050543,0.00141559926007632,0.00148472334593358,0.00155462231617105,0.00162522635776969,0.00169657224734207,0.00176858246215739,0.00184128292924902,0.00191459755274841,0.00198855971401751,0.00206309700201124,0.00213823884606977,0.00221390072824991,0.00229011116259731,0.00236680486462171,0.00244400829749066,0.00252162755198125,0.00259969797970943,0.00267815865187951,0.00275699238707736,0.00283615837128198,0.00291565612145279,0.00299541841485039,0.00307546448692997,0.00315571496135086,0.00323618817147553,0.00331680810872238,0.00339759728367135,0.00347847095160714,0.00355945302176549,0.00364045801355713,0.00372150918167257,0.00380252443146705,0.00388352332653804,0.00396441939294586,0.00404523950814305,0.00412589779852082,0.00420640941352418,0.00428669196063306,0.00436677183715796,0.00444655832219488,0.00452607142285969,0.00460523299812716,0.00468405241174943,0.00476245700439543,0.00484045896527067,0.00491797738602819,0.00499502806288758,0.00507153449172557,0.00514750647263878,0.00522286825169297,0.00529763742455926,0.00537172869163968,0.00544516922139953,0.00551787118353883,0.00558985293524430,0.00566103527103497,0.00573144265123906,0.00580098564991597,0.00586968765119188,0.00593746838994351,0.00600434983803940,0.00607024842906037,0.00613519087025290,0.00619909431264470,0.00626198447081014,0.00632378326795535,0.00638450943149585,0.00644408895040196,0.00650254590547329,0.00655980744711141,0.00661588627447452,0.00667072799918385,0.00672432819037122,0.00677663947138492,0.00682768236374154,0.00687736905905883,0.00692574000472932,0.00697272478704275,0.00701834618360901,0.00706252998982802,0.00710531700968379,0.00714663590495778,0.00718651793431919,0.00722488904769741,0.00726178404355198,0.00729713536227356,0.00733098161601232,0.00736324911867787,0.00739397283653819,0.00742308861567995,0.00745063588114617,0.00747655471033473,0.00750088351933540,0.00752358594716239,0.00754458981333914,0.00756409696501951,0.00758179825183088,0.00759784259849564,0.00761229933332533,0.00762504302649761,0.00763608902775804,0.00764540147060972,0.00765305616151820,0.00765901092419252,0.00766330901046047,0.00766587753821448,0.00766675273805942,0.00766587753821448,0.00766330901046047,0.00765901092419252,0.00765305616151820,0.00764540147060972,0.00763608902775804,0.00762504302649761,0.00761229933332533,0.00759784259849564,0.00758179825183088,0.00756409696501951,0.00754458981333914,0.00752358594716239,0.00750088351933540,0.00747655471033473,0.00745063588114617,0.00742308861567995,0.00739397283653819,0.00736324911867787,0.00733098161601232,0.00729713536227356,0.00726178404355198,0.00722488904769741,0.00718651793431919,0.00714663590495778,0.00710531700968379,0.00706252998982802,0.00701834618360901,0.00697272478704275,0.00692574000472932,0.00687736905905883,0.00682768236374154,0.00677663947138492,0.00672432819037122,0.00667072799918385,0.00661588627447452,0.00655980744711141,0.00650254590547329,0.00644408895040196,0.00638450943149585,0.00632378326795535,0.00626198447081014,0.00619909431264470,0.00613519087025290,0.00607024842906037,0.00600434983803940,0.00593746838994351,0.00586968765119188,0.00580098564991597,0.00573144265123906,0.00566103527103497,0.00558985293524430,0.00551787118353883,0.00544516922139953,0.00537172869163968,0.00529763742455926,0.00522286825169297,0.00514750647263878,0.00507153449172557,0.00499502806288758,0.00491797738602819,0.00484045896527067,0.00476245700439543,0.00468405241174943,0.00460523299812716,0.00452607142285969,0.00444655832219488,0.00436677183715796,0.00428669196063306,0.00420640941352418,0.00412589779852082,0.00404523950814305,0.00396441939294586,0.00388352332653804,0.00380252443146705,0.00372150918167257,0.00364045801355713,0.00355945302176549,0.00347847095160714,0.00339759728367135,0.00331680810872238,0.00323618817147553,0.00315571496135086,0.00307546448692997,0.00299541841485039,0.00291565612145279,0.00283615837128198,0.00275699238707736,0.00267815865187951,0.00259969797970943,0.00252162755198125,0.00244400829749066,0.00236680486462171,0.00229011116259731,0.00221390072824991,0.00213823884606977,0.00206309700201124,0.00198855971401751,0.00191459755274841,0.00184128292924902,0.00176858246215739,0.00169657224734207,0.00162522635776969,0.00155462231617105,0.00148472334593358,0.00141559926007632,0.00134721656050543,0.00127964433977491,0.00121285144909704,0.00114690729131758,0.00108179886173261,0.00101744016505917,0.000954285968006834,0.000891610805735191,0.000829970522307472,0.000769388541859999,0.000709697566461656,0.000650913251078771,0.000593018858295846,0.000536114684986784,0.000480188155160603,0.000425307568836569,0.000371418819152364,0.000318565611876379,0.000266689624197460,0.000215841768100563,0.000165983944139365,0.000117191616338152,6.94323009513472e-05,2.27625335533632e-05,-2.28794438080612e-05,-6.74507400684875e-05,-0.000110976175697161,-0.000153368898219367,-0.000194716827169965,-0.000235052758243438,-0.000274245563333016,-0.000312408368692861,-0.000349502134594296,-0.000385501049354140,-0.000420448886202128,-0.000454297513830191,-0.000487098981962924,-0.000518807513898672,-0.000549475774246958,-0.000579052925211444,-0.000607599524831393,-0.000635071485390282,-0.000661526837259463,-0.000686917624821226,-0.000711299284505064,-0.000734626928834411,-0.000756965390748311,-0.000778267181962160,-0.000798585263593298,-0.000817881480076714,-0.000836221973888521,-0.000853551942177660,-0.000869930983136534,-0.000885336927490748,-0.000899806758135586,-0.000913320924087787,-0.000925928377171722,-0.000937596065739797,-0.000948385818031927,-0.000958264039244553,-0.000967288065601140,-0.000975422972123854,-0.000982732992102768,-0.000989177746920166,-0.000994825353485726,-0.000999635478971382,-0.00100367346039378,-0.00100690493368262,-0.00100939628756502,-0.00101110505554398,-0.00101210187328732,-0.00101235351098358,-0.00101192337158328,-0.00101077273655536,-0.00100897383676334,-0.00100648840624336,-0.00100338311657393,-0.000999628829097016,-0.000995284735363989,-0.000990324037377398,-0.000984810266900323,-0.000978711973604225,-0.000972089771447915,-0.000964921530360863,-0.000957261986294197,-0.000949086060018919,-0.000940458537303633,-0.000931338658585429,-0.000921802153901256,-0.000911814674348623,-0.000901438137545911,-0.000890637435221489,-0.000879485936148014,-0.000867944860467514,-0.000856081097708123,-0.000843861392426555,-0.000831355249152882,-0.000818525790246944,-0.000805441953162336,-0.000792065298356889,-0.000778464178550299,-0.000764607158999140,-0.000750555856424403,-0.000736275102735415,-0.000721827244014612,-0.000707187922078236,-0.000692398308266890,-0.000677469185517700,-0.000662409644078374,-0.000647217512299709,-0.000631961341646842,-0.000616582149430934,-0.000601147785670214,-0.000585633180585017,-0.000570099611785535,-0.000554498478923738,-0.000538894259811931,-0.000523252423509115,-0.000507642405972521,-0.000492024291725858,-0.000476456042465615,-0.000460893486056534,-0.000445402498835314,-0.000429943114160157,-0.000414577645869965,-0.000399266350183937,-0.000384067514638573,-0.000368940000057705,-0.000353941141697300,-0.000339293525390544,-0.000323659926009783,-0.000309919414694471,-0.000295511902559795,-0.000280950954741698,-0.000266585641496571,-0.000252470845359277,-0.000238655847061092,-0.000225058353634012,-0.000211686103776358,-0.000198461680990827,-0.000185423796157919,-0.000172526661615128,-0.000159835502885448,-0.000147323340028758,-0.000135054315559792,-0.000122985841755254,-0.000111172538004551,-9.95753499456616e-05,-8.82595020734168e-05,-7.71842389048310e-05,-6.63771677159197e-05,-5.57696549960879e-05,-4.54400540476134e-05,-3.53685936727762e-05,-2.55059336623784e-05,-1.59122093433135e-05,-6.56931307847807e-06,2.54292500196256e-06,1.13821152641660e-05,1.99914518409441e-05,2.83209585120617e-05,3.64139273452416e-05,4.42231681996236e-05,5.17992381686803e-05,5.90889832100480e-05,6.61402958740868e-05,7.29040532996320e-05,7.94321200124609e-05,8.56795752837880e-05,9.16972189395753e-05,9.74338509893395e-05,0.000102945185568431,0.000108188827741880,0.000113210405380578,0.000117960899953838,0.000122502135281101,0.000126785208378037,0.000130850035109437,0.000134668974387507,0.000138276287207714,0.000141639738605470,0.000144804508496339,0.000147727128052281,0.000150454334301420,0.000152947120543503,0.000155251833704162,0.000157324536310168,0.000159218674345988,0.000160885927205586,0.000162382142477805,0.000163662215298123,0.000164776375709852,0.000165681305169779,0.000166435386587554,0.000166989668474572,0.000167395706802922,0.000167612685064898,0.000167695463616470,0.000167596209337449,0.000167371723214612,0.000166977599482095,0.000166463912498038,0.000165793466249417,0.000165012261291061,0.000164080891749098,0.000163049900688652,0.000161881818578263,0.000160617329536024,0.000159227356866869,0.000157751952845765,0.000156153810147281,0.000154490393439512,0.000152708088942437,0.000150860785731569,0.000148909574313128,0.000146905389952095,0.000144799075700251,0.000142648282511158,0.000140409372311920,0.000138133925206034,0.000135776345148193,0.000133391093408533,0.000130933154839445,0.000128458153510528,0.000125918478995377,0.000123362723034365,0.000120750142158453,0.000118132310510059,0.000115468005236082,0.000112800194579151,0.000110108789624883,0.000107392581745282,0.000104680883398715,0.000101967709425591,9.92159696337733e-05,9.65017643446230e-05,9.37690208083129e-05,9.10537173189225e-05,8.83159101463416e-05,8.56184309621784e-05,8.29159682127781e-05,8.02551921214269e-05,7.75832509210173e-05,7.49510402638714e-05,7.23168342051082e-05,6.97343785522134e-05,6.71538849242156e-05,6.46252421778573e-05,6.21002944473329e-05,5.96287108106926e-05,5.71675693367720e-05,5.47674639762573e-05,5.23815850259164e-05,5.00586518032772e-05,0.000524469736335486]; */
    i6 = b_num->size[0] * b_num->size[1];
    b_num->size[0] = 1;
    b_num->size[1] = 723;
    emxEnsureCapacity_real_T(b_num, i6);
    for (i6 = 0; i6 < 723; i6++) {
      b_num->data[b_num->size[0] * i6] = dv2[i6];
    }

    /* 'dynIdenf:561' den = 1; */
    /*  transient width 3 */
  } else if ((0.7 <= fpass) && (fpass < 1.41)) {
    /* 'dynIdenf:562' elseif fdisc(3) <= fpass && fpass < fdisc(4) */
    /* 'dynIdenf:563' num = [0.000442188562486940,-0.000110657522425279,-9.98455797241236e-05,-9.14897478869632e-05,-8.52361666827068e-05,-8.06959517311162e-05,-7.76411712637086e-05,-7.57899038105713e-05,-7.49994662714565e-05,-7.50603085922735e-05,-7.58901010430049e-05,-7.73273137334636e-05,-7.93246119168430e-05,-8.17496582065717e-05,-8.45883226012555e-05,-8.77382876114558e-05,-9.11922203705935e-05,-9.48664849214488e-05,-9.87761589957704e-05,-0.000102844175988080,-0.000107090784794005,-0.000111440909688352,-0.000115928813200131,-0.000120477160220656,-0.000125117715117357,-0.000129794392580865,-0.000134536155367815,-0.000139275765199167,-0.000144062787810614,-0.000148824518200912,-0.000153595780710076,-0.000158325434209622,-0.000163036019066201,-0.000167679678011588,-0.000172283561123639,-0.000176796250544662,-0.000181250958338245,-0.000185583835789150,-0.000189833189429288,-0.000193936404147956,-0.000197929142837483,-0.000201752638548829,-0.000205448168113863,-0.000208952400045264,-0.000212298526996463,-0.000215431951730521,-0.000218370006534548,-0.000221087779149871,-0.000223568899355575,-0.000225816461603436,-0.000227823347060940,-0.000229527953926879,-0.000230995790176611,-0.000232140202409209,-0.000233000793400310,-0.000233539688331506,-0.000233786335275040,-0.000233659708473604,-0.000233206053130778,-0.000232378805458336,-0.000231219358034977,-0.000229656785942769,-0.000227727315820203,-0.000225380645567152,-0.000222658595449340,-0.000219509062812890,-0.000215967665766597,-0.000211968448558799,-0.000207561615470521,-0.000202704140518313,-0.000197438201725542,-0.000191686800015090,-0.000185498343418373,-0.000178850781994576,-0.000171780962719214,-0.000164197976718070,-0.000156190026970194,-0.000147717289209701,-0.000138762707983870,-0.000129375123024395,-0.000119521441536058,-0.000109202912399607,-9.84484818706821e-05,-8.72242658070373e-05,-7.55832297779589e-05,-6.34806820680050e-05,-5.09741936199240e-05,-3.80211853750399e-05,-2.46764348809454e-05,-1.09007929801398e-05,3.24223471168255e-06,1.77916103388651e-05,3.26898421180593e-05,4.79656255396375e-05,6.35660022709180e-05,7.95144815924272e-05,9.57553359401000e-05,0.000112311171021065,0.000129115507243385,0.000146208523407364,0.000163510711781697,0.000181055990914089,0.000198779793426531,0.000216700488791655,0.000234745799216152,0.000252948215081320,0.000271225228480540,0.000289602307914910,0.000308006636583080,0.000326458055947502,0.000344882743063236,0.000363296018691483,0.000381627256894602,0.000399890129212652,0.000418002629144219,0.000435987929285981,0.000453767065420231,0.000471355470721825,0.000488669514310265,0.000505726649129089,0.000522444131355595,0.000538838194008873,0.000554824556931667,0.000570381891315664,0.000585743639518650,0.000599712098309406,0.000614748659898585,0.000628041748771862,0.000640655606610034,0.000653112783043952,0.000665075491543359,0.000676258910887898,0.000686503321973215,0.000695927732178446,0.000704573751514583,0.000712526700146433,0.000719690675126741,0.000726032690090996,0.000731413762698865,0.000735829824257724,0.000739201036584629,0.000741561040118516,0.000742864552642016,0.000743163355418230,0.000742393457652286,0.000740580422933247,0.000737638852671666,0.000733592725282046,0.000728345882414530,0.000721918005128956,0.000714242597601836,0.000705348970537505,0.000695168304813553,0.000683764815943311,0.000671073326074807,0.000657138705447006,0.000641920386361616,0.000625443291086012,0.000607663334750444,0.000588610658995366,0.000568235949272277,0.000546583747008114,0.000523600383315921,0.000499344029128501,0.000473767964040842,0.000446928977698078,0.000418784086169086,0.000389399293243541,0.000358735991847034,0.000326863787781629,0.000293764207557915,0.000259496360454215,0.000224062935335397,0.000187495405025190,0.000149812495177780,0.000111067243459504,7.12362989078123e-05,3.04355904575503e-05,-1.13692218873472e-05,-5.40934650222404e-05,-9.77379655364106e-05,-0.000142218358284181,-0.000187553808130795,-0.000233629127780396,-0.000280437394148562,-0.000327882387234898,-0.000375979311793093,-0.000424620534910132,-0.000473783966984720,-0.000523354431241435,-0.000573322387778651,-0.000623584537328208,-0.000674136790271356,-0.000724845720752713,-0.000775679872534808,-0.000826525753613857,-0.000877391824269090,-0.000928142377229690,-0.000978713698804689,-0.00102899532653899,-0.00107902169607718,-0.00112860084065069,-0.00117766930819020,-0.00122622838046955,-0.00127408587839297,-0.00132123592440593,-0.00136757653727756,-0.00141299760787155,-0.00145745599673916,-0.00150081806841540,-0.00154305312456678,-0.00158402096786208,-0.00162368591796622,-0.00166191133409511,-0.00169866257719637,-0.00173379937250004,-0.00176729056400442,-0.00179900352340403,-0.00182889587715962,-0.00185683950506910,-0.00188279319885088,-0.00190663066128216,-0.00192831551866106,-0.00194771439987576,-0.00196480771427863,-0.00197945287606819,-0.00199162656239649,-0.00200120800579529,-0.00200816746672251,-0.00201238080748018,-0.00201383794301468,-0.00201240896567734,-0.00200808023856785,-0.00200073871634229,-0.00199036976433981,-0.00197686642391883,-0.00196021582547589,-0.00194031922937830,-0.00191717081627876,-0.00189066708440931,-0.00186081656382531,-0.00182752251490315,-0.00179079178272426,-0.00175053388265721,-0.00170676987426391,-0.00165941766709796,-0.00160850114906796,-0.00155394758862404,-0.00149574927976897,-0.00143399143098721,-0.00136831299612247,-0.00129932561355133,-0.00122639037844494,-0.00114972371991373,-0.00106955335669744,-0.000985669737727666,-0.000898040252148599,-0.000806691221734975,-0.000711774350863481,-0.000613254881865684,-0.000511158423926788,-0.000405424015373930,-0.000296148158273841,-0.000183347113779057,-6.71406343458254e-05,5.24827105073237e-05,0.000175453268935852,0.000301784361577890,0.000431367973228352,0.000564188884829856,0.000700117036511363,0.000839124467772486,0.000981073903589326,0.00112595807659575,0.00127365802791523,0.00142414162003024,0.00157727829340564,0.00173303853158450,0.00189124867543978,0.00205186680827992,0.00221474221786150,0.00237979773402877,0.00254690057435000,0.00271597700409001,0.00288688153935119,0.00305954623686013,0.00323380501568201,0.00340958986361004,0.00358671627162477,0.00376510452014456,0.00394456924668566,0.00412502877604286,0.00430629431747686,0.00448828951085397,0.00467083359079199,0.00485383705373487,0.00503712657977312,0.00522058688840571,0.00540405069335647,0.00558738694575325,0.00577040966215496,0.00595303526933757,0.00613504619542937,0.00631636903698379,0.00649681040951683,0.00667625932295169,0.00685451323038891,0.00703148439372521,0.00720696612105708,0.00738087270312673,0.00755301874920811,0.00772331650642570,0.00789155222299408,0.00805761456478647,0.00822131875322062,0.00838259356390575,0.00854125907927656,0.00869722344317690,0.00885027046921507,0.00900031283405329,0.00914719655546697,0.00929088009437544,0.00943115585333465,0.00956789639043661,0.00970096259291177,0.00983038335230680,0.00995589630104736,0.0100773480018851,0.0101948546792788,0.0103080638144166,0.0104170283034853,0.0105216018540367,0.0106216692956179,0.0107171859212136,0.0108080160980720,0.0108941318936983,0.0109753996703859,0.0110517981558641,0.0111232094631204,0.0111896227087287,0.0112509221691577,0.0113071058431389,0.0113580739117312,0.0114038233874121,0.0114442687045384,0.0114794193958816,0.0115091971892624,0.0115336216751115,0.0115526155187787,0.0115662222486294,0.0115743673395176,0.0115770980653248,0.0115743673395176,0.0115662222486294,0.0115526155187787,0.0115336216751115,0.0115091971892624,0.0114794193958816,0.0114442687045384,0.0114038233874121,0.0113580739117312,0.0113071058431389,0.0112509221691577,0.0111896227087287,0.0111232094631204,0.0110517981558641,0.0109753996703859,0.0108941318936983,0.0108080160980720,0.0107171859212136,0.0106216692956179,0.0105216018540367,0.0104170283034853,0.0103080638144166,0.0101948546792788,0.0100773480018851,0.00995589630104736,0.00983038335230680,0.00970096259291177,0.00956789639043661,0.00943115585333465,0.00929088009437544,0.00914719655546697,0.00900031283405329,0.00885027046921507,0.00869722344317690,0.00854125907927656,0.00838259356390575,0.00822131875322062,0.00805761456478647,0.00789155222299408,0.00772331650642570,0.00755301874920811,0.00738087270312673,0.00720696612105708,0.00703148439372521,0.00685451323038891,0.00667625932295169,0.00649681040951683,0.00631636903698379,0.00613504619542937,0.00595303526933757,0.00577040966215496,0.00558738694575325,0.00540405069335647,0.00522058688840571,0.00503712657977312,0.00485383705373487,0.00467083359079199,0.00448828951085397,0.00430629431747686,0.00412502877604286,0.00394456924668566,0.00376510452014456,0.00358671627162477,0.00340958986361004,0.00323380501568201,0.00305954623686013,0.00288688153935119,0.00271597700409001,0.00254690057435000,0.00237979773402877,0.00221474221786150,0.00205186680827992,0.00189124867543978,0.00173303853158450,0.00157727829340564,0.00142414162003024,0.00127365802791523,0.00112595807659575,0.000981073903589326,0.000839124467772486,0.000700117036511363,0.000564188884829856,0.000431367973228352,0.000301784361577890,0.000175453268935852,5.24827105073237e-05,-6.71406343458254e-05,-0.000183347113779057,-0.000296148158273841,-0.000405424015373930,-0.000511158423926788,-0.000613254881865684,-0.000711774350863481,-0.000806691221734975,-0.000898040252148599,-0.000985669737727666,-0.00106955335669744,-0.00114972371991373,-0.00122639037844494,-0.00129932561355133,-0.00136831299612247,-0.00143399143098721,-0.00149574927976897,-0.00155394758862404,-0.00160850114906796,-0.00165941766709796,-0.00170676987426391,-0.00175053388265721,-0.00179079178272426,-0.00182752251490315,-0.00186081656382531,-0.00189066708440931,-0.00191717081627876,-0.00194031922937830,-0.00196021582547589,-0.00197686642391883,-0.00199036976433981,-0.00200073871634229,-0.00200808023856785,-0.00201240896567734,-0.00201383794301468,-0.00201238080748018,-0.00200816746672251,-0.00200120800579529,-0.00199162656239649,-0.00197945287606819,-0.00196480771427863,-0.00194771439987576,-0.00192831551866106,-0.00190663066128216,-0.00188279319885088,-0.00185683950506910,-0.00182889587715962,-0.00179900352340403,-0.00176729056400442,-0.00173379937250004,-0.00169866257719637,-0.00166191133409511,-0.00162368591796622,-0.00158402096786208,-0.00154305312456678,-0.00150081806841540,-0.00145745599673916,-0.00141299760787155,-0.00136757653727756,-0.00132123592440593,-0.00127408587839297,-0.00122622838046955,-0.00117766930819020,-0.00112860084065069,-0.00107902169607718,-0.00102899532653899,-0.000978713698804689,-0.000928142377229690,-0.000877391824269090,-0.000826525753613857,-0.000775679872534808,-0.000724845720752713,-0.000674136790271356,-0.000623584537328208,-0.000573322387778651,-0.000523354431241435,-0.000473783966984720,-0.000424620534910132,-0.000375979311793093,-0.000327882387234898,-0.000280437394148562,-0.000233629127780396,-0.000187553808130795,-0.000142218358284181,-9.77379655364106e-05,-5.40934650222404e-05,-1.13692218873472e-05,3.04355904575503e-05,7.12362989078123e-05,0.000111067243459504,0.000149812495177780,0.000187495405025190,0.000224062935335397,0.000259496360454215,0.000293764207557915,0.000326863787781629,0.000358735991847034,0.000389399293243541,0.000418784086169086,0.000446928977698078,0.000473767964040842,0.000499344029128501,0.000523600383315921,0.000546583747008114,0.000568235949272277,0.000588610658995366,0.000607663334750444,0.000625443291086012,0.000641920386361616,0.000657138705447006,0.000671073326074807,0.000683764815943311,0.000695168304813553,0.000705348970537505,0.000714242597601836,0.000721918005128956,0.000728345882414530,0.000733592725282046,0.000737638852671666,0.000740580422933247,0.000742393457652286,0.000743163355418230,0.000742864552642016,0.000741561040118516,0.000739201036584629,0.000735829824257724,0.000731413762698865,0.000726032690090996,0.000719690675126741,0.000712526700146433,0.000704573751514583,0.000695927732178446,0.000686503321973215,0.000676258910887898,0.000665075491543359,0.000653112783043952,0.000640655606610034,0.000628041748771862,0.000614748659898585,0.000599712098309406,0.000585743639518650,0.000570381891315664,0.000554824556931667,0.000538838194008873,0.000522444131355595,0.000505726649129089,0.000488669514310265,0.000471355470721825,0.000453767065420231,0.000435987929285981,0.000418002629144219,0.000399890129212652,0.000381627256894602,0.000363296018691483,0.000344882743063236,0.000326458055947502,0.000308006636583080,0.000289602307914910,0.000271225228480540,0.000252948215081320,0.000234745799216152,0.000216700488791655,0.000198779793426531,0.000181055990914089,0.000163510711781697,0.000146208523407364,0.000129115507243385,0.000112311171021065,9.57553359401000e-05,7.95144815924272e-05,6.35660022709180e-05,4.79656255396375e-05,3.26898421180593e-05,1.77916103388651e-05,3.24223471168255e-06,-1.09007929801398e-05,-2.46764348809454e-05,-3.80211853750399e-05,-5.09741936199240e-05,-6.34806820680050e-05,-7.55832297779589e-05,-8.72242658070373e-05,-9.84484818706821e-05,-0.000109202912399607,-0.000119521441536058,-0.000129375123024395,-0.000138762707983870,-0.000147717289209701,-0.000156190026970194,-0.000164197976718070,-0.000171780962719214,-0.000178850781994576,-0.000185498343418373,-0.000191686800015090,-0.000197438201725542,-0.000202704140518313,-0.000207561615470521,-0.000211968448558799,-0.000215967665766597,-0.000219509062812890,-0.000222658595449340,-0.000225380645567152,-0.000227727315820203,-0.000229656785942769,-0.000231219358034977,-0.000232378805458336,-0.000233206053130778,-0.000233659708473604,-0.000233786335275040,-0.000233539688331506,-0.000233000793400310,-0.000232140202409209,-0.000230995790176611,-0.000229527953926879,-0.000227823347060940,-0.000225816461603436,-0.000223568899355575,-0.000221087779149871,-0.000218370006534548,-0.000215431951730521,-0.000212298526996463,-0.000208952400045264,-0.000205448168113863,-0.000201752638548829,-0.000197929142837483,-0.000193936404147956,-0.000189833189429288,-0.000185583835789150,-0.000181250958338245,-0.000176796250544662,-0.000172283561123639,-0.000167679678011588,-0.000163036019066201,-0.000158325434209622,-0.000153595780710076,-0.000148824518200912,-0.000144062787810614,-0.000139275765199167,-0.000134536155367815,-0.000129794392580865,-0.000125117715117357,-0.000120477160220656,-0.000115928813200131,-0.000111440909688352,-0.000107090784794005,-0.000102844175988080,-9.87761589957704e-05,-9.48664849214488e-05,-9.11922203705935e-05,-8.77382876114558e-05,-8.45883226012555e-05,-8.17496582065717e-05,-7.93246119168430e-05,-7.73273137334636e-05,-7.58901010430049e-05,-7.50603085922735e-05,-7.49994662714565e-05,-7.57899038105713e-05,-7.76411712637086e-05,-8.06959517311162e-05,-8.52361666827068e-05,-9.14897478869632e-05,-9.98455797241236e-05,-0.000110657522425279,0.000442188562486940]; */
    i6 = b_num->size[0] * b_num->size[1];
    b_num->size[0] = 1;
    b_num->size[1] = 703;
    emxEnsureCapacity_real_T(b_num, i6);
    for (i6 = 0; i6 < 703; i6++) {
      b_num->data[b_num->size[0] * i6] = dv3[i6];
    }

    /* 'dynIdenf:564' den = 1; */
    /*  transient width 2 */
  } else if ((1.41 <= fpass) && (fpass < 2.82)) {
    /* 'dynIdenf:565' elseif fdisc(4) <= fpass && fpass < fdisc(5) */
    /* 'dynIdenf:566' num = [0.000532723527805931,6.90845510945624e-05,7.32098567809705e-05,7.72828817351650e-05,8.12279651616150e-05,8.50730286220366e-05,8.87432006630032e-05,9.22632831417983e-05,9.55562989757991e-05,9.86445413311909e-05,0.000101452080284138,0.000104006770866074,0.000106234517995692,0.000108156839530409,0.000109697674089670,0.000110875128450003,0.000111621547646505,0.000111970687838987,0.000111849884291550,0.000111282242197444,0.000110192908129537,0.000108623081333228,0.000106515908723708,0.000103892652709569,0.000100677354878430,9.69389345979157e-05,9.26226107885545e-05,8.77215523829069e-05,8.22521094557147e-05,7.62066042336095e-05,6.95679622187701e-05,6.23782262138511e-05,5.45955644700710e-05,4.62793390590956e-05,3.73910309205438e-05,2.79965806474123e-05,1.80639059740427e-05,7.66317268571757e-06,-3.23421059797320e-06,-1.45519362022193e-05,-2.63058607512367e-05,-3.84220186694711e-05,-5.09117413208545e-05,-6.36905034503717e-05,-7.67661277548565e-05,-9.00437851570947e-05,-0.000103535649615391,-0.000117140359067039,-0.000130853687142025,-0.000144579619052287,-0.000158316330340026,-0.000171948811459607,-0.000185476850131979,-0.000198794093871137,-0.000211878668856251,-0.000224633417558677,-0.000237034781155883,-0.000248978615566231,-0.000260448698693465,-0.000271332401327684,-0.000281620445417402,-0.000291202506939356,-0.000300068275069311,-0.000308108733299146,-0.000315319515874973,-0.000321597411331245,-0.000326934577519601,-0.000331241670082228,-0.000334505681801338,-0.000336650130867220,-0.000337671660988197,-0.000337481977392415,-0.000336106485795569,-0.000333455311220054,-0.000329550165379097,-0.000324325703452802,-0.000317811953425090,-0.000309939057334787,-0.000300754241979969,-0.000290205543465510,-0.000278344626220354,-0.000265132115339844,-0.000250628169585483,-0.000234812182294999,-0.000217759861217467,-0.000199454100181876,-0.000179986263437099,-0.000159353716962562,-0.000137659970668377,-0.000114912452474471,-9.12241339784792e-05,-6.66119812977303e-05,-4.11936497765334e-05,-1.50167581043864e-05,1.18423295007193e-05,3.92565201720357e-05,6.71914198631196e-05,9.55422954605711e-05,0.000124142970867157,0.000152973890653300,0.000181870653206300,0.000210742122882789,0.000239413971555361,0.000267824245622387,0.000295802383646803,0.000323268902375579,0.000350044687543034,0.000376064507731952,0.000401179043921861,0.000425332695313922,0.000448370380215483,0.000470221483334989,0.000490702014082546,0.000509715635440249,0.000527055147364876,0.000542619714948337,0.000556249759298917,0.000567969650086771,0.000577807639663091,0.000585912834301797,0.000591969798291343,0.000594305253361897,0.000595860696550544,0.000594399332766709,0.000590571980802667,0.000584290296204713,0.000575463781540990,0.000564118894623614,0.000550184412602726,0.000533706380462012,0.000514638074566508,0.000493050047060811,0.000468913400443583,0.000442313253832272,0.000413241328784440,0.000381803274818962,0.000348023547052433,0.000312034702218758,0.000273870106753618,0.000233681198056428,0.000191528676530613,0.000147596411187783,0.000101960287388390,5.48016986829699e-05,6.22557969660877e-06,-4.35362854074230e-05,-9.43882319573520e-05,-0.000146110116296165,-0.000198526522947237,-0.000251444850832230,-0.000304678802873628,-0.000357993619188933,-0.000411225582184539,-0.000464109785038842,-0.000516482314140015,-0.000568068863470274,-0.000618700274194817,-0.000668100604982624,-0.000716098863454397,-0.000762415967824261,-0.000806880267618045,-0.000849226128650888,-0.000889282065659916,-0.000926788086023849,-0.000961585620209408,-0.000993422610739849,-0.00102216328155141,-0.00104756127527116,-0.00106949414293076,-0.00108774785205145,-0.00110221455827652,-0.00111269268098981,-0.00111911216097625,-0.00112129581941446,-0.00111918690073059,-0.00111265273925094,-0.00110165907175774,-0.00108610324359315,-0.00106598733591590,-0.00104123445629234,-0.00101189035018557,-0.000977909727087050,-0.000939371894413249,-0.000896270612403707,-0.000848725929165545,-0.000796768602540623,-0.000740552848045237,-0.000680156633333393,-0.000615760383245651,-0.000547494664223412,-0.000475569076402148,-0.000400142017380367,-0.000321485002025148,-0.000239764850790327,-0.000155292447891851,-6.82857751072828e-05,2.09289179787751e-05,0.000112109297622250,0.000204883420977705,0.000298981574344912,0.000394014750820599,0.000489685968326800,0.000585583218912135,0.000681382407825788,0.000776660024773496,0.000871085500740466,0.000964216924345655,0.00105571593690070,0.00114513675352477,0.00123214019124730,0.00131628101314346,0.00139722080626986,0.00147451648209021,0.00154782628278128,0.00161675648683660,0.00168089004027324,0.00173997814549448,0.00179355111402787,0.00184129196834712,0.00188297511181281,0.00191820056831365,0.00194671777139435,0.00196823189611401,0.00198256648146320,0.00198943348273961,0.00198866068745110,0.00198000124960408,0.00196335543997715,0.00193855339959633,0.00190555323215255,0.00186424362293629,0.00181463879312254,0.00175667524340481,0.00169041467832458,0.00161583185530486,0.00153304327567335,0.00144212179286118,0.00134331731841221,0.00123681071188208,0.00112283496099542,0.00100142115160551,0.000872866227640752,0.000738010318540991,0.000596478847215025,0.000449101772719016,0.000296158967640106,0.000138068317705179,-2.46513410261154e-05,-0.000191559202568725,-0.000362075082993075,-0.000535698168265624,-0.000711795700978104,-0.000889822901721817,-0.00106910560605327,-0.00124905648482067,-0.00142896209429875,-0.00160819122016243,-0.00178600345365117,-0.00196174870426848,-0.00213465415341023,-0.00230404356324101,-0.00246912606313768,-0.00262922616543297,-0.00278355295102357,-0.00293140840277156,-0.00307200120028402,-0.00320467215722319,-0.00332862998939508,-0.00344320039522590,-0.00354766165289471,-0.00364133959867959,-0.00372353791259458,-0.00379364301321045,-0.00385097568858530,-0.00389498262148901,-0.00392502914591490,-0.00394061946397722,-0.00394117996376844,-0.00392627630440225,-0.00389539996021668,-0.00384818977002379,-0.00378421825818613,-0.00370319322404993,-0.00360477304508528,-0.00348875245919541,-0.00335487603222865,-0.00320303723850247,-0.00303306169439496,-0.00284494561585071,-0.00263862046510879,-0.00241416824095201,-0.00217162262420810,-0.00191117876624268,-0.00163296066456587,-0.00133726271460791,-0.00102432543856643,-0.000694533535606485,-0.000348233195879792,1.40895738088142e-05,0.000391991534902159,0.000784879339972465,0.00119222116184219,0.00161332794293357,0.00204757197540031,0.00249417437950656,0.00295242425888869,0.00342146101503524,0.00390048954278040,0.00438858491491148,0.00488486434194883,0.00538835584449848,0.00589810685900003,0.00641307176552732,0.00693228058897330,0.00745460951227381,0.00797904473983032,0.00850445114532105,0.00902977164682324,0.00955383211944169,0.0100755711107632,0.0105938092222330,0.0111074713522684,0.0116153757129336,0.0121164515317209,0.0126095392988345,0.0130935850412763,0.0135674431488346,0.0140300933926742,0.0144804281527954,0.0149174721511036,0.0153401670861071,0.0157475906337756,0.0161387441993727,0.0165127735255988,0.0168687716270319,0.0172058854778318,0.0175234683086044,0.0178204989017597,0.0180965132230385,0.0183508338529866,0.0185827142905703,0.0187916926310988,0.0189772718951223,0.0191390675550726,0.0192766147002736,0.0193896183087324,0.0194777543833272,0.0195408769319307,0.0195787891742262,0.0195914531025123,0.0195787891742262,0.0195408769319307,0.0194777543833272,0.0193896183087324,0.0192766147002736,0.0191390675550726,0.0189772718951223,0.0187916926310988,0.0185827142905703,0.0183508338529866,0.0180965132230385,0.0178204989017597,0.0175234683086044,0.0172058854778318,0.0168687716270319,0.0165127735255988,0.0161387441993727,0.0157475906337756,0.0153401670861071,0.0149174721511036,0.0144804281527954,0.0140300933926742,0.0135674431488346,0.0130935850412763,0.0126095392988345,0.0121164515317209,0.0116153757129336,0.0111074713522684,0.0105938092222330,0.0100755711107632,0.00955383211944169,0.00902977164682324,0.00850445114532105,0.00797904473983032,0.00745460951227381,0.00693228058897330,0.00641307176552732,0.00589810685900003,0.00538835584449848,0.00488486434194883,0.00438858491491148,0.00390048954278040,0.00342146101503524,0.00295242425888869,0.00249417437950656,0.00204757197540031,0.00161332794293357,0.00119222116184219,0.000784879339972465,0.000391991534902159,1.40895738088142e-05,-0.000348233195879792,-0.000694533535606485,-0.00102432543856643,-0.00133726271460791,-0.00163296066456587,-0.00191117876624268,-0.00217162262420810,-0.00241416824095201,-0.00263862046510879,-0.00284494561585071,-0.00303306169439496,-0.00320303723850247,-0.00335487603222865,-0.00348875245919541,-0.00360477304508528,-0.00370319322404993,-0.00378421825818613,-0.00384818977002379,-0.00389539996021668,-0.00392627630440225,-0.00394117996376844,-0.00394061946397722,-0.00392502914591490,-0.00389498262148901,-0.00385097568858530,-0.00379364301321045,-0.00372353791259458,-0.00364133959867959,-0.00354766165289471,-0.00344320039522590,-0.00332862998939508,-0.00320467215722319,-0.00307200120028402,-0.00293140840277156,-0.00278355295102357,-0.00262922616543297,-0.00246912606313768,-0.00230404356324101,-0.00213465415341023,-0.00196174870426848,-0.00178600345365117,-0.00160819122016243,-0.00142896209429875,-0.00124905648482067,-0.00106910560605327,-0.000889822901721817,-0.000711795700978104,-0.000535698168265624,-0.000362075082993075,-0.000191559202568725,-2.46513410261154e-05,0.000138068317705179,0.000296158967640106,0.000449101772719016,0.000596478847215025,0.000738010318540991,0.000872866227640752,0.00100142115160551,0.00112283496099542,0.00123681071188208,0.00134331731841221,0.00144212179286118,0.00153304327567335,0.00161583185530486,0.00169041467832458,0.00175667524340481,0.00181463879312254,0.00186424362293629,0.00190555323215255,0.00193855339959633,0.00196335543997715,0.00198000124960408,0.00198866068745110,0.00198943348273961,0.00198256648146320,0.00196823189611401,0.00194671777139435,0.00191820056831365,0.00188297511181281,0.00184129196834712,0.00179355111402787,0.00173997814549448,0.00168089004027324,0.00161675648683660,0.00154782628278128,0.00147451648209021,0.00139722080626986,0.00131628101314346,0.00123214019124730,0.00114513675352477,0.00105571593690070,0.000964216924345655,0.000871085500740466,0.000776660024773496,0.000681382407825788,0.000585583218912135,0.000489685968326800,0.000394014750820599,0.000298981574344912,0.000204883420977705,0.000112109297622250,2.09289179787751e-05,-6.82857751072828e-05,-0.000155292447891851,-0.000239764850790327,-0.000321485002025148,-0.000400142017380367,-0.000475569076402148,-0.000547494664223412,-0.000615760383245651,-0.000680156633333393,-0.000740552848045237,-0.000796768602540623,-0.000848725929165545,-0.000896270612403707,-0.000939371894413249,-0.000977909727087050,-0.00101189035018557,-0.00104123445629234,-0.00106598733591590,-0.00108610324359315,-0.00110165907175774,-0.00111265273925094,-0.00111918690073059,-0.00112129581941446,-0.00111911216097625,-0.00111269268098981,-0.00110221455827652,-0.00108774785205145,-0.00106949414293076,-0.00104756127527116,-0.00102216328155141,-0.000993422610739849,-0.000961585620209408,-0.000926788086023849,-0.000889282065659916,-0.000849226128650888,-0.000806880267618045,-0.000762415967824261,-0.000716098863454397,-0.000668100604982624,-0.000618700274194817,-0.000568068863470274,-0.000516482314140015,-0.000464109785038842,-0.000411225582184539,-0.000357993619188933,-0.000304678802873628,-0.000251444850832230,-0.000198526522947237,-0.000146110116296165,-9.43882319573520e-05,-4.35362854074230e-05,6.22557969660877e-06,5.48016986829699e-05,0.000101960287388390,0.000147596411187783,0.000191528676530613,0.000233681198056428,0.000273870106753618,0.000312034702218758,0.000348023547052433,0.000381803274818962,0.000413241328784440,0.000442313253832272,0.000468913400443583,0.000493050047060811,0.000514638074566508,0.000533706380462012,0.000550184412602726,0.000564118894623614,0.000575463781540990,0.000584290296204713,0.000590571980802667,0.000594399332766709,0.000595860696550544,0.000594305253361897,0.000591969798291343,0.000585912834301797,0.000577807639663091,0.000567969650086771,0.000556249759298917,0.000542619714948337,0.000527055147364876,0.000509715635440249,0.000490702014082546,0.000470221483334989,0.000448370380215483,0.000425332695313922,0.000401179043921861,0.000376064507731952,0.000350044687543034,0.000323268902375579,0.000295802383646803,0.000267824245622387,0.000239413971555361,0.000210742122882789,0.000181870653206300,0.000152973890653300,0.000124142970867157,9.55422954605711e-05,6.71914198631196e-05,3.92565201720357e-05,1.18423295007193e-05,-1.50167581043864e-05,-4.11936497765334e-05,-6.66119812977303e-05,-9.12241339784792e-05,-0.000114912452474471,-0.000137659970668377,-0.000159353716962562,-0.000179986263437099,-0.000199454100181876,-0.000217759861217467,-0.000234812182294999,-0.000250628169585483,-0.000265132115339844,-0.000278344626220354,-0.000290205543465510,-0.000300754241979969,-0.000309939057334787,-0.000317811953425090,-0.000324325703452802,-0.000329550165379097,-0.000333455311220054,-0.000336106485795569,-0.000337481977392415,-0.000337671660988197,-0.000336650130867220,-0.000334505681801338,-0.000331241670082228,-0.000326934577519601,-0.000321597411331245,-0.000315319515874973,-0.000308108733299146,-0.000300068275069311,-0.000291202506939356,-0.000281620445417402,-0.000271332401327684,-0.000260448698693465,-0.000248978615566231,-0.000237034781155883,-0.000224633417558677,-0.000211878668856251,-0.000198794093871137,-0.000185476850131979,-0.000171948811459607,-0.000158316330340026,-0.000144579619052287,-0.000130853687142025,-0.000117140359067039,-0.000103535649615391,-9.00437851570947e-05,-7.67661277548565e-05,-6.36905034503717e-05,-5.09117413208545e-05,-3.84220186694711e-05,-2.63058607512367e-05,-1.45519362022193e-05,-3.23421059797320e-06,7.66317268571757e-06,1.80639059740427e-05,2.79965806474123e-05,3.73910309205438e-05,4.62793390590956e-05,5.45955644700710e-05,6.23782262138511e-05,6.95679622187701e-05,7.62066042336095e-05,8.22521094557147e-05,8.77215523829069e-05,9.26226107885545e-05,9.69389345979157e-05,0.000100677354878430,0.000103892652709569,0.000106515908723708,0.000108623081333228,0.000110192908129537,0.000111282242197444,0.000111849884291550,0.000111970687838987,0.000111621547646505,0.000110875128450003,0.000109697674089670,0.000108156839530409,0.000106234517995692,0.000104006770866074,0.000101452080284138,9.86445413311909e-05,9.55562989757991e-05,9.22632831417983e-05,8.87432006630032e-05,8.50730286220366e-05,8.12279651616150e-05,7.72828817351650e-05,7.32098567809705e-05,6.90845510945624e-05,0.000532723527805931]; */
    i6 = b_num->size[0] * b_num->size[1];
    b_num->size[0] = 1;
    b_num->size[1] = 699;
    emxEnsureCapacity_real_T(b_num, i6);
    for (i6 = 0; i6 < 699; i6++) {
      b_num->data[b_num->size[0] * i6] = dv4[i6];
    }

    /* 'dynIdenf:567' den = 1; */
    /*  transient width 1.5 */
  } else if ((2.82 <= fpass) && (fpass < 5.65)) {
    /* 'dynIdenf:568' elseif fdisc(5) <= fpass && fpass < fdisc(6) */
    /* 'dynIdenf:569' num = [0.000541983718811927,9.39155210670181e-05,0.000100958414377355,0.000107595463015365,0.000113688228537891,0.000119177093506800,0.000123929134733962,0.000127894687581549,0.000130945415051525,0.000133044030326802,0.000134078026482925,0.000134025874269021,0.000132791129855200,0.000130373124118671,0.000126701511593129,0.000121799681582674,0.000115621965373439,0.000108220396625604,9.95750422745516e-05,8.97652851564876e-05,7.88062325959248e-05,6.68078245004558e-05,5.38098574168716e-05,3.99522311415092e-05,2.53040869264146e-05,1.00309029859570e-05,-5.77254316079002e-06,-2.19231126541713e-05,-3.83000510956379e-05,-5.47116603161720e-05,-7.10064533753939e-05,-8.70064897013418e-05,-0.000102529979431761,-0.000117407979761088,-0.000131481469000169,-0.000144534880486384,-0.000156470017010974,-0.000167084529934335,-0.000176261149179698,-0.000183837175084160,-0.000189747215029418,-0.000193841772251799,-0.000196073541936795,-0.000196346119718796,-0.000194665923081674,-0.000190965804054552,-0.000185281259751758,-0.000177600347515248,-0.000168020586519177,-0.000156570952883738,-0.000143383749251984,-0.000128534875844589,-0.000112203089934051,-9.45098568138851e-05,-7.56833010219895e-05,-5.58768211555122e-05,-3.53418834782659e-05,-1.42714666447916e-05,7.06086733384884e-06,2.84329633486768e-05,4.95590136185014e-05,7.02071138914894e-05,9.00914671901327e-05,0.000108987167847066,0.000126612317432528,0.000142760685517128,0.000157159733488689,0.000169645517325974,0.000179985578337734,0.000188051801111294,0.000193659121415496,0.000196723579657171,0.000197109088634444,0.000194802675707469,0.000189736391522837,0.000181945744112124,0.000171419148763144,0.000158270235612336,0.000142560624355696,0.000124455756549172,0.000104074080976090,8.16637642732029e-05,5.74156605995404e-05,3.16134719413095e-05,4.49760535647356e-06,-2.35688788073929e-05,-5.23127671566118e-05,-8.13640577031261e-05,-0.000110381212543967,-0.000139001076920169,-0.000166872755363170,-0.000193622797261251,-0.000218930579462653,-0.000242428006148879,-0.000263829521034822,-0.000282792065964986,-0.000299072886162210,-0.000312368664535981,-0.000322485228460860,-0.000329171873533593,-0.000332299375483446,-0.000331700761792967,-0.000327336417669472,-0.000319136839763799,-0.000307173721793938,-0.000291487872366930,-0.000272250469485773,-0.000249588877743172,-0.000223736461708317,-0.000194860168732485,-0.000163234645663197,-0.000129096464368348,-9.28427841098224e-05,-5.48945053903947e-05,-1.58536847988518e-05,2.37618534986138e-05,6.35611116325544e-05,0.000103589939605899,0.000143704940542049,0.000181025976607809,0.000217931456555117,0.000252357113836307,0.000284257828711435,0.000313182758504491,0.000338680260437589,0.000360413706777267,0.000377998077818869,0.000391181527305812,0.000399671906379765,0.000403315939413201,0.000401925594888096,0.000395460692994347,0.000383854155585259,0.000367187655733720,0.000345519399544020,0.000319057917826926,0.000287981876933402,0.000252620897756928,0.000213277183864515,0.000170392162781004,0.000124370213198492,7.57533418379797e-05,2.50361530035493e-05,-2.71611454259490e-05,-8.02768371814556e-05,-0.000133640898306514,-0.000186647951479871,-0.000238612044359073,-0.000288902070780961,-0.000336869074109680,-0.000381865684983965,-0.000423315587749875,-0.000460629099728361,-0.000493239085464758,-0.000520728588172605,-0.000542597316576951,-0.000558500829399108,-0.000568103825078294,-0.000571215822744036,-0.000567612524697722,-0.000557245516783168,-0.000540076800221215,-0.000516232780295362,-0.000485819662611941,-0.000449113416335854,-0.000406397565443762,-0.000358124919507632,-0.000304728818471515,-0.000246799554346799,-0.000184906233867251,-0.000119764102369112,-5.20614307239963e-05,1.73789449575496e-05,8.77957765650170e-05,0.000158312370733294,0.000228118320448514,0.000296318447014287,0.000362088969726770,0.000424554418003875,0.000482925129093528,0.000536391099262198,0.000584245077633135,0.000625772979312749,0.000660385014627503,0.000687485682784789,0.000706653899194518,0.000717462680602415,0.000719672891575588,0.000713058339562528,0.000697573224095843,0.000673195852659310,0.000640110647994494,0.000598517083474895,0.000548801415044521,0.000491374900037892,0.000426844595346898,0.000355818916285970,0.000279073495233852,0.000197387044022900,0.000111712026596431,2.29692786579564e-05,-6.77874598097355e-05,-0.000159529303918018,-0.000251102810739723,-0.000341447274286012,-0.000429399162350645,-0.000513872111887188,-0.000593763701806401,-0.000668037007475755,-0.000735672733369139,-0.000795757970971713,-0.000847393382183875,-0.000889830045262367,-0.000922339566245509,-0.000944383237710939,-0.000955454995225140,-0.000955264249053625,-0.000943564158373563,-0.000920339162179491,-0.000885626875170628,-0.000839683905906734,-0.000782812749022580,-0.000715536243565728,-0.000638420749418111,-0.000552260439891737,-0.000457903260615247,-0.000356409314647111,-0.000248851915276288,-0.000136464578717182,-2.04444307204847e-05,9.78643173315946e-05,0.000217112766350226,0.000335761607475920,0.000452337885564500,0.000565409273805852,0.000673758101284352,0.000775880149805064,0.000870008572946229,0.000955579555048026,0.00103084522155815,0.00109493959888634,0.00114684740229226,0.00118568928711056,0.00121079692202438,0.00122156944464255,0.00121766042770094,0.00119880438152662,0.00116500239871139,0.00111634393692054,0.00105319209384983,0.000976000891026337,0.000885492683406536,0.000782475751035737,0.000668013805650651,0.000543234859783321,0.000409507246153594,0.000268243354692785,0.000121059821395169,-3.04152432307783e-05,-0.000184379347411724,-0.000339059179704988,-0.000492548025635855,-0.000643010218370379,-0.000788522783484884,-0.000927273197598017,-0.00105741642025622,-0.00117724186294024,-0.00128509963403905,-0.00137943108525202,-0.00145889954833536,-0.00152219868611472,-0.00156826211329350,-0.00159627209149843,-0.00160546004472087,-0.00159540022296035,-0.00156583223311964,-0.00151678445534883,-0.00144841035655361,-0.00136121045709383,-0.00125585035904654,-0.00113330657899723,-0.000994673305073961,-0.000841350754803934,-0.000674867067745653,-0.000497031431954590,-0.000309729372712078,-0.000115084270895235,8.47444798042395e-05,0.000287402789594611,0.000490536127404496,0.000691646548084861,0.000888309185380098,0.00107801310643884,0.00125837311499210,0.00142698525992880,0.00158161297276232,0.00172006918497051,0.00184038233686247,0.00194070288808355,0.00201944312775646,0.00207518637052296,0.00210681645356649,0.00211341152631996,0.00209441885249221,0.00204948439910114,0.00197866066579206,0.00188222024557260,0.00176082763131060,0.00161537857650242,0.00144717165278660,0.00125770315050660,0.00104883948423773,0.000822638973055532,0.000581505818928723,0.000327975991355641,6.48740710262127e-05,-0.000204878269979111,-0.000478120672898087,-0.000751685750156150,-0.00102224662003798,-0.00128651260628978,-0.00154110992486032,-0.00178281151318399,-0.00200836474362844,-0.00221472015192204,-0.00239892427229273,-0.00255826248660913,-0.00269021370473722,-0.00279257189999333,-0.00286337500792373,-0.00290105581462172,-0.00290433180329040,-0.00287237823233321,-0.00280467899129824,-0.00270121258812355,-0.00256228356681253,-0.00238872096189341,-0.00218169844396906,-0.00194290595192604,-0.00167436411638681,-0.00137857899476873,-0.00105833699746803,-0.000716870711910641,-0.000357662081547121,1.54069234904593e-05,0.000398269941917141,0.000786554082636692,0.00117582302896910,0.00156143613107137,0.00193873592622316,0.00230289725078119,0.00264921375954365,0.00297305133436034,0.00327001605698777,0.00353565573027526,0.00376581900389305,0.00395697189829865,0.00410539360304297,0.00420817378241546,0.00426264813419011,0.00426664918670783,0.00421857702165025,0.00411729676153594,0.00396233701324010,0.00375373166757278,0.00349221769263537,0.00317906140829903,0.00281624784751191,0.00240629391400286,0.00195242310964760,0.00145836996294566,0.000928535446452230,0.000367776971615253,-0.000218437652524706,-0.000824257773321898,-0.00144331686525717,-0.00206896141715662,-0.00269412620596608,-0.00331158028138154,-0.00391380697172656,-0.00449325707533294,-0.00504223586171052,-0.00555315588129342,-0.00601843603020415,-0.00643073558805818,-0.00678291423619245,-0.00706812188302244,-0.00728002984015002,-0.00741252860038501,-0.00746025640048564,-0.00741838763488326,-0.00728259643511458,-0.00704943210024094,-0.00671611785605616,-0.00628069922505195,-0.00574193764251256,-0.00509956676718777,-0.00435411953700510,-0.00350708574003026,-0.00256073033453729,-0.00151829096644007,-0.000383799589256310,0.000837765270517561,0.00214069365337124,0.00351840339682411,0.00496366814771220,0.00646848017572530,0.00802427702998046,0.00962183372004456,0.0112515176906817,0.0129031792166922,0.0145664200115191,0.0162304954604694,0.0178845643813347,0.0195176178574509,0.0211187156283010,0.0226769516122276,0.0241816643667645,0.0256224202753077,0.0269891976615814,0.0282723656290463,0.0294628998851187,0.0305522889367809,0.0315328063885665,0.0323973751703697,0.0331398083189740,0.0337547063219542,0.0342376726977091,0.0345851566038690,0.0347946753205550,0.0348646677825939,0.0347946753205550,0.0345851566038690,0.0342376726977091,0.0337547063219542,0.0331398083189740,0.0323973751703697,0.0315328063885665,0.0305522889367809,0.0294628998851187,0.0282723656290463,0.0269891976615814,0.0256224202753077,0.0241816643667645,0.0226769516122276,0.0211187156283010,0.0195176178574509,0.0178845643813347,0.0162304954604694,0.0145664200115191,0.0129031792166922,0.0112515176906817,0.00962183372004456,0.00802427702998046,0.00646848017572530,0.00496366814771220,0.00351840339682411,0.00214069365337124,0.000837765270517561,-0.000383799589256310,-0.00151829096644007,-0.00256073033453729,-0.00350708574003026,-0.00435411953700510,-0.00509956676718777,-0.00574193764251256,-0.00628069922505195,-0.00671611785605616,-0.00704943210024094,-0.00728259643511458,-0.00741838763488326,-0.00746025640048564,-0.00741252860038501,-0.00728002984015002,-0.00706812188302244,-0.00678291423619245,-0.00643073558805818,-0.00601843603020415,-0.00555315588129342,-0.00504223586171052,-0.00449325707533294,-0.00391380697172656,-0.00331158028138154,-0.00269412620596608,-0.00206896141715662,-0.00144331686525717,-0.000824257773321898,-0.000218437652524706,0.000367776971615253,0.000928535446452230,0.00145836996294566,0.00195242310964760,0.00240629391400286,0.00281624784751191,0.00317906140829903,0.00349221769263537,0.00375373166757278,0.00396233701324010,0.00411729676153594,0.00421857702165025,0.00426664918670783,0.00426264813419011,0.00420817378241546,0.00410539360304297,0.00395697189829865,0.00376581900389305,0.00353565573027526,0.00327001605698777,0.00297305133436034,0.00264921375954365,0.00230289725078119,0.00193873592622316,0.00156143613107137,0.00117582302896910,0.000786554082636692,0.000398269941917141,1.54069234904593e-05,-0.000357662081547121,-0.000716870711910641,-0.00105833699746803,-0.00137857899476873,-0.00167436411638681,-0.00194290595192604,-0.00218169844396906,-0.00238872096189341,-0.00256228356681253,-0.00270121258812355,-0.00280467899129824,-0.00287237823233321,-0.00290433180329040,-0.00290105581462172,-0.00286337500792373,-0.00279257189999333,-0.00269021370473722,-0.00255826248660913,-0.00239892427229273,-0.00221472015192204,-0.00200836474362844,-0.00178281151318399,-0.00154110992486032,-0.00128651260628978,-0.00102224662003798,-0.000751685750156150,-0.000478120672898087,-0.000204878269979111,6.48740710262127e-05,0.000327975991355641,0.000581505818928723,0.000822638973055532,0.00104883948423773,0.00125770315050660,0.00144717165278660,0.00161537857650242,0.00176082763131060,0.00188222024557260,0.00197866066579206,0.00204948439910114,0.00209441885249221,0.00211341152631996,0.00210681645356649,0.00207518637052296,0.00201944312775646,0.00194070288808355,0.00184038233686247,0.00172006918497051,0.00158161297276232,0.00142698525992880,0.00125837311499210,0.00107801310643884,0.000888309185380098,0.000691646548084861,0.000490536127404496,0.000287402789594611,8.47444798042395e-05,-0.000115084270895235,-0.000309729372712078,-0.000497031431954590,-0.000674867067745653,-0.000841350754803934,-0.000994673305073961,-0.00113330657899723,-0.00125585035904654,-0.00136121045709383,-0.00144841035655361,-0.00151678445534883,-0.00156583223311964,-0.00159540022296035,-0.00160546004472087,-0.00159627209149843,-0.00156826211329350,-0.00152219868611472,-0.00145889954833536,-0.00137943108525202,-0.00128509963403905,-0.00117724186294024,-0.00105741642025622,-0.000927273197598017,-0.000788522783484884,-0.000643010218370379,-0.000492548025635855,-0.000339059179704988,-0.000184379347411724,-3.04152432307783e-05,0.000121059821395169,0.000268243354692785,0.000409507246153594,0.000543234859783321,0.000668013805650651,0.000782475751035737,0.000885492683406536,0.000976000891026337,0.00105319209384983,0.00111634393692054,0.00116500239871139,0.00119880438152662,0.00121766042770094,0.00122156944464255,0.00121079692202438,0.00118568928711056,0.00114684740229226,0.00109493959888634,0.00103084522155815,0.000955579555048026,0.000870008572946229,0.000775880149805064,0.000673758101284352,0.000565409273805852,0.000452337885564500,0.000335761607475920,0.000217112766350226,9.78643173315946e-05,-2.04444307204847e-05,-0.000136464578717182,-0.000248851915276288,-0.000356409314647111,-0.000457903260615247,-0.000552260439891737,-0.000638420749418111,-0.000715536243565728,-0.000782812749022580,-0.000839683905906734,-0.000885626875170628,-0.000920339162179491,-0.000943564158373563,-0.000955264249053625,-0.000955454995225140,-0.000944383237710939,-0.000922339566245509,-0.000889830045262367,-0.000847393382183875,-0.000795757970971713,-0.000735672733369139,-0.000668037007475755,-0.000593763701806401,-0.000513872111887188,-0.000429399162350645,-0.000341447274286012,-0.000251102810739723,-0.000159529303918018,-6.77874598097355e-05,2.29692786579564e-05,0.000111712026596431,0.000197387044022900,0.000279073495233852,0.000355818916285970,0.000426844595346898,0.000491374900037892,0.000548801415044521,0.000598517083474895,0.000640110647994494,0.000673195852659310,0.000697573224095843,0.000713058339562528,0.000719672891575588,0.000717462680602415,0.000706653899194518,0.000687485682784789,0.000660385014627503,0.000625772979312749,0.000584245077633135,0.000536391099262198,0.000482925129093528,0.000424554418003875,0.000362088969726770,0.000296318447014287,0.000228118320448514,0.000158312370733294,8.77957765650170e-05,1.73789449575496e-05,-5.20614307239963e-05,-0.000119764102369112,-0.000184906233867251,-0.000246799554346799,-0.000304728818471515,-0.000358124919507632,-0.000406397565443762,-0.000449113416335854,-0.000485819662611941,-0.000516232780295362,-0.000540076800221215,-0.000557245516783168,-0.000567612524697722,-0.000571215822744036,-0.000568103825078294,-0.000558500829399108,-0.000542597316576951,-0.000520728588172605,-0.000493239085464758,-0.000460629099728361,-0.000423315587749875,-0.000381865684983965,-0.000336869074109680,-0.000288902070780961,-0.000238612044359073,-0.000186647951479871,-0.000133640898306514,-8.02768371814556e-05,-2.71611454259490e-05,2.50361530035493e-05,7.57533418379797e-05,0.000124370213198492,0.000170392162781004,0.000213277183864515,0.000252620897756928,0.000287981876933402,0.000319057917826926,0.000345519399544020,0.000367187655733720,0.000383854155585259,0.000395460692994347,0.000401925594888096,0.000403315939413201,0.000399671906379765,0.000391181527305812,0.000377998077818869,0.000360413706777267,0.000338680260437589,0.000313182758504491,0.000284257828711435,0.000252357113836307,0.000217931456555117,0.000181025976607809,0.000143704940542049,0.000103589939605899,6.35611116325544e-05,2.37618534986138e-05,-1.58536847988518e-05,-5.48945053903947e-05,-9.28427841098224e-05,-0.000129096464368348,-0.000163234645663197,-0.000194860168732485,-0.000223736461708317,-0.000249588877743172,-0.000272250469485773,-0.000291487872366930,-0.000307173721793938,-0.000319136839763799,-0.000327336417669472,-0.000331700761792967,-0.000332299375483446,-0.000329171873533593,-0.000322485228460860,-0.000312368664535981,-0.000299072886162210,-0.000282792065964986,-0.000263829521034822,-0.000242428006148879,-0.000218930579462653,-0.000193622797261251,-0.000166872755363170,-0.000139001076920169,-0.000110381212543967,-8.13640577031261e-05,-5.23127671566118e-05,-2.35688788073929e-05,4.49760535647356e-06,3.16134719413095e-05,5.74156605995404e-05,8.16637642732029e-05,0.000104074080976090,0.000124455756549172,0.000142560624355696,0.000158270235612336,0.000171419148763144,0.000181945744112124,0.000189736391522837,0.000194802675707469,0.000197109088634444,0.000196723579657171,0.000193659121415496,0.000188051801111294,0.000179985578337734,0.000169645517325974,0.000157159733488689,0.000142760685517128,0.000126612317432528,0.000108987167847066,9.00914671901327e-05,7.02071138914894e-05,4.95590136185014e-05,2.84329633486768e-05,7.06086733384884e-06,-1.42714666447916e-05,-3.53418834782659e-05,-5.58768211555122e-05,-7.56833010219895e-05,-9.45098568138851e-05,-0.000112203089934051,-0.000128534875844589,-0.000143383749251984,-0.000156570952883738,-0.000168020586519177,-0.000177600347515248,-0.000185281259751758,-0.000190965804054552,-0.000194665923081674,-0.000196346119718796,-0.000196073541936795,-0.000193841772251799,-0.000189747215029418,-0.000183837175084160,-0.000176261149179698,-0.000167084529934335,-0.000156470017010974,-0.000144534880486384,-0.000131481469000169,-0.000117407979761088,-0.000102529979431761,-8.70064897013418e-05,-7.10064533753939e-05,-5.47116603161720e-05,-3.83000510956379e-05,-2.19231126541713e-05,-5.77254316079002e-06,1.00309029859570e-05,2.53040869264146e-05,3.99522311415092e-05,5.38098574168716e-05,6.68078245004558e-05,7.88062325959248e-05,8.97652851564876e-05,9.95750422745516e-05,0.000108220396625604,0.000115621965373439,0.000121799681582674,0.000126701511593129,0.000130373124118671,0.000132791129855200,0.000134025874269021,0.000134078026482925,0.000133044030326802,0.000130945415051525,0.000127894687581549,0.000123929134733962,0.000119177093506800,0.000113688228537891,0.000107595463015365,0.000100958414377355,9.39155210670181e-05,0.000541983718811927]; */
    i6 = b_num->size[0] * b_num->size[1];
    b_num->size[0] = 1;
    b_num->size[1] = 861;
    emxEnsureCapacity_real_T(b_num, i6);
    for (i6 = 0; i6 < 861; i6++) {
      b_num->data[b_num->size[0] * i6] = dv5[i6];
    }

    /* 'dynIdenf:570' den = 1; */
    /*  transient width 1.2 */
  } else if ((5.65 <= fpass) && (fpass < 11.31)) {
    /* 'dynIdenf:571' elseif fdisc(6) <= fpass && fpass < fdisc(7) */
    /* 'dynIdenf:572' num = [0.000591089306393010,0.000203135102293136,0.000228702342098083,0.000249343617535587,0.000263506976273500,0.000269934205063887,0.000267468276220013,0.000255402479651796,0.000233273621757195,0.000201219812933285,0.000159682886147079,0.000109724050959972,5.27085192408347e-05,-9.40957449163658e-06,-7.44675389209232e-05,-0.000139841547166302,-0.000202885801604255,-0.000260742843701583,-0.000310787524703063,-0.000350412679782949,-0.000377500244119596,-0.000390236144409474,-0.000387546998560333,-0.000368850380778751,-0.000334438591051295,-0.000285181900113470,-0.000222869394888480,-0.000149852306955432,-6.92610818611087e-05,1.54578194554596e-05,0.000100328930107147,0.000181312682430688,0.000254228552368588,0.000315289227791654,0.000360922754044401,0.000388245755760958,0.000395072723587488,0.000380257274688715,0.000343456936944303,0.000285548262117333,0.000208452985533168,0.000115049129259730,9.22773619128056e-06,-0.000104352698013944,-0.000220482509693878,-0.000333510603247651,-0.000437774350204889,-0.000527652824120159,-0.000598146628842731,-0.000644809129631383,-0.000664324261454271,-0.000654418809838929,-0.000614326456707724,-0.000544556021427374,-0.000447273028712333,-0.000325932862966167,-0.000185480230425671,-3.18423145588313e-05,0.000127976430894124,0.000286536638081843,0.000435981691500021,0.000568730119468750,0.000677529452718088,0.000756169921240225,0.000799491385412860,0.000804009174058410,0.000767835944931364,0.000691169259725330,0.000576024658008797,0.000426553487272439,0.000248616157020845,4.98649714080354e-05,-0.000160884175849547,-0.000373698214029173,-0.000578284175225528,-0.000764238148896663,-0.000921730736103272,-0.00104195550586188,-0.00111743946208126,-0.00114289576019665,-0.00111495567630999,-0.00103287630883038,-0.000898679009764979,-0.000716715885260000,-0.000494043479284734,-0.000239947494484822,3.43194002594898e-05,0.000316273875969929,0.000592355522177910,0.000848893105460828,0.00107254336564001,0.00125127031009951,0.00137460646711804,0.00143453992511530,0.00142569620683721,0.00134603262510240,0.00119672309564059,0.000982528974077999,0.000711371827311320,0.000394431615214206,4.54340164485031e-05,-0.000319611898991610,-0.000683490853806374,-0.00102828762426528,-0.00133657646482561,-0.00159194769783486,-0.00178011022049491,-0.00188927817714030,-0.00191116859766644,-0.00184130095806268,-0.00167978390347217,-0.00143121554544908,-0.00110479063604924,-0.000713421120918264,-0.000273458733759465,0.000195673531830978,0.000671366376562233,0.00113021460921350,0.00155149140592462,0.00191065109022734,0.00218952949178247,0.00237100500122057,0.00244271393319439,0.00239715681481954,0.00223215536250292,0.00195142503243853,0.00156428896227672,0.00108576770776838,0.000535720184246959,-6.16152130609063e-05,-0.000678958280353827,-0.00128686237729246,-0.00185542519534156,-0.00235536537880196,-0.00275979008980575,-0.00304525170506758,-0.00319337838493626,-0.00319160458969291,-0.00303434112247834,-0.00272315976994321,-0.00226733239506472,-0.00168333400111216,-0.000994682377620397,-0.000230714554582418,0.000574271356711954,0.00138280136505085,0.00215554933784529,0.00285353542133235,0.00343968761494919,0.00388099477089373,0.00414992270721462,0.00422634082660383,0.00409845181704085,0.00376392649007728,0.00323021244990824,0.00251482837201619,0.00164466215588651,0.000655298804813174,-0.000410399793453483,-0.00150401349518631,-0.00257368112952168,-0.00356627803023731,-0.00443007902205741,-0.00511708466630530,-0.00558567080162908,-0.00580256924799627,-0.00574517437109170,-0.00540273308381670,-0.00477772903587716,-0.00388608741957690,-0.00275741920859663,-0.00143393039765822,3.05684537051032e-05,0.00157291354971703,0.00312270604093914,0.00460555918146983,0.00594595608074391,0.00707100558031540,0.00791357328678103,0.00841599967138872,0.00853287838122850,0.00823414920527687,0.00750696640123341,0.00635761364944972,0.00481202042811583,0.00291616685206864,0.000734945152006006,-0.00164911359246092,-0.00413854912587785,-0.00662376584037583,-0.00898713211466898,-0.0111068908850509,-0.0128620309563478,-0.0141368054464089,-0.0148254937011961,-0.0148372599010804,-0.0140992200430354,-0.0125617908181220,-0.0101995897012186,-0.00701448013087040,-0.00303685260159460,0.00167519934115319,0.00703647678686433,0.0129367099188529,0.0192433657583256,0.0258060198556134,0.0324606280440578,0.0390349134111635,0.0453535974623025,0.0512445768528114,0.0565444150513016,0.0611042135615666,0.0647943156695586,0.0675090341357033,0.0691699572095659,0.0697290376973575,0.0691699572095659,0.0675090341357033,0.0647943156695586,0.0611042135615666,0.0565444150513016,0.0512445768528114,0.0453535974623025,0.0390349134111635,0.0324606280440578,0.0258060198556134,0.0192433657583256,0.0129367099188529,0.00703647678686433,0.00167519934115319,-0.00303685260159460,-0.00701448013087040,-0.0101995897012186,-0.0125617908181220,-0.0140992200430354,-0.0148372599010804,-0.0148254937011961,-0.0141368054464089,-0.0128620309563478,-0.0111068908850509,-0.00898713211466898,-0.00662376584037583,-0.00413854912587785,-0.00164911359246092,0.000734945152006006,0.00291616685206864,0.00481202042811583,0.00635761364944972,0.00750696640123341,0.00823414920527687,0.00853287838122850,0.00841599967138872,0.00791357328678103,0.00707100558031540,0.00594595608074391,0.00460555918146983,0.00312270604093914,0.00157291354971703,3.05684537051032e-05,-0.00143393039765822,-0.00275741920859663,-0.00388608741957690,-0.00477772903587716,-0.00540273308381670,-0.00574517437109170,-0.00580256924799627,-0.00558567080162908,-0.00511708466630530,-0.00443007902205741,-0.00356627803023731,-0.00257368112952168,-0.00150401349518631,-0.000410399793453483,0.000655298804813174,0.00164466215588651,0.00251482837201619,0.00323021244990824,0.00376392649007728,0.00409845181704085,0.00422634082660383,0.00414992270721462,0.00388099477089373,0.00343968761494919,0.00285353542133235,0.00215554933784529,0.00138280136505085,0.000574271356711954,-0.000230714554582418,-0.000994682377620397,-0.00168333400111216,-0.00226733239506472,-0.00272315976994321,-0.00303434112247834,-0.00319160458969291,-0.00319337838493626,-0.00304525170506758,-0.00275979008980575,-0.00235536537880196,-0.00185542519534156,-0.00128686237729246,-0.000678958280353827,-6.16152130609063e-05,0.000535720184246959,0.00108576770776838,0.00156428896227672,0.00195142503243853,0.00223215536250292,0.00239715681481954,0.00244271393319439,0.00237100500122057,0.00218952949178247,0.00191065109022734,0.00155149140592462,0.00113021460921350,0.000671366376562233,0.000195673531830978,-0.000273458733759465,-0.000713421120918264,-0.00110479063604924,-0.00143121554544908,-0.00167978390347217,-0.00184130095806268,-0.00191116859766644,-0.00188927817714030,-0.00178011022049491,-0.00159194769783486,-0.00133657646482561,-0.00102828762426528,-0.000683490853806374,-0.000319611898991610,4.54340164485031e-05,0.000394431615214206,0.000711371827311320,0.000982528974077999,0.00119672309564059,0.00134603262510240,0.00142569620683721,0.00143453992511530,0.00137460646711804,0.00125127031009951,0.00107254336564001,0.000848893105460828,0.000592355522177910,0.000316273875969929,3.43194002594898e-05,-0.000239947494484822,-0.000494043479284734,-0.000716715885260000,-0.000898679009764979,-0.00103287630883038,-0.00111495567630999,-0.00114289576019665,-0.00111743946208126,-0.00104195550586188,-0.000921730736103272,-0.000764238148896663,-0.000578284175225528,-0.000373698214029173,-0.000160884175849547,4.98649714080354e-05,0.000248616157020845,0.000426553487272439,0.000576024658008797,0.000691169259725330,0.000767835944931364,0.000804009174058410,0.000799491385412860,0.000756169921240225,0.000677529452718088,0.000568730119468750,0.000435981691500021,0.000286536638081843,0.000127976430894124,-3.18423145588313e-05,-0.000185480230425671,-0.000325932862966167,-0.000447273028712333,-0.000544556021427374,-0.000614326456707724,-0.000654418809838929,-0.000664324261454271,-0.000644809129631383,-0.000598146628842731,-0.000527652824120159,-0.000437774350204889,-0.000333510603247651,-0.000220482509693878,-0.000104352698013944,9.22773619128056e-06,0.000115049129259730,0.000208452985533168,0.000285548262117333,0.000343456936944303,0.000380257274688715,0.000395072723587488,0.000388245755760958,0.000360922754044401,0.000315289227791654,0.000254228552368588,0.000181312682430688,0.000100328930107147,1.54578194554596e-05,-6.92610818611087e-05,-0.000149852306955432,-0.000222869394888480,-0.000285181900113470,-0.000334438591051295,-0.000368850380778751,-0.000387546998560333,-0.000390236144409474,-0.000377500244119596,-0.000350412679782949,-0.000310787524703063,-0.000260742843701583,-0.000202885801604255,-0.000139841547166302,-7.44675389209232e-05,-9.40957449163658e-06,5.27085192408347e-05,0.000109724050959972,0.000159682886147079,0.000201219812933285,0.000233273621757195,0.000255402479651796,0.000267468276220013,0.000269934205063887,0.000263506976273500,0.000249343617535587,0.000228702342098083,0.000203135102293136,0.000591089306393010]; */
    i6 = b_num->size[0] * b_num->size[1];
    b_num->size[0] = 1;
    b_num->size[1] = 431;
    emxEnsureCapacity_real_T(b_num, i6);
    for (i6 = 0; i6 < 431; i6++) {
      b_num->data[b_num->size[0] * i6] = dv6[i6];
    }

    /* 'dynIdenf:573' den = 1; */
    /*  transient width 1.2 */
  } else if ((11.31 <= fpass) && (fpass < 22.62)) {
    /* 'dynIdenf:574' elseif fdisc(7) <= fpass && fpass < fdisc(8) */
    /* 'dynIdenf:575' num = [0.000703874537800322,0.000467770792417846,0.000539242180785141,0.000548714935447943,0.000481522011841131,0.000335028055626118,0.000121223276456724,-0.000133359726248962,-0.000390715672681121,-0.000607205143402996,-0.000741333556204524,-0.000762141773286633,-0.000656572560307136,-0.000434098146732689,-0.000127505132900021,0.000211213864529279,0.000518476440093037,0.000731495160441613,0.000799153119402637,0.000695246620796837,0.000424316615262130,2.47090795265547e-05,-0.000435793659924257,-0.000871359045494241,-0.00119303737495923,-0.00132606332118820,-0.00122651508944986,-0.000892783882200666,-0.000369395944141702,0.000257559799410267,0.000873723084736679,0.00135682421064802,0.00160057389483701,0.00153681996160046,0.00115260727769925,0.000497054237143859,-0.000322657933790586,-0.00115812613514566,-0.00184572049236201,-0.00223771577754856,-0.00223288551008581,-0.00180018865702504,-0.000990574744246535,6.65230524788864e-05,0.00118273010625971,0.00214338417588515,0.00274778393008912,0.00284990682126778,0.00239141435582874,0.00141991192268376,8.75329572949450e-05,-0.00137083732615000,-0.00267793941285247,-0.00356533814117648,-0.00382699219718207,-0.00336439837467214,-0.00221403160456867,-0.000550206509922040,0.00133894443554213,0.00309954977047986,0.00437648940085857,0.00488292289402572,0.00446165445896134,0.00312562937150274,0.00106802814472366,-0.00136197880230862,-0.00371559282089712,-0.00552480925577379,-0.00639220079314186,-0.00607420576815403,-0.00454007720411235,-0.00199452874172680,0.00114373496967410,0.00430664656510399,0.00687535360020624,0.00829614461523189,0.00819326036702013,0.00645661505509234,0.00328519219444469,-0.000825238131840536,-0.00515198261913290,-0.00886492643018391,-0.0111763865229752,-0.0114955598573719,-0.00956066216535001,-0.00551970592973219,5.66013232584069e-05,0.00624100851967568,0.0118874959665648,0.0158228554581239,0.0170617328212118,0.0150099407000975,0.00961969466734620,0.00146523736031698,-0.00828158313067491,-0.0179786968695017,-0.0257289781554320,-0.0296559727194340,-0.0282029907072548,-0.0204040677918449,-0.00607822692183737,0.0140682470909551,0.0384820357992972,0.0649166378856776,0.0907022576471758,0.113083794771166,0.129583642532735,0.138334705982925,0.138334705982925,0.129583642532735,0.113083794771166,0.0907022576471758,0.0649166378856776,0.0384820357992972,0.0140682470909551,-0.00607822692183737,-0.0204040677918449,-0.0282029907072548,-0.0296559727194340,-0.0257289781554320,-0.0179786968695017,-0.00828158313067491,0.00146523736031698,0.00961969466734620,0.0150099407000975,0.0170617328212118,0.0158228554581239,0.0118874959665648,0.00624100851967568,5.66013232584069e-05,-0.00551970592973219,-0.00956066216535001,-0.0114955598573719,-0.0111763865229752,-0.00886492643018391,-0.00515198261913290,-0.000825238131840536,0.00328519219444469,0.00645661505509234,0.00819326036702013,0.00829614461523189,0.00687535360020624,0.00430664656510399,0.00114373496967410,-0.00199452874172680,-0.00454007720411235,-0.00607420576815403,-0.00639220079314186,-0.00552480925577379,-0.00371559282089712,-0.00136197880230862,0.00106802814472366,0.00312562937150274,0.00446165445896134,0.00488292289402572,0.00437648940085857,0.00309954977047986,0.00133894443554213,-0.000550206509922040,-0.00221403160456867,-0.00336439837467214,-0.00382699219718207,-0.00356533814117648,-0.00267793941285247,-0.00137083732615000,8.75329572949450e-05,0.00141991192268376,0.00239141435582874,0.00284990682126778,0.00274778393008912,0.00214338417588515,0.00118273010625971,6.65230524788864e-05,-0.000990574744246535,-0.00180018865702504,-0.00223288551008581,-0.00223771577754856,-0.00184572049236201,-0.00115812613514566,-0.000322657933790586,0.000497054237143859,0.00115260727769925,0.00153681996160046,0.00160057389483701,0.00135682421064802,0.000873723084736679,0.000257559799410267,-0.000369395944141702,-0.000892783882200666,-0.00122651508944986,-0.00132606332118820,-0.00119303737495923,-0.000871359045494241,-0.000435793659924257,2.47090795265547e-05,0.000424316615262130,0.000695246620796837,0.000799153119402637,0.000731495160441613,0.000518476440093037,0.000211213864529279,-0.000127505132900021,-0.000434098146732689,-0.000656572560307136,-0.000762141773286633,-0.000741333556204524,-0.000607205143402996,-0.000390715672681121,-0.000133359726248962,0.000121223276456724,0.000335028055626118,0.000481522011841131,0.000548714935447943,0.000539242180785141,0.000467770792417846,0.000703874537800322]; */
    i6 = b_num->size[0] * b_num->size[1];
    b_num->size[0] = 1;
    b_num->size[1] = 216;
    emxEnsureCapacity_real_T(b_num, i6);
    for (i6 = 0; i6 < 216; i6++) {
      b_num->data[b_num->size[0] * i6] = dv7[i6];
    }

    /* 'dynIdenf:576' den = 1; */
    /*  transient width 1.2 */
  } else if ((22.62 <= fpass) && (fpass < 45.25)) {
    /* 'dynIdenf:577' elseif fdisc(8) <= fpass && fpass < fdisc(9) */
    /* 'dynIdenf:578' num = [0.000701399010505143,0.00140839928277612,0.00111499214990723,0.000816359876385092,-0.000344933864345112,-0.00108141865642759,-0.00124645997896249,-0.000328553861170572,0.000953439235176911,0.00179667825296385,0.00135345312991716,-0.000219458433219591,-0.00194482183672888,-0.00246303978496230,-0.00115724040777242,0.00129882662458625,0.00317855269358360,0.00291164222037781,0.000326935175146764,-0.00295096520726589,-0.00448154193949628,-0.00277888085477348,0.00136069032614137,0.00510846062072828,0.00548563343107216,0.00166478513201453,-0.00404551341994001,-0.00750445686013723,-0.00567597805715282,0.000839555412816575,0.00771068709028747,0.00965911058655907,0.00440494762676344,-0.00510107364052596,-0.0121479368907136,-0.0108638160387891,-0.000888239685836671,0.0114580689600834,0.0169540627855660,0.0101044811810528,-0.00598125671977792,-0.0204435882733855,-0.0215909420459943,-0.00571315147256970,0.0184513279393411,0.0336982429435121,0.0254652063742116,-0.00656621578009933,-0.0444036515025244,-0.0593295389161810,-0.0280398177296805,0.0517667833144886,0.156163072823834,0.244442499000399,0.278942485275155,0.244442499000399,0.156163072823834,0.0517667833144886,-0.0280398177296805,-0.0593295389161810,-0.0444036515025244,-0.00656621578009933,0.0254652063742116,0.0336982429435121,0.0184513279393411,-0.00571315147256970,-0.0215909420459943,-0.0204435882733855,-0.00598125671977792,0.0101044811810528,0.0169540627855660,0.0114580689600834,-0.000888239685836671,-0.0108638160387891,-0.0121479368907136,-0.00510107364052596,0.00440494762676344,0.00965911058655907,0.00771068709028747,0.000839555412816575,-0.00567597805715282,-0.00750445686013723,-0.00404551341994001,0.00166478513201453,0.00548563343107216,0.00510846062072828,0.00136069032614137,-0.00277888085477348,-0.00448154193949628,-0.00295096520726589,0.000326935175146764,0.00291164222037781,0.00317855269358360,0.00129882662458625,-0.00115724040777242,-0.00246303978496230,-0.00194482183672888,-0.000219458433219591,0.00135345312991716,0.00179667825296385,0.000953439235176911,-0.000328553861170572,-0.00124645997896249,-0.00108141865642759,-0.000344933864345112,0.000816359876385092,0.00111499214990723,0.00140839928277612,0.000701399010505143]; */
    i6 = b_num->size[0] * b_num->size[1];
    b_num->size[0] = 1;
    b_num->size[1] = 109;
    emxEnsureCapacity_real_T(b_num, i6);
    for (i6 = 0; i6 < 109; i6++) {
      b_num->data[b_num->size[0] * i6] = dv8[i6];
    }

    /* 'dynIdenf:579' den = 1; */
    /*  transient width 1.2 */
  } else if ((45.25 <= fpass) && (fpass < 80.0)) {
    /* 'dynIdenf:580' elseif fdisc(9) <= fpass && fpass < fdisc(10) */
    /* 'dynIdenf:581' num = [0.00136160164139038,0.00267369398929317,-0.000697117089159221,-0.00234837763528271,0.00191581346362666,0.00269900377250914,-0.00389621556543998,-0.00225279038429998,0.00633783534093455,0.000605005371784838,-0.00893041710195627,0.00276678504202593,0.0109185935758989,-0.00812278495383672,-0.0112857210024778,0.0154360766314228,0.00873812519907081,-0.0242867022296137,-0.00170259008266850,0.0338768389503649,-0.0120304430378272,-0.0431322973914692,0.0369550170539418,0.0508632874364068,-0.0888408228750916,-0.0560010731682915,0.312337552871639,0.557802460840468,0.312337552871639,-0.0560010731682915,-0.0888408228750916,0.0508632874364068,0.0369550170539418,-0.0431322973914692,-0.0120304430378272,0.0338768389503649,-0.00170259008266850,-0.0242867022296137,0.00873812519907081,0.0154360766314228,-0.0112857210024778,-0.00812278495383672,0.0109185935758989,0.00276678504202593,-0.00893041710195627,0.000605005371784838,0.00633783534093455,-0.00225279038429998,-0.00389621556543998,0.00269900377250914,0.00191581346362666,-0.00234837763528271,-0.000697117089159221,0.00267369398929317,0.00136160164139038]; */
    i6 = b_num->size[0] * b_num->size[1];
    b_num->size[0] = 1;
    b_num->size[1] = 55;
    emxEnsureCapacity_real_T(b_num, i6);
    for (i6 = 0; i6 < 55; i6++) {
      b_num->data[b_num->size[0] * i6] = dv9[i6];
    }

    /* 'dynIdenf:582' den = 1; */
    /*  transient width 1.2 */
  } else if ((80.0 <= fpass) && (fpass < 100.0)) {
    /* 'dynIdenf:583' elseif fdisc(10) <= fpass && fpass < fdisc(11) */
    /* 'dynIdenf:584' num = [0.00396227519738595,-0.00432580986830259,0.00537645334783375,-0.00525174059537846,0.00329357485741764,0.000923114411272612,-0.00738284264798429,0.0154672131632766,-0.0238706513909533,0.0306418755845219,-0.0333198455184126,0.0290689001544471,-0.0146804244699826,-0.0141610568835452,0.0662458874783531,-0.172852103294901,0.623027351924843,0.623027351924843,-0.172852103294901,0.0662458874783531,-0.0141610568835452,-0.0146804244699826,0.0290689001544471,-0.0333198455184126,0.0306418755845219,-0.0238706513909533,0.0154672131632766,-0.00738284264798429,0.000923114411272612,0.00329357485741764,-0.00525174059537846,0.00537645334783375,-0.00432580986830259,0.00396227519738595]; */
    i6 = b_num->size[0] * b_num->size[1];
    b_num->size[0] = 1;
    b_num->size[1] = 34;
    emxEnsureCapacity_real_T(b_num, i6);
    for (i6 = 0; i6 < 34; i6++) {
      b_num->data[b_num->size[0] * i6] = dv10[i6];
    }

    /* 'dynIdenf:585' den = 1; */
    /*  transient width 1.2 */
  } else {
    /* 'dynIdenf:586' else */
    /* 'dynIdenf:587' num = [1 0]; */
    i6 = b_num->size[0] * b_num->size[1];
    b_num->size[0] = 1;
    b_num->size[1] = 2;
    emxEnsureCapacity_real_T(b_num, i6);
    for (i6 = 0; i6 < 2; i6++) {
      b_num->data[b_num->size[0] * i6] = 1.0 - (double)i6;
    }

    /* 'dynIdenf:588' den = 1; */
  }

  *b_den = 1.0;
}

/*
 * Arguments    : emxArray_int32_T *x
 *                int idx
 *                int xstart
 *                int xend
 *                const cell_wrap_4 cmp_tunableEnvironment[2]
 * Return Type  : void
 */
static void heapify(emxArray_int32_T *x, int idx, int xstart, int xend, const
                    cell_wrap_4 cmp_tunableEnvironment[2])
{
  bool changed;
  int extremumIdx;
  int leftIdx;
  int extremum;
  bool varargout_1;
  int cmpIdx;
  int xcmp;
  changed = true;
  extremumIdx = (idx + xstart) - 2;
  leftIdx = ((idx << 1) + xstart) - 2;
  while (changed && (leftIdx + 1 < xend)) {
    changed = false;
    extremum = x->data[extremumIdx];
    cmpIdx = leftIdx;
    xcmp = x->data[leftIdx];
    varargout_1 = ((cmp_tunableEnvironment[0].f1->data[x->data[leftIdx] - 1] <
                    cmp_tunableEnvironment[0].f1->data[x->data[leftIdx + 1] - 1])
                   || ((cmp_tunableEnvironment[0].f1->data[x->data[leftIdx] - 1]
                        == cmp_tunableEnvironment[0].f1->data[x->data[leftIdx +
                        1] - 1]) && (cmp_tunableEnvironment[1].f1->data[x->
      data[leftIdx] - 1] < cmp_tunableEnvironment[1].f1->data[x->data[leftIdx +
      1] - 1])));
    if (varargout_1) {
      cmpIdx = leftIdx + 1;
      xcmp = x->data[leftIdx + 1];
    }

    varargout_1 = ((cmp_tunableEnvironment[0].f1->data[x->data[extremumIdx] - 1]
                    < cmp_tunableEnvironment[0].f1->data[xcmp - 1]) ||
                   ((cmp_tunableEnvironment[0].f1->data[x->data[extremumIdx] - 1]
                     == cmp_tunableEnvironment[0].f1->data[xcmp - 1]) &&
                    (cmp_tunableEnvironment[1].f1->data[x->data[extremumIdx] - 1]
                     < cmp_tunableEnvironment[1].f1->data[xcmp - 1])));
    if (varargout_1) {
      x->data[extremumIdx] = xcmp;
      x->data[cmpIdx] = extremum;
      extremumIdx = cmpIdx;
      leftIdx = ((((cmpIdx - xstart) + 2) << 1) + xstart) - 2;
      changed = true;
    }
  }

  if (changed && (leftIdx + 1 <= xend)) {
    extremum = x->data[extremumIdx];
    varargout_1 = ((cmp_tunableEnvironment[0].f1->data[x->data[extremumIdx] - 1]
                    < cmp_tunableEnvironment[0].f1->data[x->data[leftIdx] - 1]) ||
                   ((cmp_tunableEnvironment[0].f1->data[x->data[extremumIdx] - 1]
                     == cmp_tunableEnvironment[0].f1->data[x->data[leftIdx] - 1])
                    && (cmp_tunableEnvironment[1].f1->data[x->data[extremumIdx]
                        - 1] < cmp_tunableEnvironment[1].f1->data[x->
                        data[leftIdx] - 1])));
    if (varargout_1) {
      x->data[extremumIdx] = x->data[leftIdx];
      x->data[leftIdx] = extremum;
    }
  }
}

/*
 * Arguments    : emxArray_int32_T *x
 *                int xstart
 *                int xend
 *                const cell_wrap_4 cmp_tunableEnvironment[2]
 * Return Type  : void
 */
static void insertionsort(emxArray_int32_T *x, int xstart, int xend, const
  cell_wrap_4 cmp_tunableEnvironment[2])
{
  int k;
  int xc;
  int idx;
  bool exitg1;
  bool varargout_1;
  for (k = xstart; k < xend; k++) {
    xc = x->data[k] - 1;
    idx = k;
    exitg1 = false;
    while ((!exitg1) && (idx >= xstart)) {
      varargout_1 = ((cmp_tunableEnvironment[0].f1->data[xc] <
                      cmp_tunableEnvironment[0].f1->data[x->data[idx - 1] - 1]) ||
                     ((cmp_tunableEnvironment[0].f1->data[xc] ==
                       cmp_tunableEnvironment[0].f1->data[x->data[idx - 1] - 1])
                      && (cmp_tunableEnvironment[1].f1->data[xc] <
                          cmp_tunableEnvironment[1].f1->data[x->data[idx - 1] -
                          1])));
      if (varargout_1) {
        x->data[idx] = x->data[idx - 1];
        idx--;
      } else {
        exitg1 = true;
      }
    }

    x->data[idx] = xc + 1;
  }
}

/*
 * Arguments    : emxArray_int32_T *x
 *                int xend
 *                const cell_wrap_4 cmp_tunableEnvironment[2]
 * Return Type  : void
 */
static void introsort(emxArray_int32_T *x, int xend, const cell_wrap_4
                      cmp_tunableEnvironment[2])
{
  int pmax;
  int pmin;
  bool exitg1;
  emxArray_struct_T *st_d;
  int p;
  int MAXDEPTH;
  int pow2p;
  struct_T frame;
  signed char unnamed_idx_0;
  bool varargout_1;
  int pivot;
  int exitg2;
  int exitg3;
  if (!(1 >= xend)) {
    if (xend <= 32) {
      insertionsort(x, 1, xend, cmp_tunableEnvironment);
    } else {
      pmax = 31;
      pmin = 0;
      exitg1 = false;
      while ((!exitg1) && (pmax - pmin > 1)) {
        p = (pmin + pmax) >> 1;
        pow2p = 1 << p;
        if (pow2p == xend) {
          pmax = p;
          exitg1 = true;
        } else if (pow2p > xend) {
          pmax = p;
        } else {
          pmin = p;
        }
      }

      emxInit_struct_T(&st_d, 1);
      MAXDEPTH = (pmax - 1) << 1;
      frame.xstart = 1;
      frame.xend = xend;
      frame.depth = 0;
      unnamed_idx_0 = (signed char)(MAXDEPTH << 1);
      pmax = st_d->size[0];
      st_d->size[0] = unnamed_idx_0;
      emxEnsureCapacity_struct_T(st_d, pmax);
      pmin = unnamed_idx_0;
      for (pmax = 0; pmax < pmin; pmax++) {
        st_d->data[pmax] = frame;
      }

      st_d->data[0] = frame;
      p = 1;
      while (p > 0) {
        frame = st_d->data[p - 1];
        p--;
        if ((frame.xend - frame.xstart) + 1 <= 32) {
          insertionsort(x, frame.xstart, frame.xend, cmp_tunableEnvironment);
        } else if (frame.depth == MAXDEPTH) {
          b_heapsort(x, frame.xstart, frame.xend, cmp_tunableEnvironment);
        } else {
          pmax = (frame.xstart + (frame.xend - frame.xstart) / 2) - 1;
          varargout_1 = ((cmp_tunableEnvironment[0].f1->data[x->data[pmax] - 1] <
                          cmp_tunableEnvironment[0].f1->data[x->
                          data[frame.xstart - 1] - 1]) ||
                         ((cmp_tunableEnvironment[0].f1->data[x->data[pmax] - 1]
                           == cmp_tunableEnvironment[0].f1->data[x->
                           data[frame.xstart - 1] - 1]) &&
                          (cmp_tunableEnvironment[1].f1->data[x->data[pmax] - 1]
                           < cmp_tunableEnvironment[1].f1->data[x->
                           data[frame.xstart - 1] - 1])));
          if (varargout_1) {
            pow2p = x->data[frame.xstart - 1];
            x->data[frame.xstart - 1] = x->data[pmax];
            x->data[pmax] = pow2p;
          }

          varargout_1 = ((cmp_tunableEnvironment[0].f1->data[x->data[frame.xend
                          - 1] - 1] < cmp_tunableEnvironment[0].f1->data[x->
                          data[frame.xstart - 1] - 1]) ||
                         ((cmp_tunableEnvironment[0].f1->data[x->data[frame.xend
                           - 1] - 1] == cmp_tunableEnvironment[0].f1->data
                           [x->data[frame.xstart - 1] - 1]) &&
                          (cmp_tunableEnvironment[1].f1->data[x->data[frame.xend
                           - 1] - 1] < cmp_tunableEnvironment[1].f1->data
                           [x->data[frame.xstart - 1] - 1])));
          if (varargout_1) {
            pow2p = x->data[frame.xstart - 1];
            x->data[frame.xstart - 1] = x->data[frame.xend - 1];
            x->data[frame.xend - 1] = pow2p;
          }

          varargout_1 = ((cmp_tunableEnvironment[0].f1->data[x->data[frame.xend
                          - 1] - 1] < cmp_tunableEnvironment[0].f1->data[x->
                          data[pmax] - 1]) || ((cmp_tunableEnvironment[0]
            .f1->data[x->data[frame.xend - 1] - 1] == cmp_tunableEnvironment[0].
            f1->data[x->data[pmax] - 1]) && (cmp_tunableEnvironment[1].f1->
            data[x->data[frame.xend - 1] - 1] < cmp_tunableEnvironment[1]
            .f1->data[x->data[pmax] - 1])));
          if (varargout_1) {
            pow2p = x->data[pmax];
            x->data[pmax] = x->data[frame.xend - 1];
            x->data[frame.xend - 1] = pow2p;
          }

          pivot = x->data[pmax] - 1;
          x->data[pmax] = x->data[frame.xend - 2];
          x->data[frame.xend - 2] = pivot + 1;
          pmax = frame.xstart - 1;
          pmin = frame.xend - 2;
          do {
            exitg2 = 0;
            pmax++;
            do {
              exitg3 = 0;
              varargout_1 = ((cmp_tunableEnvironment[0].f1->data[x->data[pmax] -
                              1] < cmp_tunableEnvironment[0].f1->data[pivot]) ||
                             ((cmp_tunableEnvironment[0].f1->data[x->data[pmax]
                               - 1] == cmp_tunableEnvironment[0].f1->data[pivot])
                              && (cmp_tunableEnvironment[1].f1->data[x->
                                  data[pmax] - 1] < cmp_tunableEnvironment[1].
                                  f1->data[pivot])));
              if (varargout_1) {
                pmax++;
              } else {
                exitg3 = 1;
              }
            } while (exitg3 == 0);

            pmin--;
            do {
              exitg3 = 0;
              varargout_1 = ((cmp_tunableEnvironment[0].f1->data[pivot] <
                              cmp_tunableEnvironment[0].f1->data[x->data[pmin] -
                              1]) || ((cmp_tunableEnvironment[0].f1->data[pivot]
                == cmp_tunableEnvironment[0].f1->data[x->data[pmin] - 1]) &&
                (cmp_tunableEnvironment[1].f1->data[pivot] <
                 cmp_tunableEnvironment[1].f1->data[x->data[pmin] - 1])));
              if (varargout_1) {
                pmin--;
              } else {
                exitg3 = 1;
              }
            } while (exitg3 == 0);

            if (pmax + 1 >= pmin + 1) {
              exitg2 = 1;
            } else {
              pow2p = x->data[pmax];
              x->data[pmax] = x->data[pmin];
              x->data[pmin] = pow2p;
            }
          } while (exitg2 == 0);

          x->data[frame.xend - 2] = x->data[pmax];
          x->data[pmax] = pivot + 1;
          if (pmax + 2 < frame.xend) {
            st_d->data[p].xstart = pmax + 2;
            st_d->data[p].xend = frame.xend;
            st_d->data[p].depth = frame.depth + 1;
            p++;
          }

          if (frame.xstart < pmax + 1) {
            st_d->data[p].xstart = frame.xstart;
            st_d->data[p].xend = pmax + 1;
            st_d->data[p].depth = frame.depth + 1;
            p++;
          }
        }
      }

      emxFree_struct_T(&st_d);
    }
  }
}

/*
 * function [phi, rk] = lsqSVD(K, tau, phi_pre, count,...
 *     phi_r, noise_err, cond_max, lambda)
 * the function solves the least-square with SVD decomposition
 *  inputs:
 *  K - the matrix calculated from the measured motion states, the D-H
 *  parameters and the gravity
 *  tau - the measured torque
 *  phi_pre - the previously solved dynamic parameters
 *  count - the counter
 * Arguments    : emxArray_real_T *K
 *                emxArray_real_T *tau
 *                const emxArray_real_T *b_phi_pre
 *                double b_count
 *                const emxArray_real_T *b_phi_r
 *                double b_noise_err
 *                double b_cond_max
 *                double b_lambda
 *                emxArray_real_T *phi
 *                double *rk
 * Return Type  : void
 */
static void lsqSVD(emxArray_real_T *K, emxArray_real_T *tau, const
                   emxArray_real_T *b_phi_pre, double b_count, const
                   emxArray_real_T *b_phi_r, double b_noise_err, double
                   b_cond_max, double b_lambda, emxArray_real_T *phi, double *rk)
{
  emxArray_real_T *b_tau;
  int n;
  int nx;
  double step_ratio;
  int i15;
  int loop_ub;
  int idx;
  int i16;
  emxArray_boolean_T *r6;
  emxArray_boolean_T *r7;
  emxArray_real_T *scv;
  emxArray_real_T *b_scv;
  emxArray_real_T *b;
  emxArray_real_T *varargin_1;
  bool empty_non_axis_sizes;
  int result;
  emxArray_real_T *c_tau;
  emxArray_real_T *b_varargin_1;
  emxArray_real_T *s;
  emxArray_real_T *V;
  double ex;
  emxArray_boolean_T *e1;
  emxArray_int32_T *ii;
  bool exitg1;
  emxArray_real_T *b_b;
  emxArray_real_T *c_b;
  emxArray_int32_T *r8;
  emxArray_real_T *d_tau;
  emxInit_real_T(&b_tau, 2);

  /*  parameters: */
  /*  phi_r - the priori values of the dynamic parameters */
  /*  noise_err - the threshold for elements in K not zero */
  /*  cond_max - the maximum condition number of K */
  /*  lambda - the factor for considering the priori */
  /*  outputs: */
  /*  phi - the solved dynamic parameters in current iteration */
  /*  rk - the rank of the matrix K */
  /*  obtain the dimensions of data */
  /* 'dynIdenf:388' [pidenf, n] = size(tau); */
  n = tau->size[1];
  nx = tau->size[0];

  /* 'dynIdenf:389' nparMinSet = length(phi_pre); */
  /*  calculate the step_ratio according to the counter */
  /* 'dynIdenf:391' step_ratio = 1 / (count + 1); */
  step_ratio = 1.0 / (b_count + 1.0);

  /*  reshape the measured torque in matrix into a column vector with row wise */
  /* 'dynIdenf:393' tau = tau'; */
  i15 = b_tau->size[0] * b_tau->size[1];
  b_tau->size[0] = tau->size[1];
  b_tau->size[1] = tau->size[0];
  emxEnsureCapacity_real_T(b_tau, i15);
  loop_ub = tau->size[0];
  for (i15 = 0; i15 < loop_ub; i15++) {
    idx = tau->size[1];
    for (i16 = 0; i16 < idx; i16++) {
      b_tau->data[i16 + b_tau->size[0] * i15] = tau->data[i15 + tau->size[0] *
        i16];
    }
  }

  i15 = tau->size[0] * tau->size[1];
  tau->size[0] = b_tau->size[0];
  tau->size[1] = b_tau->size[1];
  emxEnsureCapacity_real_T(tau, i15);
  loop_ub = b_tau->size[1];
  for (i15 = 0; i15 < loop_ub; i15++) {
    idx = b_tau->size[0];
    for (i16 = 0; i16 < idx; i16++) {
      tau->data[i16 + tau->size[0] * i15] = b_tau->data[i16 + b_tau->size[0] *
        i15];
    }
  }

  emxInit_boolean_T1(&r6, 2);

  /* 'dynIdenf:394' tau = tau(:); */
  /*  detect near-zero elements in K and set them to be zero */
  /* 'dynIdenf:396' K = K .* (abs(K)>noise_err); */
  b_abs(K, b_tau);
  i15 = r6->size[0] * r6->size[1];
  r6->size[0] = b_tau->size[0];
  r6->size[1] = b_tau->size[1];
  emxEnsureCapacity_boolean_T(r6, i15);
  loop_ub = b_tau->size[1];
  for (i15 = 0; i15 < loop_ub; i15++) {
    idx = b_tau->size[0];
    for (i16 = 0; i16 < idx; i16++) {
      r6->data[i16 + r6->size[0] * i15] = (b_tau->data[i16 + b_tau->size[0] *
        i15] > b_noise_err);
    }
  }

  i15 = K->size[0] * K->size[1];
  emxEnsureCapacity_real_T(K, i15);
  loop_ub = K->size[1];
  for (i15 = 0; i15 < loop_ub; i15++) {
    idx = K->size[0];
    for (i16 = 0; i16 < idx; i16++) {
      K->data[i16 + K->size[0] * i15] *= (double)r6->data[i16 + r6->size[0] *
        i15];
    }
  }

  emxFree_boolean_T(&r6);
  emxInit_boolean_T1(&r7, 2);
  emxInit_real_T(&scv, 2);

  /*  calculate the norm of each column of K */
  /* 'dynIdenf:398' scv = ((vecnorm(K)==0) + vecnorm(K))'; */
  vecnorm(K, scv);
  i15 = r7->size[0] * r7->size[1];
  r7->size[0] = 1;
  r7->size[1] = scv->size[1];
  emxEnsureCapacity_boolean_T(r7, i15);
  loop_ub = scv->size[1];
  for (i15 = 0; i15 < loop_ub; i15++) {
    r7->data[r7->size[0] * i15] = (scv->data[scv->size[0] * i15] == 0.0);
  }

  emxInit_real_T1(&b_scv, 1);
  vecnorm(K, scv);
  i15 = b_scv->size[0];
  b_scv->size[0] = r7->size[1];
  emxEnsureCapacity_real_T1(b_scv, i15);
  loop_ub = r7->size[1];
  for (i15 = 0; i15 < loop_ub; i15++) {
    b_scv->data[i15] = (double)r7->data[r7->size[0] * i15] + scv->data[scv->
      size[0] * i15];
  }

  emxFree_boolean_T(&r7);
  emxInit_real_T(&b, 2);

  /*  normalize each column of K and combine it with a lambda-scaled */
  /*  identity */
  /* 'dynIdenf:401' K = [K ./ repmat(scv',n*pidenf,1); lambda * eye(nparMinSet)]; */
  b_eye(b_phi_pre->size[0], b);
  i15 = b->size[0] * b->size[1];
  emxEnsureCapacity_real_T(b, i15);
  loop_ub = b->size[1];
  for (i15 = 0; i15 < loop_ub; i15++) {
    idx = b->size[0];
    for (i16 = 0; i16 < idx; i16++) {
      b->data[i16 + b->size[0] * i15] *= b_lambda;
    }
  }

  i15 = scv->size[0] * scv->size[1];
  scv->size[0] = 1;
  scv->size[1] = b_scv->size[0];
  emxEnsureCapacity_real_T(scv, i15);
  loop_ub = b_scv->size[0];
  for (i15 = 0; i15 < loop_ub; i15++) {
    scv->data[scv->size[0] * i15] = b_scv->data[i15];
  }

  emxInit_real_T(&varargin_1, 2);
  repmat(scv, (double)n * (double)nx, b_tau);
  rdivide(K, b_tau, varargin_1);
  emxFree_real_T(&scv);
  if (!((varargin_1->size[0] == 0) || (varargin_1->size[1] == 0))) {
    idx = varargin_1->size[1];
  } else if (!((b->size[0] == 0) || (b->size[1] == 0))) {
    idx = b->size[1];
  } else {
    idx = varargin_1->size[1];
    if (!(idx > 0)) {
      idx = 0;
    }

    if (b->size[1] > idx) {
      idx = b->size[1];
    }
  }

  empty_non_axis_sizes = (idx == 0);
  if (empty_non_axis_sizes || (!((varargin_1->size[0] == 0) || (varargin_1->
         size[1] == 0)))) {
    n = varargin_1->size[0];
  } else {
    n = 0;
  }

  if (empty_non_axis_sizes || (!((b->size[0] == 0) || (b->size[1] == 0)))) {
    result = b->size[0];
  } else {
    result = 0;
  }

  /*  Initialize the current dynamic parameters by the previous dynamic parameters */
  /*  and scale it by the norm vector of K */
  /* 'dynIdenf:404' phi = phi_pre .* scv; */
  i15 = phi->size[0];
  phi->size[0] = b_phi_pre->size[0];
  emxEnsureCapacity_real_T1(phi, i15);
  loop_ub = b_phi_pre->size[0];
  for (i15 = 0; i15 < loop_ub; i15++) {
    phi->data[i15] = b_phi_pre->data[i15] * b_scv->data[i15];
  }

  emxInit_real_T1(&c_tau, 1);

  /*  combine the torque matrix with the priori dynamic parameters */
  /* 'dynIdenf:406' tau = [tau; lambda * phi_r.*scv]; */
  nx = tau->size[0] * tau->size[1];
  i15 = c_tau->size[0];
  c_tau->size[0] = nx + b_phi_r->size[0];
  emxEnsureCapacity_real_T1(c_tau, i15);
  for (i15 = 0; i15 < nx; i15++) {
    c_tau->data[i15] = tau->data[i15];
  }

  loop_ub = b_phi_r->size[0];
  for (i15 = 0; i15 < loop_ub; i15++) {
    c_tau->data[i15 + nx] = b_lambda * b_phi_r->data[i15] * b_scv->data[i15];
  }

  emxInit_real_T(&b_varargin_1, 2);

  /*  SVD decomposition of K, i.e. K = U*S*V' */
  /* 'dynIdenf:408' [U,S,V] = svd(K); */
  i15 = b_varargin_1->size[0] * b_varargin_1->size[1];
  b_varargin_1->size[0] = n + result;
  b_varargin_1->size[1] = idx;
  emxEnsureCapacity_real_T(b_varargin_1, i15);
  for (i15 = 0; i15 < idx; i15++) {
    for (i16 = 0; i16 < n; i16++) {
      b_varargin_1->data[i16 + b_varargin_1->size[0] * i15] = varargin_1->
        data[i16 + n * i15];
    }
  }

  for (i15 = 0; i15 < idx; i15++) {
    for (i16 = 0; i16 < result; i16++) {
      b_varargin_1->data[(i16 + n) + b_varargin_1->size[0] * i15] = b->data[i16
        + result * i15];
    }
  }

  emxInit_real_T1(&s, 1);
  emxInit_real_T(&V, 2);
  svd(b_varargin_1, varargin_1, s, V);
  n = varargin_1->size[1];
  nx = V->size[1];
  i15 = b_tau->size[0] * b_tau->size[1];
  b_tau->size[0] = n;
  b_tau->size[1] = nx;
  emxEnsureCapacity_real_T(b_tau, i15);
  emxFree_real_T(&b_varargin_1);
  for (i15 = 0; i15 < nx; i15++) {
    for (i16 = 0; i16 < n; i16++) {
      b_tau->data[i16 + b_tau->size[0] * i15] = 0.0;
    }
  }

  for (n = 0; n < s->size[0]; n++) {
    b_tau->data[n + b_tau->size[0] * n] = s->data[n];
  }

  /*  obtain the singular values */
  /* 'dynIdenf:410' s = diag(S); */
  diag(b_tau, s);

  /*  truncate the smaller singular values to make the condition number */
  /*  below cond_max */
  /* 'dynIdenf:413' e1 = s >= max(s)/cond_max; */
  emxFree_real_T(&b_tau);
  if (s->size[0] <= 2) {
    if (s->size[0] == 0) {
      ex = 0.0;
    } else if (s->size[0] == 1) {
      ex = s->data[0];
    } else if (s->data[0] < s->data[1]) {
      ex = s->data[1];
    } else {
      ex = s->data[0];
    }
  } else {
    ex = s->data[0];
    for (n = 1; n < s->size[0]; n++) {
      if (ex < s->data[n]) {
        ex = s->data[n];
      }
    }
  }

  emxInit_boolean_T(&e1, 1);
  ex /= b_cond_max;
  i15 = e1->size[0];
  e1->size[0] = s->size[0];
  emxEnsureCapacity_boolean_T1(e1, i15);
  loop_ub = s->size[0];
  for (i15 = 0; i15 < loop_ub; i15++) {
    e1->data[i15] = (s->data[i15] >= ex);
  }

  emxInit_int32_T1(&ii, 1);

  /*  calculate the rank of the truncated matrix K */
  /* 'dynIdenf:415' rk = length(find(e1)); */
  nx = e1->size[0];
  idx = 0;
  i15 = ii->size[0];
  ii->size[0] = e1->size[0];
  emxEnsureCapacity_int32_T(ii, i15);
  n = 1;
  exitg1 = false;
  while ((!exitg1) && (n <= nx)) {
    if (e1->data[n - 1]) {
      idx++;
      ii->data[idx - 1] = n;
      if (idx >= nx) {
        exitg1 = true;
      } else {
        n++;
      }
    } else {
      n++;
    }
  }

  if (e1->size[0] == 1) {
    if (idx == 0) {
      i15 = ii->size[0];
      ii->size[0] = 0;
      emxEnsureCapacity_int32_T(ii, i15);
    }
  } else {
    i15 = ii->size[0];
    if (1 > idx) {
      ii->size[0] = 0;
    } else {
      ii->size[0] = idx;
    }

    emxEnsureCapacity_int32_T(ii, i15);
  }

  result = ii->size[0];

  /*  transform the torque by U */
  /* 'dynIdenf:417' tau = U' * tau; */
  i15 = b->size[0] * b->size[1];
  b->size[0] = varargin_1->size[1];
  b->size[1] = varargin_1->size[0];
  emxEnsureCapacity_real_T(b, i15);
  loop_ub = varargin_1->size[0];
  emxFree_int32_T(&ii);
  for (i15 = 0; i15 < loop_ub; i15++) {
    idx = varargin_1->size[1];
    for (i16 = 0; i16 < idx; i16++) {
      b->data[i16 + b->size[0] * i15] = varargin_1->data[i15 + varargin_1->size
        [0] * i16];
    }
  }

  emxFree_real_T(&varargin_1);
  emxInit_real_T1(&b_b, 1);
  i15 = b_b->size[0];
  b_b->size[0] = c_tau->size[0];
  emxEnsureCapacity_real_T1(b_b, i15);
  loop_ub = c_tau->size[0];
  for (i15 = 0; i15 < loop_ub; i15++) {
    b_b->data[i15] = c_tau->data[i15];
  }

  emxInit_real_T1(&c_b, 1);
  if ((b->size[1] == 1) || (c_tau->size[0] == 1)) {
    i15 = c_b->size[0];
    c_b->size[0] = b->size[0];
    emxEnsureCapacity_real_T1(c_b, i15);
    loop_ub = b->size[0];
    for (i15 = 0; i15 < loop_ub; i15++) {
      c_b->data[i15] = 0.0;
      idx = b->size[1];
      for (i16 = 0; i16 < idx; i16++) {
        c_b->data[i15] += b->data[i15 + b->size[0] * i16] * c_tau->data[i16];
      }
    }

    i15 = c_tau->size[0];
    c_tau->size[0] = c_b->size[0];
    emxEnsureCapacity_real_T1(c_tau, i15);
    loop_ub = c_b->size[0];
    for (i15 = 0; i15 < loop_ub; i15++) {
      c_tau->data[i15] = c_b->data[i15];
    }
  } else {
    nx = b->size[0];
    i15 = c_tau->size[0];
    c_tau->size[0] = b->size[0];
    emxEnsureCapacity_real_T1(c_tau, i15);
    for (idx = 1; idx <= nx; idx++) {
      c_tau->data[idx - 1] = 0.0;
    }

    for (n = 0; n < b->size[1]; n++) {
      if (b_b->data[n] != 0.0) {
        for (idx = 0; idx < nx; idx++) {
          c_tau->data[idx] += b_b->data[n] * b->data[idx + b->size[0] * n];
        }
      }
    }
  }

  /*  transform the dynamic parameters by V */
  /* 'dynIdenf:419' phi = V' * phi; */
  i15 = b->size[0] * b->size[1];
  b->size[0] = V->size[1];
  b->size[1] = V->size[0];
  emxEnsureCapacity_real_T(b, i15);
  loop_ub = V->size[0];
  for (i15 = 0; i15 < loop_ub; i15++) {
    idx = V->size[1];
    for (i16 = 0; i16 < idx; i16++) {
      b->data[i16 + b->size[0] * i15] = V->data[i15 + V->size[0] * i16];
    }
  }

  i15 = b_b->size[0];
  b_b->size[0] = phi->size[0];
  emxEnsureCapacity_real_T1(b_b, i15);
  loop_ub = phi->size[0];
  for (i15 = 0; i15 < loop_ub; i15++) {
    b_b->data[i15] = phi->data[i15];
  }

  if ((b->size[1] == 1) || (phi->size[0] == 1)) {
    i15 = b_b->size[0];
    b_b->size[0] = b->size[0];
    emxEnsureCapacity_real_T1(b_b, i15);
    loop_ub = b->size[0];
    for (i15 = 0; i15 < loop_ub; i15++) {
      b_b->data[i15] = 0.0;
      idx = b->size[1];
      for (i16 = 0; i16 < idx; i16++) {
        b_b->data[i15] += b->data[i15 + b->size[0] * i16] * phi->data[i16];
      }
    }

    i15 = phi->size[0];
    phi->size[0] = b_b->size[0];
    emxEnsureCapacity_real_T1(phi, i15);
    loop_ub = b_b->size[0];
    for (i15 = 0; i15 < loop_ub; i15++) {
      phi->data[i15] = b_b->data[i15];
    }
  } else {
    nx = b->size[0];
    i15 = phi->size[0];
    phi->size[0] = b->size[0];
    emxEnsureCapacity_real_T1(phi, i15);
    for (idx = 1; idx <= nx; idx++) {
      phi->data[idx - 1] = 0.0;
    }

    for (n = 0; n < b->size[1]; n++) {
      if (b_b->data[n] != 0.0) {
        for (idx = 0; idx < nx; idx++) {
          phi->data[idx] += b_b->data[n] * b->data[idx + b->size[0] * n];
        }
      }
    }
  }

  emxFree_real_T(&b);

  /*  get the solution in the transoformed form */
  /* 'dynIdenf:421' phi(e1) = tau(e1) ./ s(e1); */
  nx = e1->size[0] - 1;
  n = 0;
  for (idx = 0; idx <= nx; idx++) {
    if (e1->data[idx]) {
      n++;
    }
  }

  emxInit_int32_T1(&r8, 1);
  i15 = r8->size[0];
  r8->size[0] = n;
  emxEnsureCapacity_int32_T(r8, i15);
  n = 0;
  for (idx = 0; idx <= nx; idx++) {
    if (e1->data[idx]) {
      r8->data[n] = idx + 1;
      n++;
    }
  }

  emxInit_real_T1(&d_tau, 1);
  i15 = d_tau->size[0];
  d_tau->size[0] = r8->size[0];
  emxEnsureCapacity_real_T1(d_tau, i15);
  loop_ub = r8->size[0];
  for (i15 = 0; i15 < loop_ub; i15++) {
    d_tau->data[i15] = c_tau->data[r8->data[i15] - 1];
  }

  emxFree_real_T(&c_tau);
  i15 = c_b->size[0];
  c_b->size[0] = r8->size[0];
  emxEnsureCapacity_real_T1(c_b, i15);
  loop_ub = r8->size[0];
  for (i15 = 0; i15 < loop_ub; i15++) {
    c_b->data[i15] = s->data[r8->data[i15] - 1];
  }

  emxFree_real_T(&s);
  emxFree_int32_T(&r8);
  b_rdivide(d_tau, c_b, b_b);
  nx = e1->size[0];
  n = 0;
  emxFree_real_T(&d_tau);
  for (idx = 0; idx < nx; idx++) {
    if (e1->data[idx]) {
      phi->data[idx] = b_b->data[n];
      n++;
    }
  }

  emxFree_boolean_T(&e1);

  /*  restore the solved dynamic parameters by inverse transformation */
  /*  and scaling */
  /* 'dynIdenf:424' phi = V * phi; */
  i15 = b_b->size[0];
  b_b->size[0] = phi->size[0];
  emxEnsureCapacity_real_T1(b_b, i15);
  loop_ub = phi->size[0];
  for (i15 = 0; i15 < loop_ub; i15++) {
    b_b->data[i15] = phi->data[i15];
  }

  if ((V->size[1] == 1) || (phi->size[0] == 1)) {
    i15 = c_b->size[0];
    c_b->size[0] = V->size[0];
    emxEnsureCapacity_real_T1(c_b, i15);
    loop_ub = V->size[0];
    for (i15 = 0; i15 < loop_ub; i15++) {
      c_b->data[i15] = 0.0;
      idx = V->size[1];
      for (i16 = 0; i16 < idx; i16++) {
        c_b->data[i15] += V->data[i15 + V->size[0] * i16] * phi->data[i16];
      }
    }

    i15 = phi->size[0];
    phi->size[0] = c_b->size[0];
    emxEnsureCapacity_real_T1(phi, i15);
    loop_ub = c_b->size[0];
    for (i15 = 0; i15 < loop_ub; i15++) {
      phi->data[i15] = c_b->data[i15];
    }
  } else {
    nx = V->size[0];
    i15 = phi->size[0];
    phi->size[0] = V->size[0];
    emxEnsureCapacity_real_T1(phi, i15);
    for (idx = 1; idx <= nx; idx++) {
      phi->data[idx - 1] = 0.0;
    }

    for (n = 0; n < V->size[1]; n++) {
      if (b_b->data[n] != 0.0) {
        for (idx = 0; idx < nx; idx++) {
          phi->data[idx] += b_b->data[n] * V->data[idx + V->size[0] * n];
        }
      }
    }
  }

  emxFree_real_T(&b_b);
  emxFree_real_T(&V);

  /* 'dynIdenf:425' phi = phi ./ scv; */
  i15 = c_b->size[0];
  c_b->size[0] = phi->size[0];
  emxEnsureCapacity_real_T1(c_b, i15);
  loop_ub = phi->size[0];
  for (i15 = 0; i15 < loop_ub; i15++) {
    c_b->data[i15] = phi->data[i15];
  }

  b_rdivide(c_b, b_scv, phi);

  /*  weight the currently calculated dynamic parameters and the */
  /*  previous one */
  /* 'dynIdenf:428' phi = phi_pre + step_ratio * (phi - phi_pre); */
  i15 = phi->size[0];
  phi->size[0] = b_phi_pre->size[0];
  emxEnsureCapacity_real_T1(phi, i15);
  loop_ub = b_phi_pre->size[0];
  emxFree_real_T(&c_b);
  emxFree_real_T(&b_scv);
  for (i15 = 0; i15 < loop_ub; i15++) {
    phi->data[i15] = b_phi_pre->data[i15] + step_ratio * (phi->data[i15] -
      b_phi_pre->data[i15]);
  }

  *rk = result;
}

/*
 * Arguments    : emxArray_int32_T *idx
 *                emxArray_real_T *x
 *                int offset
 *                int np
 *                int nq
 *                emxArray_int32_T *iwork
 *                emxArray_real_T *xwork
 * Return Type  : void
 */
static void merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int np,
                  int nq, emxArray_int32_T *iwork, emxArray_real_T *xwork)
{
  int n;
  int qend;
  int p;
  int iout;
  int exitg1;
  if ((np == 0) || (nq == 0)) {
  } else {
    n = np + nq;
    for (qend = 0; qend < n; qend++) {
      iwork->data[qend] = idx->data[offset + qend];
      xwork->data[qend] = x->data[offset + qend];
    }

    p = 0;
    n = np;
    qend = np + nq;
    iout = offset - 1;
    do {
      exitg1 = 0;
      iout++;
      if (xwork->data[p] >= xwork->data[n]) {
        idx->data[iout] = iwork->data[p];
        x->data[iout] = xwork->data[p];
        if (p + 1 < np) {
          p++;
        } else {
          exitg1 = 1;
        }
      } else {
        idx->data[iout] = iwork->data[n];
        x->data[iout] = xwork->data[n];
        if (n + 1 < qend) {
          n++;
        } else {
          n = iout - p;
          while (p + 1 <= np) {
            idx->data[(n + p) + 1] = iwork->data[p];
            x->data[(n + p) + 1] = xwork->data[p];
            p++;
          }

          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}

/*
 * Arguments    : emxArray_int32_T *idx
 *                emxArray_real_T *x
 *                int offset
 *                int n
 *                int preSortLevel
 *                emxArray_int32_T *iwork
 *                emxArray_real_T *xwork
 * Return Type  : void
 */
static void merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
  int n, int preSortLevel, emxArray_int32_T *iwork, emxArray_real_T *xwork)
{
  int nPairs;
  int bLen;
  int tailOffset;
  int nTail;
  nPairs = n >> preSortLevel;
  bLen = 1 << preSortLevel;
  while (nPairs > 1) {
    if ((nPairs & 1) != 0) {
      nPairs--;
      tailOffset = bLen * nPairs;
      nTail = n - tailOffset;
      if (nTail > bLen) {
        merge(idx, x, offset + tailOffset, bLen, nTail - bLen, iwork, xwork);
      }
    }

    tailOffset = bLen << 1;
    nPairs >>= 1;
    for (nTail = 1; nTail <= nPairs; nTail++) {
      merge(idx, x, offset + (nTail - 1) * tailOffset, bLen, bLen, iwork, xwork);
    }

    bLen = tailOffset;
  }

  if (n > bLen) {
    merge(idx, x, offset, bLen, n - bLen, iwork, xwork);
  }
}

/*
 * Arguments    : const emxArray_real_T *A
 *                emxArray_real_T *B
 * Return Type  : void
 */
static void mldivide(const emxArray_real_T *A, emxArray_real_T *B)
{
  emxArray_real_T *b_A;
  emxArray_real_T *tau;
  emxArray_int32_T *jpvt;
  emxArray_real_T *b_B;
  unsigned int unnamed_idx_0;
  int mn;
  int n;
  int m;
  int loop_ub;
  int kAcol;
  double wj;
  int i;
  emxInit_real_T(&b_A, 2);
  emxInit_real_T1(&tau, 1);
  emxInit_int32_T(&jpvt, 2);
  emxInit_real_T1(&b_B, 1);
  if ((A->size[0] == 0) || (A->size[1] == 0) || (B->size[0] == 0)) {
    unnamed_idx_0 = (unsigned int)A->size[1];
    mn = B->size[0];
    B->size[0] = (int)unnamed_idx_0;
    emxEnsureCapacity_real_T1(B, mn);
    m = (int)unnamed_idx_0;
    for (mn = 0; mn < m; mn++) {
      B->data[mn] = 0.0;
    }
  } else if (A->size[0] == A->size[1]) {
    n = A->size[1];
    mn = b_A->size[0] * b_A->size[1];
    b_A->size[0] = A->size[0];
    b_A->size[1] = A->size[1];
    emxEnsureCapacity_real_T(b_A, mn);
    m = A->size[1];
    for (mn = 0; mn < m; mn++) {
      loop_ub = A->size[0];
      for (kAcol = 0; kAcol < loop_ub; kAcol++) {
        b_A->data[kAcol + b_A->size[0] * mn] = A->data[kAcol + A->size[0] * mn];
      }
    }

    xgetrf(A->size[1], A->size[1], b_A, A->size[1], jpvt, &loop_ub);
    for (loop_ub = 0; loop_ub < n - 1; loop_ub++) {
      if (jpvt->data[loop_ub] != loop_ub + 1) {
        wj = B->data[loop_ub];
        B->data[loop_ub] = B->data[jpvt->data[loop_ub] - 1];
        B->data[jpvt->data[loop_ub] - 1] = wj;
      }
    }

    if (B->size[0] != 0) {
      for (loop_ub = 0; loop_ub < n; loop_ub++) {
        kAcol = n * loop_ub;
        if (B->data[loop_ub] != 0.0) {
          for (i = loop_ub + 1; i < n; i++) {
            B->data[i] -= B->data[loop_ub] * b_A->data[i + kAcol];
          }
        }
      }
    }

    if (B->size[0] != 0) {
      for (loop_ub = A->size[1] - 1; loop_ub + 1 > 0; loop_ub--) {
        kAcol = n * loop_ub;
        if (B->data[loop_ub] != 0.0) {
          B->data[loop_ub] /= b_A->data[loop_ub + kAcol];
          for (i = 0; i < loop_ub; i++) {
            B->data[i] -= B->data[loop_ub] * b_A->data[i + kAcol];
          }
        }
      }
    }
  } else {
    mn = b_A->size[0] * b_A->size[1];
    b_A->size[0] = A->size[0];
    b_A->size[1] = A->size[1];
    emxEnsureCapacity_real_T(b_A, mn);
    m = A->size[1];
    for (mn = 0; mn < m; mn++) {
      loop_ub = A->size[0];
      for (kAcol = 0; kAcol < loop_ub; kAcol++) {
        b_A->data[kAcol + b_A->size[0] * mn] = A->data[kAcol + A->size[0] * mn];
      }
    }

    xgeqp3(b_A, tau, jpvt);
    n = rankFromQR(b_A);
    mn = b_B->size[0];
    b_B->size[0] = B->size[0];
    emxEnsureCapacity_real_T1(b_B, mn);
    m = B->size[0];
    for (mn = 0; mn < m; mn++) {
      b_B->data[mn] = B->data[mn];
    }

    loop_ub = b_A->size[1];
    mn = B->size[0];
    B->size[0] = loop_ub;
    emxEnsureCapacity_real_T1(B, mn);
    for (mn = 0; mn < loop_ub; mn++) {
      B->data[mn] = 0.0;
    }

    m = b_A->size[0];
    loop_ub = b_A->size[0];
    mn = b_A->size[1];
    if (loop_ub < mn) {
      mn = loop_ub;
    }

    for (kAcol = 0; kAcol < mn; kAcol++) {
      if (tau->data[kAcol] != 0.0) {
        wj = b_B->data[kAcol];
        for (i = kAcol + 1; i < m; i++) {
          wj += b_A->data[i + b_A->size[0] * kAcol] * b_B->data[i];
        }

        wj *= tau->data[kAcol];
        if (wj != 0.0) {
          b_B->data[kAcol] -= wj;
          for (i = kAcol + 1; i < m; i++) {
            b_B->data[i] -= b_A->data[i + b_A->size[0] * kAcol] * wj;
          }
        }
      }
    }

    for (i = 0; i < n; i++) {
      B->data[jpvt->data[i] - 1] = b_B->data[i];
    }

    for (kAcol = n - 1; kAcol + 1 > 0; kAcol--) {
      B->data[jpvt->data[kAcol] - 1] /= b_A->data[kAcol + b_A->size[0] * kAcol];
      for (i = 0; i < kAcol; i++) {
        B->data[jpvt->data[i] - 1] -= B->data[jpvt->data[kAcol] - 1] * b_A->
          data[i + b_A->size[0] * kAcol];
      }
    }
  }

  emxFree_real_T(&b_B);
  emxFree_int32_T(&jpvt);
  emxFree_real_T(&tau);
  emxFree_real_T(&b_A);
}

/*
 * Arguments    : emxArray_real_T *A
 *                const emxArray_real_T *B
 * Return Type  : void
 */
static void mrdivide(emxArray_real_T *A, const emxArray_real_T *B)
{
  emxArray_real_T *Y;
  emxArray_real_T *b_B;
  emxArray_real_T *b_A;
  emxArray_real_T *tau;
  emxArray_int32_T *jpvt;
  unsigned int unnamed_idx_0;
  int mn;
  int n;
  unsigned int unnamed_idx_1;
  int m;
  int jBcol;
  int j;
  int nb;
  int k;
  int b_nb;
  double wj;
  int i;
  emxInit_real_T(&Y, 2);
  emxInit_real_T(&b_B, 2);
  emxInit_real_T(&b_A, 2);
  emxInit_real_T1(&tau, 1);
  emxInit_int32_T(&jpvt, 2);
  if ((A->size[0] == 0) || (A->size[1] == 0) || ((B->size[0] == 0) || (B->size[1]
        == 0))) {
    unnamed_idx_0 = (unsigned int)A->size[0];
    unnamed_idx_1 = (unsigned int)B->size[0];
    mn = A->size[0] * A->size[1];
    A->size[0] = (int)unnamed_idx_0;
    A->size[1] = (int)unnamed_idx_1;
    emxEnsureCapacity_real_T(A, mn);
    m = (int)unnamed_idx_1;
    for (mn = 0; mn < m; mn++) {
      jBcol = (int)unnamed_idx_0;
      for (j = 0; j < jBcol; j++) {
        A->data[j + A->size[0] * mn] = 0.0;
      }
    }
  } else if (B->size[0] == B->size[1]) {
    n = B->size[1];
    mn = b_A->size[0] * b_A->size[1];
    b_A->size[0] = B->size[0];
    b_A->size[1] = B->size[1];
    emxEnsureCapacity_real_T(b_A, mn);
    m = B->size[1];
    for (mn = 0; mn < m; mn++) {
      jBcol = B->size[0];
      for (j = 0; j < jBcol; j++) {
        b_A->data[j + b_A->size[0] * mn] = B->data[j + B->size[0] * mn];
      }
    }

    xgetrf(B->size[1], B->size[1], b_A, B->size[1], jpvt, &m);
    nb = A->size[0];
    if (!((A->size[0] == 0) || (A->size[1] == 0))) {
      for (j = 0; j < n; j++) {
        jBcol = nb * j;
        m = n * j;
        for (k = 1; k <= j; k++) {
          b_nb = nb * (k - 1);
          if (b_A->data[(k + m) - 1] != 0.0) {
            for (i = 0; i < nb; i++) {
              A->data[i + jBcol] -= b_A->data[(k + m) - 1] * A->data[i + b_nb];
            }
          }
        }

        wj = 1.0 / b_A->data[j + m];
        for (i = 0; i < nb; i++) {
          A->data[i + jBcol] *= wj;
        }
      }
    }

    if (!((A->size[0] == 0) || (A->size[1] == 0))) {
      for (j = B->size[1]; j > 0; j--) {
        jBcol = nb * (j - 1);
        m = n * (j - 1);
        for (k = j; k < n; k++) {
          b_nb = nb * k;
          if (b_A->data[k + m] != 0.0) {
            for (i = 0; i < nb; i++) {
              A->data[i + jBcol] -= b_A->data[k + m] * A->data[i + b_nb];
            }
          }
        }
      }
    }

    for (m = B->size[1] - 2; m + 1 > 0; m--) {
      if (jpvt->data[m] != m + 1) {
        jBcol = jpvt->data[m] - 1;
        for (b_nb = 0; b_nb < nb; b_nb++) {
          wj = A->data[b_nb + A->size[0] * m];
          A->data[b_nb + A->size[0] * m] = A->data[b_nb + A->size[0] * jBcol];
          A->data[b_nb + A->size[0] * jBcol] = wj;
        }
      }
    }
  } else {
    mn = b_B->size[0] * b_B->size[1];
    b_B->size[0] = A->size[1];
    b_B->size[1] = A->size[0];
    emxEnsureCapacity_real_T(b_B, mn);
    m = A->size[0];
    for (mn = 0; mn < m; mn++) {
      jBcol = A->size[1];
      for (j = 0; j < jBcol; j++) {
        b_B->data[j + b_B->size[0] * mn] = A->data[mn + A->size[0] * j];
      }
    }

    mn = b_A->size[0] * b_A->size[1];
    b_A->size[0] = B->size[1];
    b_A->size[1] = B->size[0];
    emxEnsureCapacity_real_T(b_A, mn);
    m = B->size[0];
    for (mn = 0; mn < m; mn++) {
      jBcol = B->size[1];
      for (j = 0; j < jBcol; j++) {
        b_A->data[j + b_A->size[0] * mn] = B->data[mn + B->size[0] * j];
      }
    }

    xgeqp3(b_A, tau, jpvt);
    n = rankFromQR(b_A);
    nb = b_B->size[1];
    jBcol = b_A->size[1];
    b_nb = b_B->size[1];
    mn = Y->size[0] * Y->size[1];
    Y->size[0] = jBcol;
    Y->size[1] = b_nb;
    emxEnsureCapacity_real_T(Y, mn);
    for (mn = 0; mn < b_nb; mn++) {
      for (j = 0; j < jBcol; j++) {
        Y->data[j + Y->size[0] * mn] = 0.0;
      }
    }

    m = b_A->size[0];
    b_nb = b_B->size[1];
    jBcol = b_A->size[0];
    mn = b_A->size[1];
    if (jBcol < mn) {
      mn = jBcol;
    }

    for (j = 0; j < mn; j++) {
      if (tau->data[j] != 0.0) {
        for (k = 0; k < b_nb; k++) {
          wj = b_B->data[j + b_B->size[0] * k];
          for (i = j + 1; i < m; i++) {
            wj += b_A->data[i + b_A->size[0] * j] * b_B->data[i + b_B->size[0] *
              k];
          }

          wj *= tau->data[j];
          if (wj != 0.0) {
            b_B->data[j + b_B->size[0] * k] -= wj;
            for (i = j + 1; i < m; i++) {
              b_B->data[i + b_B->size[0] * k] -= b_A->data[i + b_A->size[0] * j]
                * wj;
            }
          }
        }
      }
    }

    for (k = 0; k < nb; k++) {
      for (i = 0; i < n; i++) {
        Y->data[(jpvt->data[i] + Y->size[0] * k) - 1] = b_B->data[i + b_B->size
          [0] * k];
      }

      for (j = n - 1; j + 1 > 0; j--) {
        Y->data[(jpvt->data[j] + Y->size[0] * k) - 1] /= b_A->data[j + b_A->
          size[0] * j];
        for (i = 0; i < j; i++) {
          Y->data[(jpvt->data[i] + Y->size[0] * k) - 1] -= Y->data[(jpvt->data[j]
            + Y->size[0] * k) - 1] * b_A->data[i + b_A->size[0] * j];
        }
      }
    }

    mn = A->size[0] * A->size[1];
    A->size[0] = Y->size[1];
    A->size[1] = Y->size[0];
    emxEnsureCapacity_real_T(A, mn);
    m = Y->size[0];
    for (mn = 0; mn < m; mn++) {
      jBcol = Y->size[1];
      for (j = 0; j < jBcol; j++) {
        A->data[j + A->size[0] * mn] = Y->data[mn + Y->size[0] * j];
      }
    }
  }

  emxFree_int32_T(&jpvt);
  emxFree_real_T(&tau);
  emxFree_real_T(&b_A);
  emxFree_real_T(&b_B);
  emxFree_real_T(&Y);
}

/*
 * function y_out = myfiltfilt(b,a,x)
 * the function filter the data x with the filter specified by coefficients
 *  b and a
 * Arguments    : emxArray_real_T *b
 *                double b_a
 *                const emxArray_real_T *x
 *                emxArray_real_T *y_out
 * Return Type  : void
 */
static void myfiltfilt(emxArray_real_T *b, double b_a, const emxArray_real_T *x,
  emxArray_real_T *y_out)
{
  emxArray_real_T *b_b;
  int c_b;
  int i7;
  int loop_ub;
  int nfilt;
  double nfact;
  emxArray_real_T *a1;
  emxArray_real_T *y;
  emxArray_real_T *b_y;
  emxArray_int32_T *c_y;
  emxArray_int32_T *d_y;
  emxArray_int32_T *e_y;
  int i8;
  emxArray_real_T *f_y;
  int i9;
  emxArray_real_T *sp_d;
  coder_internal_sparse expl_temp;
  emxArray_int32_T *sp_colidx;
  emxArray_int32_T *sp_rowidx;
  int sp_n;
  emxArray_real_T *zi;
  emxArray_real_T *r4;
  double d_b;
  emxArray_real_T *g_y;
  double d1;
  emxArray_real_T *h_y;
  emxInit_real_T(&b_b, 2);

  /*  inputs: */
  /*  x - the raw data */
  /*  parameters */
  /*  b - the numerator of the transfer function of the filter */
  /*  a - the denominator of the transfer function of the filter */
  /*  outputs */
  /*  y - the filtered data */
  /*  length of the data */
  /* 'dynIdenf:446' len = size(x,1); */
  /*  convert the coefficients into row vector */
  /* 'dynIdenf:448' b = b(:).'; */
  c_b = b->size[1];
  i7 = b_b->size[0] * b_b->size[1];
  b_b->size[0] = 1;
  b_b->size[1] = c_b;
  emxEnsureCapacity_real_T(b_b, i7);
  for (i7 = 0; i7 < c_b; i7++) {
    b_b->data[b_b->size[0] * i7] = b->data[i7];
  }

  i7 = b->size[0] * b->size[1];
  b->size[0] = 1;
  b->size[1] = b_b->size[1];
  emxEnsureCapacity_real_T(b, i7);
  loop_ub = b_b->size[1];
  for (i7 = 0; i7 < loop_ub; i7++) {
    b->data[b->size[0] * i7] = b_b->data[b_b->size[0] * i7];
  }

  /* 'dynIdenf:449' a = a(:).'; */
  /*  number of elements in the numerator and denominator */
  /* 'dynIdenf:451' nb = length(b); */
  /* 'dynIdenf:452' na = length(a); */
  /*  the maximum number of elements */
  /* 'dynIdenf:454' nfilt = max(nb,na); */
  nfilt = (int)fmax(b->size[1], 1.0);

  /*  length of edge transients */
  /* 'dynIdenf:456' nfact = 3*(nfilt-1); */
  nfact = 3.0 * ((double)nfilt - 1.0);

  /*  set up filter's initial conditions to remove dc offset problems at the */
  /*  beginning and end of the sequence */
  /* 'dynIdenf:459' if nb < nfilt */
  if (b->size[1] < nfilt) {
    /* 'dynIdenf:460' b = [b, zeros(1,nfilt-nb)]; */
    c_b = b->size[1];
    c_b = nfilt - c_b;
    i7 = b_b->size[0] * b_b->size[1];
    b_b->size[0] = 1;
    b_b->size[1] = b->size[1] + c_b;
    emxEnsureCapacity_real_T(b_b, i7);
    loop_ub = b->size[1];
    for (i7 = 0; i7 < loop_ub; i7++) {
      b_b->data[b_b->size[0] * i7] = b->data[b->size[0] * i7];
    }

    for (i7 = 0; i7 < c_b; i7++) {
      b_b->data[b_b->size[0] * (i7 + b->size[1])] = 0.0;
    }

    i7 = b->size[0] * b->size[1];
    b->size[0] = 1;
    b->size[1] = b_b->size[1];
    emxEnsureCapacity_real_T(b, i7);
    loop_ub = b_b->size[1];
    for (i7 = 0; i7 < loop_ub; i7++) {
      b->data[b->size[0] * i7] = b_b->data[b_b->size[0] * i7];
    }
  }

  /*  zero-pad if necessary */
  /* 'dynIdenf:463' if na < nfilt */
  emxInit_real_T(&a1, 2);
  if (1 < nfilt) {
    /* 'dynIdenf:464' a1 = [a, zeros(1,nfilt-na)]; */
    i7 = a1->size[0] * a1->size[1];
    a1->size[0] = 1;
    a1->size[1] = nfilt;
    emxEnsureCapacity_real_T(a1, i7);
    a1->data[0] = b_a;
    loop_ub = nfilt - 1;
    for (i7 = 0; i7 < loop_ub; i7++) {
      a1->data[a1->size[0] * (i7 + 1)] = 0.0;
    }
  } else {
    /* 'dynIdenf:465' else */
    /* 'dynIdenf:466' a1 = a; */
    i7 = a1->size[0] * a1->size[1];
    a1->size[0] = 1;
    a1->size[1] = 1;
    emxEnsureCapacity_real_T(a1, i7);
    a1->data[0] = b_a;
  }

  /*  use sparse matrix to solve system of linear equations for initial conditions */
  /* 'dynIdenf:469' rows = [1:nfilt-1  2:nfilt-1  1:nfilt-2]; */
  emxInit_real_T(&y, 2);
  if (nfilt - 1 < 1) {
    i7 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = 0;
    emxEnsureCapacity_real_T(y, i7);
  } else {
    i7 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = nfilt - 1;
    emxEnsureCapacity_real_T(y, i7);
    loop_ub = nfilt - 2;
    for (i7 = 0; i7 <= loop_ub; i7++) {
      y->data[y->size[0] * i7] = 1.0 + (double)i7;
    }
  }

  emxInit_real_T(&b_y, 2);
  if (nfilt - 1 < 2) {
    i7 = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    b_y->size[1] = 0;
    emxEnsureCapacity_real_T(b_y, i7);
  } else {
    i7 = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    b_y->size[1] = nfilt - 2;
    emxEnsureCapacity_real_T(b_y, i7);
    loop_ub = nfilt - 3;
    for (i7 = 0; i7 <= loop_ub; i7++) {
      b_y->data[b_y->size[0] * i7] = 2.0 + (double)i7;
    }
  }

  emxInit_int32_T(&c_y, 2);
  if (nfilt - 2 < 1) {
    i7 = c_y->size[0] * c_y->size[1];
    c_y->size[0] = 1;
    c_y->size[1] = 0;
    emxEnsureCapacity_int32_T1(c_y, i7);
  } else {
    loop_ub = nfilt - 3;
    i7 = c_y->size[0] * c_y->size[1];
    c_y->size[0] = 1;
    c_y->size[1] = nfilt - 2;
    emxEnsureCapacity_int32_T1(c_y, i7);
    for (i7 = 0; i7 <= loop_ub; i7++) {
      c_y->data[c_y->size[0] * i7] = 1 + i7;
    }
  }

  /* 'dynIdenf:470' cols = [ones(1,nfilt-1) 2:nfilt-1  2:nfilt-1]; */
  emxInit_int32_T(&d_y, 2);
  if (nfilt - 1 < 2) {
    i7 = d_y->size[0] * d_y->size[1];
    d_y->size[0] = 1;
    d_y->size[1] = 0;
    emxEnsureCapacity_int32_T1(d_y, i7);
  } else {
    loop_ub = nfilt - 3;
    i7 = d_y->size[0] * d_y->size[1];
    d_y->size[0] = 1;
    d_y->size[1] = nfilt - 2;
    emxEnsureCapacity_int32_T1(d_y, i7);
    for (i7 = 0; i7 <= loop_ub; i7++) {
      d_y->data[d_y->size[0] * i7] = 2 + i7;
    }
  }

  emxInit_int32_T(&e_y, 2);
  if (nfilt - 1 < 2) {
    i7 = e_y->size[0] * e_y->size[1];
    e_y->size[0] = 1;
    e_y->size[1] = 0;
    emxEnsureCapacity_int32_T1(e_y, i7);
  } else {
    loop_ub = nfilt - 3;
    i7 = e_y->size[0] * e_y->size[1];
    e_y->size[0] = 1;
    e_y->size[1] = nfilt - 2;
    emxEnsureCapacity_int32_T1(e_y, i7);
    for (i7 = 0; i7 <= loop_ub; i7++) {
      e_y->data[e_y->size[0] * i7] = 2 + i7;
    }
  }

  /* 'dynIdenf:471' data = [1+a1(2) a1(3:nfilt) ones(1,nfilt-2)  -ones(1,nfilt-2)]; */
  if (3 > nfilt) {
    i7 = -1;
    i8 = 0;
  } else {
    i7 = 1;
    i8 = nfilt;
  }

  emxInit_real_T(&f_y, 2);

  /* 'dynIdenf:472' sp = sparse(rows,cols,data); */
  i9 = f_y->size[0] * f_y->size[1];
  f_y->size[0] = 1;
  f_y->size[1] = (y->size[1] + b_y->size[1]) + c_y->size[1];
  emxEnsureCapacity_real_T(f_y, i9);
  loop_ub = y->size[1];
  for (i9 = 0; i9 < loop_ub; i9++) {
    f_y->data[f_y->size[0] * i9] = y->data[y->size[0] * i9];
  }

  loop_ub = b_y->size[1];
  for (i9 = 0; i9 < loop_ub; i9++) {
    f_y->data[f_y->size[0] * (i9 + y->size[1])] = b_y->data[b_y->size[0] * i9];
  }

  loop_ub = c_y->size[1];
  for (i9 = 0; i9 < loop_ub; i9++) {
    f_y->data[f_y->size[0] * ((i9 + y->size[1]) + b_y->size[1])] = c_y->data
      [c_y->size[0] * i9];
  }

  emxFree_int32_T(&c_y);
  i9 = b_b->size[0] * b_b->size[1];
  b_b->size[0] = 1;
  b_b->size[1] = ((nfilt + d_y->size[1]) + e_y->size[1]) - 1;
  emxEnsureCapacity_real_T(b_b, i9);
  loop_ub = nfilt - 1;
  for (i9 = 0; i9 < loop_ub; i9++) {
    b_b->data[b_b->size[0] * i9] = 1.0;
  }

  loop_ub = d_y->size[1];
  for (i9 = 0; i9 < loop_ub; i9++) {
    b_b->data[b_b->size[0] * ((i9 + nfilt) - 1)] = d_y->data[d_y->size[0] * i9];
  }

  loop_ub = e_y->size[1];
  for (i9 = 0; i9 < loop_ub; i9++) {
    b_b->data[b_b->size[0] * (((i9 + nfilt) + d_y->size[1]) - 1)] = e_y->
      data[e_y->size[0] * i9];
  }

  emxFree_int32_T(&e_y);
  emxFree_int32_T(&d_y);
  i9 = b_y->size[0] * b_y->size[1];
  b_y->size[0] = 1;
  b_y->size[1] = (((i8 - i7) + nfilt) + nfilt) - 4;
  emxEnsureCapacity_real_T(b_y, i9);
  b_y->data[0] = 1.0 + a1->data[1];
  loop_ub = i8 - i7;
  for (i9 = 0; i9 <= loop_ub - 2; i9++) {
    b_y->data[b_y->size[0] * (i9 + 1)] = a1->data[(i7 + i9) + 1];
  }

  loop_ub = nfilt - 2;
  for (i9 = 0; i9 < loop_ub; i9++) {
    b_y->data[b_y->size[0] * ((i9 + i8) - i7)] = 1.0;
  }

  loop_ub = nfilt - 2;
  for (i9 = 0; i9 < loop_ub; i9++) {
    b_y->data[b_y->size[0] * ((((i9 + i8) - i7) + nfilt) - 2)] = -1.0;
  }

  emxInit_real_T1(&sp_d, 1);
  c_emxInitStruct_coder_internal_(&expl_temp);
  sparse(f_y, b_b, b_y, &expl_temp);
  i7 = sp_d->size[0];
  sp_d->size[0] = expl_temp.d->size[0];
  emxEnsureCapacity_real_T1(sp_d, i7);
  loop_ub = expl_temp.d->size[0];
  emxFree_real_T(&f_y);
  for (i7 = 0; i7 < loop_ub; i7++) {
    sp_d->data[i7] = expl_temp.d->data[i7];
  }

  emxInit_int32_T1(&sp_colidx, 1);
  i7 = sp_colidx->size[0];
  sp_colidx->size[0] = expl_temp.colidx->size[0];
  emxEnsureCapacity_int32_T(sp_colidx, i7);
  loop_ub = expl_temp.colidx->size[0];
  for (i7 = 0; i7 < loop_ub; i7++) {
    sp_colidx->data[i7] = expl_temp.colidx->data[i7];
  }

  emxInit_int32_T1(&sp_rowidx, 1);
  i7 = sp_rowidx->size[0];
  sp_rowidx->size[0] = expl_temp.rowidx->size[0];
  emxEnsureCapacity_int32_T(sp_rowidx, i7);
  loop_ub = expl_temp.rowidx->size[0];
  for (i7 = 0; i7 < loop_ub; i7++) {
    sp_rowidx->data[i7] = expl_temp.rowidx->data[i7];
  }

  c_b = expl_temp.m;
  sp_n = expl_temp.n;

  /*  zi are the steady-state states of the filter b(z)/a(z) in the state-space */
  /*  non-sparse: */
  /*  zi = ( eye(nfilt-1) - [-a(2:nfilt).' [eye(nfilt-2); zeros(1,nfilt-2)]] ) \ ... */
  /*       ( b(2:nfilt).' - a(2:nfilt).'*b(1) ); */
  /* 'dynIdenf:477' zi = full(sp) \ ( b(2:nfilt).' - a1(2:nfilt).'*b(1) ); */
  c_emxFreeStruct_coder_internal_(&expl_temp);
  if (2 > nfilt) {
    i7 = 0;
    i8 = 0;
  } else {
    i7 = 1;
    i8 = nfilt;
  }

  i9 = !(2 > nfilt);
  emxInit_real_T1(&zi, 1);
  emxInit_real_T(&r4, 2);
  sparse_full(sp_d, sp_colidx, sp_rowidx, c_b, sp_n, r4);
  d_b = b->data[0];
  c_b = zi->size[0];
  zi->size[0] = i8 - i7;
  emxEnsureCapacity_real_T1(zi, c_b);
  loop_ub = i8 - i7;
  emxFree_int32_T(&sp_rowidx);
  for (i8 = 0; i8 < loop_ub; i8++) {
    zi->data[i8] = b->data[i7 + i8] - a1->data[i9 + i8] * d_b;
  }

  mldivide(r4, zi);

  /*  Extrapolate beginning and end of data sequence using a "reflection */
  /*  method".  Slopes of original and extrapolated sequences match at */
  /*  the end points. This reduces end effects. */
  /* 'dynIdenf:481' y = [2*x(1)-x((nfact+1):-1:2);x;2*x(len)-x((len-1):-1:len-nfact)]; */
  emxFree_real_T(&r4);
  if (2.0 > nfact + 1.0) {
    i7 = 1;
    i8 = 1;
    i9 = 0;
  } else {
    i7 = (int)(nfact + 1.0);
    i8 = -1;
    i9 = 2;
  }

  sp_n = x->size[0];
  c_b = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = (int)-((0.0 - nfact) - -1.0) + 1;
  emxEnsureCapacity_real_T(y, c_b);
  loop_ub = (int)-((0.0 - nfact) - -1.0);
  for (c_b = 0; c_b <= loop_ub; c_b++) {
    y->data[y->size[0] * c_b] = (double)sp_n + (-1.0 - (double)c_b);
  }

  emxInit_real_T1(&g_y, 1);
  d_b = 2.0 * x->data[0];
  d1 = 2.0 * x->data[x->size[0] - 1];
  c_b = g_y->size[0];
  g_y->size[0] = ((div_s32_floor(i9 - i7, i8) + x->size[0]) + y->size[1]) + 1;
  emxEnsureCapacity_real_T1(g_y, c_b);
  loop_ub = div_s32_floor(i9 - i7, i8);
  for (c_b = 0; c_b <= loop_ub; c_b++) {
    g_y->data[c_b] = d_b - x->data[(i7 + i8 * c_b) - 1];
  }

  loop_ub = x->size[0];
  for (c_b = 0; c_b < loop_ub; c_b++) {
    g_y->data[(c_b + div_s32_floor(i9 - i7, i8)) + 1] = x->data[c_b];
  }

  loop_ub = y->size[1];
  for (c_b = 0; c_b < loop_ub; c_b++) {
    g_y->data[((c_b + div_s32_floor(i9 - i7, i8)) + x->size[0]) + 1] = d1 -
      x->data[(int)y->data[y->size[0] * c_b] - 1];
  }

  /*  filter, reverse data, filter again, and reverse data again */
  /* 'dynIdenf:483' y = filter(b,a1,y,zi*y(1)); */
  i7 = b_b->size[0] * b_b->size[1];
  b_b->size[0] = 1;
  b_b->size[1] = b->size[1];
  emxEnsureCapacity_real_T(b_b, i7);
  loop_ub = b->size[1];
  for (i7 = 0; i7 < loop_ub; i7++) {
    b_b->data[b_b->size[0] * i7] = b->data[b->size[0] * i7];
  }

  i7 = b_y->size[0] * b_y->size[1];
  b_y->size[0] = 1;
  b_y->size[1] = a1->size[1];
  emxEnsureCapacity_real_T(b_y, i7);
  loop_ub = a1->size[1];
  for (i7 = 0; i7 < loop_ub; i7++) {
    b_y->data[b_y->size[0] * i7] = a1->data[a1->size[0] * i7];
  }

  d_b = g_y->data[0];
  i7 = sp_d->size[0];
  sp_d->size[0] = zi->size[0];
  emxEnsureCapacity_real_T1(sp_d, i7);
  loop_ub = zi->size[0];
  for (i7 = 0; i7 < loop_ub; i7++) {
    sp_d->data[i7] = zi->data[i7] * d_b;
  }

  emxInit_real_T1(&h_y, 1);
  i7 = h_y->size[0];
  h_y->size[0] = g_y->size[0];
  emxEnsureCapacity_real_T1(h_y, i7);
  loop_ub = g_y->size[0];
  for (i7 = 0; i7 < loop_ub; i7++) {
    h_y->data[i7] = g_y->data[i7];
  }

  filter(b_b, b_y, h_y, sp_d, g_y);

  /* 'dynIdenf:484' y = y(length(y):-1:1); */
  if (1 > g_y->size[0]) {
    i7 = 1;
    i8 = 1;
    i9 = 0;
  } else {
    i7 = g_y->size[0];
    i8 = -1;
    i9 = 1;
  }

  c_b = h_y->size[0];
  h_y->size[0] = div_s32_floor(i9 - i7, i8) + 1;
  emxEnsureCapacity_real_T1(h_y, c_b);
  loop_ub = div_s32_floor(i9 - i7, i8);
  for (i9 = 0; i9 <= loop_ub; i9++) {
    h_y->data[i9] = g_y->data[(i7 + i8 * i9) - 1];
  }

  i7 = g_y->size[0];
  g_y->size[0] = h_y->size[0];
  emxEnsureCapacity_real_T1(g_y, i7);
  loop_ub = h_y->size[0];
  for (i7 = 0; i7 < loop_ub; i7++) {
    g_y->data[i7] = h_y->data[i7];
  }

  /* 'dynIdenf:485' y = filter(b,a1,y,zi*y(1)); */
  i7 = b_b->size[0] * b_b->size[1];
  b_b->size[0] = 1;
  b_b->size[1] = b->size[1];
  emxEnsureCapacity_real_T(b_b, i7);
  loop_ub = b->size[1];
  for (i7 = 0; i7 < loop_ub; i7++) {
    b_b->data[b_b->size[0] * i7] = b->data[b->size[0] * i7];
  }

  d_b = g_y->data[0];
  i7 = sp_d->size[0];
  sp_d->size[0] = zi->size[0];
  emxEnsureCapacity_real_T1(sp_d, i7);
  loop_ub = zi->size[0];
  for (i7 = 0; i7 < loop_ub; i7++) {
    sp_d->data[i7] = zi->data[i7] * d_b;
  }

  emxFree_real_T(&zi);
  i7 = h_y->size[0];
  h_y->size[0] = g_y->size[0];
  emxEnsureCapacity_real_T1(h_y, i7);
  loop_ub = g_y->size[0];
  for (i7 = 0; i7 < loop_ub; i7++) {
    h_y->data[i7] = g_y->data[i7];
  }

  filter(b_b, a1, h_y, sp_d, g_y);

  /* 'dynIdenf:486' y = y(length(y):-1:1); */
  emxFree_real_T(&b_b);
  emxFree_real_T(&sp_d);
  emxFree_real_T(&a1);
  if (1 > g_y->size[0]) {
    i7 = 1;
    i8 = 1;
    i9 = 0;
  } else {
    i7 = g_y->size[0];
    i8 = -1;
    i9 = 1;
  }

  c_b = h_y->size[0];
  h_y->size[0] = div_s32_floor(i9 - i7, i8) + 1;
  emxEnsureCapacity_real_T1(h_y, c_b);
  loop_ub = div_s32_floor(i9 - i7, i8);
  for (i9 = 0; i9 <= loop_ub; i9++) {
    h_y->data[i9] = g_y->data[(i7 + i8 * i9) - 1];
  }

  i7 = g_y->size[0];
  g_y->size[0] = h_y->size[0];
  emxEnsureCapacity_real_T1(g_y, i7);
  loop_ub = h_y->size[0];
  for (i7 = 0; i7 < loop_ub; i7++) {
    g_y->data[i7] = h_y->data[i7];
  }

  emxFree_real_T(&h_y);

  /*  remove extrapolated pieces of y */
  /* 'dynIdenf:488' y_out = y([nfact+1:len+nfact len+2*nfact+1:length(y)]); */
  d_b = (double)x->size[0] + nfact;
  if (d_b < nfact + 1.0) {
    i7 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = 0;
    emxEnsureCapacity_real_T(y, i7);
  } else {
    i7 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = (int)(d_b - (nfact + 1.0)) + 1;
    emxEnsureCapacity_real_T(y, i7);
    loop_ub = (int)(d_b - (nfact + 1.0));
    for (i7 = 0; i7 <= loop_ub; i7++) {
      y->data[y->size[0] * i7] = (nfact + 1.0) + (double)i7;
    }
  }

  d_b = ((double)x->size[0] + 2.0 * nfact) + 1.0;
  if (g_y->size[0] < d_b) {
    i7 = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    b_y->size[1] = 0;
    emxEnsureCapacity_real_T(b_y, i7);
  } else {
    i7 = g_y->size[0];
    i8 = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    b_y->size[1] = (int)((double)i7 - d_b) + 1;
    emxEnsureCapacity_real_T(b_y, i8);
    loop_ub = (int)((double)i7 - d_b);
    for (i7 = 0; i7 <= loop_ub; i7++) {
      b_y->data[b_y->size[0] * i7] = d_b + (double)i7;
    }
  }

  i7 = sp_colidx->size[0];
  sp_colidx->size[0] = y->size[1] + b_y->size[1];
  emxEnsureCapacity_int32_T(sp_colidx, i7);
  loop_ub = y->size[1];
  for (i7 = 0; i7 < loop_ub; i7++) {
    sp_colidx->data[i7] = (int)y->data[y->size[0] * i7] - 1;
  }

  loop_ub = b_y->size[1];
  for (i7 = 0; i7 < loop_ub; i7++) {
    sp_colidx->data[i7 + y->size[1]] = (int)b_y->data[b_y->size[0] * i7] - 1;
  }

  i7 = y_out->size[0];
  y_out->size[0] = y->size[1] + b_y->size[1];
  emxEnsureCapacity_real_T1(y_out, i7);
  loop_ub = y->size[1] + b_y->size[1];
  emxFree_real_T(&b_y);
  emxFree_real_T(&y);
  for (i7 = 0; i7 < loop_ub; i7++) {
    y_out->data[i7] = g_y->data[sp_colidx->data[i7]];
  }

  emxFree_real_T(&g_y);
  emxFree_int32_T(&sp_colidx);
}

/*
 * Arguments    : const emxArray_int32_T *idx
 *                emxArray_int32_T *y
 * Return Type  : void
 */
static void permuteVector(const emxArray_int32_T *idx, emxArray_int32_T *y)
{
  emxArray_int32_T *t;
  int ny;
  int k;
  int loop_ub;
  emxInit_int32_T1(&t, 1);
  ny = y->size[0];
  k = t->size[0];
  t->size[0] = y->size[0];
  emxEnsureCapacity_int32_T(t, k);
  loop_ub = y->size[0];
  for (k = 0; k < loop_ub; k++) {
    t->data[k] = y->data[k];
  }

  for (k = 0; k < ny; k++) {
    y->data[k] = t->data[idx->data[k] - 1];
  }

  emxFree_int32_T(&t);
}

/*
 * Arguments    : const emxArray_real_T *A
 * Return Type  : int
 */
static int rankFromQR(const emxArray_real_T *A)
{
  int r;
  int minmn;
  int maxmn;
  double tol;
  r = 0;
  if (A->size[0] < A->size[1]) {
    minmn = A->size[0];
    maxmn = A->size[1];
  } else {
    minmn = A->size[1];
    maxmn = A->size[0];
  }

  if (minmn > 0) {
    tol = (double)maxmn * fabs(A->data[0]) * 2.2204460492503131E-16;
    while ((r < minmn) && (!(fabs(A->data[r + A->size[0] * r]) <= tol))) {
      r++;
    }
  }

  return r;
}

/*
 * Arguments    : const emxArray_real_T *x
 *                const emxArray_real_T *y
 *                emxArray_real_T *z
 * Return Type  : void
 */
static void rdivide(const emxArray_real_T *x, const emxArray_real_T *y,
                    emxArray_real_T *z)
{
  int i20;
  int loop_ub;
  int b_loop_ub;
  int i21;
  i20 = z->size[0] * z->size[1];
  z->size[0] = x->size[0];
  z->size[1] = x->size[1];
  emxEnsureCapacity_real_T(z, i20);
  loop_ub = x->size[1];
  for (i20 = 0; i20 < loop_ub; i20++) {
    b_loop_ub = x->size[0];
    for (i21 = 0; i21 < b_loop_ub; i21++) {
      z->data[i21 + z->size[0] * i20] = x->data[i21 + x->size[0] * i20] /
        y->data[i21 + y->size[0] * i20];
    }
  }
}

/*
 * function [Rpk, Rphi] = regroup(isLevel, nparJoint, nparMinSet, a, alpha, d)
 * the function calculates the matrix that regroups the dynamic
 *  parameters into a minimum set
 *  inputs:
 * Arguments    : bool isLevel
 *                double nparJoint
 *                double nparMinSet
 *                const emxArray_real_T *b_a
 *                const emxArray_real_T *b_alpha
 *                const emxArray_real_T *b_d
 *                emxArray_real_T *b_Rpk
 *                emxArray_real_T *b_Rphi
 * Return Type  : void
 */
static void regroup(bool isLevel, double nparJoint, double nparMinSet, const
                    emxArray_real_T *b_a, const emxArray_real_T *b_alpha, const
                    emxArray_real_T *b_d, emxArray_real_T *b_Rpk,
                    emxArray_real_T *b_Rphi)
{
  emxArray_real_T *Rp1;
  int n;
  int i1;
  int j;
  emxArray_int32_T *ia;
  int fnop;
  emxArray_real_T *fnops;
  static const signed char iv1[7] = { 1, 4, 5, 6, 7, 8, 9 };

  emxArray_real_T *nopset;
  double x;
  double b_x;
  double c_x;
  double d_x;
  double e_x;
  double f_x;
  int unnamed_idx_1;
  double c_d;
  double d0;
  int i2;
  int i3;
  int loop_ub;
  emxArray_int32_T *r3;
  double dv0[3][12];
  int iv2[3];
  static const signed char nops[3] = { 1, 4, 8 };

  static const signed char iv3[3] = { -1, 0, 0 };

  emxArray_real_T *pset;
  emxArray_int32_T *ib;
  emxArray_real_T *invRp1;
  emxInit_real_T(&Rp1, 2);

  /*  parameters */
  /*  nparJoint - the number dynamic parameters per joint */
  /*  nparMinSet - the number of dynamic parameters in minimum set */
  /*  a, alpha, d - the D-H parameters */
  /*  outputs: */
  /*  Rpk - the matrix for regrouping the matrix K */
  /*  Rphi - the matrix for regrouping the dynamics parameters phi */
  /*  tau = K*phi = (K*Rpk)*(Rphi*phi) */
  /*  obtain the number of joints */
  /* 'dynIdenf:206' n = size(a,2)-1; */
  n = b_a->size[1] - 2;

  /*  initialize Rp1, Rp1 is the complete form of Rpk */
  /* 'dynIdenf:208' Rp1 = eye(nparJoint*n, nparJoint*n); */
  eye(nparJoint * ((double)b_a->size[1] - 1.0), nparJoint * ((double)b_a->size[1]
       - 1.0), Rp1);

  /*  the index of parameters to be truncated in each joint or link */
  /*  [m, mcx, mcy, mcz, Ioxx, Ioxy, Ioxz, Ioyy, Ioyz, Iozz, cf, cv] */
  /* 'dynIdenf:211' nops = [1 4 8]; */
  /*  the loop indicates the number of joints */
  /* 'dynIdenf:213' for j = n:-1:2 */
  i1 = (int)((2.0 + (-1.0 - ((double)b_a->size[1] - 1.0))) / -1.0);
  j = 0;
  emxInit_int32_T1(&ia, 1);
  while (j <= i1 - 1) {
    fnop = n - j;

    /*  the revised elements in the matrix Rp1 */
    /* 'dynIdenf:215' temp = [-1                                     0                                       0 */
    /* 'dynIdenf:216'         -a(j)                                  0                                       0 */
    /* 'dynIdenf:217'         d(j+1)*sin(alpha(j))                   sin(alpha(j))                           0 */
    /* 'dynIdenf:218'         -d(j+1)*cos(alpha(j))                  -cos(alpha(j))                          0 */
    /* 'dynIdenf:219'         -d(j+1)^2                              -2*d(j+1)                               -1 */
    /* 'dynIdenf:220'         -a(j)*d(j+1)*sin(alpha(j))             -a(j)*sin(alpha(j))                     0 */
    /* 'dynIdenf:221'         a(j)*d(j+1)*cos(alpha(j))              a(j)*cos(alpha(j))                      0 */
    /* 'dynIdenf:222'         -(a(j)^2+d(j+1)^2*cos(alpha(j))^2)     -2*d(j+1)*cos(alpha(j))^2               -cos(alpha(j))^2 */
    /* 'dynIdenf:223'         -d(j+1)^2*cos(alpha(j))*sin(alpha(j))  -2*d(j+1)*cos(alpha(j))*sin(alpha(j))   -cos(alpha(j))*sin(alpha(j)) */
    /* 'dynIdenf:224'         -(a(j)^2+d(j+1)^2*sin(alpha(j))^2)     -2*d(j+1)*sin(alpha(j))^2               -sin(alpha(j))^2 */
    /* 'dynIdenf:225'         0                                      0                                       0 */
    /* 'dynIdenf:226'         0                                      0                                       0]; */
    x = cos(b_alpha->data[fnop]);
    b_x = cos(b_alpha->data[fnop]);
    c_x = cos(b_alpha->data[fnop]);
    d_x = sin(b_alpha->data[fnop]);
    e_x = sin(b_alpha->data[fnop]);
    f_x = sin(b_alpha->data[fnop]);

    /*  revise the elements in the matrix Rp1 */
    /* 'dynIdenf:228' Rp1(nparJoint*(j-2)+1:nparJoint*(j-1), nparJoint*(j-1)+nops) = temp; */
    c_d = nparJoint * ((double)(fnop + 1) - 2.0) + 1.0;
    d0 = nparJoint * ((double)(fnop + 1) - 1.0);
    if (c_d > d0) {
      i2 = 0;
      i3 = 0;
    } else {
      i2 = (int)c_d - 1;
      i3 = (int)d0;
    }

    unnamed_idx_1 = ia->size[0];
    ia->size[0] = i3 - i2;
    emxEnsureCapacity_int32_T(ia, unnamed_idx_1);
    loop_ub = i3 - i2;
    for (i3 = 0; i3 < loop_ub; i3++) {
      ia->data[i3] = i2 + i3;
    }

    c_d = nparJoint * ((double)(fnop + 1) - 1.0);
    for (i2 = 0; i2 < 3; i2++) {
      iv2[i2] = (int)(c_d + (double)nops[i2]) - 1;
      dv0[i2][0] = iv3[i2];
    }

    dv0[0][1] = -b_a->data[fnop];
    dv0[1][1] = 0.0;
    dv0[2][1] = 0.0;
    dv0[0][2] = b_d->data[fnop + 1] * sin(b_alpha->data[fnop]);
    dv0[1][2] = sin(b_alpha->data[fnop]);
    dv0[2][2] = 0.0;
    dv0[0][3] = -b_d->data[fnop + 1] * cos(b_alpha->data[fnop]);
    dv0[1][3] = -cos(b_alpha->data[fnop]);
    dv0[2][3] = 0.0;
    dv0[0][4] = -(b_d->data[fnop + 1] * b_d->data[fnop + 1]);
    dv0[1][4] = -2.0 * b_d->data[fnop + 1];
    dv0[2][4] = -1.0;
    dv0[0][5] = -b_a->data[fnop] * b_d->data[fnop + 1] * sin(b_alpha->data[fnop]);
    dv0[1][5] = -b_a->data[fnop] * sin(b_alpha->data[fnop]);
    dv0[2][5] = 0.0;
    dv0[0][6] = b_a->data[fnop] * b_d->data[fnop + 1] * cos(b_alpha->data[fnop]);
    dv0[1][6] = b_a->data[fnop] * cos(b_alpha->data[fnop]);
    dv0[2][6] = 0.0;
    dv0[0][7] = -(b_a->data[fnop] * b_a->data[fnop] + b_d->data[fnop + 1] *
                  b_d->data[fnop + 1] * (x * x));
    dv0[1][7] = -2.0 * b_d->data[fnop + 1] * (b_x * b_x);
    dv0[2][7] = -(c_x * c_x);
    dv0[0][8] = -(b_d->data[fnop + 1] * b_d->data[fnop + 1]) * cos(b_alpha->
      data[fnop]) * sin(b_alpha->data[fnop]);
    dv0[1][8] = -2.0 * b_d->data[fnop + 1] * cos(b_alpha->data[fnop]) * sin
      (b_alpha->data[fnop]);
    dv0[2][8] = -cos(b_alpha->data[fnop]) * sin(b_alpha->data[fnop]);
    dv0[0][9] = -(b_a->data[fnop] * b_a->data[fnop] + b_d->data[fnop + 1] *
                  b_d->data[fnop + 1] * (d_x * d_x));
    dv0[1][9] = -2.0 * b_d->data[fnop + 1] * (e_x * e_x);
    dv0[2][9] = -(f_x * f_x);
    for (i2 = 0; i2 < 3; i2++) {
      dv0[i2][10] = 0.0;
      dv0[i2][11] = 0.0;
    }

    unnamed_idx_1 = ia->size[0];
    for (i2 = 0; i2 < 3; i2++) {
      for (i3 = 0; i3 < unnamed_idx_1; i3++) {
        Rp1->data[ia->data[i3] + Rp1->size[0] * iv2[i2]] = (&dv0[0][0])[i3 +
          unnamed_idx_1 * i2];
      }
    }

    /* 'dynIdenf:229' Rp1(nparJoint*(j-1)+5, nparJoint*(j-1)+8) = 1; */
    Rp1->data[((int)(nparJoint * ((double)(fnop + 1) - 1.0) + 5.0) + Rp1->size[0]
               * ((int)(nparJoint * ((double)(fnop + 1) - 1.0) + 8.0) - 1)) - 1]
      = 1.0;
    j++;
  }

  /*  if the base of the robot is horizontal, the second and third */
  /*  dynamic parameters, i.e. mcx, mcy can be truncated */
  /* 'dynIdenf:233' if isLevel ~= 0 */
  emxInit_real_T(&fnops, 2);
  if (isLevel) {
    /* 'dynIdenf:234' fnops = [1 2 3 4 5 6 7 8 9]; */
    i1 = fnops->size[0] * fnops->size[1];
    fnops->size[0] = 1;
    fnops->size[1] = 9;
    emxEnsureCapacity_real_T(fnops, i1);
    for (i1 = 0; i1 < 9; i1++) {
      fnops->data[fnops->size[0] * i1] = 1.0 + (double)i1;
    }

    /*  otherwise, mcx, mcy must be considered */
  } else {
    /* 'dynIdenf:236' else */
    /* 'dynIdenf:237' fnops = [1 4 5 6 7 8 9]; */
    i1 = fnops->size[0] * fnops->size[1];
    fnops->size[0] = 1;
    fnops->size[1] = 7;
    emxEnsureCapacity_real_T(fnops, i1);
    for (i1 = 0; i1 < 7; i1++) {
      fnops->data[fnops->size[0] * i1] = iv1[i1];
    }
  }

  emxInit_real_T(&nopset, 2);

  /*  the number of truncated parameters in the first joint or link */
  /* 'dynIdenf:240' nop = numel(nops); */
  /*  the number of truncated parameters in each of other joint */
  /* 'dynIdenf:242' fnop = numel(fnops); */
  fnop = fnops->size[1];

  /*  initialize the index set of all truncated parameters */
  /* 'dynIdenf:244' nopset = zeros(1, fnop+nop*(n-1)); */
  unnamed_idx_1 = (int)((double)fnops->size[1] + 3.0 * (((double)b_a->size[1] -
    1.0) - 1.0));
  i1 = nopset->size[0] * nopset->size[1];
  nopset->size[0] = 1;
  nopset->size[1] = unnamed_idx_1;
  emxEnsureCapacity_real_T(nopset, i1);
  for (i1 = 0; i1 < unnamed_idx_1; i1++) {
    nopset->data[nopset->size[0] * i1] = 0.0;
  }

  /*  assign the index set of all truncated parameters */
  /* 'dynIdenf:246' nopset(1:fnop) = fnops; */
  loop_ub = fnops->size[1];
  for (i1 = 0; i1 < loop_ub; i1++) {
    nopset->data[i1] = fnops->data[fnops->size[0] * i1];
  }

  /* 'dynIdenf:247' for j = 2:n */
  j = 0;
  emxInit_int32_T(&r3, 2);
  while (j <= n - 1) {
    /* 'dynIdenf:248' nopset(fnop+nop*(j-2)+1:fnop+nop*(j-1)) = nparJoint*(j-1)+nops; */
    c_d = ((double)fnop + 3.0 * ((2.0 + (double)j) - 2.0)) + 1.0;
    d0 = (double)fnop + 3.0 * ((2.0 + (double)j) - 1.0);
    if (c_d > d0) {
      i1 = 0;
      i2 = 0;
    } else {
      i1 = (int)c_d - 1;
      i2 = (int)d0;
    }

    i3 = r3->size[0] * r3->size[1];
    r3->size[0] = 1;
    r3->size[1] = i2 - i1;
    emxEnsureCapacity_int32_T1(r3, i3);
    loop_ub = i2 - i1;
    for (i2 = 0; i2 < loop_ub; i2++) {
      r3->data[r3->size[0] * i2] = i1 + i2;
    }

    c_d = nparJoint * ((2.0 + (double)j) - 1.0);
    loop_ub = r3->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      nopset->data[r3->data[r3->size[0] * i1]] = c_d + (double)nops[(*(int (*)[2])
        r3->size)[0] * i1];
    }

    j++;
  }

  emxFree_int32_T(&r3);

  /*  the index of the minimum set of dynamic parameters */
  /* 'dynIdenf:251' pset = setdiff(1:nparJoint*n, nopset); */
  c_d = nparJoint * ((double)b_a->size[1] - 1.0);
  if (c_d < 1.0) {
    i1 = fnops->size[0] * fnops->size[1];
    fnops->size[0] = 1;
    fnops->size[1] = 0;
    emxEnsureCapacity_real_T(fnops, i1);
  } else {
    i1 = fnops->size[0] * fnops->size[1];
    fnops->size[0] = 1;
    fnops->size[1] = (int)floor(c_d - 1.0) + 1;
    emxEnsureCapacity_real_T(fnops, i1);
    loop_ub = (int)floor(c_d - 1.0);
    for (i1 = 0; i1 <= loop_ub; i1++) {
      fnops->data[fnops->size[0] * i1] = 1.0 + (double)i1;
    }
  }

  emxInit_real_T(&pset, 2);
  emxInit_int32_T1(&ib, 1);
  do_vectors(fnops, nopset, pset, ia, ib);

  /*  extract the pset columns of Rp1 and assign them to Rpk */
  /* 'dynIdenf:253' Rpk = zeros(nparJoint*n,nparMinSet); */
  unnamed_idx_1 = (int)(nparJoint * ((double)b_a->size[1] - 1.0));
  i1 = b_Rpk->size[0] * b_Rpk->size[1];
  b_Rpk->size[0] = unnamed_idx_1;
  b_Rpk->size[1] = (int)nparMinSet;
  emxEnsureCapacity_real_T(b_Rpk, i1);
  loop_ub = (int)nparMinSet;
  emxFree_int32_T(&ib);
  emxFree_int32_T(&ia);
  emxFree_real_T(&nopset);
  emxFree_real_T(&fnops);
  for (i1 = 0; i1 < loop_ub; i1++) {
    for (i2 = 0; i2 < unnamed_idx_1; i2++) {
      b_Rpk->data[i2 + b_Rpk->size[0] * i1] = 0.0;
    }
  }

  emxInit_real_T(&invRp1, 2);

  /* 'dynIdenf:254' Rpk(:,:) = Rp1(:,pset); */
  loop_ub = Rp1->size[0];
  i1 = invRp1->size[0] * invRp1->size[1];
  invRp1->size[0] = loop_ub;
  invRp1->size[1] = pset->size[1];
  emxEnsureCapacity_real_T(invRp1, i1);
  unnamed_idx_1 = pset->size[1];
  for (i1 = 0; i1 < unnamed_idx_1; i1++) {
    for (i2 = 0; i2 < loop_ub; i2++) {
      invRp1->data[i2 + invRp1->size[0] * i1] = Rp1->data[i2 + Rp1->size[0] *
        ((int)pset->data[pset->size[0] * i1] - 1)];
    }
  }

  loop_ub = invRp1->size[1];
  for (i1 = 0; i1 < loop_ub; i1++) {
    unnamed_idx_1 = invRp1->size[0];
    for (i2 = 0; i2 < unnamed_idx_1; i2++) {
      b_Rpk->data[i2 + b_Rpk->size[0] * i1] = invRp1->data[i2 + invRp1->size[0] *
        i1];
    }
  }

  /*  inverse of Rp1 */
  /* 'dynIdenf:256' invRp1 = eye(nparJoint*n, nparJoint*n) / Rp1; */
  eye(nparJoint * ((double)b_a->size[1] - 1.0), nparJoint * ((double)b_a->size[1]
       - 1.0), invRp1);
  mrdivide(invRp1, Rp1);

  /*  extract the pset rows of invRp1 and assign them to Rphi */
  /* 'dynIdenf:258' Rphi = invRp1(pset,:); */
  loop_ub = invRp1->size[1];
  i1 = b_Rphi->size[0] * b_Rphi->size[1];
  b_Rphi->size[0] = pset->size[1];
  b_Rphi->size[1] = loop_ub;
  emxEnsureCapacity_real_T(b_Rphi, i1);
  emxFree_real_T(&Rp1);
  for (i1 = 0; i1 < loop_ub; i1++) {
    unnamed_idx_1 = pset->size[1];
    for (i2 = 0; i2 < unnamed_idx_1; i2++) {
      b_Rphi->data[i2 + b_Rphi->size[0] * i1] = invRp1->data[((int)pset->
        data[pset->size[0] * i2] + invRp1->size[0] * i1) - 1];
    }
  }

  emxFree_real_T(&invRp1);
  emxFree_real_T(&pset);
}

/*
 * Arguments    : const emxArray_real_T *b_a
 *                double varargin_1
 *                emxArray_real_T *b
 * Return Type  : void
 */
static void repmat(const emxArray_real_T *b_a, double varargin_1,
                   emxArray_real_T *b)
{
  int outsize_idx_1;
  int t;
  outsize_idx_1 = b_a->size[1];
  t = b->size[0] * b->size[1];
  b->size[0] = (int)varargin_1;
  b->size[1] = outsize_idx_1;
  emxEnsureCapacity_real_T(b, t);
  if (!((b->size[0] == 0) || (b->size[1] == 0))) {
    for (outsize_idx_1 = 0; outsize_idx_1 < b_a->size[1]; outsize_idx_1++) {
      for (t = 0; t < (int)varargin_1; t++) {
        b->data[t + b->size[0] * outsize_idx_1] = b_a->data[b_a->size[0] *
          outsize_idx_1];
      }
    }
  }
}

/*
 * Arguments    : const emxArray_real_T *x
 * Return Type  : double
 */
static double rms(const emxArray_real_T *x)
{
  double y;
  emxArray_real_T *b_x;
  int i11;
  int loop_ub;
  emxInit_real_T1(&b_x, 1);
  i11 = b_x->size[0];
  b_x->size[0] = x->size[0];
  emxEnsureCapacity_real_T1(b_x, i11);
  loop_ub = x->size[0];
  for (i11 = 0; i11 < loop_ub; i11++) {
    b_x->data[i11] = x->data[i11] * x->data[i11];
  }

  y = sqrt(combineVectorElements(b_x) / (double)b_x->size[0]);
  emxFree_real_T(&b_x);
  return y;
}

/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_hypotd(double u0, double u1)
{
  double y;
  double b_a;
  double b;
  b_a = fabs(u0);
  b = fabs(u1);
  if (b_a < b) {
    b_a /= b;
    y = b * sqrt(b_a * b_a + 1.0);
  } else if (b_a > b) {
    b /= b_a;
    y = b_a * sqrt(b * b + 1.0);
  } else {
    y = b_a * 1.4142135623730951;
  }

  return y;
}

/*
 * Arguments    : int *k
 *                const emxArray_real_T *x
 * Return Type  : double
 */
static double skip_to_last_equal_value(int *k, const emxArray_real_T *x)
{
  double xk;
  bool exitg1;
  double absxk;
  int exponent;
  bool p;
  xk = x->data[*k - 1];
  exitg1 = false;
  while ((!exitg1) && (*k < x->size[1])) {
    absxk = fabs(xk / 2.0);
    if (absxk <= 2.2250738585072014E-308) {
      absxk = 4.94065645841247E-324;
    } else {
      frexp(absxk, &exponent);
      absxk = ldexp(1.0, exponent - 53);
    }

    p = (fabs(xk - x->data[*k]) < absxk);
    if (p) {
      (*k)++;
    } else {
      exitg1 = true;
    }
  }

  return xk;
}

/*
 * Arguments    : emxArray_real_T *x
 * Return Type  : void
 */
static void sort(emxArray_real_T *x)
{
  int dim;
  emxArray_real_T *vwork;
  int i31;
  int pagesize;
  int vstride;
  int k;
  int npages;
  emxArray_int32_T *b_vwork;
  int pageoffset;
  int j;
  int idx0;
  dim = 1;
  if (x->size[0] != 1) {
    dim = 0;
  }

  emxInit_real_T1(&vwork, 1);
  i31 = x->size[dim];
  pagesize = vwork->size[0];
  vwork->size[0] = i31;
  emxEnsureCapacity_real_T1(vwork, pagesize);
  vstride = 1;
  k = 1;
  while (k <= dim) {
    pagesize = x->size[0];
    vstride *= pagesize;
    k = 2;
  }

  npages = 1;
  k = dim + 2;
  while (k < 3) {
    pagesize = x->size[1];
    npages *= pagesize;
    k = 3;
  }

  pagesize = i31 * vstride;
  dim = 1;
  emxInit_int32_T1(&b_vwork, 1);
  while (dim <= npages) {
    pageoffset = (dim - 1) * pagesize;
    for (j = 0; j < vstride; j++) {
      idx0 = pageoffset + j;
      for (k = 0; k < i31; k++) {
        vwork->data[k] = x->data[idx0 + k * vstride];
      }

      sortIdx(vwork, b_vwork);
      for (k = 0; k < i31; k++) {
        x->data[idx0 + k * vstride] = vwork->data[k];
      }
    }

    dim++;
  }

  emxFree_int32_T(&b_vwork);
  emxFree_real_T(&vwork);
}

/*
 * Arguments    : emxArray_real_T *x
 *                emxArray_int32_T *idx
 * Return Type  : void
 */
static void sortIdx(emxArray_real_T *x, emxArray_int32_T *idx)
{
  int b_x[2];
  int i32;
  unsigned int uv2[2];
  int i;
  emxArray_int32_T *b_idx;
  emxArray_real_T *c_x;
  int nBlocks;
  emxArray_int32_T *iwork;
  double x4[4];
  int ib;
  int idx4[4];
  int q;
  emxArray_real_T *xwork;
  int k;
  signed char perm[4];
  int b;
  int offset;
  int b_b;
  int p;
  int bLen;
  int bLen2;
  int nPairs;
  int b_iwork[256];
  double b_xwork[256];
  int exitg1;
  b_x[0] = x->size[0];
  b_x[1] = 1;
  for (i32 = 0; i32 < 2; i32++) {
    uv2[i32] = (unsigned int)b_x[i32];
  }

  i32 = idx->size[0];
  idx->size[0] = (int)uv2[0];
  emxEnsureCapacity_int32_T(idx, i32);
  i = (int)uv2[0];
  for (i32 = 0; i32 < i; i32++) {
    idx->data[i32] = 0;
  }

  if (x->size[0] != 0) {
    emxInit_int32_T1(&b_idx, 1);
    i = (int)uv2[0];
    i32 = b_idx->size[0];
    b_idx->size[0] = (int)uv2[0];
    emxEnsureCapacity_int32_T(b_idx, i32);
    for (i32 = 0; i32 < i; i32++) {
      b_idx->data[i32] = 0;
    }

    emxInit_real_T1(&c_x, 1);
    i = x->size[0];
    i32 = c_x->size[0];
    c_x->size[0] = i;
    emxEnsureCapacity_real_T1(c_x, i32);
    for (i32 = 0; i32 < i; i32++) {
      c_x->data[i32] = x->data[i32];
    }

    i32 = x->size[0];
    nBlocks = x->size[0];
    for (i = 0; i < 4; i++) {
      x4[i] = 0.0;
      idx4[i] = 0;
    }

    emxInit_int32_T1(&iwork, 1);
    ib = b_idx->size[0];
    q = iwork->size[0];
    iwork->size[0] = ib;
    emxEnsureCapacity_int32_T(iwork, q);
    for (q = 0; q < ib; q++) {
      iwork->data[q] = 0;
    }

    emxInit_real_T1(&xwork, 1);
    i = x->size[0];
    q = xwork->size[0];
    xwork->size[0] = i;
    emxEnsureCapacity_real_T1(xwork, q);
    for (q = 0; q < i; q++) {
      xwork->data[q] = 0.0;
    }

    ib = 0;
    for (k = 0; k < nBlocks; k++) {
      ib++;
      idx4[ib - 1] = k + 1;
      x4[ib - 1] = c_x->data[k];
      if (ib == 4) {
        if (x4[0] >= x4[1]) {
          i = 1;
          q = 2;
        } else {
          i = 2;
          q = 1;
        }

        if (x4[2] >= x4[3]) {
          ib = 3;
          p = 4;
        } else {
          ib = 4;
          p = 3;
        }

        if (x4[i - 1] >= x4[ib - 1]) {
          if (x4[q - 1] >= x4[ib - 1]) {
            perm[0] = (signed char)i;
            perm[1] = (signed char)q;
            perm[2] = (signed char)ib;
            perm[3] = (signed char)p;
          } else if (x4[q - 1] >= x4[p - 1]) {
            perm[0] = (signed char)i;
            perm[1] = (signed char)ib;
            perm[2] = (signed char)q;
            perm[3] = (signed char)p;
          } else {
            perm[0] = (signed char)i;
            perm[1] = (signed char)ib;
            perm[2] = (signed char)p;
            perm[3] = (signed char)q;
          }
        } else if (x4[i - 1] >= x4[p - 1]) {
          if (x4[q - 1] >= x4[p - 1]) {
            perm[0] = (signed char)ib;
            perm[1] = (signed char)i;
            perm[2] = (signed char)q;
            perm[3] = (signed char)p;
          } else {
            perm[0] = (signed char)ib;
            perm[1] = (signed char)i;
            perm[2] = (signed char)p;
            perm[3] = (signed char)q;
          }
        } else {
          perm[0] = (signed char)ib;
          perm[1] = (signed char)p;
          perm[2] = (signed char)i;
          perm[3] = (signed char)q;
        }

        b_idx->data[k - 3] = idx4[perm[0] - 1];
        b_idx->data[k - 2] = idx4[perm[1] - 1];
        b_idx->data[k - 1] = idx4[perm[2] - 1];
        b_idx->data[k] = idx4[perm[3] - 1];
        c_x->data[k - 3] = x4[perm[0] - 1];
        c_x->data[k - 2] = x4[perm[1] - 1];
        c_x->data[k - 1] = x4[perm[2] - 1];
        c_x->data[k] = x4[perm[3] - 1];
        ib = 0;
      }
    }

    if (ib > 0) {
      for (i = 0; i < 4; i++) {
        perm[i] = 0;
      }

      if (ib == 1) {
        perm[0] = 1;
      } else if (ib == 2) {
        if (x4[0] >= x4[1]) {
          perm[0] = 1;
          perm[1] = 2;
        } else {
          perm[0] = 2;
          perm[1] = 1;
        }
      } else if (x4[0] >= x4[1]) {
        if (x4[1] >= x4[2]) {
          perm[0] = 1;
          perm[1] = 2;
          perm[2] = 3;
        } else if (x4[0] >= x4[2]) {
          perm[0] = 1;
          perm[1] = 3;
          perm[2] = 2;
        } else {
          perm[0] = 3;
          perm[1] = 1;
          perm[2] = 2;
        }
      } else if (x4[0] >= x4[2]) {
        perm[0] = 2;
        perm[1] = 1;
        perm[2] = 3;
      } else if (x4[1] >= x4[2]) {
        perm[0] = 2;
        perm[1] = 3;
        perm[2] = 1;
      } else {
        perm[0] = 3;
        perm[1] = 2;
        perm[2] = 1;
      }

      for (k = 1; k <= ib; k++) {
        b_idx->data[((nBlocks - ib) + k) - 1] = idx4[perm[k - 1] - 1];
        c_x->data[((nBlocks - ib) + k) - 1] = x4[perm[k - 1] - 1];
      }
    }

    i = 2;
    if (i32 > 1) {
      if (i32 >= 256) {
        nBlocks = i32 >> 8;
        for (b = 1; b <= nBlocks; b++) {
          offset = (b - 1) << 8;
          for (b_b = 0; b_b < 6; b_b++) {
            bLen = 1 << (b_b + 2);
            bLen2 = bLen << 1;
            nPairs = 256 >> (b_b + 3);
            for (k = 1; k <= nPairs; k++) {
              ib = offset + (k - 1) * bLen2;
              for (i = 0; i < bLen2; i++) {
                b_iwork[i] = b_idx->data[ib + i];
                b_xwork[i] = c_x->data[ib + i];
              }

              p = 0;
              q = bLen;
              i = ib - 1;
              do {
                exitg1 = 0;
                i++;
                if (b_xwork[p] >= b_xwork[q]) {
                  b_idx->data[i] = b_iwork[p];
                  c_x->data[i] = b_xwork[p];
                  if (p + 1 < bLen) {
                    p++;
                  } else {
                    exitg1 = 1;
                  }
                } else {
                  b_idx->data[i] = b_iwork[q];
                  c_x->data[i] = b_xwork[q];
                  if (q + 1 < bLen2) {
                    q++;
                  } else {
                    i = (i - p) + 1;
                    while (p + 1 <= bLen) {
                      b_idx->data[i + p] = b_iwork[p];
                      c_x->data[i + p] = b_xwork[p];
                      p++;
                    }

                    exitg1 = 1;
                  }
                }
              } while (exitg1 == 0);
            }
          }
        }

        i = nBlocks << 8;
        ib = i32 - i;
        if (ib > 0) {
          merge_block(b_idx, c_x, i, ib, 2, iwork, xwork);
        }

        i = 8;
      }

      merge_block(b_idx, c_x, 0, i32, i, iwork, xwork);
    }

    emxFree_real_T(&xwork);
    emxFree_int32_T(&iwork);
    i = b_idx->size[0];
    for (i32 = 0; i32 < i; i32++) {
      idx->data[i32] = b_idx->data[i32];
    }

    emxFree_int32_T(&b_idx);
    i = c_x->size[0];
    for (i32 = 0; i32 < i; i32++) {
      x->data[i32] = c_x->data[i32];
    }

    emxFree_real_T(&c_x);
  }
}

/*
 * Arguments    : const emxArray_real_T *varargin_1
 *                const emxArray_real_T *varargin_2
 *                const emxArray_real_T *varargin_3
 *                coder_internal_sparse *y
 * Return Type  : void
 */
static void sparse(const emxArray_real_T *varargin_1, const emxArray_real_T
                   *varargin_2, const emxArray_real_T *varargin_3,
                   coder_internal_sparse *y)
{
  emxArray_int32_T *ridxInt;
  emxArray_int32_T *cidxInt;
  emxArray_int32_T *sortedIndices;
  int nc;
  int i10;
  int numalloc;
  cell_wrap_4 tunableEnvironment[2];
  cell_wrap_4 this_tunableEnvironment[2];
  int thism;
  int thisn;
  int c;
  double val;
  emxInit_int32_T1(&ridxInt, 1);
  emxInit_int32_T1(&cidxInt, 1);
  emxInit_int32_T1(&sortedIndices, 1);
  nc = varargin_2->size[1];
  assertValidIndexArg(varargin_1, ridxInt);
  assertValidIndexArg(varargin_2, cidxInt);
  i10 = sortedIndices->size[0];
  sortedIndices->size[0] = varargin_2->size[1];
  emxEnsureCapacity_int32_T(sortedIndices, i10);
  for (numalloc = 1; numalloc <= nc; numalloc++) {
    sortedIndices->data[numalloc - 1] = numalloc;
  }

  emxInitMatrix_cell_wrap_4(tunableEnvironment);
  i10 = tunableEnvironment[0].f1->size[0];
  tunableEnvironment[0].f1->size[0] = cidxInt->size[0];
  emxEnsureCapacity_int32_T(tunableEnvironment[0].f1, i10);
  numalloc = cidxInt->size[0];
  for (i10 = 0; i10 < numalloc; i10++) {
    tunableEnvironment[0].f1->data[i10] = cidxInt->data[i10];
  }

  i10 = tunableEnvironment[1].f1->size[0];
  tunableEnvironment[1].f1->size[0] = ridxInt->size[0];
  emxEnsureCapacity_int32_T(tunableEnvironment[1].f1, i10);
  numalloc = ridxInt->size[0];
  for (i10 = 0; i10 < numalloc; i10++) {
    tunableEnvironment[1].f1->data[i10] = ridxInt->data[i10];
  }

  emxInitMatrix_cell_wrap_4(this_tunableEnvironment);
  for (i10 = 0; i10 < 2; i10++) {
    emxCopyStruct_cell_wrap_4(&this_tunableEnvironment[i10],
      &tunableEnvironment[i10]);
  }

  emxFreeMatrix_cell_wrap_4(tunableEnvironment);
  introsort(sortedIndices, cidxInt->size[0], this_tunableEnvironment);
  permuteVector(sortedIndices, cidxInt);
  permuteVector(sortedIndices, ridxInt);
  emxFreeMatrix_cell_wrap_4(this_tunableEnvironment);
  if ((ridxInt->size[0] == 0) || (cidxInt->size[0] == 0)) {
    thism = 0;
    thisn = 0;
  } else {
    thism = ridxInt->data[0];
    for (numalloc = 1; numalloc < ridxInt->size[0]; numalloc++) {
      if (thism < ridxInt->data[numalloc]) {
        thism = ridxInt->data[numalloc];
      }
    }

    thisn = cidxInt->data[cidxInt->size[0] - 1];
  }

  y->m = thism;
  y->n = thisn;
  if (varargin_2->size[1] >= 1) {
    numalloc = varargin_2->size[1];
  } else {
    numalloc = 1;
  }

  i10 = y->d->size[0];
  y->d->size[0] = numalloc;
  emxEnsureCapacity_real_T1(y->d, i10);
  for (i10 = 0; i10 < numalloc; i10++) {
    y->d->data[i10] = 0.0;
  }

  i10 = y->colidx->size[0];
  y->colidx->size[0] = thisn + 1;
  emxEnsureCapacity_int32_T(y->colidx, i10);
  y->colidx->data[0] = 1;
  i10 = y->rowidx->size[0];
  y->rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(y->rowidx, i10);
  for (i10 = 0; i10 < numalloc; i10++) {
    y->rowidx->data[i10] = 0;
  }

  thism = 0;
  for (c = 1; c <= thisn; c++) {
    while ((thism + 1 <= nc) && (cidxInt->data[thism] == c)) {
      y->rowidx->data[thism] = ridxInt->data[thism];
      thism++;
    }

    y->colidx->data[c] = thism + 1;
  }

  emxFree_int32_T(&cidxInt);
  emxFree_int32_T(&ridxInt);
  for (numalloc = 0; numalloc < nc; numalloc++) {
    y->d->data[numalloc] = varargin_3->data[sortedIndices->data[numalloc] - 1];
  }

  emxFree_int32_T(&sortedIndices);
  thism = 1;
  i10 = y->colidx->size[0] - 1;
  for (c = 1; c <= i10; c++) {
    numalloc = y->colidx->data[c - 1];
    y->colidx->data[c - 1] = thism;
    while (numalloc < y->colidx->data[c]) {
      val = 0.0;
      thisn = y->rowidx->data[numalloc - 1];
      while ((numalloc < y->colidx->data[c]) && (y->rowidx->data[numalloc - 1] ==
              thisn)) {
        val += y->d->data[numalloc - 1];
        numalloc++;
      }

      if (val != 0.0) {
        y->d->data[thism - 1] = val;
        y->rowidx->data[thism - 1] = thisn;
        thism++;
      }
    }
  }

  y->colidx->data[y->colidx->size[0] - 1] = thism;
}

/*
 * Arguments    : const emxArray_real_T *this_d
 *                const emxArray_int32_T *this_colidx
 *                const emxArray_int32_T *this_rowidx
 *                int this_m
 *                int this_n
 *                emxArray_real_T *y
 * Return Type  : void
 */
static void sparse_full(const emxArray_real_T *this_d, const emxArray_int32_T
  *this_colidx, const emxArray_int32_T *this_rowidx, int this_m, int this_n,
  emxArray_real_T *y)
{
  int c;
  int idx;
  c = y->size[0] * y->size[1];
  y->size[0] = this_m;
  y->size[1] = this_n;
  emxEnsureCapacity_real_T(y, c);
  for (c = 0; c < this_n; c++) {
    for (idx = 0; idx < this_m; idx++) {
      y->data[idx + y->size[0] * c] = 0.0;
    }
  }

  for (c = 1; c <= this_n; c++) {
    for (idx = this_colidx->data[c - 1]; idx < this_colidx->data[c]; idx++) {
      y->data[(this_rowidx->data[idx - 1] + y->size[0] * (c - 1)) - 1] =
        this_d->data[idx - 1];
    }
  }
}

/*
 * Arguments    : const emxArray_real_T *A
 *                emxArray_real_T *U
 *                emxArray_real_T *s
 *                emxArray_real_T *V
 * Return Type  : void
 */
static void svd(const emxArray_real_T *A, emxArray_real_T *U, emxArray_real_T *s,
                emxArray_real_T *V)
{
  emxArray_real_T *b_A;
  int m;
  int ns;
  int n;
  int qs;
  int p;
  int iter;
  int minnp;
  emxArray_real_T *b_s;
  emxArray_real_T *e;
  emxArray_real_T *work;
  int nrt;
  int nct;
  int q;
  int nmq;
  bool apply_transform;
  double ztest0;
  int mm;
  double ztest;
  double snorm;
  bool exitg1;
  double f;
  double scale;
  double sqds;
  double b;
  emxInit_real_T(&b_A, 2);
  m = b_A->size[0] * b_A->size[1];
  b_A->size[0] = A->size[0];
  b_A->size[1] = A->size[1];
  emxEnsureCapacity_real_T(b_A, m);
  ns = A->size[1];
  for (m = 0; m < ns; m++) {
    qs = A->size[0];
    for (iter = 0; iter < qs; iter++) {
      b_A->data[iter + b_A->size[0] * m] = A->data[iter + A->size[0] * m];
    }
  }

  n = A->size[0];
  p = A->size[1];
  qs = A->size[0] + 1;
  ns = A->size[1];
  if (qs < ns) {
    ns = qs;
  }

  qs = A->size[0];
  minnp = A->size[1];
  if (qs < minnp) {
    minnp = qs;
  }

  emxInit_real_T1(&b_s, 1);
  m = b_s->size[0];
  b_s->size[0] = ns;
  emxEnsureCapacity_real_T1(b_s, m);
  for (m = 0; m < ns; m++) {
    b_s->data[m] = 0.0;
  }

  emxInit_real_T1(&e, 1);
  qs = A->size[1];
  m = e->size[0];
  e->size[0] = qs;
  emxEnsureCapacity_real_T1(e, m);
  for (m = 0; m < qs; m++) {
    e->data[m] = 0.0;
  }

  emxInit_real_T1(&work, 1);
  qs = A->size[0];
  m = work->size[0];
  work->size[0] = qs;
  emxEnsureCapacity_real_T1(work, m);
  for (m = 0; m < qs; m++) {
    work->data[m] = 0.0;
  }

  qs = A->size[0];
  ns = A->size[0];
  m = U->size[0] * U->size[1];
  U->size[0] = qs;
  U->size[1] = ns;
  emxEnsureCapacity_real_T(U, m);
  for (m = 0; m < ns; m++) {
    for (iter = 0; iter < qs; iter++) {
      U->data[iter + U->size[0] * m] = 0.0;
    }
  }

  qs = A->size[1];
  ns = A->size[1];
  m = V->size[0] * V->size[1];
  V->size[0] = qs;
  V->size[1] = ns;
  emxEnsureCapacity_real_T(V, m);
  for (m = 0; m < ns; m++) {
    for (iter = 0; iter < qs; iter++) {
      V->data[iter + V->size[0] * m] = 0.0;
    }
  }

  if ((A->size[0] == 0) || (A->size[1] == 0)) {
    m = A->size[0];
    for (ns = 0; ns < m; ns++) {
      U->data[ns + U->size[0] * ns] = 1.0;
    }

    m = A->size[1];
    for (ns = 0; ns < m; ns++) {
      V->data[ns + V->size[0] * ns] = 1.0;
    }
  } else {
    if (A->size[1] > 2) {
      qs = A->size[1] - 2;
    } else {
      qs = 0;
    }

    nrt = A->size[0];
    if (qs < nrt) {
      nrt = qs;
    }

    if (A->size[0] > 1) {
      qs = A->size[0] - 1;
    } else {
      qs = 0;
    }

    nct = A->size[1];
    if (qs < nct) {
      nct = qs;
    }

    if (nct > nrt) {
      m = nct;
    } else {
      m = nrt;
    }

    for (q = 0; q < m; q++) {
      qs = q + n * q;
      nmq = n - q;
      apply_transform = false;
      if (q + 1 <= nct) {
        ztest0 = xnrm2(nmq, b_A, qs + 1);
        if (ztest0 > 0.0) {
          apply_transform = true;
          if (b_A->data[qs] < 0.0) {
            ztest0 = -ztest0;
          }

          b_s->data[q] = ztest0;
          if (fabs(b_s->data[q]) >= 1.0020841800044864E-292) {
            xscal(nmq, 1.0 / b_s->data[q], b_A, qs + 1);
          } else {
            iter = qs + nmq;
            for (ns = qs; ns < iter; ns++) {
              b_A->data[ns] /= b_s->data[q];
            }
          }

          b_A->data[qs]++;
          b_s->data[q] = -b_s->data[q];
        } else {
          b_s->data[q] = 0.0;
        }
      }

      for (mm = q + 1; mm < p; mm++) {
        ns = q + n * mm;
        if (apply_transform) {
          ztest0 = -(xdotc(nmq, b_A, qs + 1, b_A, ns + 1) / b_A->data[q +
                     b_A->size[0] * q]);
          xaxpy(nmq, ztest0, qs + 1, b_A, ns + 1);
        }

        e->data[mm] = b_A->data[ns];
      }

      if (q + 1 <= nct) {
        for (ns = q; ns < n; ns++) {
          U->data[ns + U->size[0] * q] = b_A->data[ns + b_A->size[0] * q];
        }
      }

      if (q + 1 <= nrt) {
        iter = p - q;
        ztest0 = b_xnrm2(iter - 1, e, q + 2);
        if (ztest0 == 0.0) {
          e->data[q] = 0.0;
        } else {
          if (e->data[q + 1] < 0.0) {
            ztest0 = -ztest0;
          }

          e->data[q] = ztest0;
          ztest0 = e->data[q];
          if (fabs(e->data[q]) >= 1.0020841800044864E-292) {
            ztest0 = 1.0 / e->data[q];
            iter += q;
            for (ns = q + 1; ns < iter; ns++) {
              e->data[ns] *= ztest0;
            }
          } else {
            iter += q;
            for (ns = q + 1; ns < iter; ns++) {
              e->data[ns] /= ztest0;
            }
          }

          e->data[q + 1]++;
          e->data[q] = -e->data[q];
          if (q + 2 <= n) {
            for (ns = q + 1; ns < n; ns++) {
              work->data[ns] = 0.0;
            }

            for (mm = q + 1; mm < p; mm++) {
              b_xaxpy(nmq - 1, e->data[mm], b_A, (q + n * mm) + 2, work, q + 2);
            }

            for (mm = q + 1; mm < p; mm++) {
              c_xaxpy(nmq - 1, -e->data[mm] / e->data[q + 1], work, q + 2, b_A,
                      (q + n * mm) + 2);
            }
          }
        }

        for (ns = q + 1; ns < p; ns++) {
          V->data[ns + V->size[0] * q] = e->data[ns];
        }
      }
    }

    qs = A->size[1];
    m = A->size[0] + 1;
    if (qs < m) {
      m = qs;
    }

    if (nct < A->size[1]) {
      b_s->data[nct] = b_A->data[nct + b_A->size[0] * nct];
    }

    if (A->size[0] < m) {
      b_s->data[m - 1] = 0.0;
    }

    if (nrt + 1 < m) {
      e->data[nrt] = b_A->data[nrt + b_A->size[0] * (m - 1)];
    }

    e->data[m - 1] = 0.0;
    if (nct + 1 <= A->size[0]) {
      for (mm = nct; mm < n; mm++) {
        for (ns = 1; ns <= n; ns++) {
          U->data[(ns + U->size[0] * mm) - 1] = 0.0;
        }

        U->data[mm + U->size[0] * mm] = 1.0;
      }
    }

    for (q = nct - 1; q + 1 > 0; q--) {
      nmq = n - q;
      qs = q + n * q;
      if (b_s->data[q] != 0.0) {
        for (mm = q + 1; mm < n; mm++) {
          ns = (q + n * mm) + 1;
          ztest0 = -(xdotc(nmq, U, qs + 1, U, ns) / U->data[qs]);
          xaxpy(nmq, ztest0, qs + 1, U, ns);
        }

        for (ns = q; ns < n; ns++) {
          U->data[ns + U->size[0] * q] = -U->data[ns + U->size[0] * q];
        }

        U->data[qs]++;
        for (ns = 1; ns <= q; ns++) {
          U->data[(ns + U->size[0] * q) - 1] = 0.0;
        }
      } else {
        for (ns = 1; ns <= n; ns++) {
          U->data[(ns + U->size[0] * q) - 1] = 0.0;
        }

        U->data[qs] = 1.0;
      }
    }

    for (q = A->size[1] - 1; q + 1 > 0; q--) {
      if ((q + 1 <= nrt) && (e->data[q] != 0.0)) {
        iter = (p - q) - 1;
        ns = (q + p * q) + 2;
        for (mm = q + 1; mm < p; mm++) {
          qs = (q + p * mm) + 2;
          ztest0 = -(xdotc(iter, V, ns, V, qs) / V->data[ns - 1]);
          xaxpy(iter, ztest0, ns, V, qs);
        }
      }

      for (ns = 1; ns <= p; ns++) {
        V->data[(ns + V->size[0] * q) - 1] = 0.0;
      }

      V->data[q + V->size[0] * q] = 1.0;
    }

    for (q = 0; q < m; q++) {
      if (b_s->data[q] != 0.0) {
        ztest = fabs(b_s->data[q]);
        ztest0 = b_s->data[q] / ztest;
        b_s->data[q] = ztest;
        if (q + 1 < m) {
          e->data[q] /= ztest0;
        }

        if (q + 1 <= n) {
          xscal(n, ztest0, U, 1 + n * q);
        }
      }

      if ((q + 1 < m) && (e->data[q] != 0.0)) {
        ztest = fabs(e->data[q]);
        ztest0 = ztest / e->data[q];
        e->data[q] = ztest;
        b_s->data[q + 1] *= ztest0;
        xscal(p, ztest0, V, 1 + p * (q + 1));
      }
    }

    mm = m;
    iter = 0;
    snorm = 0.0;
    for (ns = 0; ns < m; ns++) {
      snorm = fmax(snorm, fmax(fabs(b_s->data[ns]), fabs(e->data[ns])));
    }

    while ((m > 0) && (!(iter >= 75))) {
      q = m - 1;
      exitg1 = false;
      while (!(exitg1 || (q == 0))) {
        ztest0 = fabs(e->data[q - 1]);
        if ((ztest0 <= 2.2204460492503131E-16 * (fabs(b_s->data[q - 1]) + fabs
              (b_s->data[q]))) || (ztest0 <= 1.0020841800044864E-292) || ((iter >
              20) && (ztest0 <= 2.2204460492503131E-16 * snorm))) {
          e->data[q - 1] = 0.0;
          exitg1 = true;
        } else {
          q--;
        }
      }

      if (q == m - 1) {
        ns = 4;
      } else {
        qs = m;
        ns = m;
        exitg1 = false;
        while ((!exitg1) && (ns >= q)) {
          qs = ns;
          if (ns == q) {
            exitg1 = true;
          } else {
            ztest0 = 0.0;
            if (ns < m) {
              ztest0 = fabs(e->data[ns - 1]);
            }

            if (ns > q + 1) {
              ztest0 += fabs(e->data[ns - 2]);
            }

            ztest = fabs(b_s->data[ns - 1]);
            if ((ztest <= 2.2204460492503131E-16 * ztest0) || (ztest <=
                 1.0020841800044864E-292)) {
              b_s->data[ns - 1] = 0.0;
              exitg1 = true;
            } else {
              ns--;
            }
          }
        }

        if (qs == q) {
          ns = 3;
        } else if (qs == m) {
          ns = 1;
        } else {
          ns = 2;
          q = qs;
        }
      }

      switch (ns) {
       case 1:
        f = e->data[m - 2];
        e->data[m - 2] = 0.0;
        for (ns = m - 3; ns + 2 >= q + 1; ns--) {
          xrotg(&b_s->data[ns + 1], &f, &ztest0, &ztest);
          if (ns + 2 > q + 1) {
            f = -ztest * e->data[ns];
            e->data[ns] *= ztest0;
          }

          xrot(p, V, 1 + p * (ns + 1), 1 + p * (m - 1), ztest0, ztest);
        }
        break;

       case 2:
        f = e->data[q - 1];
        e->data[q - 1] = 0.0;
        for (ns = q; ns < m; ns++) {
          xrotg(&b_s->data[ns], &f, &ztest0, &ztest);
          f = -ztest * e->data[ns];
          e->data[ns] *= ztest0;
          xrot(n, U, 1 + n * ns, 1 + n * (q - 1), ztest0, ztest);
        }
        break;

       case 3:
        scale = fmax(fmax(fmax(fmax(fabs(b_s->data[m - 1]), fabs(b_s->data[m - 2])),
          fabs(e->data[m - 2])), fabs(b_s->data[q])), fabs(e->data[q]));
        f = b_s->data[m - 1] / scale;
        ztest0 = b_s->data[m - 2] / scale;
        ztest = e->data[m - 2] / scale;
        sqds = b_s->data[q] / scale;
        b = ((ztest0 + f) * (ztest0 - f) + ztest * ztest) / 2.0;
        ztest0 = f * ztest;
        ztest0 *= ztest0;
        if ((b != 0.0) || (ztest0 != 0.0)) {
          ztest = b * b + ztest0;
          b_sqrt(&ztest);
          if (b < 0.0) {
            ztest = -ztest;
          }

          ztest = ztest0 / (b + ztest);
        } else {
          ztest = 0.0;
        }

        f = (sqds + f) * (sqds - f) + ztest;
        b = sqds * (e->data[q] / scale);
        for (ns = q + 1; ns < m; ns++) {
          xrotg(&f, &b, &ztest0, &ztest);
          if (ns > q + 1) {
            e->data[ns - 2] = f;
          }

          f = ztest0 * b_s->data[ns - 1] + ztest * e->data[ns - 1];
          e->data[ns - 1] = ztest0 * e->data[ns - 1] - ztest * b_s->data[ns - 1];
          b = ztest * b_s->data[ns];
          b_s->data[ns] *= ztest0;
          xrot(p, V, 1 + p * (ns - 1), 1 + p * ns, ztest0, ztest);
          b_s->data[ns - 1] = f;
          xrotg(&b_s->data[ns - 1], &b, &ztest0, &ztest);
          f = ztest0 * e->data[ns - 1] + ztest * b_s->data[ns];
          b_s->data[ns] = -ztest * e->data[ns - 1] + ztest0 * b_s->data[ns];
          b = ztest * e->data[ns];
          e->data[ns] *= ztest0;
          if (ns < n) {
            xrot(n, U, 1 + n * (ns - 1), 1 + n * ns, ztest0, ztest);
          }
        }

        e->data[m - 2] = f;
        iter++;
        break;

       default:
        if (b_s->data[q] < 0.0) {
          b_s->data[q] = -b_s->data[q];
          xscal(p, -1.0, V, 1 + p * q);
        }

        qs = q + 1;
        while ((q + 1 < mm) && (b_s->data[q] < b_s->data[qs])) {
          ztest = b_s->data[q];
          b_s->data[q] = b_s->data[qs];
          b_s->data[qs] = ztest;
          if (q + 1 < p) {
            xswap(p, V, 1 + p * q, 1 + p * (q + 1));
          }

          if (q + 1 < n) {
            xswap(n, U, 1 + n * q, 1 + n * (q + 1));
          }

          q = qs;
          qs++;
        }

        iter = 0;
        m--;
        break;
      }
    }
  }

  emxFree_real_T(&work);
  emxFree_real_T(&e);
  emxFree_real_T(&b_A);
  m = s->size[0];
  s->size[0] = minnp;
  emxEnsureCapacity_real_T1(s, m);
  for (ns = 0; ns < minnp; ns++) {
    s->data[ns] = b_s->data[ns];
  }

  emxFree_real_T(&b_s);
}

/*
 * Arguments    : const emxArray_real_T *x
 *                emxArray_real_T *y
 * Return Type  : void
 */
static void vecnorm(const emxArray_real_T *x, emxArray_real_T *y)
{
  int i19;
  unsigned int szy[2];
  int outsize_idx_0;
  int k;
  emxArray_real_T *xv;
  double yv;
  double scale;
  double absxk;
  double t;
  for (i19 = 0; i19 < 2; i19++) {
    szy[i19] = (unsigned int)x->size[i19];
  }

  i19 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = (int)szy[1];
  emxEnsureCapacity_real_T(y, i19);
  outsize_idx_0 = (int)szy[1];
  for (i19 = 0; i19 < outsize_idx_0; i19++) {
    y->data[y->size[0] * i19] = 0.0;
  }

  k = 1;
  emxInit_real_T1(&xv, 1);
  while (k <= x->size[1]) {
    outsize_idx_0 = x->size[0];
    i19 = xv->size[0];
    xv->size[0] = outsize_idx_0;
    emxEnsureCapacity_real_T1(xv, i19);
    for (i19 = 0; i19 < outsize_idx_0; i19++) {
      xv->data[i19] = 0.0;
    }

    for (outsize_idx_0 = 0; outsize_idx_0 < x->size[0]; outsize_idx_0++) {
      xv->data[outsize_idx_0] = x->data[outsize_idx_0 + x->size[0] * (k - 1)];
    }

    if (xv->size[0] == 0) {
      yv = 0.0;
    } else {
      yv = 0.0;
      if (xv->size[0] == 1) {
        yv = fabs(xv->data[0]);
      } else {
        scale = 3.3121686421112381E-170;
        for (outsize_idx_0 = 1; outsize_idx_0 <= xv->size[0]; outsize_idx_0++) {
          absxk = fabs(xv->data[outsize_idx_0 - 1]);
          if (absxk > scale) {
            t = scale / absxk;
            yv = 1.0 + yv * t * t;
            scale = absxk;
          } else {
            t = absxk / scale;
            yv += t * t;
          }
        }

        yv = scale * sqrt(yv);
      }
    }

    y->data[y->size[0] * (k - 1)] = yv;
    k++;
  }

  emxFree_real_T(&xv);
}

/*
 * Arguments    : int n
 *                double b_a
 *                int ix0
 *                emxArray_real_T *y
 *                int iy0
 * Return Type  : void
 */
static void xaxpy(int n, double b_a, int ix0, emxArray_real_T *y, int iy0)
{
  int ix;
  int iy;
  int k;
  if ((n < 1) || (b_a == 0.0)) {
  } else {
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k < n; k++) {
      y->data[iy] += b_a * y->data[ix];
      ix++;
      iy++;
    }
  }
}

/*
 * Arguments    : int n
 *                const emxArray_real_T *x
 *                int ix0
 *                const emxArray_real_T *y
 *                int iy0
 * Return Type  : double
 */
static double xdotc(int n, const emxArray_real_T *x, int ix0, const
                    emxArray_real_T *y, int iy0)
{
  double b_d;
  int ix;
  int iy;
  int k;
  b_d = 0.0;
  if (!(n < 1)) {
    ix = ix0;
    iy = iy0;
    for (k = 1; k <= n; k++) {
      b_d += x->data[ix - 1] * y->data[iy - 1];
      ix++;
      iy++;
    }
  }

  return b_d;
}

/*
 * Arguments    : emxArray_real_T *A
 *                emxArray_real_T *tau
 *                emxArray_int32_T *jpvt
 * Return Type  : void
 */
static void xgeqp3(emxArray_real_T *A, emxArray_real_T *tau, emxArray_int32_T
                   *jpvt)
{
  int m;
  int n;
  int pvt;
  int mn;
  int i29;
  emxArray_real_T *work;
  emxArray_real_T *vn1;
  emxArray_real_T *vn2;
  int k;
  int nmi;
  int i;
  int i_i;
  int mmi;
  int ix;
  double smax;
  double atmp;
  double d2;
  double s;
  int i_ip1;
  int lastv;
  int lastc;
  bool exitg2;
  int ia;
  int exitg1;
  m = A->size[0];
  n = A->size[1];
  pvt = A->size[0];
  mn = A->size[1];
  if (pvt < mn) {
    mn = pvt;
  }

  i29 = tau->size[0];
  tau->size[0] = mn;
  emxEnsureCapacity_real_T1(tau, i29);
  eml_signed_integer_colon(A->size[1], jpvt);
  if (!((A->size[0] == 0) || (A->size[1] == 0))) {
    emxInit_real_T1(&work, 1);
    pvt = A->size[1];
    i29 = work->size[0];
    work->size[0] = pvt;
    emxEnsureCapacity_real_T1(work, i29);
    for (i29 = 0; i29 < pvt; i29++) {
      work->data[i29] = 0.0;
    }

    emxInit_real_T1(&vn1, 1);
    emxInit_real_T1(&vn2, 1);
    pvt = A->size[1];
    i29 = vn1->size[0];
    vn1->size[0] = pvt;
    emxEnsureCapacity_real_T1(vn1, i29);
    i29 = vn2->size[0];
    vn2->size[0] = vn1->size[0];
    emxEnsureCapacity_real_T1(vn2, i29);
    k = 1;
    for (nmi = 0; nmi < n; nmi++) {
      vn1->data[nmi] = xnrm2(m, A, k);
      vn2->data[nmi] = vn1->data[nmi];
      k += m;
    }

    for (i = 0; i < mn; i++) {
      i_i = i + i * m;
      nmi = n - i;
      mmi = (m - i) - 1;
      if (nmi < 1) {
        pvt = 0;
      } else {
        pvt = 1;
        if (nmi > 1) {
          ix = i;
          smax = fabs(vn1->data[i]);
          for (k = 2; k <= nmi; k++) {
            ix++;
            s = fabs(vn1->data[ix]);
            if (s > smax) {
              pvt = k;
              smax = s;
            }
          }
        }
      }

      pvt = (i + pvt) - 1;
      if (pvt + 1 != i + 1) {
        xswap(m, A, 1 + m * pvt, 1 + m * i);
        k = jpvt->data[pvt];
        jpvt->data[pvt] = jpvt->data[i];
        jpvt->data[i] = k;
        vn1->data[pvt] = vn1->data[i];
        vn2->data[pvt] = vn2->data[i];
      }

      if (i + 1 < m) {
        atmp = A->data[i_i];
        d2 = 0.0;
        if (!(1 + mmi <= 0)) {
          smax = xnrm2(mmi, A, i_i + 2);
          if (smax != 0.0) {
            s = rt_hypotd(A->data[i_i], smax);
            if (A->data[i_i] >= 0.0) {
              s = -s;
            }

            if (fabs(s) < 1.0020841800044864E-292) {
              pvt = 0;
              do {
                pvt++;
                xscal(mmi, 9.9792015476736E+291, A, i_i + 2);
                s *= 9.9792015476736E+291;
                atmp *= 9.9792015476736E+291;
              } while (!(fabs(s) >= 1.0020841800044864E-292));

              s = rt_hypotd(atmp, xnrm2(mmi, A, i_i + 2));
              if (atmp >= 0.0) {
                s = -s;
              }

              d2 = (s - atmp) / s;
              xscal(mmi, 1.0 / (atmp - s), A, i_i + 2);
              for (k = 1; k <= pvt; k++) {
                s *= 1.0020841800044864E-292;
              }

              atmp = s;
            } else {
              d2 = (s - A->data[i_i]) / s;
              smax = 1.0 / (A->data[i_i] - s);
              xscal(mmi, smax, A, i_i + 2);
              atmp = s;
            }
          }
        }

        tau->data[i] = d2;
        A->data[i_i] = atmp;
      } else {
        tau->data[i] = 0.0;
      }

      if (i + 1 < n) {
        atmp = A->data[i_i];
        A->data[i_i] = 1.0;
        i_ip1 = (i + (i + 1) * m) + 1;
        if (tau->data[i] != 0.0) {
          lastv = mmi;
          pvt = i_i + mmi;
          while ((lastv + 1 > 0) && (A->data[pvt] == 0.0)) {
            lastv--;
            pvt--;
          }

          lastc = nmi - 1;
          exitg2 = false;
          while ((!exitg2) && (lastc > 0)) {
            pvt = i_ip1 + (lastc - 1) * m;
            ia = pvt;
            do {
              exitg1 = 0;
              if (ia <= pvt + lastv) {
                if (A->data[ia - 1] != 0.0) {
                  exitg1 = 1;
                } else {
                  ia++;
                }
              } else {
                lastc--;
                exitg1 = 2;
              }
            } while (exitg1 == 0);

            if (exitg1 == 1) {
              exitg2 = true;
            }
          }
        } else {
          lastv = -1;
          lastc = 0;
        }

        if (lastv + 1 > 0) {
          if (lastc != 0) {
            for (pvt = 1; pvt <= lastc; pvt++) {
              work->data[pvt - 1] = 0.0;
            }

            pvt = 0;
            i29 = i_ip1 + m * (lastc - 1);
            k = i_ip1;
            while (((m > 0) && (k <= i29)) || ((m < 0) && (k >= i29))) {
              ix = i_i;
              smax = 0.0;
              nmi = k + lastv;
              for (ia = k; ia <= nmi; ia++) {
                smax += A->data[ia - 1] * A->data[ix];
                ix++;
              }

              work->data[pvt] += smax;
              pvt++;
              k += m;
            }
          }

          if (!(-tau->data[i] == 0.0)) {
            pvt = 0;
            for (nmi = 1; nmi <= lastc; nmi++) {
              if (work->data[pvt] != 0.0) {
                smax = work->data[pvt] * -tau->data[i];
                ix = i_i;
                i29 = lastv + i_ip1;
                for (k = i_ip1; k <= i29; k++) {
                  A->data[k - 1] += A->data[ix] * smax;
                  ix++;
                }
              }

              pvt++;
              i_ip1 += m;
            }
          }
        }

        A->data[i_i] = atmp;
      }

      for (nmi = i + 1; nmi < n; nmi++) {
        if (vn1->data[nmi] != 0.0) {
          smax = fabs(A->data[i + A->size[0] * nmi]) / vn1->data[nmi];
          smax = 1.0 - smax * smax;
          if (smax < 0.0) {
            smax = 0.0;
          }

          s = vn1->data[nmi] / vn2->data[nmi];
          s = smax * (s * s);
          if (s <= 1.4901161193847656E-8) {
            if (i + 1 < m) {
              vn1->data[nmi] = xnrm2(mmi, A, (i + m * nmi) + 2);
              vn2->data[nmi] = vn1->data[nmi];
            } else {
              vn1->data[nmi] = 0.0;
              vn2->data[nmi] = 0.0;
            }
          } else {
            vn1->data[nmi] *= sqrt(smax);
          }
        }
      }
    }

    emxFree_real_T(&vn2);
    emxFree_real_T(&vn1);
    emxFree_real_T(&work);
  }
}

/*
 * Arguments    : int m
 *                int n
 *                emxArray_real_T *A
 *                int lda
 *                emxArray_int32_T *ipiv
 *                int *info
 * Return Type  : void
 */
static void xgetrf(int m, int n, emxArray_real_T *A, int lda, emxArray_int32_T
                   *ipiv, int *info)
{
  int iy;
  int u0;
  int j;
  int mmj;
  int c;
  int ix;
  double smax;
  int jy;
  int i28;
  int ixinc;
  double s;
  int iyinc;
  int b_j;
  int ijA;
  if (m < n) {
    iy = m;
  } else {
    iy = n;
  }

  eml_signed_integer_colon(iy, ipiv);
  *info = 0;
  if ((m < 1) || (n < 1)) {
  } else {
    u0 = m - 1;
    if (!(u0 < n)) {
      u0 = n;
    }

    for (j = 0; j < u0; j++) {
      mmj = m - j;
      c = j * (lda + 1);
      if (mmj < 1) {
        iy = -1;
      } else {
        iy = 0;
        if (mmj > 1) {
          ix = c;
          smax = fabs(A->data[c]);
          for (jy = 1; jy < mmj; jy++) {
            ix++;
            s = fabs(A->data[ix]);
            if (s > smax) {
              iy = jy;
              smax = s;
            }
          }
        }
      }

      if (A->data[c + iy] != 0.0) {
        if (iy != 0) {
          ipiv->data[j] = (j + iy) + 1;
          ix = j;
          iy += j;
          if (lda < 0) {
            ixinc = -lda;
            iyinc = -lda;
          } else {
            ixinc = lda;
            iyinc = lda;
          }

          for (jy = 1; jy <= n; jy++) {
            smax = A->data[ix];
            A->data[ix] = A->data[iy];
            A->data[iy] = smax;
            if (lda < 0) {
              ix -= ixinc;
              iy -= iyinc;
            } else {
              ix += ixinc;
              iy += iyinc;
            }
          }
        }

        i28 = c + mmj;
        for (iy = c + 1; iy < i28; iy++) {
          A->data[iy] /= A->data[c];
        }
      } else {
        *info = j + 1;
      }

      iy = (n - j) - 1;
      if (lda < 0) {
        ixinc = -lda;
      } else {
        ixinc = lda;
      }

      iyinc = c + lda;
      jy = c + lda;
      for (b_j = 1; b_j <= iy; b_j++) {
        smax = A->data[jy];
        if (A->data[jy] != 0.0) {
          ix = c + 1;
          i28 = mmj + iyinc;
          for (ijA = 1 + iyinc; ijA < i28; ijA++) {
            A->data[ijA] += A->data[ix] * -smax;
            ix++;
          }
        }

        if (lda < 0) {
          jy -= ixinc;
        } else {
          jy += ixinc;
        }

        iyinc += lda;
      }
    }

    if ((*info == 0) && (m <= n) && (!(A->data[(m + A->size[0] * (m - 1)) - 1]
          != 0.0))) {
      *info = m;
    }
  }
}

/*
 * Arguments    : int n
 *                const emxArray_real_T *x
 *                int ix0
 * Return Type  : double
 */
static double xnrm2(int n, const emxArray_real_T *x, int ix0)
{
  double y;
  double scale;
  int kend;
  int k;
  double absxk;
  double t;
  y = 0.0;
  if (!(n < 1)) {
    if (n == 1) {
      y = fabs(x->data[ix0 - 1]);
    } else {
      scale = 3.3121686421112381E-170;
      kend = (ix0 + n) - 1;
      for (k = ix0; k <= kend; k++) {
        absxk = fabs(x->data[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          y = 1.0 + y * t * t;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * sqrt(y);
    }
  }

  return y;
}

/*
 * Arguments    : int n
 *                emxArray_real_T *x
 *                int ix0
 *                int iy0
 *                double c
 *                double s
 * Return Type  : void
 */
static void xrot(int n, emxArray_real_T *x, int ix0, int iy0, double c, double s)
{
  int ix;
  int iy;
  int k;
  double temp;
  if (!(n < 1)) {
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 1; k <= n; k++) {
      temp = c * x->data[ix] + s * x->data[iy];
      x->data[iy] = c * x->data[iy] - s * x->data[ix];
      x->data[ix] = temp;
      iy++;
      ix++;
    }
  }
}

/*
 * Arguments    : double *b_a
 *                double *b
 *                double *c
 *                double *s
 * Return Type  : void
 */
static void xrotg(double *b_a, double *b, double *c, double *s)
{
  double roe;
  double absa;
  double absb;
  double scale;
  double ads;
  double bds;
  roe = *b;
  absa = fabs(*b_a);
  absb = fabs(*b);
  if (absa > absb) {
    roe = *b_a;
  }

  scale = absa + absb;
  if (scale == 0.0) {
    *s = 0.0;
    *c = 1.0;
    scale = 0.0;
    *b = 0.0;
  } else {
    ads = absa / scale;
    bds = absb / scale;
    scale *= sqrt(ads * ads + bds * bds);
    if (roe < 0.0) {
      scale = -scale;
    }

    *c = *b_a / scale;
    *s = *b / scale;
    if (absa > absb) {
      *b = *s;
    } else if (*c != 0.0) {
      *b = 1.0 / *c;
    } else {
      *b = 1.0;
    }
  }

  *b_a = scale;
}

/*
 * Arguments    : int n
 *                double b_a
 *                emxArray_real_T *x
 *                int ix0
 * Return Type  : void
 */
static void xscal(int n, double b_a, emxArray_real_T *x, int ix0)
{
  int i30;
  int k;
  i30 = (ix0 + n) - 1;
  for (k = ix0; k <= i30; k++) {
    x->data[k - 1] *= b_a;
  }
}

/*
 * Arguments    : int n
 *                emxArray_real_T *x
 *                int ix0
 *                int iy0
 * Return Type  : void
 */
static void xswap(int n, emxArray_real_T *x, int ix0, int iy0)
{
  int ix;
  int iy;
  int k;
  double temp;
  ix = ix0 - 1;
  iy = iy0 - 1;
  for (k = 1; k <= n; k++) {
    temp = x->data[ix];
    x->data[ix] = x->data[iy];
    x->data[iy] = temp;
    ix++;
    iy++;
  }
}

/*
 * function [phi, tau_pos, tau_pre, tau_filt, errs] = dynIdenf(theta, theta_dot, theta_ddot, tau, pars)
 * the function identify the dynamic parameters
 * Arguments    : const emxArray_real_T *theta
 *                const emxArray_real_T *theta_dot
 *                const emxArray_real_T *theta_ddot
 *                const emxArray_real_T *tau
 *                const cell_0 *pars
 *                emxArray_real_T *phi
 *                emxArray_real_T *tau_pos
 *                emxArray_real_T *tau_pre
 *                emxArray_real_T *tau_filt
 *                cell_1 *errs
 * Return Type  : void
 */
void dynIdenf(const emxArray_real_T *theta, const emxArray_real_T *theta_dot,
              const emxArray_real_T *theta_ddot, const emxArray_real_T *tau,
              const cell_0 *pars, emxArray_real_T *phi, emxArray_real_T *tau_pos,
              emxArray_real_T *tau_pre, emxArray_real_T *tau_filt, cell_1 *errs)
{
  emxArray_real_T *filtErr;
  emxArray_real_T *x;
  emxArray_real_T *r0;
  int i0;
  int loop_ub;
  emxArray_real_T *b_theta_dot;
  int i;
  int k;
  emxArray_real_T *b_theta_ddot;
  emxArray_real_T *b_tau;
  double ndbl;
  bool isLevel;
  double b_x;
  double b_d;
  emxArray_real_T *setFilt;
  double apnd;
  double cdiff;
  double absa;
  double absb;
  emxArray_real_T *tauf;
  emxArray_real_T *r1;
  emxArray_real_T *c_x;
  int m;
  emxArray_real_T *d_x;
  emxArray_real_T *c_theta_dot;
  emxArray_real_T *r2;
  signed char iv0[2];
  emxArray_real_T *c_theta_ddot;
  emxArray_real_T *K;

  /*  inputs: */
  /*  theta - the raw angular position of joints */
  /*  theta_dot - the raw angular velocity of joints */
  /*  theta_ddot - the raw angular acceleration of joints */
  /*  tau - the raw torque of joints */
  /*  outputs: */
  /*  phi - the identified dynamic parameters */
  /*  tau_pos - the torque calculated by currently solved dynamic parameters */
  /*  tau_pre - the torque calculated by previously solved dynamic parameters */
  /*  errs - the errs of torque, errs = {distErr, filtErr, convErr} */
  /*  tau_filt - the measured, filtered and resampled torque of joints */
  /*  distErr - the percentage distribution of error in segErr, between tau_pos */
  /*  and tau_filt */
  /*  filtErr - the error of torque induced by the filtering, between tau_filt and tau */
  /*  convErr - convergence error of the iterations, between the tau_pre and */
  /*  tau_pos */
  /*  rk - the rank of the matrix K, representing the number of dynamic */
  /*  parameters that are updated in current iteration */
  /*  pars - the parameters */
  /*  pars = {1:a, 2:alpha, 3:d, 4:g, 5:phi_re, 6:pidenf, 7:peval, 8:noise_err, 9:cond_max, 10:lambda, 11:fpass} */
  /*  default values */
  /*  a = [0, 0, 0.416, 0.4208, 0, 0, 0]; */
  /*  alpha = [0, pi/2, 0, 0, -pi/2, pi/2, 0]; */
  /*  d = [0, 0.1181, 0, 0, 0.1301, 0.1021, 0.0568]; */
  /*  g = [0; 0; -9.80665]; */
  /*  pfilt = 5001; */
  /*  pidenf = [53, 100]; */
  /*  peval = [53, 100]; */
  /*  noise_err = 1e-6; */
  /*  cond_max = 100; */
  /*  lambda = 0; */
  /*  fpass = 2.0; */
  /*  segErr = [0.1 0.2 1.0]; */
  /*  phi_r0 = [2,-0.0165019647020000,-0.0256289655860000,-0.0456037864700000,0.00388830591570456,7.63638706167896e-05,0.000340020397992622,0.00333118131630822,0.000215211872926742,0.00262696792835856,0,0; */
  /*      3.42000000000000,0.426439873502640,-0.0124271511264000,0.386835887459880,0.0475012505609916,3.00901079796083e-06,-0.0528168520144169,0.131182004246096,0.00172637529085978,0.0864328868443469,0,0; */
  /*      1.26000000000000,0.139078355556420,-1.07045177400000e-05,0.0411503703017400,0.00223362215511204,5.03356734903604e-07,-0.00375135880026783,0.0264742256951927,3.69281505903818e-07,0.0248796360944246,0,0; */
  /*      0.800000000000000,1.39460800000000e-07,0.000225010768000000,-0.00478865628240000,0.000772782642110361,5.04623525623476e-09,8.09638107508407e-09,0.000537499391466433,-6.01459794336313e-06,0.000649992384815521,0,0; */
  /*      0.800000000000000,-2.63220000000000e-06,-0.000320608410400000,-0.00460317585600000,0.000770846701360872,-6.16295700192826e-08,-4.41585931650860e-08,0.000531210606854337,6.03332741029044e-06,0.000653231083089124,0,0; */
  /*      0.350000000000000,8.29500000000000e-11,-3.71160097000000e-05,-0.00683040118320000,0.000261196607049744,-4.06640243700716e-12,1.37724047043447e-11,0.000261849371147053,-4.02812899656085e-07,0.000175309188669306,0,0]'; */
  /*  the counter to count the iterations, and store the dynamic parameters in */
  /*  initialize the parameters */
  /* 'dynIdenf:50' if isempty(count) */
  emxInit_real_T(&filtErr, 2);
  emxInit_real_T(&x, 2);
  emxInit_real_T(&r0, 2);
  if (!count_not_empty) {
    /*  set the numbers and dimensions */
    /*  n - the number of joints */
    /*  pfilt - the number of remained points after filtering and truncation */
    /*  pidenf - the number of resampled points for parameter identification */
    /*  per cycle, [points, interval] */
    /*  peval - the number of resampled points for torque evaluation per */
    /*  cycle, [points, interval] */
    /*  setResample - the index of sampled points for parameter */
    /*  identification */
    /*  setEvalTorq - the index of sampled points for torque evaluation */
    /* 'dynIdenf:61' n = size(theta,2); */
    /* 'dynIdenf:62' pfilt = pars{6}; */
    pfilt = pars->f6;

    /* 'dynIdenf:63' pidenf = pars{7}; */
    /* 'dynIdenf:64' peval = pars{8}; */
    /* 'dynIdenf:65' setResample = ((0:pidenf(1)-1)*pidenf(2))+1; */
    if (pars->f7[0] - 1.0 < 0.0) {
      i0 = setResample->size[0] * setResample->size[1];
      setResample->size[0] = 1;
      setResample->size[1] = 0;
      emxEnsureCapacity_real_T(setResample, i0);
    } else {
      i0 = setResample->size[0] * setResample->size[1];
      setResample->size[0] = 1;
      setResample->size[1] = (int)floor(pars->f7[0] - 1.0) + 1;
      emxEnsureCapacity_real_T(setResample, i0);
      loop_ub = (int)floor(pars->f7[0] - 1.0);
      for (i0 = 0; i0 <= loop_ub; i0++) {
        setResample->data[setResample->size[0] * i0] = i0;
      }
    }

    i0 = setResample->size[0] * setResample->size[1];
    setResample->size[0] = 1;
    emxEnsureCapacity_real_T(setResample, i0);
    loop_ub = setResample->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      setResample->data[setResample->size[0] * i0] = setResample->
        data[setResample->size[0] * i0] * pars->f7[1] + 1.0;
    }

    /* 'dynIdenf:66' setEvalTorq = ((0:peval(1)-1)*peval(2))+1; */
    if (pars->f8[0] - 1.0 < 0.0) {
      i0 = setEvalTorq->size[0] * setEvalTorq->size[1];
      setEvalTorq->size[0] = 1;
      setEvalTorq->size[1] = 0;
      emxEnsureCapacity_real_T(setEvalTorq, i0);
    } else {
      i0 = setEvalTorq->size[0] * setEvalTorq->size[1];
      setEvalTorq->size[0] = 1;
      setEvalTorq->size[1] = (int)floor(pars->f8[0] - 1.0) + 1;
      emxEnsureCapacity_real_T(setEvalTorq, i0);
      loop_ub = (int)floor(pars->f8[0] - 1.0);
      for (i0 = 0; i0 <= loop_ub; i0++) {
        setEvalTorq->data[setEvalTorq->size[0] * i0] = i0;
      }
    }

    i0 = setEvalTorq->size[0] * setEvalTorq->size[1];
    setEvalTorq->size[0] = 1;
    emxEnsureCapacity_real_T(setEvalTorq, i0);
    loop_ub = setEvalTorq->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      setEvalTorq->data[setEvalTorq->size[0] * i0] = setEvalTorq->
        data[setEvalTorq->size[0] * i0] * pars->f8[1] + 1.0;
    }

    /*  calculate the number of dynamic parameters, considering the gravity */
    /*  g - the gravity */
    /*  m - the spatial dimensions */
    /*  isLevel - whether the base is horizontal or declined */
    /*  nparJoint - the number of parameters per joint; */
    /*  nparMinSet - the number of minimum set of parameters, 48(horizontal base) */
    /*  or 50 (declined base) */
    /* 'dynIdenf:75' g = pars{4}; */
    for (i = 0; i < 3; i++) {
      g[i] = pars->f4[i];
    }

    /* 'dynIdenf:76' m = length(g); */
    /* 'dynIdenf:77' isLevel = g(1)==0 && g(2)==0; */
    if ((g[0] == 0.0) && (g[1] == 0.0)) {
      isLevel = true;
    } else {
      isLevel = false;
    }

    /* 'dynIdenf:78' nparJoint = (m^2+3*m+6)/2; */
    /* 'dynIdenf:79' nparMinSet = [nparJoint*n - 3*n - 4, nparJoint*n - 3*n - 6]; */
    /*  the D-H parameters and Regroup matrix */
    /* 'dynIdenf:82' a = pars{1}; */
    i0 = a->size[0] * a->size[1];
    a->size[0] = 1;
    a->size[1] = pars->f1->size[1];
    emxEnsureCapacity_real_T(a, i0);
    loop_ub = pars->f1->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      a->data[a->size[0] * i0] = pars->f1->data[pars->f1->size[0] * i0];
    }

    /* 'dynIdenf:83' alpha = pars{2}; */
    i0 = alpha->size[0] * alpha->size[1];
    alpha->size[0] = 1;
    alpha->size[1] = pars->f2->size[1];
    emxEnsureCapacity_real_T(alpha, i0);
    loop_ub = pars->f2->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      alpha->data[alpha->size[0] * i0] = pars->f2->data[pars->f2->size[0] * i0];
    }

    /* 'dynIdenf:84' d = pars{3}; */
    i0 = d->size[0] * d->size[1];
    d->size[0] = 1;
    d->size[1] = pars->f3->size[1];
    emxEnsureCapacity_real_T(d, i0);
    loop_ub = pars->f3->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      d->data[d->size[0] * i0] = pars->f3->data[pars->f3->size[0] * i0];
    }

    /* 'dynIdenf:85' [Rpk, Rphi] = regroup(isLevel, nparJoint, nparMinSet(isLevel+1), a(1:n+1), alpha(1:n+1), d(1:n+1)); */
    loop_ub = theta->size[1] + 1;
    i = theta->size[1] + 1;
    m = theta->size[1] + 1;
    i0 = r0->size[0] * r0->size[1];
    r0->size[0] = 1;
    r0->size[1] = loop_ub;
    emxEnsureCapacity_real_T(r0, i0);
    for (i0 = 0; i0 < loop_ub; i0++) {
      r0->data[r0->size[0] * i0] = a->data[i0];
    }

    i0 = filtErr->size[0] * filtErr->size[1];
    filtErr->size[0] = 1;
    filtErr->size[1] = i;
    emxEnsureCapacity_real_T(filtErr, i0);
    for (i0 = 0; i0 < i; i0++) {
      filtErr->data[filtErr->size[0] * i0] = alpha->data[i0];
    }

    emxInit_real_T(&r2, 2);
    i0 = r2->size[0] * r2->size[1];
    r2->size[0] = 1;
    r2->size[1] = m;
    emxEnsureCapacity_real_T(r2, i0);
    for (i0 = 0; i0 < m; i0++) {
      r2->data[r2->size[0] * i0] = d->data[i0];
    }

    iv0[0] = (signed char)((12 * theta->size[1] - 3 * theta->size[1]) - 4);
    iv0[1] = (signed char)((12 * theta->size[1] - 3 * theta->size[1]) - 6);
    regroup(isLevel, 12.0, iv0[(int)isLevel], r0, filtErr, r2, Rpk, Rphi);

    /*  about the least-square resolution */
    /*  noise_err - the threshold for elements in K not zero */
    /*  cond_max - the maximum condition number of K */
    /*  lambda - the factor for considering the priori */
    /* 'dynIdenf:91' noise_err = pars{9}; */
    noise_err = pars->f9;

    /* 'dynIdenf:92' cond_max = pars{10}; */
    cond_max = pars->f10;

    /* 'dynIdenf:93' lambda = pars{11}; */
    lambda = pars->f11;

    /*  set parameters of the filter */
    /*  fpass - pass frequency of the filter */
    /*  num - the numerator of the transfer function of the filter */
    /*  den - the denominator of the transfer function of the filter */
    /* 'dynIdenf:99' fpass = pars{12}; */
    /* 'dynIdenf:100' [num, den] = getCoeffs(fpass); */
    getCoeffs(pars->f12, num, &ndbl);
    den = ndbl;

    /*  set the priori values of the parameters */
    /*  phi_r0 the full form of the priori dynamic parameters */
    /*  phi_r - the regrouped form of the priori dynamic parameters */
    /* 'dynIdenf:105' phi_r0 = pars{5}; */
    i0 = phi_r0->size[0] * phi_r0->size[1];
    phi_r0->size[0] = 12;
    phi_r0->size[1] = pars->f5->size[1];
    emxEnsureCapacity_real_T(phi_r0, i0);
    loop_ub = pars->f5->size[1];
    emxFree_real_T(&r2);
    for (i0 = 0; i0 < loop_ub; i0++) {
      for (k = 0; k < 12; k++) {
        phi_r0->data[k + phi_r0->size[0] * i0] = pars->f5->data[k + pars->
          f5->size[0] * i0];
      }
    }

    /* 'dynIdenf:106' phi_r = Rphi * reshape(phi_r0(1:nparJoint,1:n),n*nparJoint,1); */
    if (1 > theta->size[1]) {
      loop_ub = 0;
    } else {
      loop_ub = theta->size[1];
    }

    i0 = x->size[0] * x->size[1];
    x->size[0] = 12;
    x->size[1] = loop_ub;
    emxEnsureCapacity_real_T(x, i0);
    for (i0 = 0; i0 < loop_ub; i0++) {
      for (k = 0; k < 12; k++) {
        x->data[k + x->size[0] * i0] = phi_r0->data[k + phi_r0->size[0] * i0];
      }
    }

    m = theta->size[1] * 12;
    if ((Rphi->size[1] == 1) || (m == 1)) {
      i0 = phi_r->size[0];
      phi_r->size[0] = Rphi->size[0];
      emxEnsureCapacity_real_T1(phi_r, i0);
      loop_ub = Rphi->size[0];
      for (i0 = 0; i0 < loop_ub; i0++) {
        phi_r->data[i0] = 0.0;
        i = Rphi->size[1];
        for (k = 0; k < i; k++) {
          phi_r->data[i0] += Rphi->data[i0 + Rphi->size[0] * k] * x->data[k];
        }
      }
    } else {
      m = Rphi->size[0];
      i0 = phi_r->size[0];
      phi_r->size[0] = Rphi->size[0];
      emxEnsureCapacity_real_T1(phi_r, i0);
      for (i = 1; i <= m; i++) {
        phi_r->data[i - 1] = 0.0;
      }

      for (k = 0; k < Rphi->size[1]; k++) {
        if (x->data[k] != 0.0) {
          for (i = 0; i < m; i++) {
            phi_r->data[i] += x->data[k] * Rphi->data[i + Rphi->size[0] * k];
          }
        }
      }
    }

    /*  set the error segments */
    /*  segErr - segments of the error */
    /* 'dynIdenf:110' segErr = pars{13}; */
    i0 = segErr->size[0] * segErr->size[1];
    segErr->size[0] = 1;
    segErr->size[1] = pars->f13->size[1];
    emxEnsureCapacity_real_T(segErr, i0);
    loop_ub = pars->f13->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      segErr->data[segErr->size[0] * i0] = pars->f13->data[pars->f13->size[0] *
        i0];
    }

    /*  the counter */
    /* 'dynIdenf:113' count = 0; */
    count = 0.0;
    count_not_empty = true;

    /*  the previous parameters */
    /* 'dynIdenf:115' phi_pre = phi_r; */
    i0 = phi_pre->size[0];
    phi_pre->size[0] = phi_r->size[0];
    emxEnsureCapacity_real_T1(phi_pre, i0);
    loop_ub = phi_r->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      phi_pre->data[i0] = phi_r->data[i0];
    }
  }

  /*  filter the data and evaluate the filtering error of the torque */
  /* 'dynIdenf:119' [theta, theta_dot, theta_ddot, tau, filtErr] = filtTrimData(theta, theta_dot, theta_ddot, tau, pfilt, num, den); */
  i0 = x->size[0] * x->size[1];
  x->size[0] = theta->size[0];
  x->size[1] = theta->size[1];
  emxEnsureCapacity_real_T(x, i0);
  loop_ub = theta->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    i = theta->size[0];
    for (k = 0; k < i; k++) {
      x->data[k + x->size[0] * i0] = theta->data[k + theta->size[0] * i0];
    }
  }

  emxInit_real_T(&b_theta_dot, 2);
  i0 = b_theta_dot->size[0] * b_theta_dot->size[1];
  b_theta_dot->size[0] = theta_dot->size[0];
  b_theta_dot->size[1] = theta_dot->size[1];
  emxEnsureCapacity_real_T(b_theta_dot, i0);
  loop_ub = theta_dot->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    i = theta_dot->size[0];
    for (k = 0; k < i; k++) {
      b_theta_dot->data[k + b_theta_dot->size[0] * i0] = theta_dot->data[k +
        theta_dot->size[0] * i0];
    }
  }

  emxInit_real_T(&b_theta_ddot, 2);
  i0 = b_theta_ddot->size[0] * b_theta_ddot->size[1];
  b_theta_ddot->size[0] = theta_ddot->size[0];
  b_theta_ddot->size[1] = theta_ddot->size[1];
  emxEnsureCapacity_real_T(b_theta_ddot, i0);
  loop_ub = theta_ddot->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    i = theta_ddot->size[0];
    for (k = 0; k < i; k++) {
      b_theta_ddot->data[k + b_theta_ddot->size[0] * i0] = theta_ddot->data[k +
        theta_ddot->size[0] * i0];
    }
  }

  emxInit_real_T(&b_tau, 2);
  i0 = b_tau->size[0] * b_tau->size[1];
  b_tau->size[0] = tau->size[0];
  b_tau->size[1] = tau->size[1];
  emxEnsureCapacity_real_T(b_tau, i0);
  loop_ub = tau->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    i = tau->size[0];
    for (k = 0; k < i; k++) {
      b_tau->data[k + b_tau->size[0] * i0] = tau->data[k + tau->size[0] * i0];
    }
  }

  /*  this function filter the data */
  /*  inputs */
  /*  theta - the raw angular position of joints */
  /*  theta_dot - the raw angular velocity of joints */
  /*  theta_ddot - the raw angular acceleration of joints */
  /*  tau - the raw torque of joints */
  /*  parameters */
  /*  pfilt - the number of remained points after filtering and truncation */
  /*  num - numerator of the transfer function of the filter */
  /*  den - denominator of the transfer function of the filter */
  /*  outputs */
  /*  theta - the filtered, resampled angular position of joints */
  /*  theta_dot - the filtered, resampled angular velocity of joints */
  /*  theta_ddot - the filtered, resampled angular acceleration of joints */
  /*  tau- the filtered, resampled torque of joints */
  /*  filtErr - the error of torque induced by error */
  /*  obtain the dimensions of data */
  /* 'dynIdenf:168' [p, n] = size(theta); */
  /*  obtain the index set of the remained points after filtering */
  /* 'dynIdenf:170' setFilt = floor((p-pfilt)/2)+1:floor((p-pfilt)/2)+pfilt; */
  ndbl = ((double)theta->size[0] - pfilt) / 2.0;
  b_x = floor(ndbl);
  b_d = floor(((double)theta->size[0] - pfilt) / 2.0) + pfilt;
  emxInit_real_T(&setFilt, 2);
  if (b_d < b_x + 1.0) {
    i0 = setFilt->size[0] * setFilt->size[1];
    setFilt->size[0] = 1;
    setFilt->size[1] = 0;
    emxEnsureCapacity_real_T(setFilt, i0);
  } else if (floor(ndbl) + 1.0 == b_x + 1.0) {
    i0 = setFilt->size[0] * setFilt->size[1];
    setFilt->size[0] = 1;
    setFilt->size[1] = (int)floor(b_d - (b_x + 1.0)) + 1;
    emxEnsureCapacity_real_T(setFilt, i0);
    loop_ub = (int)floor(b_d - (b_x + 1.0));
    for (i0 = 0; i0 <= loop_ub; i0++) {
      setFilt->data[setFilt->size[0] * i0] = (b_x + 1.0) + (double)i0;
    }
  } else {
    ndbl = floor((b_d - (b_x + 1.0)) + 0.5);
    apnd = (b_x + 1.0) + ndbl;
    cdiff = apnd - b_d;
    absa = fabs(b_x + 1.0);
    absb = fabs(b_d);
    if (fabs(cdiff) < 4.4408920985006262E-16 * fmax(absa, absb)) {
      ndbl++;
      apnd = b_d;
    } else if (cdiff > 0.0) {
      apnd = (b_x + 1.0) + (ndbl - 1.0);
    } else {
      ndbl++;
    }

    if (ndbl >= 0.0) {
      i = (int)ndbl;
    } else {
      i = 0;
    }

    i0 = setFilt->size[0] * setFilt->size[1];
    setFilt->size[0] = 1;
    setFilt->size[1] = i;
    emxEnsureCapacity_real_T(setFilt, i0);
    if (i > 0) {
      setFilt->data[0] = b_x + 1.0;
      if (i > 1) {
        setFilt->data[i - 1] = apnd;
        m = (i - 1) / 2;
        for (k = 1; k < m; k++) {
          setFilt->data[k] = (b_x + 1.0) + (double)k;
          setFilt->data[(i - k) - 1] = apnd - (double)k;
        }

        if (m << 1 == i - 1) {
          setFilt->data[m] = ((b_x + 1.0) + apnd) / 2.0;
        } else {
          setFilt->data[m] = (b_x + 1.0) + (double)m;
          setFilt->data[m + 1] = apnd - (double)m;
        }
      }
    }
  }

  /*  initialize the filtering error */
  /* 'dynIdenf:172' filtErr = zeros(1,n); */
  /*  filter the data */
  /* 'dynIdenf:174' for i = 1:n */
  i0 = filtErr->size[0] * filtErr->size[1];
  filtErr->size[0] = 1;
  filtErr->size[1] = theta->size[1];
  emxEnsureCapacity_real_T(filtErr, i0);
  i = 0;
  emxInit_real_T1(&tauf, 1);
  emxInit_real_T1(&r1, 1);
  emxInit_real_T1(&c_x, 1);
  while (i <= theta->size[1] - 1) {
    /* 'dynIdenf:175' theta(:,i) = myfiltfilt(num,den,theta(:,i)); */
    i0 = r0->size[0] * r0->size[1];
    r0->size[0] = 1;
    r0->size[1] = num->size[1];
    emxEnsureCapacity_real_T(r0, i0);
    loop_ub = num->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      r0->data[r0->size[0] * i0] = num->data[num->size[0] * i0];
    }

    loop_ub = x->size[0];
    i0 = c_x->size[0];
    c_x->size[0] = loop_ub;
    emxEnsureCapacity_real_T1(c_x, i0);
    for (i0 = 0; i0 < loop_ub; i0++) {
      c_x->data[i0] = x->data[i0 + x->size[0] * i];
    }

    myfiltfilt(r0, den, c_x, tauf);
    loop_ub = tauf->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      x->data[i0 + x->size[0] * i] = tauf->data[i0];
    }

    /* 'dynIdenf:176' theta_dot(:,i) = myfiltfilt(num,den,theta_dot(:,i)); */
    i0 = r0->size[0] * r0->size[1];
    r0->size[0] = 1;
    r0->size[1] = num->size[1];
    emxEnsureCapacity_real_T(r0, i0);
    loop_ub = num->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      r0->data[r0->size[0] * i0] = num->data[num->size[0] * i0];
    }

    loop_ub = b_theta_dot->size[0];
    i0 = c_x->size[0];
    c_x->size[0] = loop_ub;
    emxEnsureCapacity_real_T1(c_x, i0);
    for (i0 = 0; i0 < loop_ub; i0++) {
      c_x->data[i0] = b_theta_dot->data[i0 + b_theta_dot->size[0] * i];
    }

    myfiltfilt(r0, den, c_x, tauf);
    loop_ub = tauf->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_theta_dot->data[i0 + b_theta_dot->size[0] * i] = tauf->data[i0];
    }

    /* 'dynIdenf:177' theta_ddot(:,i) = myfiltfilt(num,den,theta_ddot(:,i)); */
    i0 = r0->size[0] * r0->size[1];
    r0->size[0] = 1;
    r0->size[1] = num->size[1];
    emxEnsureCapacity_real_T(r0, i0);
    loop_ub = num->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      r0->data[r0->size[0] * i0] = num->data[num->size[0] * i0];
    }

    loop_ub = b_theta_ddot->size[0];
    i0 = c_x->size[0];
    c_x->size[0] = loop_ub;
    emxEnsureCapacity_real_T1(c_x, i0);
    for (i0 = 0; i0 < loop_ub; i0++) {
      c_x->data[i0] = b_theta_ddot->data[i0 + b_theta_ddot->size[0] * i];
    }

    myfiltfilt(r0, den, c_x, tauf);
    loop_ub = tauf->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_theta_ddot->data[i0 + b_theta_ddot->size[0] * i] = tauf->data[i0];
    }

    /* 'dynIdenf:178' tauf = myfiltfilt(num,den,tau(:,i)); */
    i0 = r0->size[0] * r0->size[1];
    r0->size[0] = 1;
    r0->size[1] = num->size[1];
    emxEnsureCapacity_real_T(r0, i0);
    loop_ub = num->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      r0->data[r0->size[0] * i0] = num->data[num->size[0] * i0];
    }

    loop_ub = b_tau->size[0];
    i0 = c_x->size[0];
    c_x->size[0] = loop_ub;
    emxEnsureCapacity_real_T1(c_x, i0);
    for (i0 = 0; i0 < loop_ub; i0++) {
      c_x->data[i0] = b_tau->data[i0 + b_tau->size[0] * i];
    }

    myfiltfilt(r0, den, c_x, tauf);

    /*  calculate the filtering error of the torque */
    /* 'dynIdenf:180' filtErr(i) = rms(tauf(setFilt)-tau(setFilt,i))/rms(tau(setFilt,i)); */
    i0 = r1->size[0];
    r1->size[0] = setFilt->size[1];
    emxEnsureCapacity_real_T1(r1, i0);
    loop_ub = setFilt->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      r1->data[i0] = b_tau->data[((int)setFilt->data[setFilt->size[0] * i0] +
        b_tau->size[0] * i) - 1];
    }

    i0 = c_x->size[0];
    c_x->size[0] = setFilt->size[1];
    emxEnsureCapacity_real_T1(c_x, i0);
    loop_ub = setFilt->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      c_x->data[i0] = b_tau->data[((int)setFilt->data[setFilt->size[0] * i0] +
        b_tau->size[0] * i) - 1];
    }

    ndbl = rms(c_x);
    i0 = c_x->size[0];
    c_x->size[0] = setFilt->size[1];
    emxEnsureCapacity_real_T1(c_x, i0);
    loop_ub = setFilt->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      c_x->data[i0] = tauf->data[(int)setFilt->data[setFilt->size[0] * i0] - 1]
        - r1->data[i0];
    }

    filtErr->data[i] = rms(c_x) / ndbl;

    /* 'dynIdenf:181' tau(:,i) = tauf; */
    loop_ub = tauf->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_tau->data[i0 + b_tau->size[0] * i] = tauf->data[i0];
    }

    i++;
  }

  emxFree_real_T(&c_x);
  emxFree_real_T(&r0);
  emxFree_real_T(&r1);
  emxFree_real_T(&tauf);
  emxInit_real_T(&d_x, 2);

  /*  trancate the data */
  /* 'dynIdenf:184' theta = theta(setFilt,:); */
  m = x->size[1];
  i0 = d_x->size[0] * d_x->size[1];
  d_x->size[0] = setFilt->size[1];
  d_x->size[1] = m;
  emxEnsureCapacity_real_T(d_x, i0);
  for (i0 = 0; i0 < m; i0++) {
    loop_ub = setFilt->size[1];
    for (k = 0; k < loop_ub; k++) {
      d_x->data[k + d_x->size[0] * i0] = x->data[((int)setFilt->data
        [setFilt->size[0] * k] + x->size[0] * i0) - 1];
    }
  }

  i0 = x->size[0] * x->size[1];
  x->size[0] = d_x->size[0];
  x->size[1] = d_x->size[1];
  emxEnsureCapacity_real_T(x, i0);
  loop_ub = d_x->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    i = d_x->size[0];
    for (k = 0; k < i; k++) {
      x->data[k + x->size[0] * i0] = d_x->data[k + d_x->size[0] * i0];
    }
  }

  emxInit_real_T(&c_theta_dot, 2);

  /* 'dynIdenf:185' theta_dot = theta_dot(setFilt,:); */
  m = b_theta_dot->size[1];
  i0 = c_theta_dot->size[0] * c_theta_dot->size[1];
  c_theta_dot->size[0] = setFilt->size[1];
  c_theta_dot->size[1] = m;
  emxEnsureCapacity_real_T(c_theta_dot, i0);
  for (i0 = 0; i0 < m; i0++) {
    loop_ub = setFilt->size[1];
    for (k = 0; k < loop_ub; k++) {
      c_theta_dot->data[k + c_theta_dot->size[0] * i0] = b_theta_dot->data[((int)
        setFilt->data[setFilt->size[0] * k] + b_theta_dot->size[0] * i0) - 1];
    }
  }

  i0 = b_theta_dot->size[0] * b_theta_dot->size[1];
  b_theta_dot->size[0] = c_theta_dot->size[0];
  b_theta_dot->size[1] = c_theta_dot->size[1];
  emxEnsureCapacity_real_T(b_theta_dot, i0);
  loop_ub = c_theta_dot->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    i = c_theta_dot->size[0];
    for (k = 0; k < i; k++) {
      b_theta_dot->data[k + b_theta_dot->size[0] * i0] = c_theta_dot->data[k +
        c_theta_dot->size[0] * i0];
    }
  }

  emxInit_real_T(&c_theta_ddot, 2);

  /* 'dynIdenf:186' theta_ddot = theta_ddot(setFilt,:); */
  m = b_theta_ddot->size[1];
  i0 = c_theta_ddot->size[0] * c_theta_ddot->size[1];
  c_theta_ddot->size[0] = setFilt->size[1];
  c_theta_ddot->size[1] = m;
  emxEnsureCapacity_real_T(c_theta_ddot, i0);
  for (i0 = 0; i0 < m; i0++) {
    loop_ub = setFilt->size[1];
    for (k = 0; k < loop_ub; k++) {
      c_theta_ddot->data[k + c_theta_ddot->size[0] * i0] = b_theta_ddot->data
        [((int)setFilt->data[setFilt->size[0] * k] + b_theta_ddot->size[0] * i0)
        - 1];
    }
  }

  i0 = b_theta_ddot->size[0] * b_theta_ddot->size[1];
  b_theta_ddot->size[0] = c_theta_ddot->size[0];
  b_theta_ddot->size[1] = c_theta_ddot->size[1];
  emxEnsureCapacity_real_T(b_theta_ddot, i0);
  loop_ub = c_theta_ddot->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    i = c_theta_ddot->size[0];
    for (k = 0; k < i; k++) {
      b_theta_ddot->data[k + b_theta_ddot->size[0] * i0] = c_theta_ddot->data[k
        + c_theta_ddot->size[0] * i0];
    }
  }

  /* 'dynIdenf:187' tau = tau(setFilt,:); */
  m = b_tau->size[1];
  i0 = d_x->size[0] * d_x->size[1];
  d_x->size[0] = setFilt->size[1];
  d_x->size[1] = m;
  emxEnsureCapacity_real_T(d_x, i0);
  for (i0 = 0; i0 < m; i0++) {
    loop_ub = setFilt->size[1];
    for (k = 0; k < loop_ub; k++) {
      d_x->data[k + d_x->size[0] * i0] = b_tau->data[((int)setFilt->data
        [setFilt->size[0] * k] + b_tau->size[0] * i0) - 1];
    }
  }

  i0 = b_tau->size[0] * b_tau->size[1];
  b_tau->size[0] = d_x->size[0];
  b_tau->size[1] = d_x->size[1];
  emxEnsureCapacity_real_T(b_tau, i0);
  loop_ub = d_x->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    i = d_x->size[0];
    for (k = 0; k < i; k++) {
      b_tau->data[k + b_tau->size[0] * i0] = d_x->data[k + d_x->size[0] * i0];
    }
  }

  /*  parameter identification by fordKinematics() and lsqSVD(). */
  /*  calculate the K according to forward kinematics */
  /* 'dynIdenf:123' K = fordKinematics(theta(setResample,:), theta_dot(setResample,:), theta_ddot(setResample,:), g, a, alpha, d, Rpk); */
  loop_ub = x->size[1];
  i0 = d_x->size[0] * d_x->size[1];
  d_x->size[0] = setResample->size[1];
  d_x->size[1] = loop_ub;
  emxEnsureCapacity_real_T(d_x, i0);
  for (i0 = 0; i0 < loop_ub; i0++) {
    i = setResample->size[1];
    for (k = 0; k < i; k++) {
      d_x->data[k + d_x->size[0] * i0] = x->data[((int)setResample->
        data[setResample->size[0] * k] + x->size[0] * i0) - 1];
    }
  }

  loop_ub = b_theta_dot->size[1];
  i0 = c_theta_dot->size[0] * c_theta_dot->size[1];
  c_theta_dot->size[0] = setResample->size[1];
  c_theta_dot->size[1] = loop_ub;
  emxEnsureCapacity_real_T(c_theta_dot, i0);
  for (i0 = 0; i0 < loop_ub; i0++) {
    i = setResample->size[1];
    for (k = 0; k < i; k++) {
      c_theta_dot->data[k + c_theta_dot->size[0] * i0] = b_theta_dot->data[((int)
        setResample->data[setResample->size[0] * k] + b_theta_dot->size[0] * i0)
        - 1];
    }
  }

  loop_ub = b_theta_ddot->size[1];
  i0 = c_theta_ddot->size[0] * c_theta_ddot->size[1];
  c_theta_ddot->size[0] = setResample->size[1];
  c_theta_ddot->size[1] = loop_ub;
  emxEnsureCapacity_real_T(c_theta_ddot, i0);
  for (i0 = 0; i0 < loop_ub; i0++) {
    i = setResample->size[1];
    for (k = 0; k < i; k++) {
      c_theta_ddot->data[k + c_theta_ddot->size[0] * i0] = b_theta_ddot->data
        [((int)setResample->data[setResample->size[0] * k] + b_theta_ddot->size
          [0] * i0) - 1];
    }
  }

  emxInit_real_T(&K, 2);
  fordKinematics(d_x, c_theta_dot, c_theta_ddot, g, a, alpha, d, Rpk, K);

  /*  resolve the dynamic parameters by least-square method with SVD */
  /*  decomposition, and record the rank of matrix K */
  /* 'dynIdenf:127' [phi, rk] = lsqSVD(K, tau(setResample,:), phi_pre, count, phi_r, noise_err, cond_max, lambda); */
  loop_ub = b_tau->size[1];
  i0 = d_x->size[0] * d_x->size[1];
  d_x->size[0] = setResample->size[1];
  d_x->size[1] = loop_ub;
  emxEnsureCapacity_real_T(d_x, i0);
  for (i0 = 0; i0 < loop_ub; i0++) {
    i = setResample->size[1];
    for (k = 0; k < i; k++) {
      d_x->data[k + d_x->size[0] * i0] = b_tau->data[((int)setResample->
        data[setResample->size[0] * k] + b_tau->size[0] * i0) - 1];
    }
  }

  lsqSVD(K, d_x, phi_pre, count, phi_r, noise_err, cond_max, lambda, phi, &ndbl);

  /*  evaluate the torque by fordKinematics() and evalTorque() */
  /*  calculate the K according to forward kinematics */
  /* 'dynIdenf:131' K = fordKinematics(theta(setEvalTorq,:), theta_dot(setEvalTorq,:), theta_ddot(setEvalTorq,:), g, a, alpha, d, Rpk); */
  loop_ub = x->size[1];
  i0 = d_x->size[0] * d_x->size[1];
  d_x->size[0] = setEvalTorq->size[1];
  d_x->size[1] = loop_ub;
  emxEnsureCapacity_real_T(d_x, i0);
  for (i0 = 0; i0 < loop_ub; i0++) {
    i = setEvalTorq->size[1];
    for (k = 0; k < i; k++) {
      d_x->data[k + d_x->size[0] * i0] = x->data[((int)setEvalTorq->
        data[setEvalTorq->size[0] * k] + x->size[0] * i0) - 1];
    }
  }

  loop_ub = b_theta_dot->size[1];
  i0 = c_theta_dot->size[0] * c_theta_dot->size[1];
  c_theta_dot->size[0] = setEvalTorq->size[1];
  c_theta_dot->size[1] = loop_ub;
  emxEnsureCapacity_real_T(c_theta_dot, i0);
  for (i0 = 0; i0 < loop_ub; i0++) {
    i = setEvalTorq->size[1];
    for (k = 0; k < i; k++) {
      c_theta_dot->data[k + c_theta_dot->size[0] * i0] = b_theta_dot->data[((int)
        setEvalTorq->data[setEvalTorq->size[0] * k] + b_theta_dot->size[0] * i0)
        - 1];
    }
  }

  emxFree_real_T(&b_theta_dot);
  loop_ub = b_theta_ddot->size[1];
  i0 = c_theta_ddot->size[0] * c_theta_ddot->size[1];
  c_theta_ddot->size[0] = setEvalTorq->size[1];
  c_theta_ddot->size[1] = loop_ub;
  emxEnsureCapacity_real_T(c_theta_ddot, i0);
  for (i0 = 0; i0 < loop_ub; i0++) {
    i = setEvalTorq->size[1];
    for (k = 0; k < i; k++) {
      c_theta_ddot->data[k + c_theta_ddot->size[0] * i0] = b_theta_ddot->data
        [((int)setEvalTorq->data[setEvalTorq->size[0] * k] + b_theta_ddot->size
          [0] * i0) - 1];
    }
  }

  emxFree_real_T(&b_theta_ddot);
  fordKinematics(d_x, c_theta_dot, c_theta_ddot, g, a, alpha, d, Rpk, K);

  /*  evaluate the torques and evaluate the error of the torque */
  /* 'dynIdenf:134' tau_filt = tau(setEvalTorq,:); */
  loop_ub = b_tau->size[1];
  i0 = tau_filt->size[0] * tau_filt->size[1];
  tau_filt->size[0] = setEvalTorq->size[1];
  tau_filt->size[1] = loop_ub;
  emxEnsureCapacity_real_T(tau_filt, i0);
  emxFree_real_T(&c_theta_ddot);
  emxFree_real_T(&c_theta_dot);
  emxFree_real_T(&d_x);
  for (i0 = 0; i0 < loop_ub; i0++) {
    i = setEvalTorq->size[1];
    for (k = 0; k < i; k++) {
      tau_filt->data[k + tau_filt->size[0] * i0] = b_tau->data[((int)
        setEvalTorq->data[setEvalTorq->size[0] * k] + b_tau->size[0] * i0) - 1];
    }
  }

  emxFree_real_T(&b_tau);

  /* 'dynIdenf:135' [tau_pos, tau_pre, distErr, convErr] = evalTorque(K, phi, phi_pre, tau_filt, segErr); */
  evalTorque(K, phi, phi_pre, tau_filt, segErr, tau_pos, tau_pre, x, setFilt);

  /*  pack and output the errors */
  /* 'dynIdenf:138' errs = {distErr, filtErr, convErr, rk}; */
  i0 = errs->f1->size[0] * errs->f1->size[1];
  errs->f1->size[0] = x->size[0];
  errs->f1->size[1] = x->size[1];
  emxEnsureCapacity_real_T(errs->f1, i0);
  loop_ub = x->size[1];
  emxFree_real_T(&K);
  for (i0 = 0; i0 < loop_ub; i0++) {
    i = x->size[0];
    for (k = 0; k < i; k++) {
      errs->f1->data[k + errs->f1->size[0] * i0] = x->data[k + x->size[0] * i0];
    }
  }

  emxFree_real_T(&x);
  i0 = errs->f2->size[0] * errs->f2->size[1];
  errs->f2->size[0] = 1;
  errs->f2->size[1] = filtErr->size[1];
  emxEnsureCapacity_real_T(errs->f2, i0);
  loop_ub = filtErr->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    errs->f2->data[errs->f2->size[0] * i0] = filtErr->data[filtErr->size[0] * i0];
  }

  emxFree_real_T(&filtErr);
  i0 = errs->f3->size[0] * errs->f3->size[1];
  errs->f3->size[0] = 1;
  errs->f3->size[1] = setFilt->size[1];
  emxEnsureCapacity_real_T(errs->f3, i0);
  loop_ub = setFilt->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    errs->f3->data[errs->f3->size[0] * i0] = setFilt->data[setFilt->size[0] * i0];
  }

  emxFree_real_T(&setFilt);
  errs->f4 = ndbl;

  /*  increment counter and update the previous dynamic parameters */
  /* 'dynIdenf:141' count = count + 1; */
  count++;

  /* 'dynIdenf:142' phi_pre = phi; */
  i0 = phi_pre->size[0];
  phi_pre->size[0] = phi->size[0];
  emxEnsureCapacity_real_T1(phi_pre, i0);
  loop_ub = phi->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    phi_pre->data[i0] = phi->data[i0];
  }
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void dynIdenf_initialize(void)
{
  count_not_empty = false;
  dynIdenf_init();
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void dynIdenf_terminate(void)
{
  dynIdenf_free();
}

/*
 * File trailer for dynIdenf.c
 *
 * [EOF]
 */
