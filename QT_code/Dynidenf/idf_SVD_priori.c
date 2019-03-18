/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: idf_SVD_priori.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 24-Jul-2018 20:21:10
 */

/* Include Files */
#include <math.h>
#include "rt_nonfinite.h"
#include "idf_SVD_priori.h"
#include "matlab_emxutil.h"

/* Function Declarations */
static void b_abs(const emxArray_real_T *x, emxArray_real_T *y);
static void b_rdivide(const emxArray_real_T *x, const emxArray_real_T *y,
                      emxArray_real_T *z);
static void b_sqrt(double *x);
static void b_svd(const emxArray_real_T *A, emxArray_real_T *U, emxArray_real_T *
                  s, emxArray_real_T *V);
static void b_xaxpy(int n, double a, const emxArray_real_T *x, int ix0,
                    emxArray_real_T *y, int iy0);
static double b_xnrm2(int n, const emxArray_real_T *x, int ix0);
static void c_xaxpy(int n, double a, const emxArray_real_T *x, int ix0,
                    emxArray_real_T *y, int iy0);
static void diag(const emxArray_real_T *v, emxArray_real_T *d);
static void eye(double varargin_1, emxArray_real_T *I);
static void rdivide(const emxArray_real_T *x, const emxArray_real_T *y,
                    emxArray_real_T *z);
static void repmat(const emxArray_real_T *a, double varargin_1, emxArray_real_T *
                   b);
static void svd(const emxArray_real_T *A, emxArray_real_T *U, emxArray_real_T *S,
                emxArray_real_T *V);
static void vecnorm(const emxArray_real_T *x, emxArray_real_T *y);
static void xaxpy(int n, double a, int ix0, emxArray_real_T *y, int iy0);
static double xdotc(int n, const emxArray_real_T *x, int ix0, const
                    emxArray_real_T *y, int iy0);
static double xnrm2(int n, const emxArray_real_T *x, int ix0);
static void xrot(int n, emxArray_real_T *x, int ix0, int iy0, double c, double s);
static void xrotg(double *a, double *b, double *c, double *s);
static void xscal(int n, double a, emxArray_real_T *x, int ix0);
static void xswap(int n, emxArray_real_T *x, int ix0, int iy0);

/* Function Definitions */

/*
 * Arguments    : const emxArray_real_T *x
 *                emxArray_real_T *y
 * Return Type  : void
 */
static void b_abs(const emxArray_real_T *x, emxArray_real_T *y)
{
  int nx;
  int k;
  unsigned int uv0[2];
  nx = x->size[0] * x->size[1];
  for (k = 0; k < 2; k++) {
    uv0[k] = (unsigned int)x->size[k];
  }

  k = y->size[0] * y->size[1];
  y->size[0] = (int)uv0[0];
  y->size[1] = (int)uv0[1];
  emxEnsureCapacity_real_T(y, k);
  for (k = 0; k < nx; k++) {
    y->data[k] = fabs(x->data[k]);
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
  int i3;
  int loop_ub;
  i3 = z->size[0];
  z->size[0] = x->size[0];
  emxEnsureCapacity_real_T1(z, i3);
  loop_ub = x->size[0];
  for (i3 = 0; i3 < loop_ub; i3++) {
    z->data[i3] = x->data[i3] / y->data[i3];
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
 * Arguments    : const emxArray_real_T *A
 *                emxArray_real_T *U
 *                emxArray_real_T *s
 *                emxArray_real_T *V
 * Return Type  : void
 */
static void b_svd(const emxArray_real_T *A, emxArray_real_T *U, emxArray_real_T *
                  s, emxArray_real_T *V)
{
  emxArray_real_T *b_A;
  int m;
  int ns;
  int n;
  int p;
  int qs;
  int minnp;
  emxArray_real_T *b_s;
  emxArray_real_T *e;
  emxArray_real_T *work;
  int nrt;
  int nct;
  int q;
  int iter;
  int nmq;
  boolean_T apply_transform;
  double ztest0;
  int mm;
  double ztest;
  double snorm;
  boolean_T exitg1;
  double f;
  double scale;
  double sqds;
  double b;
  emxInit_real_T1(&b_A, 2);
  m = b_A->size[0] * b_A->size[1];
  b_A->size[0] = A->size[0];
  b_A->size[1] = A->size[1];
  emxEnsureCapacity_real_T(b_A, m);
  ns = A->size[0] * A->size[1];
  for (m = 0; m < ns; m++) {
    b_A->data[m] = A->data[m];
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

  emxInit_real_T(&b_s, 1);
  m = b_s->size[0];
  b_s->size[0] = ns;
  emxEnsureCapacity_real_T1(b_s, m);
  for (m = 0; m < ns; m++) {
    b_s->data[m] = 0.0;
  }

  emxInit_real_T(&e, 1);
  ns = A->size[1];
  m = e->size[0];
  e->size[0] = ns;
  emxEnsureCapacity_real_T1(e, m);
  for (m = 0; m < ns; m++) {
    e->data[m] = 0.0;
  }

  emxInit_real_T(&work, 1);
  ns = A->size[0];
  m = work->size[0];
  work->size[0] = ns;
  emxEnsureCapacity_real_T1(work, m);
  for (m = 0; m < ns; m++) {
    work->data[m] = 0.0;
  }

  ns = A->size[0];
  qs = A->size[0];
  m = U->size[0] * U->size[1];
  U->size[0] = ns;
  U->size[1] = qs;
  emxEnsureCapacity_real_T(U, m);
  ns *= qs;
  for (m = 0; m < ns; m++) {
    U->data[m] = 0.0;
  }

  ns = A->size[1];
  qs = A->size[1];
  m = V->size[0] * V->size[1];
  V->size[0] = ns;
  V->size[1] = qs;
  emxEnsureCapacity_real_T(V, m);
  ns *= qs;
  for (m = 0; m < ns; m++) {
    V->data[m] = 0.0;
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
      ns = A->size[1] - 2;
    } else {
      ns = 0;
    }

    nrt = A->size[0];
    if (ns < nrt) {
      nrt = ns;
    }

    if (A->size[0] > 1) {
      ns = A->size[0] - 1;
    } else {
      ns = 0;
    }

    nct = A->size[1];
    if (ns < nct) {
      nct = ns;
    }

    if (nct > nrt) {
      m = nct;
    } else {
      m = nrt;
    }

    for (q = 0; q < m; q++) {
      iter = q + n * q;
      nmq = n - q;
      apply_transform = false;
      if (q + 1 <= nct) {
        ztest0 = xnrm2(nmq, b_A, iter + 1);
        if (ztest0 > 0.0) {
          apply_transform = true;
          if (b_A->data[iter] < 0.0) {
            ztest0 = -ztest0;
          }

          b_s->data[q] = ztest0;
          if (fabs(b_s->data[q]) >= 1.0020841800044864E-292) {
            xscal(nmq, 1.0 / b_s->data[q], b_A, iter + 1);
          } else {
            ns = iter + nmq;
            for (qs = iter; qs < ns; qs++) {
              b_A->data[qs] /= b_s->data[q];
            }
          }

          b_A->data[iter]++;
          b_s->data[q] = -b_s->data[q];
        } else {
          b_s->data[q] = 0.0;
        }
      }

      for (mm = q + 1; mm < p; mm++) {
        ns = q + n * mm;
        if (apply_transform) {
          ztest0 = -(xdotc(nmq, b_A, iter + 1, b_A, ns + 1) / b_A->data[q +
                     b_A->size[0] * q]);
          xaxpy(nmq, ztest0, iter + 1, b_A, ns + 1);
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
            ns = q + iter;
            for (qs = q + 1; qs < ns; qs++) {
              e->data[qs] *= ztest0;
            }
          } else {
            ns = q + iter;
            for (qs = q + 1; qs < ns; qs++) {
              e->data[qs] /= ztest0;
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
      iter = q + n * q;
      if (b_s->data[q] != 0.0) {
        for (mm = q + 1; mm < n; mm++) {
          ns = (q + n * mm) + 1;
          ztest0 = -(xdotc(nmq, U, iter + 1, U, ns) / U->data[iter]);
          xaxpy(nmq, ztest0, iter + 1, U, ns);
        }

        for (ns = q; ns < n; ns++) {
          U->data[ns + U->size[0] * q] = -U->data[ns + U->size[0] * q];
        }

        U->data[iter]++;
        for (ns = 1; ns <= q; ns++) {
          U->data[(ns + U->size[0] * q) - 1] = 0.0;
        }
      } else {
        for (ns = 1; ns <= n; ns++) {
          U->data[(ns + U->size[0] * q) - 1] = 0.0;
        }

        U->data[iter] = 1.0;
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
      ztest0 = fabs(b_s->data[ns]);
      ztest = fabs(e->data[ns]);
      if ((ztest0 > ztest) || rtIsNaN(ztest)) {
      } else {
        ztest0 = ztest;
      }

      if (!((snorm > ztest0) || rtIsNaN(ztest0))) {
        snorm = ztest0;
      }
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
        for (qs = m - 3; qs + 2 >= q + 1; qs--) {
          xrotg(&b_s->data[qs + 1], &f, &ztest0, &ztest);
          if (qs + 2 > q + 1) {
            f = -ztest * e->data[qs];
            e->data[qs] *= ztest0;
          }

          xrot(p, V, 1 + p * (qs + 1), 1 + p * (m - 1), ztest0, ztest);
        }
        break;

       case 2:
        f = e->data[q - 1];
        e->data[q - 1] = 0.0;
        for (qs = q; qs < m; qs++) {
          xrotg(&b_s->data[qs], &f, &ztest0, &ztest);
          f = -ztest * e->data[qs];
          e->data[qs] *= ztest0;
          xrot(n, U, 1 + n * qs, 1 + n * (q - 1), ztest0, ztest);
        }
        break;

       case 3:
        scale = fabs(b_s->data[m - 1]);
        ztest = fabs(b_s->data[m - 2]);
        if (!((scale > ztest) || rtIsNaN(ztest))) {
          scale = ztest;
        }

        ztest = fabs(e->data[m - 2]);
        if (!((scale > ztest) || rtIsNaN(ztest))) {
          scale = ztest;
        }

        ztest = fabs(b_s->data[q]);
        if (!((scale > ztest) || rtIsNaN(ztest))) {
          scale = ztest;
        }

        ztest = fabs(e->data[q]);
        if (!((scale > ztest) || rtIsNaN(ztest))) {
          scale = ztest;
        }

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
        for (qs = q + 1; qs < m; qs++) {
          xrotg(&f, &b, &ztest0, &ztest);
          if (qs > q + 1) {
            e->data[qs - 2] = f;
          }

          f = ztest0 * b_s->data[qs - 1] + ztest * e->data[qs - 1];
          e->data[qs - 1] = ztest0 * e->data[qs - 1] - ztest * b_s->data[qs - 1];
          b = ztest * b_s->data[qs];
          b_s->data[qs] *= ztest0;
          xrot(p, V, 1 + p * (qs - 1), 1 + p * qs, ztest0, ztest);
          b_s->data[qs - 1] = f;
          xrotg(&b_s->data[qs - 1], &b, &ztest0, &ztest);
          f = ztest0 * e->data[qs - 1] + ztest * b_s->data[qs];
          b_s->data[qs] = -ztest * e->data[qs - 1] + ztest0 * b_s->data[qs];
          b = ztest * e->data[qs];
          e->data[qs] *= ztest0;
          if (qs < n) {
            xrot(n, U, 1 + n * (qs - 1), 1 + n * qs, ztest0, ztest);
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

        ns = q + 1;
        while ((q + 1 < mm) && (b_s->data[q] < b_s->data[ns])) {
          ztest = b_s->data[q];
          b_s->data[q] = b_s->data[ns];
          b_s->data[ns] = ztest;
          if (q + 1 < p) {
            xswap(p, V, 1 + p * q, 1 + p * (q + 1));
          }

          if (q + 1 < n) {
            xswap(n, U, 1 + n * q, 1 + n * (q + 1));
          }

          q = ns;
          ns++;
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
  for (qs = 0; qs < minnp; qs++) {
    s->data[qs] = b_s->data[qs];
  }

  emxFree_real_T(&b_s);
}

/*
 * Arguments    : int n
 *                double a
 *                const emxArray_real_T *x
 *                int ix0
 *                emxArray_real_T *y
 *                int iy0
 * Return Type  : void
 */
static void b_xaxpy(int n, double a, const emxArray_real_T *x, int ix0,
                    emxArray_real_T *y, int iy0)
{
  int ix;
  int iy;
  int k;
  if ((n < 1) || (a == 0.0)) {
  } else {
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k < n; k++) {
      y->data[iy] += a * x->data[ix];
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
 * Arguments    : int n
 *                double a
 *                const emxArray_real_T *x
 *                int ix0
 *                emxArray_real_T *y
 *                int iy0
 * Return Type  : void
 */
static void c_xaxpy(int n, double a, const emxArray_real_T *x, int ix0,
                    emxArray_real_T *y, int iy0)
{
  int ix;
  int iy;
  int k;
  if ((n < 1) || (a == 0.0)) {
  } else {
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k < n; k++) {
      y->data[iy] += a * x->data[ix];
      ix++;
      iy++;
    }
  }
}

/*
 * Arguments    : const emxArray_real_T *v
 *                emxArray_real_T *d
 * Return Type  : void
 */
static void diag(const emxArray_real_T *v, emxArray_real_T *d)
{
  int u0;
  int u1;
  if ((v->size[0] == 1) && (v->size[1] == 1)) {
    u0 = d->size[0];
    d->size[0] = 1;
    emxEnsureCapacity_real_T1(d, u0);
    d->data[0] = v->data[0];
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

    u0 = d->size[0];
    d->size[0] = u1;
    emxEnsureCapacity_real_T1(d, u0);
    for (u0 = 0; u0 < u1; u0++) {
      d->data[u0] = v->data[u0 + v->size[0] * u0];
    }
  }
}

/*
 * Arguments    : double varargin_1
 *                emxArray_real_T *I
 * Return Type  : void
 */
static void eye(double varargin_1, emxArray_real_T *I)
{
  double t;
  int k;
  int loop_ub;
  if (varargin_1 < 0.0) {
    t = 0.0;
  } else {
    t = varargin_1;
  }

  k = I->size[0] * I->size[1];
  I->size[0] = (int)t;
  I->size[1] = (int)t;
  emxEnsureCapacity_real_T(I, k);
  loop_ub = (int)t * (int)t;
  for (k = 0; k < loop_ub; k++) {
    I->data[k] = 0.0;
  }

  if ((int)t > 0) {
    for (k = 0; k < (int)t; k++) {
      I->data[k + I->size[0] * k] = 1.0;
    }
  }
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
  int i2;
  int loop_ub;
  i2 = z->size[0] * z->size[1];
  z->size[0] = x->size[0];
  z->size[1] = x->size[1];
  emxEnsureCapacity_real_T(z, i2);
  loop_ub = x->size[0] * x->size[1];
  for (i2 = 0; i2 < loop_ub; i2++) {
    z->data[i2] = x->data[i2] / y->data[i2];
  }
}

/*
 * Arguments    : const emxArray_real_T *a
 *                double varargin_1
 *                emxArray_real_T *b
 * Return Type  : void
 */
static void repmat(const emxArray_real_T *a, double varargin_1, emxArray_real_T *
                   b)
{
  int outsize_idx_1;
  int ibmat;
  int itilerow;
  outsize_idx_1 = a->size[1];
  ibmat = b->size[0] * b->size[1];
  b->size[0] = (int)varargin_1;
  b->size[1] = outsize_idx_1;
  emxEnsureCapacity_real_T(b, ibmat);
  for (outsize_idx_1 = 0; outsize_idx_1 < a->size[1]; outsize_idx_1++) {
    ibmat = outsize_idx_1 * (int)varargin_1;
    for (itilerow = 1; itilerow <= (int)varargin_1; itilerow++) {
      b->data[(ibmat + itilerow) - 1] = a->data[outsize_idx_1];
    }
  }
}

/*
 * Arguments    : const emxArray_real_T *A
 *                emxArray_real_T *U
 *                emxArray_real_T *S
 *                emxArray_real_T *V
 * Return Type  : void
 */
static void svd(const emxArray_real_T *A, emxArray_real_T *U, emxArray_real_T *S,
                emxArray_real_T *V)
{
  int nx;
  boolean_T p;
  int k;
  emxArray_real_T *s;
  emxArray_real_T *r4;
  unsigned int uv1[2];
  emxArray_real_T *U1;
  emxArray_real_T *V1;
  nx = A->size[0] * A->size[1];
  p = true;
  for (k = 0; k < nx; k++) {
    if (p && ((!rtIsInf(A->data[k])) && (!rtIsNaN(A->data[k])))) {
      p = true;
    } else {
      p = false;
    }
  }

  emxInit_real_T(&s, 1);
  if (p) {
    b_svd(A, U, s, V);
  } else {
    for (nx = 0; nx < 2; nx++) {
      uv1[nx] = (unsigned int)A->size[nx];
    }

    emxInit_real_T1(&r4, 2);
    nx = r4->size[0] * r4->size[1];
    r4->size[0] = (int)uv1[0];
    r4->size[1] = (int)uv1[1];
    emxEnsureCapacity_real_T(r4, nx);
    k = (int)uv1[0] * (int)uv1[1];
    for (nx = 0; nx < k; nx++) {
      r4->data[nx] = 0.0;
    }

    emxInit_real_T1(&U1, 2);
    emxInit_real_T1(&V1, 2);
    b_svd(r4, U1, s, V1);
    nx = U->size[0] * U->size[1];
    U->size[0] = U1->size[0];
    U->size[1] = U1->size[1];
    emxEnsureCapacity_real_T(U, nx);
    k = U1->size[0] * U1->size[1];
    emxFree_real_T(&r4);
    emxFree_real_T(&U1);
    for (nx = 0; nx < k; nx++) {
      U->data[nx] = rtNaN;
    }

    k = s->size[0];
    nx = s->size[0];
    s->size[0] = k;
    emxEnsureCapacity_real_T1(s, nx);
    for (nx = 0; nx < k; nx++) {
      s->data[nx] = rtNaN;
    }

    nx = V->size[0] * V->size[1];
    V->size[0] = V1->size[0];
    V->size[1] = V1->size[1];
    emxEnsureCapacity_real_T(V, nx);
    k = V1->size[0] * V1->size[1];
    emxFree_real_T(&V1);
    for (nx = 0; nx < k; nx++) {
      V->data[nx] = rtNaN;
    }
  }

  nx = S->size[0] * S->size[1];
  S->size[0] = U->size[1];
  S->size[1] = V->size[1];
  emxEnsureCapacity_real_T(S, nx);
  k = U->size[1] * V->size[1];
  for (nx = 0; nx < k; nx++) {
    S->data[nx] = 0.0;
  }

  for (k = 0; k < s->size[0]; k++) {
    S->data[k + S->size[0] * k] = s->data[k];
  }

  emxFree_real_T(&s);
}

/*
 * Arguments    : const emxArray_real_T *x
 *                emxArray_real_T *y
 * Return Type  : void
 */
static void vecnorm(const emxArray_real_T *x, emxArray_real_T *y)
{
  int nrows;
  int j;
  unsigned int sizey[2];
  nrows = x->size[0];
  for (j = 0; j < 2; j++) {
    sizey[j] = (unsigned int)x->size[j];
  }

  j = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = (int)sizey[1];
  emxEnsureCapacity_real_T(y, j);
  for (j = 0; j < x->size[1]; j++) {
    y->data[j] = xnrm2(nrows, x, 1 + j * nrows);
  }
}

/*
 * Arguments    : int n
 *                double a
 *                int ix0
 *                emxArray_real_T *y
 *                int iy0
 * Return Type  : void
 */
static void xaxpy(int n, double a, int ix0, emxArray_real_T *y, int iy0)
{
  int ix;
  int iy;
  int k;
  if ((n < 1) || (a == 0.0)) {
  } else {
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k < n; k++) {
      y->data[iy] += a * y->data[ix];
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
  double d;
  int ix;
  int iy;
  int k;
  d = 0.0;
  if (!(n < 1)) {
    ix = ix0;
    iy = iy0;
    for (k = 1; k <= n; k++) {
      d += x->data[ix - 1] * y->data[iy - 1];
      ix++;
      iy++;
    }
  }

  return d;
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
 * Arguments    : double *a
 *                double *b
 *                double *c
 *                double *s
 * Return Type  : void
 */
static void xrotg(double *a, double *b, double *c, double *s)
{
  double roe;
  double absa;
  double absb;
  double scale;
  double ads;
  double bds;
  roe = *b;
  absa = fabs(*a);
  absb = fabs(*b);
  if (absa > absb) {
    roe = *a;
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

    *c = *a / scale;
    *s = *b / scale;
    if (absa > absb) {
      *b = *s;
    } else if (*c != 0.0) {
      *b = 1.0 / *c;
    } else {
      *b = 1.0;
    }
  }

  *a = scale;
}

/*
 * Arguments    : int n
 *                double a
 *                emxArray_real_T *x
 *                int ix0
 * Return Type  : void
 */
static void xscal(int n, double a, emxArray_real_T *x, int ix0)
{
  int i4;
  int k;
  i4 = (ix0 + n) - 1;
  for (k = ix0; k <= i4; k++) {
    x->data[k - 1] *= a;
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
 * Arguments    : const emxArray_real_T *phi_pre
 *                const emxArray_real_T *phi_r
 *                const emxArray_real_T *tau
 *                const emxArray_real_T *K
 *                const double opts[4]
 *                emxArray_real_T *phi_pos
 *                double rank_cond[2]
 * Return Type  : void
 */
void idf_SVD_priori(const emxArray_real_T *phi_pre, const emxArray_real_T *phi_r,
                    const emxArray_real_T *tau, const emxArray_real_T *K, const
                    double opts[4], emxArray_real_T *phi_pos, double rank_cond[2])
{ 
  rt_InitInfAndNaN((size_t)0);
  emxArray_real_T *S;
  emxArray_boolean_T *r0;
  int i0;
  int loop_ub;
  emxArray_real_T *b_K;
  emxArray_boolean_T *r1;
  emxArray_real_T *scv;
  emxArray_real_T *b_scv;
  emxArray_real_T *b;
  emxArray_real_T *varargin_1;
  int idx;
  boolean_T empty_non_axis_sizes;
  int m;
  int result;
  emxArray_real_T *b_tau;
  emxArray_real_T *b_varargin_1;
  int i1;
  emxArray_real_T *s;
  double ex;
  boolean_T exitg1;
  emxArray_boolean_T *e1;
  double y;
  emxArray_real_T *b_b;
  emxArray_real_T *b_S;
  emxArray_int32_T *r2;
  emxArray_real_T *c_tau;
  emxArray_int32_T *r3;
  double b_ex;
  emxInit_real_T1(&S, 2);
  emxInit_boolean_T1(&r0, 2);

  /*  noise_err = 1e-6; */
  /*  cond_max = 100; */
  /*  lambda = 1; */
  /*  step_ratio = 1; */
  b_abs(K, S);
  i0 = r0->size[0] * r0->size[1];
  r0->size[0] = S->size[0];
  r0->size[1] = S->size[1];
  emxEnsureCapacity_boolean_T(r0, i0);
  loop_ub = S->size[0] * S->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    r0->data[i0] = (S->data[i0] > opts[0]);
  }

  emxInit_real_T1(&b_K, 2);
  i0 = b_K->size[0] * b_K->size[1];
  b_K->size[0] = K->size[0];
  b_K->size[1] = K->size[1];
  emxEnsureCapacity_real_T(b_K, i0);
  loop_ub = K->size[0] * K->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_K->data[i0] = K->data[i0] * (double)r0->data[i0];
  }

  emxFree_boolean_T(&r0);
  emxInit_boolean_T1(&r1, 2);
  emxInit_real_T1(&scv, 2);
  vecnorm(b_K, scv);
  i0 = r1->size[0] * r1->size[1];
  r1->size[0] = 1;
  r1->size[1] = scv->size[1];
  emxEnsureCapacity_boolean_T(r1, i0);
  loop_ub = scv->size[0] * scv->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    r1->data[i0] = (scv->data[i0] == 0.0);
  }

  emxInit_real_T(&b_scv, 1);
  vecnorm(b_K, scv);
  i0 = b_scv->size[0];
  b_scv->size[0] = r1->size[1];
  emxEnsureCapacity_real_T1(b_scv, i0);
  loop_ub = r1->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_scv->data[i0] = (double)r1->data[r1->size[0] * i0] + scv->data[scv->size[0]
      * i0];
  }

  emxFree_boolean_T(&r1);
  emxInit_real_T1(&b, 2);
  eye(K->size[1], b);
  loop_ub = b->size[0] * b->size[1] - 1;
  i0 = b->size[0] * b->size[1];
  emxEnsureCapacity_real_T(b, i0);
  for (i0 = 0; i0 <= loop_ub; i0++) {
    b->data[i0] *= opts[2];
  }

  i0 = scv->size[0] * scv->size[1];
  scv->size[0] = 1;
  scv->size[1] = b_scv->size[0];
  emxEnsureCapacity_real_T(scv, i0);
  loop_ub = b_scv->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    scv->data[scv->size[0] * i0] = b_scv->data[i0];
  }

  emxInit_real_T1(&varargin_1, 2);
  repmat(scv, K->size[0], S);
  rdivide(b_K, S, varargin_1);
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
    m = varargin_1->size[0];
  } else {
    m = 0;
  }

  if (empty_non_axis_sizes || (!((b->size[0] == 0) || (b->size[1] == 0)))) {
    result = b->size[0];
  } else {
    result = 0;
  }

  i0 = phi_pos->size[0];
  phi_pos->size[0] = phi_pre->size[0];
  emxEnsureCapacity_real_T1(phi_pos, i0);
  loop_ub = phi_pre->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    phi_pos->data[i0] = phi_pre->data[i0] * b_scv->data[i0];
  }

  emxInit_real_T(&b_tau, 1);
  i0 = b_tau->size[0];
  b_tau->size[0] = tau->size[0] + phi_r->size[0];
  emxEnsureCapacity_real_T1(b_tau, i0);
  loop_ub = tau->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_tau->data[i0] = tau->data[i0];
  }

  loop_ub = phi_r->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_tau->data[i0 + tau->size[0]] = opts[2] * phi_r->data[i0] * b_scv->data[i0];
  }

  emxInit_real_T1(&b_varargin_1, 2);
  i0 = b_varargin_1->size[0] * b_varargin_1->size[1];
  b_varargin_1->size[0] = m + result;
  b_varargin_1->size[1] = idx;
  emxEnsureCapacity_real_T(b_varargin_1, i0);
  for (i0 = 0; i0 < idx; i0++) {
    for (i1 = 0; i1 < m; i1++) {
      b_varargin_1->data[i1 + b_varargin_1->size[0] * i0] = varargin_1->data[i1
        + m * i0];
    }
  }

  emxFree_real_T(&varargin_1);
  for (i0 = 0; i0 < idx; i0++) {
    for (i1 = 0; i1 < result; i1++) {
      b_varargin_1->data[(i1 + m) + b_varargin_1->size[0] * i0] = b->data[i1 +
        result * i0];
    }
  }

  emxInit_real_T(&s, 1);
  svd(b_varargin_1, b_K, S, b);
  diag(S, s);
  emxFree_real_T(&b_varargin_1);
  if (s->size[0] <= 2) {
    if (s->size[0] == 0) {
      ex = rtNaN;
    } else if (s->size[0] == 1) {
      ex = s->data[0];
    } else if ((s->data[0] < s->data[1]) || (rtIsNaN(s->data[0]) && (!rtIsNaN
                 (s->data[1])))) {
      ex = s->data[1];
    } else {
      ex = s->data[0];
    }
  } else {
    if (!rtIsNaN(s->data[0])) {
      idx = 1;
    } else {
      idx = 0;
      result = 2;
      exitg1 = false;
      while ((!exitg1) && (result <= s->size[0])) {
        if (!rtIsNaN(s->data[result - 1])) {
          idx = result;
          exitg1 = true;
        } else {
          result++;
        }
      }
    }

    if (idx == 0) {
      ex = s->data[0];
    } else {
      ex = s->data[idx - 1];
      while (idx + 1 <= s->size[0]) {
        if (ex < s->data[idx]) {
          ex = s->data[idx];
        }

        idx++;
      }
    }
  }

  emxInit_boolean_T(&e1, 1);
  y = ex / opts[1];
  i0 = e1->size[0];
  e1->size[0] = s->size[0];
  emxEnsureCapacity_boolean_T1(e1, i0);
  loop_ub = s->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    e1->data[i0] = (s->data[i0] >= y);
  }

  i0 = S->size[0] * S->size[1];
  S->size[0] = b_K->size[1];
  S->size[1] = b_K->size[0];
  emxEnsureCapacity_real_T(S, i0);
  loop_ub = b_K->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    idx = b_K->size[1];
    for (i1 = 0; i1 < idx; i1++) {
      S->data[i1 + S->size[0] * i0] = b_K->data[i0 + b_K->size[0] * i1];
    }
  }

  emxFree_real_T(&b_K);
  emxInit_real_T(&b_b, 1);
  i0 = b_b->size[0];
  b_b->size[0] = b_tau->size[0];
  emxEnsureCapacity_real_T1(b_b, i0);
  loop_ub = b_tau->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_b->data[i0] = b_tau->data[i0];
  }

  emxInit_real_T(&b_S, 1);
  if ((S->size[1] == 1) || (b_tau->size[0] == 1)) {
    i0 = b_S->size[0];
    b_S->size[0] = S->size[0];
    emxEnsureCapacity_real_T1(b_S, i0);
    loop_ub = S->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_S->data[i0] = 0.0;
      idx = S->size[1];
      for (i1 = 0; i1 < idx; i1++) {
        b_S->data[i0] += S->data[i0 + S->size[0] * i1] * b_tau->data[i1];
      }
    }

    i0 = b_tau->size[0];
    b_tau->size[0] = b_S->size[0];
    emxEnsureCapacity_real_T1(b_tau, i0);
    loop_ub = b_S->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_tau->data[i0] = b_S->data[i0];
    }
  } else {
    m = S->size[0];
    i0 = b_tau->size[0];
    b_tau->size[0] = S->size[0];
    emxEnsureCapacity_real_T1(b_tau, i0);
    for (loop_ub = 1; loop_ub <= m; loop_ub++) {
      b_tau->data[loop_ub - 1] = 0.0;
    }

    for (result = 0; result < S->size[1]; result++) {
      if (b_b->data[result] != 0.0) {
        idx = result * m;
        for (loop_ub = 0; loop_ub < m; loop_ub++) {
          b_tau->data[loop_ub] += b_b->data[result] * S->data[idx + loop_ub];
        }
      }
    }
  }

  i0 = S->size[0] * S->size[1];
  S->size[0] = b->size[1];
  S->size[1] = b->size[0];
  emxEnsureCapacity_real_T(S, i0);
  loop_ub = b->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    idx = b->size[1];
    for (i1 = 0; i1 < idx; i1++) {
      S->data[i1 + S->size[0] * i0] = b->data[i0 + b->size[0] * i1];
    }
  }

  i0 = b_b->size[0];
  b_b->size[0] = phi_pos->size[0];
  emxEnsureCapacity_real_T1(b_b, i0);
  loop_ub = phi_pos->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_b->data[i0] = phi_pos->data[i0];
  }

  if ((S->size[1] == 1) || (phi_pos->size[0] == 1)) {
    i0 = b_S->size[0];
    b_S->size[0] = S->size[0];
    emxEnsureCapacity_real_T1(b_S, i0);
    loop_ub = S->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_S->data[i0] = 0.0;
      idx = S->size[1];
      for (i1 = 0; i1 < idx; i1++) {
        b_S->data[i0] += S->data[i0 + S->size[0] * i1] * phi_pos->data[i1];
      }
    }

    i0 = phi_pos->size[0];
    phi_pos->size[0] = b_S->size[0];
    emxEnsureCapacity_real_T1(phi_pos, i0);
    loop_ub = b_S->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      phi_pos->data[i0] = b_S->data[i0];
    }
  } else {
    m = S->size[0];
    i0 = phi_pos->size[0];
    phi_pos->size[0] = S->size[0];
    emxEnsureCapacity_real_T1(phi_pos, i0);
    for (loop_ub = 1; loop_ub <= m; loop_ub++) {
      phi_pos->data[loop_ub - 1] = 0.0;
    }

    for (result = 0; result < S->size[1]; result++) {
      if (b_b->data[result] != 0.0) {
        idx = result * m;
        for (loop_ub = 0; loop_ub < m; loop_ub++) {
          phi_pos->data[loop_ub] += b_b->data[result] * S->data[idx + loop_ub];
        }
      }
    }
  }

  emxFree_real_T(&S);
  m = e1->size[0] - 1;
  idx = 0;
  for (loop_ub = 0; loop_ub <= m; loop_ub++) {
    if (e1->data[loop_ub]) {
      idx++;
    }
  }

  emxInit_int32_T(&r2, 1);
  i0 = r2->size[0];
  r2->size[0] = idx;
  emxEnsureCapacity_int32_T(r2, i0);
  idx = 0;
  for (loop_ub = 0; loop_ub <= m; loop_ub++) {
    if (e1->data[loop_ub]) {
      r2->data[idx] = loop_ub + 1;
      idx++;
    }
  }

  emxInit_real_T(&c_tau, 1);
  i0 = c_tau->size[0];
  c_tau->size[0] = r2->size[0];
  emxEnsureCapacity_real_T1(c_tau, i0);
  loop_ub = r2->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    c_tau->data[i0] = b_tau->data[r2->data[i0] - 1];
  }

  emxFree_real_T(&b_tau);
  i0 = b_S->size[0];
  b_S->size[0] = r2->size[0];
  emxEnsureCapacity_real_T1(b_S, i0);
  loop_ub = r2->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_S->data[i0] = s->data[r2->data[i0] - 1];
  }

  emxFree_int32_T(&r2);
  b_rdivide(c_tau, b_S, b_b);
  m = e1->size[0];
  idx = 0;
  emxFree_real_T(&c_tau);
  for (loop_ub = 0; loop_ub < m; loop_ub++) {
    if (e1->data[loop_ub]) {
      phi_pos->data[loop_ub] = b_b->data[idx];
      idx++;
    }
  }

  i0 = b_b->size[0];
  b_b->size[0] = phi_pos->size[0];
  emxEnsureCapacity_real_T1(b_b, i0);
  loop_ub = phi_pos->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_b->data[i0] = phi_pos->data[i0];
  }

  if ((b->size[1] == 1) || (phi_pos->size[0] == 1)) {
    i0 = b_b->size[0];
    b_b->size[0] = b->size[0];
    emxEnsureCapacity_real_T1(b_b, i0);
    loop_ub = b->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_b->data[i0] = 0.0;
      idx = b->size[1];
      for (i1 = 0; i1 < idx; i1++) {
        b_b->data[i0] += b->data[i0 + b->size[0] * i1] * phi_pos->data[i1];
      }
    }

    i0 = phi_pos->size[0];
    phi_pos->size[0] = b_b->size[0];
    emxEnsureCapacity_real_T1(phi_pos, i0);
    loop_ub = b_b->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      phi_pos->data[i0] = b_b->data[i0];
    }
  } else {
    m = b->size[0];
    i0 = phi_pos->size[0];
    phi_pos->size[0] = b->size[0];
    emxEnsureCapacity_real_T1(phi_pos, i0);
    for (loop_ub = 1; loop_ub <= m; loop_ub++) {
      phi_pos->data[loop_ub - 1] = 0.0;
    }

    for (result = 0; result < b->size[1]; result++) {
      if (b_b->data[result] != 0.0) {
        idx = result * m;
        for (loop_ub = 0; loop_ub < m; loop_ub++) {
          phi_pos->data[loop_ub] += b_b->data[result] * b->data[idx + loop_ub];
        }
      }
    }
  }

  emxFree_real_T(&b_b);
  emxFree_real_T(&b);
  i0 = b_S->size[0];
  b_S->size[0] = phi_pos->size[0];
  emxEnsureCapacity_real_T1(b_S, i0);
  loop_ub = phi_pos->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_S->data[i0] = phi_pos->data[i0];
  }

  b_rdivide(b_S, b_scv, phi_pos);
  i0 = phi_pos->size[0];
  phi_pos->size[0] = phi_pre->size[0];
  emxEnsureCapacity_real_T1(phi_pos, i0);
  loop_ub = phi_pre->size[0];
  emxFree_real_T(&b_S);
  emxFree_real_T(&b_scv);
  for (i0 = 0; i0 < loop_ub; i0++) {
    phi_pos->data[i0] = phi_pre->data[i0] + opts[3] * (phi_pos->data[i0] -
      phi_pre->data[i0]);
  }

  m = e1->size[0] - 1;
  idx = 0;
  for (loop_ub = 0; loop_ub <= m; loop_ub++) {
    if (e1->data[loop_ub]) {
      idx++;
    }
  }

  emxInit_int32_T(&r3, 1);
  i0 = r3->size[0];
  r3->size[0] = idx;
  emxEnsureCapacity_int32_T(r3, i0);
  idx = 0;
  for (loop_ub = 0; loop_ub <= m; loop_ub++) {
    if (e1->data[loop_ub]) {
      r3->data[idx] = loop_ub + 1;
      idx++;
    }
  }

  if (e1->size[0] == 0) {
    y = 0.0;
  } else {
    y = e1->data[0];
    for (result = 2; result <= e1->size[0]; result++) {
      y += (double)e1->data[result - 1];
    }
  }

  emxFree_boolean_T(&e1);
  m = r3->size[0];
  if (r3->size[0] <= 2) {
    if (r3->size[0] == 0) {
      ex = rtNaN;
    } else if (r3->size[0] == 1) {
      ex = s->data[r3->data[0] - 1];
    } else if ((s->data[r3->data[0] - 1] < s->data[r3->data[1] - 1]) || (rtIsNaN
                (s->data[r3->data[0] - 1]) && (!rtIsNaN(s->data[r3->data[1] - 1]))))
    {
      ex = s->data[r3->data[1] - 1];
    } else {
      ex = s->data[r3->data[0] - 1];
    }
  } else {
    if (!rtIsNaN(s->data[r3->data[0] - 1])) {
      idx = 1;
    } else {
      idx = 0;
      result = 2;
      exitg1 = false;
      while ((!exitg1) && (result <= r3->size[0])) {
        if (!rtIsNaN(s->data[r3->data[result - 1] - 1])) {
          idx = result;
          exitg1 = true;
        } else {
          result++;
        }
      }
    }

    if (idx == 0) {
      ex = s->data[r3->data[0] - 1];
    } else {
      ex = s->data[r3->data[idx - 1] - 1];
      while (idx + 1 <= m) {
        if (ex < s->data[r3->data[idx] - 1]) {
          ex = s->data[r3->data[idx] - 1];
        }

        idx++;
      }
    }
  }

  m = r3->size[0];
  if (r3->size[0] <= 2) {
    if (r3->size[0] == 0) {
      b_ex = rtNaN;
    } else if (r3->size[0] == 1) {
      b_ex = s->data[r3->data[0] - 1];
    } else if ((s->data[r3->data[0] - 1] > s->data[r3->data[1] - 1]) || (rtIsNaN
                (s->data[r3->data[0] - 1]) && (!rtIsNaN(s->data[r3->data[1] - 1]))))
    {
      b_ex = s->data[r3->data[1] - 1];
    } else {
      b_ex = s->data[r3->data[0] - 1];
    }
  } else {
    if (!rtIsNaN(s->data[r3->data[0] - 1])) {
      idx = 1;
    } else {
      idx = 0;
      result = 2;
      exitg1 = false;
      while ((!exitg1) && (result <= r3->size[0])) {
        if (!rtIsNaN(s->data[r3->data[result - 1] - 1])) {
          idx = result;
          exitg1 = true;
        } else {
          result++;
        }
      }
    }

    if (idx == 0) {
      b_ex = s->data[r3->data[0] - 1];
    } else {
      b_ex = s->data[r3->data[idx - 1] - 1];
      while (idx + 1 <= m) {
        if (b_ex > s->data[r3->data[idx] - 1]) {
          b_ex = s->data[r3->data[idx] - 1];
        }

        idx++;
      }
    }
  }

  emxFree_int32_T(&r3);
  emxFree_real_T(&s);
  rank_cond[0] = y;
  rank_cond[1] = ex / b_ex;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void idf_SVD_priori_initialize(void)
{
  rt_InitInfAndNaN(8U);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void idf_SVD_priori_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for idf_SVD_priori.c
 *
 * [EOF]
 */
