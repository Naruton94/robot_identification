/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: Errestim.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 31-Aug-2018 13:44:40
 */

/* Include Files */
#include <math.h>
#include "rt_nonfinite.h"
#include "Errestim.h"
#include "matlab_emxutil.h"

/* Function Declarations */
static double rms(const emxArray_real_T *x);

/* Function Definitions */

/*
 * Arguments    : const emxArray_real_T *x
 * Return Type  : double
 */
static double rms(const emxArray_real_T *x)
{
  double y;
  emxArray_real_T *b_x;
  int k;
  int loop_ub;
  emxInit_real_T(&b_x, 1);
  k = b_x->size[0];
  b_x->size[0] = x->size[0];
  emxEnsureCapacity_real_T(b_x, k);
  loop_ub = x->size[0];
  for (k = 0; k < loop_ub; k++) {
    b_x->data[k] = x->data[k] * x->data[k];
  }

  if (b_x->size[0] == 0) {
    y = 0.0;
  } else {
    y = b_x->data[0];
    for (k = 2; k <= b_x->size[0]; k++) {
      y += b_x->data[k - 1];
    }
  }

  y = sqrt(y / (double)b_x->size[0]);
  emxFree_real_T(&b_x);
  return y;
}

/*
 * Arguments    : const emxArray_real_T *tau
 *                const emxArray_real_T *K
 *                const emxArray_real_T *phi_pre
 *                emxArray_real_T *tau_pre
 *                double err_rms_rel[2]
 * Return Type  : void
 */
void Errestim(const emxArray_real_T *tau, const emxArray_real_T *K, const
              emxArray_real_T *phi_pre, emxArray_real_T *tau_pre, double
              err_rms_rel[2])
{
  int m;
  int i0;
  int i;
  int k;
  emxArray_real_T *b_tau_pre;
  int aoffset;
  double rms_err;
  double d0;
  if ((K->size[1] == 1) || (phi_pre->size[0] == 1)) {
    i0 = tau_pre->size[0];
    tau_pre->size[0] = K->size[0];
    emxEnsureCapacity_real_T(tau_pre, i0);
    i = K->size[0];
    for (i0 = 0; i0 < i; i0++) {
      tau_pre->data[i0] = 0.0;
      k = K->size[1];
      for (aoffset = 0; aoffset < k; aoffset++) {
        tau_pre->data[i0] += K->data[i0 + K->size[0] * aoffset] * phi_pre->
          data[aoffset];
      }
    }
  } else {
    m = K->size[0];
    i0 = tau_pre->size[0];
    tau_pre->size[0] = K->size[0];
    emxEnsureCapacity_real_T(tau_pre, i0);
    for (i = 1; i <= m; i++) {
      tau_pre->data[i - 1] = 0.0;
    }

    for (k = 0; k < K->size[1]; k++) {
      if (phi_pre->data[k] != 0.0) {
        aoffset = k * m;
        for (i = 0; i < m; i++) {
          tau_pre->data[i] += phi_pre->data[k] * K->data[aoffset + i];
        }
      }
    }
  }

  emxInit_real_T(&b_tau_pre, 1);
  i0 = b_tau_pre->size[0];
  b_tau_pre->size[0] = tau_pre->size[0];
  emxEnsureCapacity_real_T(b_tau_pre, i0);
  i = tau_pre->size[0];
  for (i0 = 0; i0 < i; i0++) {
    b_tau_pre->data[i0] = tau_pre->data[i0] - tau->data[i0];
  }

  rms_err = rms(b_tau_pre);
  d0 = rms(tau);
  err_rms_rel[0] = rms_err;
  err_rms_rel[1] = rms_err / d0;
  emxFree_real_T(&b_tau_pre);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void Errestim_initialize(void)
{
  rt_InitInfAndNaN(8U);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void Errestim_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for Errestim.c
 *
 * [EOF]
 */
