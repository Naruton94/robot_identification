/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: fk.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 24-Jul-2018 17:05:24
 */

/* Include Files */
#include <math.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "fk.h"
#include "matlab_emxutil.h"

/* Function Definitions */

/*
 * Arguments    : const emxArray_real_T *theta
 *                const emxArray_real_T *theta_dot
 *                const emxArray_real_T *theta_ddot
 *                const emxArray_real_T *a_alpha_d
 *                const double g[3]
 *                emxArray_real_T *K
 * Return Type  : void
 */
void fk(const emxArray_real_T *theta, const emxArray_real_T *theta_dot, const
        emxArray_real_T *theta_ddot, const emxArray_real_T *a_alpha_d, const
        double g[3], emxArray_real_T *K)
{
  rt_InitInfAndNaN((size_t)0);
  emxArray_real_T *o_X_i;
  int i0;
  int loop_ub;
  double i_X_o[36];
  int i;
  double i_V_i[6];
  double i_A_i[6];
  emxArray_int8_T *r0;
  double i1_R_i[9];
  double i1_p_i[3];
  double dv0[9];
  int i1;
  double dv1[9];
  int i2;
  double i1_X_i[36];
  double i_V_i1_i[6];
  double i_X_i1[36];
  double b_i1_R_i[9];
  double i_A_i1_i[6];
  int b_i;
  double b_o_X_i[36];
  double i_V_i1[6];
  double b_i_X_i1[36];
  double b_i_V_i[18];
  double c_i_X_i1[6];
  double dv2[6];
  double b_i_A_i[18];
  double dv3[9];
  double i_Am_i[60];
  double x;
  double b_i_X_o[3];
  double dv4[18];
  int j;
  double c_i_X_o[12];
  double d_i_X_o[60];
  int i3;
  emxInit_real_T(&o_X_i, 3);
  i0 = o_X_i->size[0] * o_X_i->size[1] * o_X_i->size[2];
  o_X_i->size[0] = 6;
  o_X_i->size[1] = 6;
  o_X_i->size[2] = theta->size[1];
  emxEnsureCapacity_real_T(o_X_i, i0);
  loop_ub = 36 * theta->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    o_X_i->data[i0] = 0.0;
  }

  memset(&i_X_o[0], 0, 36U * sizeof(double));
  for (i = 0; i < 6; i++) {
    i_V_i[i] = 0.0;
    i_A_i[i] = 0.0;
  }

  i0 = K->size[0] * K->size[1];
  K->size[0] = theta->size[1];
  K->size[1] = 12 * theta->size[1];
  emxEnsureCapacity_real_T1(K, i0);
  loop_ub = theta->size[1] * (12 * theta->size[1]);
  for (i0 = 0; i0 < loop_ub; i0++) {
    K->data[i0] = 0.0;
  }

  i = 0;
  emxInit_int8_T(&r0, 1);
  while (i <= theta->size[1] - 1) {
    i1_R_i[0] = cos(theta->data[i]);
    i1_R_i[3] = -sin(theta->data[i]);
    i1_R_i[6] = 0.0;
    i1_R_i[1] = sin(theta->data[i]) * cos(a_alpha_d->data[1 + a_alpha_d->size[0]
      * i]);
    i1_R_i[4] = cos(theta->data[i]) * cos(a_alpha_d->data[1 + a_alpha_d->size[0]
      * i]);
    i1_R_i[7] = -sin(a_alpha_d->data[1 + a_alpha_d->size[0] * i]);
    i1_R_i[2] = sin(theta->data[i]) * sin(a_alpha_d->data[1 + a_alpha_d->size[0]
      * i]);
    i1_R_i[5] = cos(theta->data[i]) * sin(a_alpha_d->data[1 + a_alpha_d->size[0]
      * i]);
    i1_R_i[8] = cos(a_alpha_d->data[1 + a_alpha_d->size[0] * i]);
    i1_p_i[0] = a_alpha_d->data[a_alpha_d->size[0] * i];
    i1_p_i[1] = -a_alpha_d->data[2 + a_alpha_d->size[0] * (i + 1)] * sin
      (a_alpha_d->data[1 + a_alpha_d->size[0] * i]);
    i1_p_i[2] = a_alpha_d->data[2 + a_alpha_d->size[0] * (i + 1)] * cos
      (a_alpha_d->data[1 + a_alpha_d->size[0] * i]);
    dv0[0] = 0.0;
    dv0[3] = -i1_p_i[2];
    dv0[6] = i1_p_i[1];
    dv0[1] = i1_p_i[2];
    dv0[4] = 0.0;
    dv0[7] = -i1_p_i[0];
    dv0[2] = -i1_p_i[1];
    dv0[5] = i1_p_i[0];
    dv0[8] = 0.0;
    for (i0 = 0; i0 < 3; i0++) {
      for (i1 = 0; i1 < 3; i1++) {
        dv1[i0 + 3 * i1] = 0.0;
        for (i2 = 0; i2 < 3; i2++) {
          dv1[i0 + 3 * i1] += dv0[i0 + 3 * i2] * i1_R_i[i2 + 3 * i1];
        }

        i1_X_i[i1 + 6 * i0] = i1_R_i[i1 + 3 * i0];
        i1_X_i[i1 + 6 * (i0 + 3)] = 0.0;
      }
    }

    dv0[0] = 0.0;
    dv0[1] = -i1_p_i[2];
    dv0[2] = i1_p_i[1];
    dv0[3] = i1_p_i[2];
    dv0[4] = 0.0;
    dv0[5] = -i1_p_i[0];
    dv0[6] = -i1_p_i[1];
    dv0[7] = i1_p_i[0];
    dv0[8] = 0.0;
    for (i0 = 0; i0 < 3; i0++) {
      for (i1 = 0; i1 < 3; i1++) {
        i1_X_i[(i1 + 6 * i0) + 3] = dv1[i1 + 3 * i0];
        i1_X_i[(i1 + 6 * (i0 + 3)) + 3] = i1_R_i[i1 + 3 * i0];
        b_i1_R_i[i0 + 3 * i1] = 0.0;
        for (i2 = 0; i2 < 3; i2++) {
          b_i1_R_i[i0 + 3 * i1] += i1_R_i[i2 + 3 * i0] * dv0[i2 + 3 * i1];
        }

        i_X_i1[i1 + 6 * i0] = i1_R_i[i0 + 3 * i1];
        i_X_i1[i1 + 6 * (i0 + 3)] = 0.0;
      }
    }

    for (i0 = 0; i0 < 3; i0++) {
      for (i1 = 0; i1 < 3; i1++) {
        i_X_i1[(i1 + 6 * i0) + 3] = b_i1_R_i[i1 + 3 * i0];
        i_X_i1[(i1 + 6 * (i0 + 3)) + 3] = i1_R_i[i0 + 3 * i1];
      }
    }

    for (i0 = 0; i0 < 2; i0++) {
      i_V_i1_i[i0] = 0.0;
    }

    i_V_i1_i[2] = theta_dot->data[i];
    for (i0 = 0; i0 < 3; i0++) {
      i_V_i1_i[i0 + 3] = 0.0;
    }

    for (i0 = 0; i0 < 2; i0++) {
      i_A_i1_i[i0] = 0.0;
    }

    i_A_i1_i[2] = theta_ddot->data[i];
    for (i0 = 0; i0 < 3; i0++) {
      i_A_i1_i[i0 + 3] = 0.0;
    }

    if (1 + i == 1) {
      for (i0 = 0; i0 < 6; i0++) {
        for (i1 = 0; i1 < 6; i1++) {
          o_X_i->data[(i1 + o_X_i->size[0] * i0) + o_X_i->size[0] * o_X_i->size
            [1] * i] = i1_X_i[i1 + 6 * i0];
        }
      }

      memcpy(&i_X_o[0], &i_X_i1[0], 36U * sizeof(double));
      for (b_i = 0; b_i < 6; b_i++) {
        i_V_i[b_i] = i_V_i1_i[b_i];
        i_A_i[b_i] = i_A_i1_i[b_i];
      }
    } else {
      for (i0 = 0; i0 < 6; i0++) {
        for (i1 = 0; i1 < 6; i1++) {
          b_o_X_i[i0 + 6 * i1] = 0.0;
          for (i2 = 0; i2 < 6; i2++) {
            b_o_X_i[i0 + 6 * i1] += o_X_i->data[(i0 + o_X_i->size[0] * i2) +
              o_X_i->size[0] * o_X_i->size[1] * (i - 1)] * i1_X_i[i2 + 6 * i1];
          }
        }
      }

      for (i0 = 0; i0 < 6; i0++) {
        for (i1 = 0; i1 < 6; i1++) {
          o_X_i->data[(i1 + o_X_i->size[0] * i0) + o_X_i->size[0] * o_X_i->size
            [1] * i] = b_o_X_i[i1 + 6 * i0];
          b_i_X_i1[i0 + 6 * i1] = 0.0;
          for (i2 = 0; i2 < 6; i2++) {
            b_i_X_i1[i0 + 6 * i1] += i_X_i1[i0 + 6 * i2] * i_X_o[i2 + 6 * i1];
          }
        }
      }

      for (i0 = 0; i0 < 6; i0++) {
        i_V_i1[i0] = 0.0;
        for (i1 = 0; i1 < 6; i1++) {
          i_X_o[i1 + 6 * i0] = b_i_X_i1[i1 + 6 * i0];
          i_V_i1[i0] += i_X_i1[i0 + 6 * i1] * i_V_i[i1];
        }
      }

      for (b_i = 0; b_i < 6; b_i++) {
        i_V_i[b_i] = i_V_i1[b_i] + i_V_i1_i[b_i];
      }

      b_o_X_i[0] = 0.0;
      b_o_X_i[6] = -i_V_i1_i[2];
      b_o_X_i[12] = 0.0;
      b_o_X_i[1] = i_V_i1_i[2];
      b_o_X_i[7] = 0.0;
      b_o_X_i[13] = -0.0;
      b_o_X_i[2] = -0.0;
      b_o_X_i[8] = 0.0;
      b_o_X_i[14] = 0.0;
      for (i0 = 0; i0 < 3; i0++) {
        for (i1 = 0; i1 < 3; i1++) {
          b_o_X_i[i1 + 6 * (i0 + 3)] = 0.0;
        }
      }

      b_o_X_i[3] = 0.0;
      b_o_X_i[9] = -0.0;
      b_o_X_i[15] = 0.0;
      b_o_X_i[4] = 0.0;
      b_o_X_i[10] = 0.0;
      b_o_X_i[16] = -0.0;
      b_o_X_i[5] = -0.0;
      b_o_X_i[11] = 0.0;
      b_o_X_i[17] = 0.0;
      b_o_X_i[21] = 0.0;
      b_o_X_i[27] = -i_V_i1_i[2];
      b_o_X_i[33] = 0.0;
      b_o_X_i[22] = i_V_i1_i[2];
      b_o_X_i[28] = 0.0;
      b_o_X_i[34] = -0.0;
      b_o_X_i[23] = -0.0;
      b_o_X_i[29] = 0.0;
      b_o_X_i[35] = 0.0;
      for (i0 = 0; i0 < 6; i0++) {
        c_i_X_i1[i0] = 0.0;
        dv2[i0] = 0.0;
        for (i1 = 0; i1 < 6; i1++) {
          c_i_X_i1[i0] += i_X_i1[i0 + 6 * i1] * i_A_i[i1];
          dv2[i0] += b_o_X_i[i0 + 6 * i1] * i_V_i1[i1];
        }
      }

      for (i0 = 0; i0 < 6; i0++) {
        i_A_i[i0] = (c_i_X_i1[i0] - dv2[i0]) + i_A_i1_i[i0];
      }
    }

    dv0[0] = 0.0;
    dv0[3] = -i_V_i[2];
    dv0[6] = i_V_i[1];
    dv0[1] = i_V_i[2];
    dv0[4] = 0.0;
    dv0[7] = -i_V_i[0];
    dv0[2] = -i_V_i[1];
    dv0[5] = i_V_i[0];
    dv0[8] = 0.0;
    i1_R_i[0] = 0.0;
    i1_R_i[3] = -i_V_i[2];
    i1_R_i[6] = i_V_i[1];
    i1_R_i[1] = i_V_i[2];
    i1_R_i[4] = 0.0;
    i1_R_i[7] = -i_V_i[0];
    i1_R_i[2] = -i_V_i[1];
    i1_R_i[5] = i_V_i[0];
    i1_R_i[8] = 0.0;
    dv1[0] = 0.0;
    dv1[3] = -i_V_i[2];
    dv1[6] = i_V_i[1];
    dv1[1] = i_V_i[2];
    dv1[4] = 0.0;
    dv1[7] = -i_V_i[0];
    dv1[2] = -i_V_i[1];
    dv1[5] = i_V_i[0];
    dv1[8] = 0.0;
    b_i_V_i[0] = i_V_i[0];
    b_i_V_i[3] = i_V_i[1];
    b_i_V_i[6] = i_V_i[2];
    b_i_V_i[9] = 0.0;
    b_i_V_i[12] = 0.0;
    b_i_V_i[15] = 0.0;
    b_i_V_i[1] = 0.0;
    b_i_V_i[4] = i_V_i[0];
    b_i_V_i[7] = 0.0;
    b_i_V_i[10] = i_V_i[1];
    b_i_V_i[13] = i_V_i[2];
    b_i_V_i[16] = 0.0;
    b_i_V_i[2] = 0.0;
    b_i_V_i[5] = 0.0;
    b_i_V_i[8] = i_V_i[0];
    b_i_V_i[11] = 0.0;
    b_i_V_i[14] = i_V_i[1];
    b_i_V_i[17] = i_V_i[2];
    b_i_A_i[0] = i_A_i[0];
    b_i_A_i[3] = i_A_i[1];
    b_i_A_i[6] = i_A_i[2];
    b_i_A_i[9] = 0.0;
    b_i_A_i[12] = 0.0;
    b_i_A_i[15] = 0.0;
    b_i_A_i[1] = 0.0;
    b_i_A_i[4] = i_A_i[0];
    b_i_A_i[7] = 0.0;
    b_i_A_i[10] = i_A_i[1];
    b_i_A_i[13] = i_A_i[2];
    b_i_A_i[16] = 0.0;
    b_i_A_i[2] = 0.0;
    b_i_A_i[5] = 0.0;
    b_i_A_i[8] = i_A_i[0];
    b_i_A_i[11] = 0.0;
    b_i_A_i[14] = i_A_i[1];
    b_i_A_i[17] = i_A_i[2];
    dv3[0] = 0.0;
    dv3[3] = -i_A_i[2];
    dv3[6] = i_A_i[1];
    dv3[1] = i_A_i[2];
    dv3[4] = 0.0;
    dv3[7] = -i_A_i[0];
    dv3[2] = -i_A_i[1];
    dv3[5] = i_A_i[0];
    dv3[8] = 0.0;
    for (i0 = 0; i0 < 3; i0++) {
      x = 0.0;
      for (i1 = 0; i1 < 3; i1++) {
        x += dv0[i0 + 3 * i1] * i_V_i[3 + i1];
      }

      b_i_X_o[i0] = 0.0;
      for (i1 = 0; i1 < 3; i1++) {
        b_i_X_o[i0] += i_X_o[i0 + 6 * i1] * g[i1];
      }

      i1_p_i[i0] = (i_A_i[3 + i0] + x) - b_i_X_o[i0];
      for (i1 = 0; i1 < 6; i1++) {
        dv4[i0 + 3 * i1] = 0.0;
        for (i2 = 0; i2 < 3; i2++) {
          dv4[i0 + 3 * i1] += dv1[i0 + 3 * i2] * b_i_V_i[i2 + 3 * i1];
        }
      }

      for (i1 = 0; i1 < 3; i1++) {
        b_i1_R_i[i0 + 3 * i1] = 0.0;
        for (i2 = 0; i2 < 3; i2++) {
          b_i1_R_i[i0 + 3 * i1] += i1_R_i[i0 + 3 * i2] * i1_R_i[i2 + 3 * i1];
        }
      }

      i_Am_i[i0] = 0.0;
    }

    i_Am_i[6] = -0.0;
    i_Am_i[12] = -(-i1_p_i[2]);
    i_Am_i[18] = -i1_p_i[1];
    i_Am_i[7] = -i1_p_i[2];
    i_Am_i[13] = -0.0;
    i_Am_i[19] = -(-i1_p_i[0]);
    i_Am_i[8] = -(-i1_p_i[1]);
    i_Am_i[14] = -i1_p_i[0];
    i_Am_i[20] = -0.0;
    for (i0 = 0; i0 < 6; i0++) {
      for (i1 = 0; i1 < 3; i1++) {
        i_Am_i[i1 + 6 * (i0 + 4)] = b_i_A_i[i1 + 3 * i0] + dv4[i1 + 3 * i0];
      }
    }

    for (i0 = 0; i0 < 3; i0++) {
      i_Am_i[i0 + 3] = i1_p_i[i0];
      for (i1 = 0; i1 < 3; i1++) {
        i_Am_i[(i1 + 6 * (i0 + 1)) + 3] = dv3[i1 + 3 * i0] + b_i1_R_i[i1 + 3 *
          i0];
      }
    }

    for (i0 = 0; i0 < 6; i0++) {
      for (i1 = 0; i1 < 3; i1++) {
        i_Am_i[(i1 + 6 * (i0 + 4)) + 3] = 0.0;
      }
    }

    i0 = 12 * i + 1;
    i1 = 12 * (1 + i);
    if (i0 > i1) {
      i0 = 1;
      i1 = 0;
    }

    loop_ub = (signed char)i1 - (signed char)i0;
    for (j = 0; j <= i; j++) {
      i2 = r0->size[0];
      r0->size[0] = ((signed char)i1 - (signed char)i0) + 1;
      emxEnsureCapacity_int8_T(r0, i2);
      for (i2 = 0; i2 <= loop_ub; i2++) {
        r0->data[i2] = (signed char)((signed char)((signed char)i0 + (signed
          char)i2) - 1);
      }

      x = theta_dot->data[i];
      if (theta_dot->data[i] < 0.0) {
        x = -1.0;
      } else if (theta_dot->data[i] > 0.0) {
        x = 1.0;
      } else {
        if (theta_dot->data[i] == 0.0) {
          x = 0.0;
        }
      }

      for (i2 = 0; i2 < 6; i2++) {
        for (b_i = 0; b_i < 6; b_i++) {
          b_o_X_i[i2 + 6 * b_i] = 0.0;
          for (i3 = 0; i3 < 6; i3++) {
            b_o_X_i[i2 + 6 * b_i] += i_X_o[b_i + 6 * i3] * o_X_i->data[(i3 +
              o_X_i->size[0] * i2) + o_X_i->size[0] * o_X_i->size[1] * j];
          }
        }

        for (b_i = 0; b_i < 10; b_i++) {
          d_i_X_o[i2 + 6 * b_i] = 0.0;
          for (i3 = 0; i3 < 6; i3++) {
            d_i_X_o[i2 + 6 * b_i] += b_o_X_i[i2 + 6 * i3] * i_Am_i[i3 + 6 * b_i];
          }
        }
      }

      for (i2 = 0; i2 < 10; i2++) {
        c_i_X_o[i2] = d_i_X_o[2 + 6 * i2];
      }

      c_i_X_o[10] = (double)(1 + i == 1 + j) * x;
      c_i_X_o[11] = (double)(1 + i == 1 + j) * theta_dot->data[i];
      b_i = r0->size[0];
      for (i2 = 0; i2 < b_i; i2++) {
        K->data[j + K->size[0] * r0->data[i2]] = c_i_X_o[i2];
      }
    }

    i++;
  }

  emxFree_int8_T(&r0);
  emxFree_real_T(&o_X_i);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void fk_initialize(void)
{
  rt_InitInfAndNaN(8U);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void fk_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for fk.c
 *
 * [EOF]
 */
