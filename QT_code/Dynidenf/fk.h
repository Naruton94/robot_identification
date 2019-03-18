/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: fk.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 24-Jul-2018 17:05:24
 */

#ifndef FK_H
#define FK_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "matlab_types.h"

/* Function Declarations */
extern void fk(const emxArray_real_T *theta, const emxArray_real_T *theta_dot,
               const emxArray_real_T *theta_ddot, const emxArray_real_T
               *a_alpha_d, const double g[3], emxArray_real_T *K);
extern void fk_initialize(void);
extern void fk_terminate(void);

#endif

/*
 * File trailer for fk.h
 *
 * [EOF]
 */
