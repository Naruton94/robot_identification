/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: dynIdenf.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 24-Oct-2018 21:35:46
 */

#ifndef DYNIDENF_H
#define DYNIDENF_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "dynIdenf_types.h"

/* Function Declarations */
extern void dynIdenf(const emxArray_real_T *theta, const emxArray_real_T
                     *theta_dot, const emxArray_real_T *theta_ddot, const
                     emxArray_real_T *tau, const cell_0 *pars, emxArray_real_T
                     *phi, emxArray_real_T *tau_pos, emxArray_real_T *tau_pre,
                     emxArray_real_T *tau_filt, cell_1 *errs);
extern void dynIdenf_initialize(void);
extern void dynIdenf_terminate(void);

#endif

/*
 * File trailer for dynIdenf.h
 *
 * [EOF]
 */
