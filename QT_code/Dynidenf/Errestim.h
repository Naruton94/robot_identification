/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: Errestim.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 31-Aug-2018 13:44:40
 */

#ifndef ERRESTIM_H
#define ERRESTIM_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "matlab_types.h"

/* Function Declarations */
extern void Errestim(const emxArray_real_T *tau, const emxArray_real_T *K, const
                     emxArray_real_T *phi_pre, emxArray_real_T *tau_pre, double
                     err_rms_rel[2]);
extern void Errestim_initialize(void);
extern void Errestim_terminate(void);

#endif

/*
 * File trailer for Errestim.h
 *
 * [EOF]
 */
