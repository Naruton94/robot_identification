/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: idf_SVD_priori.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 24-Jul-2018 20:21:10
 */

#ifndef IDF_SVD_PRIORI_H
#define IDF_SVD_PRIORI_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "matlab_types.h"

/* Function Declarations */
extern void idf_SVD_priori(const emxArray_real_T *phi_pre, const emxArray_real_T
  *phi_r, const emxArray_real_T *tau, const emxArray_real_T *K, const double
  opts[4], emxArray_real_T *phi_pos, double rank_cond[2]);
extern void idf_SVD_priori_initialize(void);
extern void idf_SVD_priori_terminate(void);

#endif

/*
 * File trailer for idf_SVD_priori.h
 *
 * [EOF]
 */
