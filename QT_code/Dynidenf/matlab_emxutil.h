/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: matlab_emxutil.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 24-Jul-2018 20:21:10
 */

#ifndef MATLAB_EMXUTIL_H
#define MATLAB_EMXUTIL_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "matlab_types.h"

/* Function Declarations */
extern void emxEnsureCapacity_boolean_T(emxArray_boolean_T *emxArray, int
  oldNumel);
extern void emxEnsureCapacity_boolean_T1(emxArray_boolean_T *emxArray, int
  oldNumel);
extern void emxEnsureCapacity_int32_T(emxArray_int32_T *emxArray, int oldNumel);
extern void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int oldNumel);
extern void emxEnsureCapacity_real_T1(emxArray_real_T *emxArray, int oldNumel);
extern void emxFree_boolean_T(emxArray_boolean_T **pEmxArray);
extern void emxFree_int32_T(emxArray_int32_T **pEmxArray);
extern void emxFree_real_T(emxArray_real_T **pEmxArray);
extern void emxInit_boolean_T(emxArray_boolean_T **pEmxArray, int numDimensions);
extern void emxInit_boolean_T1(emxArray_boolean_T **pEmxArray, int numDimensions);
extern void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions);
extern void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);
extern void emxInit_real_T1(emxArray_real_T **pEmxArray, int numDimensions);

extern void emxEnsureCapacity_int8_T(emxArray_int8_T *emxArray, int oldNumel);
extern void emxFree_int8_T(emxArray_int8_T **pEmxArray);
extern void emxInit_int8_T(emxArray_int8_T **pEmxArray, int numDimensions);

#endif

/*
 * File trailer for matlab_emxutil.h
 *
 * [EOF]
 */
