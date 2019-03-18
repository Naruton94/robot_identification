/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: dynIdenf_emxutil.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 24-Oct-2018 21:35:46
 */

#ifndef DYNIDENF_EMXUTIL_H
#define DYNIDENF_EMXUTIL_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "dynIdenf_types.h"

/* Function Declarations */
extern void c_emxFreeStruct_coder_internal_(coder_internal_sparse *pStruct);
extern void c_emxInitStruct_coder_internal_(coder_internal_sparse *pStruct);
extern void emxCopyStruct_cell_wrap_4(cell_wrap_4 *dst, const cell_wrap_4 *src);
extern void emxEnsureCapacity_boolean_T(emxArray_boolean_T *emxArray, int
  oldNumel);
extern void emxEnsureCapacity_boolean_T1(emxArray_boolean_T *emxArray, int
  oldNumel);
extern void emxEnsureCapacity_int32_T(emxArray_int32_T *emxArray, int oldNumel);
extern void emxEnsureCapacity_int32_T1(emxArray_int32_T *emxArray, int oldNumel);
extern void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int oldNumel);
extern void emxEnsureCapacity_real_T1(emxArray_real_T *emxArray, int oldNumel);
extern void emxEnsureCapacity_real_T2(emxArray_real_T *emxArray, int oldNumel);
extern void emxEnsureCapacity_struct_T(emxArray_struct_T *emxArray, int oldNumel);
extern void emxFreeMatrix_cell_wrap_4(cell_wrap_4 pMatrix[2]);
extern void emxFreeStruct_cell_0(cell_0 *pStruct);
extern void emxFreeStruct_cell_1(cell_1 *pStruct);
extern void emxFree_boolean_T(emxArray_boolean_T **pEmxArray);
extern void emxFree_int32_T(emxArray_int32_T **pEmxArray);
extern void emxFree_real_T(emxArray_real_T **pEmxArray);
extern void emxFree_struct_T(emxArray_struct_T **pEmxArray);
extern void emxInitMatrix_cell_wrap_4(cell_wrap_4 pMatrix[2]);
extern void emxInitStruct_cell_0(cell_0 *pStruct);
extern void emxInitStruct_cell_1(cell_1 *pStruct);
extern void emxInit_boolean_T(emxArray_boolean_T **pEmxArray, int numDimensions);
extern void emxInit_boolean_T1(emxArray_boolean_T **pEmxArray, int numDimensions);
extern void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions);
extern void emxInit_int32_T1(emxArray_int32_T **pEmxArray, int numDimensions);
extern void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);
extern void emxInit_real_T1(emxArray_real_T **pEmxArray, int numDimensions);
extern void emxInit_real_T2(emxArray_real_T **pEmxArray, int numDimensions);
extern void emxInit_struct_T(emxArray_struct_T **pEmxArray, int numDimensions);

#endif

/*
 * File trailer for dynIdenf_emxutil.h
 *
 * [EOF]
 */
