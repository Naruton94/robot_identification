/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: dynIdenf_emxAPI.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 24-Oct-2018 21:35:46
 */

#ifndef DYNIDENF_EMXAPI_H
#define DYNIDENF_EMXAPI_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "dynIdenf_types.h"

/* Function Declarations */
extern emxArray_real_T *emxCreateND_real_T(int numDimensions, int *size);
extern emxArray_real_T *emxCreateWrapperND_real_T(double *data, int
  numDimensions, int *size);
extern emxArray_real_T *emxCreateWrapper_real_T(double *data, int rows, int cols);
extern emxArray_real_T *emxCreate_real_T(int rows, int cols);
extern void emxDestroyArray_real_T(emxArray_real_T *emxArray);
extern void emxDestroy_cell_0(cell_0 emxArray);
extern void emxDestroy_cell_1(cell_1 emxArray);
extern void emxInitArray_real_T(emxArray_real_T **pEmxArray, int numDimensions);
extern void emxInit_cell_0(cell_0 *pStruct);
extern void emxInit_cell_1(cell_1 *pStruct);

#endif

/*
 * File trailer for dynIdenf_emxAPI.h
 *
 * [EOF]
 */
