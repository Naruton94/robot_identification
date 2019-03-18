/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: dynIdenf_emxAPI.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 24-Oct-2018 21:35:46
 */

/* Include Files */
#include <stdlib.h>
#include "dynIdenf.h"
#include "dynIdenf_emxAPI.h"
#include "dynIdenf_emxutil.h"

/* Function Definitions */

/*
 * Arguments    : int numDimensions
 *                int *size
 * Return Type  : emxArray_real_T *
 */
emxArray_real_T *emxCreateND_real_T(int numDimensions, int *size)
{
  emxArray_real_T *emx;
  int numEl;
  int i;
  emxInit_real_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (double *)calloc((unsigned int)numEl, sizeof(double));
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  return emx;
}

/*
 * Arguments    : double *data
 *                int numDimensions
 *                int *size
 * Return Type  : emxArray_real_T *
 */
emxArray_real_T *emxCreateWrapperND_real_T(double *data, int numDimensions, int *
  size)
{
  emxArray_real_T *emx;
  int numEl;
  int i;
  emxInit_real_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  emx->canFreeData = false;
  return emx;
}

/*
 * Arguments    : double *data
 *                int rows
 *                int cols
 * Return Type  : emxArray_real_T *
 */
emxArray_real_T *emxCreateWrapper_real_T(double *data, int rows, int cols)
{
  emxArray_real_T *emx;
  int size[2];
  int numEl;
  int i;
  size[0] = rows;
  size[1] = cols;
  emxInit_real_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  emx->canFreeData = false;
  return emx;
}

/*
 * Arguments    : int rows
 *                int cols
 * Return Type  : emxArray_real_T *
 */
emxArray_real_T *emxCreate_real_T(int rows, int cols)
{
  emxArray_real_T *emx;
  int size[2];
  int numEl;
  int i;
  size[0] = rows;
  size[1] = cols;
  emxInit_real_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (double *)calloc((unsigned int)numEl, sizeof(double));
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  return emx;
}

/*
 * Arguments    : emxArray_real_T *emxArray
 * Return Type  : void
 */
void emxDestroyArray_real_T(emxArray_real_T *emxArray)
{
  emxFree_real_T(&emxArray);
}

/*
 * Arguments    : cell_0 emxArray
 * Return Type  : void
 */
void emxDestroy_cell_0(cell_0 emxArray)
{
  emxFreeStruct_cell_0(&emxArray);
}

/*
 * Arguments    : cell_1 emxArray
 * Return Type  : void
 */
void emxDestroy_cell_1(cell_1 emxArray)
{
  emxFreeStruct_cell_1(&emxArray);
}

/*
 * Arguments    : emxArray_real_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
void emxInitArray_real_T(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxInit_real_T(pEmxArray, numDimensions);
}

/*
 * Arguments    : cell_0 *pStruct
 * Return Type  : void
 */
void emxInit_cell_0(cell_0 *pStruct)
{
  emxInitStruct_cell_0(pStruct);
}

/*
 * Arguments    : cell_1 *pStruct
 * Return Type  : void
 */
void emxInit_cell_1(cell_1 *pStruct)
{
  emxInitStruct_cell_1(pStruct);
}

/*
 * File trailer for dynIdenf_emxAPI.c
 *
 * [EOF]
 */
