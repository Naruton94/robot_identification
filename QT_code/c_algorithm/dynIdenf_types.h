/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: dynIdenf_types.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 24-Oct-2018 21:35:46
 */

#ifndef DYNIDENF_TYPES_H
#define DYNIDENF_TYPES_H

/* Include Files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  double *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  bool canFreeData;
};

#endif                                 /*struct_emxArray_real_T*/

#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T

typedef struct emxArray_real_T emxArray_real_T;

#endif                                 /*typedef_emxArray_real_T*/

#ifndef typedef_cell_0
#define typedef_cell_0

typedef struct {
  emxArray_real_T *f1;
  emxArray_real_T *f2;
  emxArray_real_T *f3;
  double f4[3];
  emxArray_real_T *f5;
  double f6;
  double f7[2];
  double f8[2];
  double f9;
  double f10;
  double f11;
  double f12;
  emxArray_real_T *f13;
} cell_0;

#endif                                 /*typedef_cell_0*/

#ifndef typedef_cell_1
#define typedef_cell_1

typedef struct {
  emxArray_real_T *f1;
  emxArray_real_T *f2;
  emxArray_real_T *f3;
  double f4;
} cell_1;

#endif                                 /*typedef_cell_1*/

#ifndef struct_emxArray_int32_T
#define struct_emxArray_int32_T

struct emxArray_int32_T
{
  int *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  bool canFreeData;
};

#endif                                 /*struct_emxArray_int32_T*/

#ifndef typedef_emxArray_int32_T
#define typedef_emxArray_int32_T

typedef struct emxArray_int32_T emxArray_int32_T;

#endif                                 /*typedef_emxArray_int32_T*/

#ifndef typedef_cell_wrap_4
#define typedef_cell_wrap_4

typedef struct {
  emxArray_int32_T *f1;
} cell_wrap_4;

#endif                                 /*typedef_cell_wrap_4*/

#ifndef typedef_coder_internal_sparse
#define typedef_coder_internal_sparse

typedef struct {
  emxArray_real_T *d;
  emxArray_int32_T *colidx;
  emxArray_int32_T *rowidx;
  int m;
  bool __m_AssignmentSentinel;
  int n;
  bool __n_AssignmentSentinel;
} coder_internal_sparse;

#endif                                 /*typedef_coder_internal_sparse*/

#ifndef struct_emxArray_boolean_T
#define struct_emxArray_boolean_T

struct emxArray_boolean_T
{
  bool *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  bool canFreeData;
};

#endif                                 /*struct_emxArray_boolean_T*/

#ifndef typedef_emxArray_boolean_T
#define typedef_emxArray_boolean_T

typedef struct emxArray_boolean_T emxArray_boolean_T;

#endif                                 /*typedef_emxArray_boolean_T*/

#ifndef typedef_struct_T
#define typedef_struct_T

typedef struct {
  int xstart;
  int xend;
  int depth;
} struct_T;

#endif                                 /*typedef_struct_T*/

#ifndef typedef_emxArray_struct_T
#define typedef_emxArray_struct_T

typedef struct {
  struct_T *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  bool canFreeData;
} emxArray_struct_T;

#endif                                 /*typedef_emxArray_struct_T*/
#endif

/*
 * File trailer for dynIdenf_types.h
 *
 * [EOF]
 */
