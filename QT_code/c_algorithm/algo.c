#include <stdio.h>
#include <math.h>

#include "dynIdenf_emxAPI.h"
#include "dynIdenf_emxutil.h"
#include "dynIdenf.h"
#define p 15001
#define n 6
#define m 3
#define nparJoint ((m*m+3*m+6)/2)  // 12
#define nparMinSet (n*nparJoint - 3*n - 6)    // 48
#define nSeg 2
#define M_PI 3.14159265358979323846

double theta_data[p*n];
double theta_dot_data[p*n];
double theta_ddot_data[p*n];
double tau_data[p*n];

extern int loadfile(double *theta_data, double *theta_dot_data, double *theta_ddot_data, double *tau_data, int n1, int p1);
extern int write_results(emxArray_real_T *phi, emxArray_real_T *tau_pos, emxArray_real_T *tau_pre,
                         emxArray_real_T *tau_filt, int n1, int p1, int nparMinSet1);
extern int write_error(cell_1 *errs1, int n1, int nSeg1);

int algo(){

  //set input parameters

  double a_data[n+1] = {0, 0, 0.416, 0.4208, 0, 0, 0};
  double alpha_data[n+1] = {0, M_PI/2, 0, 0, -M_PI/2, M_PI/2, 0};
  double d_data[n+1] = {0, 0.1181, 0, 0, 0.1301, 0.1021, 0.0568};

  double g[m] = {0, 0, -9.80665};
  double phi_r0_data[n*nparJoint] = {2,-0.0165019647020000,-0.0256289655860000,-0.0456037864700000,0.00388830591570456,
  7.63638706167896e-05,0.000340020397992622,0.00333118131630822,0.000215211872926742,0.00262696792835856,0,0,
  3.42000000000000,0.426439873502640,-0.0124271511264000,0.386835887459880,0.0475012505609916,3.00901079796083e-06,
  -0.0528168520144169,0.131182004246096,0.00172637529085978,0.0864328868443469,0,0, 1.26000000000000,0.139078355556420,
  -1.07045177400000e-05,0.0411503703017400,0.00223362215511204,5.03356734903604e-07,-0.00375135880026783,0.0264742256951927,
  3.69281505903818e-07,0.0248796360944246,0,0, 0.800000000000000,1.39460800000000e-07,0.000225010768000000,-0.00478865628240000,
  0.000772782642110361,5.04623525623476e-09,8.09638107508407e-09,0.000537499391466433,-6.01459794336313e-06,0.000649992384815521,0,0,
   0.800000000000000,-2.63220000000000e-06,-0.000320608410400000,-0.00460317585600000,0.000770846701360872,-6.16295700192826e-08,
   -4.41585931650860e-08,0.000531210606854337,6.03332741029044e-06,0.000653231083089124,0,0, 0.350000000000000,8.29500000000000e-11,
   -3.71160097000000e-05,-0.00683040118320000,0.000261196607049744,-4.06640243700716e-12,1.37724047043447e-11,0.000261849371147053,
   -4.02812899656085e-07,0.000175309188669306,0,0};
  double pfilt = 5001;
  double pidenf[2] = {51,100};
  double peval[2] = {5001,1};
  double noise_err = 1e-6;
  double cond_max = 100.0;
  double lambda = 0.0;
  double fpass = 2.0;
  double segErr_data[2] = {0.2,1.0};
  emxArray_real_T *a = emxCreateWrapper_real_T(a_data,1,n+1);
  emxArray_real_T *alpha = emxCreateWrapper_real_T(alpha_data,1,n+1);
  emxArray_real_T *d = emxCreateWrapper_real_T(d_data,1,n+1);
  emxArray_real_T *phi_r0 = emxCreateWrapper_real_T(phi_r0_data,nparJoint,n);
  emxArray_real_T *segErr = emxCreateWrapper_real_T(segErr_data,1,nSeg);
  cell_0 *pars = (cell_0 *)malloc(sizeof(cell_0));
  emxInit_cell_0(pars);
  pars->f1 = a;
  pars->f2 = alpha;
  pars->f3 = d;
  for(int i=0;i<m;i++){
    pars->f4[i] = g[i];
  }
  pars->f5 = phi_r0;
  pars->f6 = pfilt;
  for(int i=0;i<2;i++){
    pars->f7[i] = pidenf[i];
    pars->f8[i] = peval[i];
  }
  pars->f9 = noise_err;
  pars->f10 = cond_max;
  pars->f11 = lambda;
  pars->f12 = fpass;
  pars->f13 = segErr;

  // set the output parameters
  cell_1 *errs = (cell_1 *)malloc(sizeof(cell_1));
  emxInit_cell_1(errs);
  errs->f1 = emxCreate_real_T(nSeg,n);
  errs->f2 = emxCreate_real_T(1,n);
  errs->f3 = emxCreate_real_T(1,n);
  errs->f4 = 0;

  // set the input variables

  emxArray_real_T *theta = emxCreateWrapper_real_T(theta_data,p,n);
  emxArray_real_T *theta_dot = emxCreateWrapper_real_T(theta_dot_data,p,n);
  emxArray_real_T *theta_ddot = emxCreateWrapper_real_T(theta_ddot_data,p,n);
  emxArray_real_T *tau = emxCreateWrapper_real_T(tau_data,p,n);

  // set the output variables
  emxArray_real_T *phi = emxCreate_real_T(nparMinSet,1);
  emxArray_real_T *tau_pos = emxCreate_real_T((int)peval[0],n);
  emxArray_real_T *tau_pre = emxCreate_real_T((int)peval[0],n);
  emxArray_real_T *tau_filt = emxCreate_real_T((int)peval[0],n);

  loadfile(theta_data, theta_dot_data, theta_ddot_data, tau_data, n, p);


  //printf("%d\t%d\n",errs->f3->size[0],errs->f3->size[1]);
  // calculate
  dynIdenf_initialize();
  dynIdenf(theta, theta_dot, theta_ddot, tau, pars, phi, tau_pos, tau_pre, tau_filt, errs);
  dynIdenf_terminate();
  write_results(phi, tau_pos, tau_pre, tau_filt , n, peval[0], nparMinSet);
    /*
  FILE *ftau = fopen("results/tau_data.txt", "w");
  for(int i=0; i<n*p;++i)
    fprintf(ftau,"%f\n", tau_data[i]);
  fclose(ftau);

  // print results
  for(int i=0;i<(int)peval[0];i++){
    for(int j=0;j<n;j++)
      printf("%lf\t",theta->data[j*(int)peval[0]+i]);
    printf("\n");
  }
  //for(int i=0;i<nparMinSet;i++)
    //printf("%f\n", phi->data[i]);
        */

  write_error( errs,  n,  nSeg);
  /*
  for(int i=0;i<nSeg;i++){
    for(int j=0;j<n;j++)
      printf("%lf\t",errs->f1->data[j*nSeg+i]);
    printf("\n");
    }
    */
    return 0;

}


