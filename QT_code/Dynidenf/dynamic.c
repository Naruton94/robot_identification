#include <stdio.h>
#include <math.h>
#include "fk.h"
#include "idf_SVD_priori.h"
#include "Errestim.h"
#include "matlab_emxAPI.h"
#include "matlab_emxutil.h"
#define PI 3.141592653589793
#define s 100
#define m 3
#define n 6

int dynamic(void){
  
	FILE *fptheta_data = fopen("theta_data.txt", "w");
	FILE *fptheta_dot_data = fopen("theta_dot_data.txt", "w");
	FILE *fptheta_ddot_data = fopen("theta_ddot_data.txt", "w");
	FILE *fptau_data = fopen("tau_data.txt", "w");
	FILE *fpresult = fopen("result.txt", "w");
    FILE *fptau_pos = fopen("tau_pos.txt", "w");

	if (fptheta_data == NULL || fptheta_dot_data == NULL ||
        fptheta_ddot_data == NULL || fptau_data == NULL || fpresult == NULL)
    {
		return 0;
	}


  double theta_data[s*n];   
  double theta_dot_data[s*n];  
  double theta_ddot_data[s*n];  
  double tau_data[s*n];  
  double phi_pre_data[(3*m+3)*n];
  double phi_r_data[(3*m+3)*n];
  double a_alpha_d_data[3*(n+1)];
  double g[m]; g[0] = 0; g[1] = 0; g[2] = -9.81;
  //
  double opts[4] = {1e-6, 100, 1, 1};

  double thetas_data[n];
  double thetas_dot_data[n];
  double thetas_ddot_data[n];
  double K_data[s*n*(3*m+3)*n];
  double Ks_data[n*(3*m+3)*n];
  double phi_pos_data[(3*m+3)*n];
  double rank_cond[2];
  double err_rms_rel[2];

  double tau_pos_data[s*n];

  // ������������
  for(int i=0; i<s*n; i++){
    theta_data[i] = (double)(rand()%1000)/1000;
    theta_dot_data[i] = (double)(rand()%1000)/1000;
    theta_ddot_data[i] = (double)(rand()%1000)/1000;
    tau_data[i] = (double)(rand()%1000)/100;


    fprintf(fptheta_data, "%.4f\n", theta_data[i]);
	// if ((i+1) % 6 == 0)
	// 	fprintf(fptheta_data, "\n");

    fprintf(fptheta_dot_data, "%.4f\n", theta_dot_data[i]);
	// if ((i + 1) % 6 == 0)
	//	fprintf(fptheta_dot_data, "\n");

    fprintf(fptheta_ddot_data, "%.4f\n", theta_ddot_data[i]);
	// if ((i + 1) % 6 == 0)
	// 	fprintf(fptheta_ddot_data, "\n");

    fprintf(fptau_data, "%.4f\n", tau_data[i]);
	// if ((i + 1) % 6 == 0)
	//	fprintf(fptau_data, "\n");
  }
  fclose(fptheta_data);
  fclose(fptheta_dot_data);
  fclose(fptheta_ddot_data);
  fclose(fptau_data);

  for(int i=0; i<(3*m+3)*n; i++){
    phi_pre_data[i] = (double)(rand()%1000)/1000;
    phi_r_data[i] = (double)(rand()%1000)/1000;
  }
 // D-H
  double a[] = {0, 0, 0.24365, 0.213, 0, 0, 0};
  double alpha[] = {0, PI/2, 0, 0, -PI/2, PI/2, 0};
  double d[] = {0, 0.1519, 0, 0, 0.0834, 0.0834, 0.0819};

  for(int i=0;i<n+1;i++){
    a_alpha_d_data[i*3] = a[i];
    a_alpha_d_data[i*3+1] = alpha[i];
    a_alpha_d_data[i*3+2] = d[i];
  }

  emxArray_real_T *tau, *phi_pre, *phi_r, *a_alpha_d;
  emxArray_real_T *thetas, *thetas_dot, *thetas_ddot, *Ks, *K, *phi_pos, *tau_pos;
  phi_pre = emxCreateWrapper_real_T(phi_pre_data,(3*m+3)*n,1);
  phi_r = emxCreateWrapper_real_T(phi_r_data,(3*m+3)*n,1);
  // trans
  tau = emxCreateWrapper_real_T(tau_data,n*s,1);
  Ks = emxCreateWrapper_real_T(Ks_data,n,(3*m+3)*n);
  phi_pos = emxCreateWrapper_real_T(phi_pos_data,(3*m+3)*n,1);
  tau_pos = emxCreateWrapper_real_T(tau_pos_data, s*n, 1);
  a_alpha_d = emxCreateWrapper_real_T(a_alpha_d_data,3,n);

  for(int i=0; i<s; i++){
    for(int j=0; j<n; j++){
      thetas_data[j] = theta_data[i*n+j];
      thetas_dot_data[j] = theta_dot_data[i*n+j];
      thetas_ddot_data[j] = theta_ddot_data[i*n+j];
    }
    thetas = emxCreateWrapper_real_T(thetas_data,1,n);
    thetas_dot = emxCreateWrapper_real_T(thetas_dot_data,1,n);
    thetas_ddot = emxCreateWrapper_real_T(thetas_ddot_data,1,n);
    fk(thetas, thetas_dot, thetas_ddot, a_alpha_d, g, Ks);
    for(int j=0; j<n; j++)
      for(int k=0; k<(3*m+3)*n; k++)
        K_data[k*n*s+i*n+j] = Ks->data[k*n+j];
  }
  K = emxCreateWrapper_real_T(K_data,s*n,(3*m+3)*n);

  idf_SVD_priori(phi_pre, phi_r, tau, K, opts, phi_pos, rank_cond);
  Errestim(tau, K, phi_pos, tau_pos ,err_rms_rel);

#if 0
  printf("rank = %f, ",rank_cond[0]);
  printf("cond = %f\n",rank_cond[1]);
  printf("rms = %f, ",err_rms_rel[0]);
  printf("rel = %f\n",err_rms_rel[1]);
#endif

  fprintf(fpresult, "%.4f\n", rank_cond[0]);
  fprintf(fpresult, "%.4f\n", rank_cond[1]);
  fprintf(fpresult, "%.4f\n", err_rms_rel[0]);
  fprintf(fpresult, "%.4f\n", err_rms_rel[1]);


  fclose(fpresult);
  for(int i=0; i<s*n; i++)
  {
      fprintf(fptau_pos, "%f\n", tau_pos->data[i]);
  }
  fclose(fptau_pos);

  return 0;
}
