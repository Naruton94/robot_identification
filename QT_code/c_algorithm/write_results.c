#include <stdio.h>
#include "dynIdenf_emxAPI.h"
#include "dynIdenf_emxutil.h"
#include "dynIdenf.h"
int write_results(emxArray_real_T *phi, emxArray_real_T *tau_pos, emxArray_real_T *tau_pre, emxArray_real_T *tau_filt,
                   int n, int p, int nparMinSet)
{
  FILE *fphi = fopen("results/phi.txt", "w");
  FILE *ftau_pos = fopen("results/tau_pos_data.txt", "w");
  FILE *ftau_pre = fopen("results/tau_pre_data.txt", "w");
  FILE *ftau_filt = fopen("results/tau_filt_data.txt", "w");
  if(fphi==NULL || ftau_pos==NULL || ftau_pre==NULL || ftau_filt==NULL)
  {
    printf("NO file");
    return 0;
  }
  for(int i=0; i<nparMinSet;++i)
    fprintf(fphi,"%.10lf\n", phi->data[i]);
  fclose(fphi);
  for(int i=0; i<p*n;++i)
    fprintf(ftau_pos,"%lf\n", tau_pos->data[i]);
  fclose(ftau_pos);
  for(int i=0; i<p*n;++i)
    fprintf(ftau_pos,"%.10lf\n", tau_pos->data[i]);
  fclose(ftau_pos);
  for(int i=0; i<p*n;++i)
    fprintf(ftau_filt,"%.10lf\n", tau_filt->data[i]);
  fclose(ftau_filt);

  for(int i=0; i<p*n;++i)
    fprintf(ftau_pre,"%.10lf\n", tau_pre->data[i]);
  fclose(ftau_pre);


}
