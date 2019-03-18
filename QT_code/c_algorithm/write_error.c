#include <stdio.h>
#include "dynIdenf_emxAPI.h"
#include "dynIdenf_emxutil.h"
#include "dynIdenf.h"
int write_error(cell_1 *errs, int n, int nSeg)
{
    FILE *ferrs = fopen("results/error.txt", "w");
    if(ferrs==NULL)
    {
      printf("NO file");
      return -1;
    }

    for(int i=0;i<nSeg;i++){
      for(int j=0;j<n;j++)
        fprintf(ferrs ,"%lf\n",errs->f1->data[j*nSeg+i]);
      }
    fclose(ferrs);
    return 0;
}
