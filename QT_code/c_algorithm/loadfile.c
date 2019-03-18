#include<stdio.h>
#include <string.h>
int loadfile(double *theta_data, double *theta_dot_data, double *theta_ddot_data, double *tau_data,int n, int p)
{

  char *path = "data\\DataOutput_";
  char *suff = ".txt";
  char *varname[] = {"theta", "velocity", "acceleration", "torque"};
  char *joint[] = {"1","2","3","4","5","6"};
  char fullname[200];
  int nums;
  double data;
  char strs[20];
  int idx;

  for(int i=0;i<4;i++){
    for(int j=0;j<n;j++){
      memset(fullname,0,200*sizeof(char));
      strcat(fullname,path);
      strcat(fullname,varname[i]);
      strcat(fullname,joint[j]);
      strcat(fullname,suff);
      FILE *fileID = fopen(fullname,"r");
      if (fileID == NULL){
        printf("No file find!");
        return 0;
      }
      for(int k=0;k<p;k++){
        fscanf(fileID,"DataOutput.OutData[%d]:=%lf%s\n",&nums,&data,strs);
        idx = j*p+k;
        if(i==0)
          theta_data[idx] = data;
        else if(i==1)
          theta_dot_data[idx] = data;
        else if(i==2)
          theta_ddot_data[idx] = data;
        else
          tau_data[idx] = data;
      }
      fclose(fileID);
      fileID = NULL;
      memset(fullname,0,200*sizeof(char));
    }

  }
  return 1;
}
