#include <stdio.h>
void ctwrite_(lvisit,visit,k)
float *lvisit;
int *visit,*k;
{
   FILE *ff;
   int i;
   ff=fopen("triolrlisting.tmp","a");
   fprintf(ff,"%f %f ",lvisit[0],lvisit[1]);
   for(i=0;i<(*k);i++)fprintf(ff,"%d ",visit[i]);
   fprintf(ff,"\n");
   fclose(ff);
}
void ctiwrite_()
{
   FILE *ff;
   ff=fopen("triolrlisting.tmp","w");
   fclose(ff);
}
