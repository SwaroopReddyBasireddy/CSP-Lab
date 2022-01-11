#include <stdio.h>
#include <complex.h>
#include <math.h>

#define PI 3.14159265


void main(){

int N=64,i,j,k,no_stages;

float complex Y[N],t[N],x[N];
int16_t Y_fp[N],t_fp[N],x_fp[N];
no_stages=round(log10(N)/log10(2));

//printf("No of stages:%d",no_stages);
for(i=0;i<N;i++)
{
    x[i]=(i+1)+ 0.0*I;
    x_fp[i]=(x[i])*(int16_t)(pow(2,11));
    Y[i]=0.0 + 0.0*I;
    Y_fp[i]=(Y[i])*(int16_t)(pow(2,11));
}
for(j=0;j<no_stages;j++)
{
   if(j==0)
       fft_stage(N,j+1,x,Y);
   else
     {
       for(i=0;i<N;i++)
         t[i]=Y[i];
       fft_stage(N,j+1,t,Y);
     }
 printf("\n stage:%d\n",j+1);
 for(k=0;k<N;k++)
   printf("%.2f %+.2fI \t",crealf(Y[k]),cimagf(Y[k]));
 }
 bit_rev(Y,N);
 printf("\nFFT of x:");
 for(k=0;k<N;k++)
   printf("%.2f %+.2fI \t",crealf(Y[k]),cimagf(Y[k]));

return(0);

}
