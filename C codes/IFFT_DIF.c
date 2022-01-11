#include <stdio.h>
#include <complex.h>
#include <math.h>

#define PI 3.14159265
void twiddle_factor(float complex *p,int n)
{
int i;
float c;
c=2*PI/n;
for(i=0;i<n;i++)
  p[i]=cos(c*i) + (sin(c*i))*I;
}

void ifft_stage(int N,int stage,float complex *x,float complex *y)
{
  int k=N/pow(2,stage-1);
  float complex w[k];
  int i,j,m;
  //printf("stage:%d\n",stage);
  twiddle_factor(w,k);
  //for(i=0;i<k;i++)
    //printf("%.3f %+.3fI \n",crealf(w[i]),cimagf(w[i]));

for(j=0;j<(int)(pow(2,stage-1));j++)
for(i=0;i<(int)(N/pow(2,stage));i++)
    {
        m=2*j*N/pow(2,stage);
            y[i+m]=x[i+m]+x[i+m+(int)(N/pow(2,stage))];
        if((int)pow(2,stage)!=N)
            y[i+m+(int)(N/pow(2,stage))]=(x[i+m]-x[i+m+(int)(N/pow(2,stage))])*w[i];
        else
            y[i+m+(int)(N/pow(2,stage))]=(x[i+m]-x[i+m+(int)(N/pow(2,stage))]);
    }
 //bit_rev(y,N);
}

void bit_rev(float complex *y,int N)
{
    float complex z[N];
    int i;
    for(i=0;i<N;i++)
    {

        if(i%2!=0 && i<N/2)
        {
        z[i]=y[i];
        y[i]=(y[i-1+N/2]);
        y[i-1+N/2]=z[i];
        }
        else
          z[i]=y[i];

    }
    for(i=0;i<N;i++)
       y[i]/=(float)N;

}

void main(){

int N=64,i,j,k,no_stages;
float complex Y[N],t[N],x[N];
no_stages=round(log10(N)/log10(2));
//printf("No of stages:%d",no_stages);
//x[N]=(1,2,2,1,2,1,1,2);
for(i=0;i<N;i++)
{
    x[i]=(i+1.0)+ 0.0*I;
    Y[i]=0.0 + 0.0*I;
}
for(j=0;j<no_stages;j++)
{
   if(j==0)
       ifft_stage(N,j+1,x,Y);
   else
     {
       for(i=0;i<N;i++)
         t[i]=Y[i];
       ifft_stage(N,j+1,t,Y);
     }
 printf("\n stage:%d\n",j+1);
 for(k=0;k<N;k++)
   printf("%.2f %+.2fI \t",crealf(Y[k]),cimagf(Y[k]));
 }
bit_rev(Y,N);
 printf("\nIFFT of x:");
 for(k=0;k<N;k++)
   printf("%.2f %+.2fI \t",crealf(Y[k]),cimagf(Y[k]));

return(0);

}
