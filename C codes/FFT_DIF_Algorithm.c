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


void fft_stage(int N,int stage,float complex *x,float complex *y)
{
  int k=N/pow(2,stage-1);
  float complex w[k];
  int i,j,m;

  twiddle_factor(w,k);

for(j=0;j<(int)pow(2,stage-1);j++)
{
    for(i=0;i<(int)(N/pow(2,stage));i++)
    {
        m=2*j*N/pow(2,stage);
            y[i+m]=x[i+m]+x[i+m+(int)(N/pow(2,stage))];
        if((int)pow(2,stage)!=N)
            y[i+m+(int)(N/pow(2,stage))]=(x[i+m]-x[i+m+(int)(N/pow(2,stage))])*w[i];
        else
            y[i+m+(int)(N/pow(2,stage))]=(x[i+m]-x[i+m+(int)(N/pow(2,stage))]);
    }
}
}


int main(){

int N=8,i,j,no_stages;
float complex z[3][N],Y[N],t[N],x[N];
no_stages=round(log(N)/log(2));
for(i=0;i<N;i++)
    x[i]=(i+1)+0*I;

for(i=0;i<no_stages;i++)
{

   if(i==0)
       fft_stage(N,i+1,x,Y);
     else
     {
       for(j=0;j<N;j++)
         t[j]=Y[j];
       fft_stage(N,i+1,t,Y);
       for(j=0;j<N;j++)
         z[i-1][j]=t[j];
     }
}
for(i=0;i<N;i++)
    printf("%.2f %+.2fI \n",crealf(z[no_stages-1][i]),cimagf(z[no_stages-1][i]));
return(0);

}
