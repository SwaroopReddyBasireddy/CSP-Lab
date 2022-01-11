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
  p[i]=cos(c*i) - (sin(c*i))*I;
}

void fft_stage(int N,int stage,float complex *x,float complex *y)
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

void bit_rev(int N,float complex *y)
{
  int i,j,n=round(log10(N)/log10(2)),bin[N][n],dec[N],t=0,p=0;
  float complex x[N];
  for(i=0;i<N;i++){
    x[i]=y[i];
    dec[i]=0;
    //y[i]=0;
    //printf("%d\t",x[i]);
  }


  for(j=0;j<N;j++){
        p=j;
    for(i=0;i<n;i++)
     {

      bin[j][i]= p%2;

        p/=2;
      //printf("%d",bin[j][i]);
        }
      //printf("\t");

      }
  printf("\n");
  for(j=0;j<N;j++){
    for(i=0;i<n;i++)
       dec[j] +=bin[j][i]*round(pow(2,n-i-1));
     //printf("%d\t",dec[j]);

  }
  //printf("\n");
for(i=0;i<N;i++)
{
    j=dec[i];
    y[i]=x[j];
    //printf("%d\t",y[j]);

    }

}
void main(){

int N=16,i,j,k,no_stages;
float complex Y[N],t[N],x[N],y[N];
no_stages=round(log10(N)/log10(2));
//printf("No of stages:%d",no_stages);
for(i=0;i<N;i++)
{
    x[i]=(i+1)+ 0.0*I;
    Y[i]=0.0 + 0.0*I;
    y[i]=0.0 + 0.0*I;
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

 }
 for(k=0;k<N;k++)
   printf("%.2f %+.2fI \t",crealf(Y[k]),cimagf(Y[k]));
 bit_rev(N,Y);
 printf("\n FFT of x:");
 for(k=0;k<N;k++)
   printf("%.2f %+.2fI \t",crealf(Y[k]),cimagf(Y[k]));

return(0);

}
