#include <stdio.h>
#include <complex.h>
#include <math.h>

void bit_rev(int N,float complex *y)
{
    float complex c;
    int i;
    while(i<N/2)//for(i=1;i<N/2;i+2)
    {
        c=*(&y[i]);
        y[i]=y[i+N/2-1];
        y[i+N/2-1]=c;
        i=i+2;
    }
}

main(){

int i,j,N=8;
float complex z[N],c=0.0+0.0*I;
for(i=0;i<N;i++)
{
  z[i]=(i+1)+(i+1)*I;
  printf("%.2f %+.2fI \t",crealf(z[i]),cimagf(z[i]));
}
//bit_reversal
for(j=0;j<N;j++)
    { if(j%2!=0)
       {
        c=z[j];
        z[j]=z[j+(int)(N/2)];
        z[j+(int)(N/2)]=c;
       }

       }

printf("\nAfter bit_reversal\n");
for(i=0;i<N;i++)
  printf("%.2f %+.2fI\t",crealf(z[i]),cimagf(z[i]));

  return(0);
}
