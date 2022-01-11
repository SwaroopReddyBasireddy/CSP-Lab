#include <stdio.h>
#include <complex.h>
#include <math.h>

void convolution(int n,int m,float complex *h, float complex *x,float complex *y)
  {
    int i,j,k;
    float complex c[n][m];
    for(i=0;i<n;i++)
        for(j=0;j<m;j++)
           c[i][j]=h[i]*x[j];

  for(k=0;k<(n+m);k++)
    for(i=0;i<n;i++)
       for(j=0;j<m;j++)
         if((i+j)==k)
           y[k]+=c[i][j];
  }

  main(){

    int n=3,m=3,i;
    float complex h[n],x[m],y[n+m-1];
    for(i=0;i<m+n;i++)
        y[i]=0.0+0.0*I;
    for(i=0;i<n;i++)
        h[i]=(i+1)+(i+1)*I;
    for(i=0;i<m;i++)
        x[i]=(i+2)+(i+2)*I;

    convolution(n,m,h,x,y);

    for(i=0;i<m+n-1;i++)
        printf("%.2f %+.2fI \n",crealf(y[i]),cimagf(y[i]));


  }
