#include <stdio.h>
#include <complex.h>
#include <math.h>

void convolution(int n,int m,float complex *h, float complex *x,float complex *y,int N)
  {
    int i,j,k;
    float complex c[n][m],z[n+m-1];
    for(i=0;i<m+n-1;i++)
        z[i]=0.0+0.0*I;
    for(i=0;i<n;i++)
        for(j=0;j<m;j++)
           c[i][j]=h[i]*x[j];

  for(k=0;k<(n+m-1);k++)
    for(i=0;i<n;i++)
       for(j=0;j<m;j++)
         if((i+j)==k)
           z[k]+=c[i][j];
for(i=0;i<n+m-1;i++)
    printf("%.2f %+.2fI \n",crealf(z[i]),cimagf(z[i]));
  if(n+m-1>N)
     for(i=0;i<n+m-1;i++)
       if(i<N)
         y[i]=z[i];

       else
          y[i-N]+=z[i];
  else
      for(i=0;i<N;i++)
          if(i<n+m-1)
            y[i]=z[i];
          else
            z[i]=0;

  }

  main(){

    int n=3,m=3,i,N=3;
    float complex h[n],x[m],y[N];
    for(i=0;i<m+n;i++)
        y[i]=0.0+0.0*I;
    for(i=0;i<n;i++)
        h[i]=(i+1)+(i+1)*I;
    for(i=0;i<m;i++)
        x[i]=(i+2)+(i+2)*I;

    convolution(n,m,h,x,y,N);

    for(i=0;i<N;i++)
        printf("\n%.2f %+.2fI \n",crealf(y[i]),cimagf(y[i]));


  }
