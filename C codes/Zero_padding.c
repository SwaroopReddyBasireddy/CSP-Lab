#include <stdio.h>
#include <complex.h>
#include <math.h>

#define PI 3.14159265

main(){
int h[5],i,j,n=3,N=64,h1[N];
for(i=0;i<5;i++)
    h[i]=i+1;
for(i=0;i<N;i++)
{
   if(i<(n+5))
     {
    if(i<n)
        h1[i]=0;
    else
        h1[i]=h[i-n];
     }
   else
     h1[i]=0;
  printf("%d\t",h1[i]);
}

}
