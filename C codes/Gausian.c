#include<stdio.h>
#include<math.h>
#include<stdlib.h>

int main(){
int N=1000,n=100,i,j,x;
int noise[N],sum=0;

for(i=0;i<N;i++)
{
noise[i]=0.0;
   for(j=0;j<n;j++)
     noise[i]+= rand() % 2;

   noise[i]-= n/2.0;
   noise[i] /= sqrt(n/12.0);
sum+=noise[i];
printf("%d\n",noise[i]);
}
printf("\t%d",sum/N);
return(0);
}
