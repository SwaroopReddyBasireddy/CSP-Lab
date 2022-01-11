#include <stdio.h>
#include <complex.h>
#include <math.h>

void add(int *p,int *q,int *r)
{
int i;
for(i=0;i<10;i++)
r[i]=p[i]+q[i];
}
main(){

int a[5],b[5],c[5],i;
printf("Enter the values of a:");
for(i=0;i<5;i++)
{
    scanf("%d",&a[i]);
}
printf("Enter the values of b:");
for(i=0;i<5;i++)
{
    scanf("%d",&b[i]);
}
add(a,b,c);
for(i=0;i<5;i++)
 printf("%d\t",c[i]);

}



for(i=0;i<N;i++){
      w[i]=twiddle_factor(N,i);
      printf("%.3f %+.3fI \n",crealf(w[i]),cimagf(w[i]));

   }

float complex twiddle_factor(N,i)
{
float c;
float complex x;
c=2*PI/N;
x=cos(c*i) + (sin(c*i))*I;

return(x);
}

for(i=0;i<N;i++){
      w[i]=twiddle_factor(N,i);
      printf("%.3f %+.3fI \n",crealf(w[i]),cimagf(w[i]));

   }
