#include <stdio.h>
#include <complex.h>
#include <math.h>

#define PI 3.14159265

void twiddle_factor(float complex *p,int N)
{
int i;
float c;
//float complex x;
c=2*PI/N;
for(i=0;i<N;i++)
  p[i]=cos(c*i) + (sin(c*i))*I;

return(0);
}

main()
{
int r,i,N=8;
float complex w[N],*p;
p=&w;
twiddle_factor(w,N);

for(i=0;i<N;i++)

   printf("%.3f %+.3fI \n",crealf(w[i]),cimagf(w[i]));


}
