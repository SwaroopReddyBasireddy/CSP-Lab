#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define NSUM 10000

float gaussrand()
{
	float x = 0;
	int i;
	for(i = 0; i < NSUM; i++)
		x += (float)rand() / RAND_MAX;

	x -= NSUM / 2.0;
	x /= sqrt(NSUM / 12.0);

	return x;
}


int main(){
long int i,N=10000;
float x1[N],x2[N];
float complex noise[N],sum=0.0+0.0*I;

for(i=0;i<N;i++){
    x1[i]=gaussrand();
    x2[i]=gaussrand();
    noise[i]=((float complex)x1[i]+(float complex)(x2[i])*I)/sqrt(2.0);
    sum+=noise[i];
    //var+=pow(x[i],2);
    printf("%.2f %+.2fI \n",crealf(noise[i]),cimagf(noise[i]));
}
printf("\t%f + %fI",creal(sum/N),cimag(sum/N));
return(0);
}

