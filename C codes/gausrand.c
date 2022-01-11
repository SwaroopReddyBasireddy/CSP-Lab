#include <stdlib.h>
#include <math.h>

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
long int sum=0,i,N=100000;
float x[N],var;

for(i=0;i<N;i++){
    x[i]=gaussrand();
    sum+=x[i];
    var+=pow(x[i],2);
    printf("%f\n",x[i]);
}
printf("\t%ld \t%f",sum/N,var/N);
return(0);
}
