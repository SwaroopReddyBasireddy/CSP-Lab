#include <math.h> // for RAND, and rand
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
float gaussrand() {
    float u = ((double) rand() / (RAND_MAX)) * 2 - 1;
    float v = ((double) rand() / (RAND_MAX)) * 2 - 1;
    float r = u * u + v * v;
    if (r == 0 || r > 1) return gaussrand();
    float c = sqrt(-2 * log(r) / r);
    return u * c;
}

main(){
int i,N=1000;
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
}
