#include <stdio.h>
#include <complex.h>

void main(){
float z1[5],z2[5];
float complex z[5];
int i;
for(i=0;i<5;i++)
{
    z1[i]=i+1;
    z2[i]=5-i;
    z[i]= (float complex)z1[i]+(float complex)(z2[i])*I;
    printf("z[%d]=%.2f %+.2fI\n",i+1,crealf(z[i]),cimagf(z[i]));
}
//z1=3.0 + 4.0*I;
//z2=3.0 - 5.0*I;
//z=z1+z2;



}
