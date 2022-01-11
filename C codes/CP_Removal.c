#include <stdio.h>
#include <math.h>


main(){

int CP=10,a[74],i,b[CP],c[64];
for(i=0;i<64+CP;i++)
    a[i]=i+1;

for(i=0;i<64;i++){
    c[i]=a[CP+i];

}
for(i=0;i<64;i++)
    printf("%d\t",c[i]);
}
