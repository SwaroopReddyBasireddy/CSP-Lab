#include <stdio.h>
#include <math.h>


main(){

int CP=10,a[64],i,c[64+CP];
for(i=0;i<64;i++)
    a[i]=i+1;

for(i=0;i<64+CP;i++){
    if(i<CP)
       c[i]=a[64-CP+i];
    else
      c[i]=a[i-CP];

}
for(i=0;i<64+CP;i++)
    printf("%d\t",c[i]);
}
