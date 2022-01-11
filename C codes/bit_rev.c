#include <stdio.h>
#include <math.h>
#include <complex.h>

void bit_rev(int N,int *y)
{
  int i,j,n=round(log10(N)/log10(2)),bin[N][n],dec[N],t=0,x[N],p;
  for(i=0;i<N;i++){
    x[i]=y[i];
    dec[i]=0;
    //y[i]=0;
    //printf("%d\t",x[i]);
  }


  for(j=0;j<N;j++){
        p=j;
    for(i=0;i<n;i++)
     {
      bin[j][i]= p%2;

        p/=2;
      //printf("%d",bin[j][i]);
        }
      //printf("\t");

      }
  printf("\n");
  for(j=0;j<N;j++){
    for(i=0;i<n;i++)
       dec[j] +=bin[j][i]*round(pow(2,n-i-1));
     //printf("%d\t",dec[j]);

  }
  //printf("\n");
for(i=0;i<N;i++)
{
    j=dec[i];
    y[j]=x[i];
    //printf("%d\t",y[j]);

    }

}

main(){
    int N=16,x[N],y[N],i;

    for(i=0;i<N;i++){
        x[i]=i;
        printf("%d\t",x[i]);
    }

    bit_rev(N,x);
    printf("\n After Bit Reversal\n");
    for(i=0;i<N;i++)
       printf("%d\t",x[i]);

}
