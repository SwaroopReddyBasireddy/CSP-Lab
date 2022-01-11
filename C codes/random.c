#include<stdio.h>
#include<math.h>
#include<stdlib.h>

int main()
{
int i,n=100000,a[n],s[n];

for(i=0;i<n;i++)
{
a[i]=rand()%2;
s[i]=2*a[i]-1;
printf("%d\t",a[i]);
}
return(0);
}

