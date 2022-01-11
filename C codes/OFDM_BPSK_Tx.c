#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>

#define PI 3.14159265
#define NSUM 100

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

void convolution(int n,int m,float complex *h, float complex *x,float complex *y,int N)
  {
    int i,j,k;
    float complex c[n][m],z[n+m-1];
    for(i=0;i<m+n-1;i++)
        z[i]=0.0+0.0*I;
    for(i=0;i<n;i++)
        for(j=0;j<m;j++)
           c[i][j]=h[i]*x[j];

  for(k=0;k<(n+m);k++)
    for(i=0;i<n;i++)
       for(j=0;j<m;j++)
         if((i+j)==k)
           z[k]+=c[i][j];

  if(n+m-1>N)
     for(i=0;i<n+m-1;i++)
       if(i<N)
         y[i]=z[i];

       else
          y[i-N]+=z[i];
  else
      for(i=0;i<N;i++)
          if(i<n+m-1)
            y[i]=z[i];
          else
            z[i]=0;
  }

void ifft_twiddle_factor(float complex *p,int n)
{
int i;
float c;
c=2*PI/n;
for(i=0;i<n;i++)
  p[i]=cos(c*i) + (sin(c*i))*I;
}

void ifft_stage(int N,int stage,float complex *x,float complex *y)
{
  int k=N/pow(2,stage-1);
  float complex w[k];
  int i,j,m;
  //printf("stage:%d\n",stage);
  ifft_twiddle_factor(w,k);
  //for(i=0;i<k;i++)
    //printf("%.3f %+.3fI \n",crealf(w[i]),cimagf(w[i]));

for(j=0;j<(int)(pow(2,stage-1));j++)
for(i=0;i<(int)(N/pow(2,stage));i++)
    {
        m=2*j*N/pow(2,stage);
            y[i+m]=x[i+m]+x[i+m+(int)(N/pow(2,stage))];
        if((int)pow(2,stage)!=N)
            y[i+m+(int)(N/pow(2,stage))]=(x[i+m]-x[i+m+(int)(N/pow(2,stage))])*w[i];
        else
            y[i+m+(int)(N/pow(2,stage))]=(x[i+m]-x[i+m+(int)(N/pow(2,stage))]);
    }
 }

 void ifft_bit_rev(float complex *y,int N)
{
    float complex z[N];
    int i;
    for(i=0;i<N;i++)
    {

        if(i%2!=0 && i<N/2)
        {
        z[i]=y[i];
        y[i]=(y[i-1+N/2]);
        y[i-1+N/2]=z[i];
        }
        else
          z[i]=y[i];

    }
    for(i=0;i<N;i++)
       y[i]/=(float)N;

}

void add_CP(int N,int CP,float complex *a,float complex *c)

    {
     int i;
     float complex b[CP];
    //for(i=0;i<CP;i++)
      //  b[i]=a[64-CP+i];

    for(i=0;i<64+CP;i++){
      if(i<CP)
        c[i]=a[64-CP+i];
      else
        c[i]=a[i-CP];

     }

}




main(){

int n_sym=1000,n_fft=64,snrlen=11,i,j,k,l,m,n,n_bits,s[n_sym*n_fft],s_in[n_sym*n_fft],no_stages,CP=10,n_tap=5,len=n_fft+CP,tau=3,tau1;
float complex x[n_fft],s_ofdm[n_fft],t[n_fft],s_ofdm_cp[n_fft+CP],h[n_tap],y[len],noise[len],r[len],pre[16];
float x1[n_tap],x2[n_tap],x3[len],x4[len];
n_bits=n_sym*n_fft;

for(i=0;i<n_bits;i++)
  {
  s[i]=rand()%2;
  s_in[i]=2*s[i]-1;
  //printf("%d\t",s_in[i]);
   }
for(i=0;i<1;i++)
 {
    for(j=i*n_fft;j<(i+1)*n_fft;j++)
     {
       x[j]=(float complex)s_in[j];
       s_ofdm[j]=0.0+0.0*I;
     }
//Serial to parallel conversions printing
    //for(k=0;k<n_fft;k++)
        //printf("%.2f %+.2fI\t",crealf(x[k]),cimagf(x[k]));
//IFFT
no_stages=round(log10(n_fft)/log10(2));
   for(k=0;k<no_stages;k++)
     {
       if(k==0)
         ifft_stage(n_fft,k+1,x,s_ofdm);
       else
        {
          for(l=0;l<n_fft;l++)
           t[l]=s_ofdm[l];
          ifft_stage(n_fft,k+1,t,s_ofdm);
        }
//printing outputs at each stage of butterfly algorithm
 //printf("\n stage:%d\n",k+1);
 //for(m=0;m<n_sym;m++)
   //printf("%.2f %+.2fI \t",crealf(s_ofdm[m]),cimagf(s_ofdm[m]));
 }
 ifft_bit_rev(s_ofdm,n_fft);
//The outputs of IFFT
 //printf("\nIFFT of s_ofdm:\n");
 //for(m=0;m<n_fft;m++)
   //printf("%.2f %+.2fI \t",crealf(s_ofdm[m]),cimagf(s_ofdm[m]));
   //printf("\n After CP\n");
  add_CP(n_fft,CP,s_ofdm,s_ofdm_cp);
  //for(m=0;m<n_fft+CP;m++)
    //printf("%.2f %+.2fI \t",crealf(s_ofdm_cp[m]),cimagf(s_ofdm_cp[m]));

//Convolve with Multi tap channel

 for(m=0;m<n_tap;m++){
    x1[m]=gaussrand();
    x2[m]=gaussrand();
    h[m]=((float complex)x1[m]+(float complex)(x2[m])*I)/sqrt(n_tap*2.0);
    }
 convolution(n_tap,n_fft+CP,h,s_ofdm_cp,y,n_fft);

//Addition of white Gausian Noise

 for(m=0;m<len;m++){
    x3[m]=gaussrand();
    x4[m]=gaussrand();
    noise[m]=((float complex)x1[m]+(float complex)(x2[m])*I)/sqrt(2.0);
    r[m]=y[m]+noise[m];
    printf("%.2f %+.2fI\t",crealf(r[m]),cimagf(r[m]));
    }


// Receiver Part
 for(m=0;m<16;m++)
     pre[m]=s_ofdm[m];

    tau1=time_sync(pre,tau,r);

}

}
