////////////////////////////////////////////////////////////////////////////
//
//	FFT routines for 3-Dimensional real data.  
//	Simple 3-pass DFT
//
//	Code for 1-D DFT (fftl) and 1-D REAL DFT (fftr)
//	are from Numerical Recipes in C.
//
//   NOTE: Input data must be stored in REAL array of size (NX+2,NY,NZ)
//	 with the original data stored in the first NX columns. Last two
//	 columns are for the sfft3du/fftr routine.
// 
//	 Input array is over-written with output. For REVERSE, last two
//	 columns are meaningless.
//
//Author:
//	Keith Doolittle
//	Routines fftr() and fftl() from Numerical Recipes in C
//
////////////////////////////////////////////////////////////////////////////
#include <math.h>
#include <TDL/TDLfft.h>

void TDLfft::pad(Array3D<float> &D)
{
   int nx = D.getXsize();
   int ny = D.getYsize();
   int nz = D.getZsize();
   Array3D<float> nD(nz,ny,nx+2);

   for(int k=0;k<nz;k++)
     for(int j=0;j<ny;j++) {
       for(int i=0;i<nx;i++)
            nD[k][j][i] = D[k][j][i];
     nD[k][j][nx]   = 0.f;
     nD[k][j][nx+1] = 0.f;
     }

   D.setDim(nz,ny,nx+2);
   D = nD;
}


void TDLfft::unpad(Array3D<float> &D)
{
   int nx = D.getXsize()-2;
   int ny = D.getYsize();
   int nz = D.getZsize();
   Array3D<float> nD(nz,ny,nx);

   for(int k=0;k<nz;k++)
     for(int j=0;j<ny;j++)
       for(int i=0;i<nx;i++)
            nD[k][j][i] = D[k][j][i];

   D.setDim(nz,ny,nx);
   D = nD;
}

void TDLfft::scale(Array3D<float> &D)
{
  int i,j,k;
  int nx = D.getXsize()-2;
  int ny = D.getYsize();
  int nz = D.getZsize();
  float fscl = 1.f/((float)(nx*ny*nz/2));

  for(k=0;k<nz;k++)
    for(j=0;j<ny;j++)
      for(i=0;i<nx;i++)
         D[k][j][i] = D[k][j][i] * fscl;
}


/***********************************************************************
   Three pass 3-D DFT. Order = XYZ (FORWARD) or ZYX (REVERSE). For
   X dimension, use fftr which takes input array of size N+2 and
   output's complex array of size (N+2)/2 (for REAL sequences). Others
   use standard fftl complex DFT.
***********************************************************************/
//
// NOTE: assumes data is already padded (nx = nx+2)
//

bool TDLfft::FFT(Array3D<float> &D, FFTtype dir)
{
  int nx = D.getXsize()-2;
  int ny = D.getYsize();
  int nz = D.getZsize();

  if((nx<1)||(ny<1)||(nz<1)||(nx>FFT_MAX_DIM)||
     (ny>FFT_MAX_DIM)||(nz>FFT_MAX_DIM))
	return(false);

  int i1,i2,i3;

  bwcomplex *x;

  x = (bwcomplex*)D.data();

  int xdim    = nx+2;
  int Halfnx  = nx/2;
  int cx      = xdim/2;
  int cxy     = cx*ny;
  int xadd    = xdim-nx;
  int pl_size = xdim*ny;

  bwcomplex *cptr;

  if(dir == TDLfft::Forward) {

  // X direction
  for(i3=0;i3<nz;i3++)
   for(i2=0;i2<ny;i2++) {
    for(i1=0;i1<Halfnx;i1++) {
	wk[i1].re = D[i3][i2][2*i1+0];
	wk[i1].im = D[i3][i2][2*i1+1];
    }
    fftr(Halfnx,dir);

    D[i3][i2][2*0+0]      = wk[0].re;
    D[i3][i2][2*0+1]      = 0.f;
    D[i3][i2][2*Halfnx+0] = wk[0].im;
    D[i3][i2][2*Halfnx+1] = 0.f;

    for(i1=1;i1<Halfnx;i1++) {
        D[i3][i2][i1*2+0] = wk[i1].re;
        D[i3][i2][i1*2+1] = wk[i1].im;
    }

    }

  // Y direction
  for(i3=0;i3<nz;i3++)
   for(i1=0;i1<=Halfnx;i1++) {
    for(i2=0;i2<ny;i2++) {
        wk[i2].re = D[i3][i2][i1*2+0];
        wk[i2].im = D[i3][i2][i1*2+1];
    }
    fftl(ny,dir);
    for(i2=0;i2<ny;i2++) {
        D[i3][i2][i1*2+0] = wk[i2].re;
        D[i3][i2][i1*2+1] = wk[i2].im;
    }
    }

  if(nz == 1) 
     return(true);

  // Z direction
  for(i1=0;i1<=Halfnx;i1++)
   for(i2=0;i2<ny;i2++) {
    for(i3=0;i3<nz;i3++) {
        wk[i3].re = D[i3][i2][i1*2+0];
        wk[i3].im = D[i3][i2][i1*2+1];
    }
    fftl(nz,dir);
    for(i3=0;i3<nz;i3++) {
        D[i3][i2][i1*2+0] = wk[i3].re;
        D[i3][i2][i1*2+1] = wk[i3].im;
    }
    }


  }
  else { /* dir == TDLfft::Reverse */

  if(nz == 1) goto SkipZ;

  // Z direction
  for(i1=0;i1<=Halfnx;i1++) 
   for(i2=0;i2<ny;i2++) {
    for(i3=0;i3<nz;i3++) {
        wk[i3].re = D[i3][i2][i1*2+0];
        wk[i3].im = D[i3][i2][i1*2+1];
    }
    fftl(nz,dir);
    for(i3=0;i3<nz;i3++) {
        D[i3][i2][i1*2+0] = wk[i3].re;
        D[i3][i2][i1*2+1] = wk[i3].im;
    }
    }

SkipZ:

  // Y direction
  for(i3=0;i3<nz;i3++)
   for(i1=0;i1<=Halfnx;i1++) {
    for(i2=0;i2<ny;i2++) {
        wk[i2].re = D[i3][i2][i1*2+0];
        wk[i2].im = D[i3][i2][i1*2+1];
    }
    fftl(ny,dir);
    for(i2=0;i2<ny;i2++) {
        D[i3][i2][i1*2+0] = wk[i2].re;
        D[i3][i2][i1*2+1] = wk[i2].im;
    }
    }

  // X direction
  for(i3=0;i3<nz;i3++)
   for(i2=0;i2<ny;i2++) {

    wk[0].re = D[i3][i2][0*2+0];
    wk[0].im = D[i3][i2][Halfnx*2+0];

    for(i1=1;i1<=Halfnx;i1++) {
        wk[i1].re = D[i3][i2][i1*2+0];
        wk[i1].im = D[i3][i2][i1*2+1];
    }
    fftr(Halfnx,dir);
    for(i1=0;i1<Halfnx;i1++) {
        D[i3][i2][i1*2+0] = wk[i1].re;
        D[i3][i2][i1*2+1] = wk[i1].im;
    }
    }
  scale(D);
 }

 return(true);
}


/* Takes FFT of real data which is 2N length */

void TDLfft::fftr(int n, FFTtype direc)
{
float *dat,*fptr;
int i,i1,i2,i3,i4,n2p3;
float c1,c2,h1r,h1i,h2r,h2i;
double wr,wi,wpr,wpi,wtemp,theta;

fptr = (float*)&wk[0];
dat  = fptr-1;

c1 = 0.5;
theta = 3.141592653589793/(double)n;
if(direc == TDLfft::Forward) {
  c2 = -0.5; fftl(n,TDLfft::Forward);
  }
else {
  c2 = 0.5; theta = -theta;
  }

wtemp = sin(0.5*theta);
wpr = -2.0*wtemp*wtemp;
wpi = sin(theta);
wr = 1.0+wpr; wi = wpi;
n2p3 = 2*n+3;
for(i=2;i<=n/2;i++) {
  i4 = 1+(i3=n2p3-(i2=1+(i1=i+i-1)));
  h1r = c1*(dat[i1]+dat[i3]);
  h1i = c1*(dat[i2]-dat[i4]);
  h2r = -c2*(dat[i2]+dat[i4]);
  h2i = c2*(dat[i1]-dat[i3]);
  dat[i1] = h1r+wr*h2r-wi*h2i;
  dat[i2] = h1i+wr*h2i+wi*h2r;
  dat[i3] = h1r-wr*h2r+wi*h2i;
  dat[i4] = -h1i+wr*h2i+wi*h2r;
  wr = (wtemp=wr)*wpr-wi*wpi+wr;
  wi = wi*wpr+wtemp*wpi+wi;
  }
if(direc == TDLfft::Forward) {
  dat[1] = (h1r=dat[1])+dat[2];
  dat[2] = h1r-dat[2];
  }
else {
  dat[1] = c1*((h1r=dat[1])+dat[2]);
  dat[2] = c1*(h1r-dat[2]);
  fftl(n,TDLfft::Reverse);
  }
}

#define FFTSWAP(a,b) { tempr = (a); (a)=(b); (b)=tempr; }

void TDLfft::fftl(int nn, FFTtype direc)
{
float *dat,*fptr;
int n,mmax,m,j,istep,i;
double wtemp,wr,wpr,wpi,wi,theta;
float tempr,tempi;

/* This routine is 1 based, so start at wk[0]-sizeof(float) */

fptr = (float*)&wk[0];
dat = fptr-1;
n = nn<<1;
j = 1;
for(i=1;i<n;i+=2) {
 if(j>i) {
  FFTSWAP(dat[j],dat[i])
  FFTSWAP(dat[j+1],dat[i+1])
  }
  m = n>>1;
  while(m >= 2 && j > m) {
   j -= m; m >>=1;
   }
  j += m;
}
mmax = 2;
while (n > mmax) {
  istep = 2*mmax;
  theta = 6.28318530717959/(direc*mmax);
  wtemp = sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi = sin(theta);
  wr = 1.0;
  wi = 0.0;
  for(m=1;m<mmax;m+=2) {
   for(i=m;i<=n;i+=istep) {
     j=i+mmax;
     tempr = wr*dat[j]-wi*dat[j+1];
     tempi = wr*dat[j+1]+wi*dat[j];
     dat[j] = dat[i]-tempr;
     dat[j+1] = dat[i+1]-tempi;
     dat[i] += tempr; dat[i+1] += tempi;
   }
   wr = (wtemp=wr)*wpr-wi*wpi+wr;
   wi = wi*wpr+wtemp*wpi+wi;
  }
  mmax=istep;
}
}


