///////////////////////////////////////////////////////////////////////////
//
// File: Convolution.C
//
// Author: Guoling Tong
//
// Purpose: Implementation of the Convolution class, which provides routines
//          for various convolutions.
// 
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: Convolution.C,v 1.2 2004/11/28 03:37:52 joeh Exp $
//
// Revision History
//
// $Log: Convolution.C,v $
// Revision 1.2  2004/11/28 03:37:52  joeh
// New Major Version will no longer compile on AIX
// but will compile on newer versions of gcc
// tested on gcc 3.2 gcc 3.3 and gcc 3.4
// It should run on Redhat 9.0/Fedora Core 1/Fedora Core 2/Fedora Core 3
//
// Revision 1.3  1999/10/08 18:41:22  sarang
// revision
//
// Revision 1.2  1999/10/08 18:16:51  sarang
// revision
//
// Revision 1.1  1999/10/01 15:29:06  kem
// Initial revision
//
// Revision 1.1  1999/09/30 17:39:06  kem
// Initial revision
//
// Revision 1.2  1999/09/29 14:44:46  tong
// Introduced file
//
// Revision 1.1  1999/07/20 16:25:34  RAZORDB
// Initial revision
//
// Revision 1.1  1999/07/15 22:52:25  rst
// Initial revision
//
///////////////////////////////////////////////////////////////////////////


#include <iostream.h>
#include <stdlib.h>
#include <math.h>

#include <ADL/Array3D.h>
#include <TDL/ITX_FFT.h>
#include <TDL/Convolution.h>

#define mmin(x,y) (((x) < (y)) ? (x) : (y))
#define mmax(x,y) (((x) > (y)) ? (x) : (y))

# define static members;
template <class T>
float Convolution<T>::_kernX[_MaxSize*_Oversample+1]; 
template <class T>
float Convolution<T>::_kernY[_MaxSize*_Oversample+1]; 
template <class T>
float Convolution<T>::_kernZ[_MaxSize*_Oversample+1]; 

template <class T>
float Convolution<T>::_RatioX; 
template <class T>
float Convolution<T>::_RatioY; 
template <class T>
float Convolution<T>::_RatioZ; 

template <class T>
int Convolution<T>::_kernLengthX; 
template <class T>
int Convolution<T>::_kernLengthY; 
template <class T>
int Convolution<T>::_kernLengthZ; 

template <class T>
void Convolution<T>:: initialize(ConvolutionType flag)
{
  if(flag == LINEAR)
    {
      initX(linear);
      initY(linear);
      initZ(linear);
    }
  else if(flag == CUBIC)
    {
      initX(cubic);
      initY(cubic);
      initZ(cubic);
    }
  else if(flag == GAUSS)
    {
      initX(gauss);
      initY(gauss);
      initZ(gauss);
    }
  else if (flag == FFT)
    ;
  else if (flag == TRI)
    ;
  else
    {
      cerr<<"Invalid resampling type!"<<endl;
      exit(1);
    }
}


template <class T>
void Convolution<T>:: convolution(ConvolutionType flag,T const *volIn,int x1,int y1,int z1,
			    T *volOut,int x2,int y2,int z2)
{

  if(x2<x1) 
    _RatioX=(float)x1/x2;
  else 
    _RatioX=(float)x2/x1;
  if(y2<y1)
    _RatioY=(float)y1/y2;
  else
    _RatioY=(float)y2/y1;
  if(z2<z1)
    _RatioZ=(float)z1/z2;
  else
    _RatioZ=(float)z2/z1;


  if(flag==LINEAR)
    {
      _kernLengthX =(int) ceil(_RatioX)*2+1;
      _kernLengthY =(int) ceil(_RatioY)*2+1;
      _kernLengthZ =(int) ceil(_RatioZ)*2+1; 
    }
  else if(flag==GAUSS)
    {
      _kernLengthX =(int) ceil(3*_RatioX)*2+1;
      _kernLengthY =(int) ceil(3*_RatioY)*2+1;
      _kernLengthZ =(int) ceil(3*_RatioZ)*2+1; 
    }
   else if(flag==CUBIC)
    {
      _kernLengthX =(int) ceil(2*_RatioX)*2+1;
      _kernLengthY =(int) ceil(2*_RatioY)*2+1;
      _kernLengthZ =(int) ceil(2*_RatioZ)*2+1; 
    }
  
  initialize(flag);

  if(flag == FFT)
    fft(volIn,x1,y1,z1,volOut,x2,y2,z2);

  else if(flag == TRI)
    trilinear(volIn,x1,y1,z1,volOut,x2,y2,z2);

  else
    coeffBin(volIn,x1,y1,z1,volOut,x2,y2,z2);  
}

template <class T>
void Convolution<T>:: convolution(ConvolutionType flag,float para,T const *volIn,
				  int x1,int y1,int z1,T *volOut)
{
  if(flag==GAUSS)
    {
      _kernLengthX =(int) ceil(3*para)*2+1; // the "3" here is because we use 3 for the
      _kernLengthY =(int) ceil(3*para)*2+1; // sigma value in Gauss kernel function.
      _kernLengthZ =(int) ceil(3*para)*2+1; 
    }
   else if(flag==CUBIC)
    {
      _kernLengthX =(int) ceil(2*para)*2+1;
      _kernLengthY =(int) ceil(2*para)*2+1;
      _kernLengthZ =(int) ceil(2*para)*2+1; 
    }
  else
    {
      _kernLengthX =(int) ceil(para)*2+1;
      _kernLengthY =(int) ceil(para)*2+1;
      _kernLengthZ =(int) ceil(para)*2+1; 
    }
  initialize(flag);

  coeffBin(volIn,x1,y1,z1,volOut,x1,y1,z1);  
}



// FFT resampling
template <class T>
void Convolution<T>::fft(T const *I,int ix,int iy,int iz,T *O,int ox,int oy,int oz)
{
  /*************************************************************************
   * If necessary, make (oz) to be a multiple of (iz), i.e, make (oz) > (iz).
   * Also, make (oy) to be a multiple of (iy), i.e, make (oy) > (iy),
   * and make (ox) to be a multiple of (ix), i.e make (ox) > (ix).
   * This enables us to do decimation by downsampling in the spatial domain.
   *************************************************************************/
 

   int zinc = 1; 
   if (oz <= iz) 
     {
       zinc = (iz / oz) + 1; 
       oz *= zinc; 		// New (oz).
     }   
   
   int yinc = 1; 
   if (oy <= iy) 
     {
       yinc = (iy / oy) + 1; 
       oy *= yinc; 		// New (oy).
     }   

   int xinc = 1; 
   if (ox <= ix) 
     {
       xinc = (ix / ox) + 1; 
       ox *= xinc; 		// New (ox).
     }   

/***************************************************************
* Work space array. The array is made big enough to hold the 
* +half frequency half of the FFT of the real output sequence.
***************************************************************/

   int oxc = 2 * ((ox + 2) / 2);

   Array3D<float> W(oz, oy, oxc); 

/******************
* Clear work space.
******************/

   W = 0.0f; 

/****************************
* Make a copy of input array.
****************************/

   int i, j, k; 
   //   int xstride = oxc - ix; 
   //   int ystride = oxc * (oy - iy); 
   int Xstride = 2*((ix+2)/2)-ix;

   T const *iptr = I; 

   FFT_real *wptr=new FFT_real[iz*iy*(ix+Xstride)];
   FFT_real *temp=wptr;
   for (i = iz; i; i--)
      for (j = iy; j; j--, wptr += Xstride)
	 for (k = ix; k; k--)
	    *wptr++ = *iptr++; 
   wptr=temp;

/*****************************
* Forward real to complex FFT.
*****************************/

   ITX_FFT fftF;
   fftF.initialize(3, ITX_FFT::R2C,ix,iy,iz);
   fftF.fftForward(wptr);

   temp=(FFT_real *)W.data();
   FFT_real *temp2=wptr;
   for (i = 0; i<iz; i++)
     for (j = 0; j<iy; j++)
       for (k=0; k<ix+Xstride; k++)
	 {
	   W[i][j][k] = *wptr++;
	 }
   wptr=temp2;
   delete [] wptr;

/******************************
* Pad with zeros in the middle.
******************************/

   int znyqst = iz / 2; 		
   int ynyqst = iy / 2; 		

// Move 3 blocks of the FFT around and insert zeros.
   int ip, jp; 
   int ixc = 2 * ((ix + 2) / 2); 
   for (i = 0; i <= znyqst; i++) {
      for (j = iy - 1, jp = oy - 1; j > ynyqst; j--, jp--) {
	 for (k = 0; k < ixc; k++) {
	    W[i][jp][k] = W[i][j][k]; 
	    W[i][j][k] = 0.0f; 
	 }
      }
   }
   for (i = iz - 1, ip = oz - 1; i > znyqst; i--, ip--) {
      for (j = iy - 1, jp = oy - 1; j > ynyqst; j--, jp--) {
	 for (k = 0; k < ixc; k++) {
	    W[ip][jp][k] = W[i][j][k]; 
	    W[i][j][k] = 0.0f; 
	 }
      }
      for (j = 0; j <= ynyqst; j++) {
	 for (k = 0; k < ixc; k++) {
	    W[ip][j][k] = W[i][j][k]; 
	    W[i][j][k] = 0.0f; 
	 }
      }
   }

/*******************************************************
* For even number of elements in the z and y dimensions.
*******************************************************/

// Z = Znyqst plane.
   if (!(iz % 2)) {
      ip = oz - znyqst; 
      for (j = 0; j < iy; j++) {
	 for (k = 0; k < ixc; k++) {
	    W[znyqst][j][k] /= 2.0f; 
	    W[ip][j][k] = W[znyqst][j][k]; 
	 }
      }
   }

// Y = Ynyqst plane.
   if (!(iy % 2)) {
      jp = oy - ynyqst; 
      for (i = 0; i< iz; i++) {
	 for (k = 0; k < ixc; k++) {
	    W[i][ynyqst][k] /= 2.0f; 
	    W[i][jp][k] = W[i][ynyqst][k];
	 }
      }
   }

/*****************************
* Complex to real inverse FFT.
*****************************/
   ITX_FFT fftI;
   fftI.initialize(3,ITX_FFT::C2R,ox,oy,oz);
   fftI.fftInverse((FFT_complex *)W.data());


/**********************************************
* Copy into output Array. Skip extra points.
* Also, normalize by length of input sequence.
**********************************************/
   T *optr = O; 
   float alpha = 1.0f / (iz * iy * ix); 
   for (k = 0; k < oz; k += zinc)
      for (j = 0; j < oy; j += yinc)
	 for (i = 0; i < ox; i += xinc)
	   *optr++ = (T) ((W[k][j][i] > 0.0f) ? 
			  (alpha * W[k][j][i]) : 0.0f); 

}

// trilinear sampling
template <class T>
void Convolution<T>::trilinear(T const *I,int ix,int iy,int iz,T *O,int ox,int oy,int oz)
{
  double ioRatioX = (double) ix/ox;
  double ioRatioY = (double) iy/oy;
  double ioRatioZ = (double) iz/oz;
  
  for (int outZ=0; outZ<oz; outZ++) 
    {
      double z = ioRatioZ * outZ;
      int z0 = (int) z;
      double dz = z-z0;
      int z1 = (z0 < iz-1) ? z0 + 1 : iz-1;
      
      for (int outY=0; outY<oy; outY++) 
	{
	  double y = ioRatioY * outY;
	  int y0 = (int) y;
	  double dy = y-y0;
	  int y1 = (y0 < iy-1) ? y0 + 1 : iy-1;

	  for (int outX=0; outX<ox; outX++) 
	    {
	      double x = ioRatioX * outX;
	      int x0 = (int) x;
	      double dx = x-x0;
	      int x1 = (x0 < ix-1) ? x0 + 1 : ix-1;

	      T in000 = I[x0+ix*(y0+iy*z0)];
	      T in001 = I[x1+ix*(y0+iy*z0)];
	      T in010 = I[x0+ix*(y1+iy*z0)];
	      T in011 = I[x1+ix*(y1+iy*z0)];
	      T in100 = I[x0+ix*(y0+iy*z1)];
	      T in101 = I[x1+ix*(y0+iy*z1)];
	      T in110 = I[x0+ix*(y1+iy*z1)];
	      T in111 = I[x1+ix*(y1+iy*z1)];

	      double in00x = in000 + dx*(in001-in000);
	      double in01x = in010 + dx*(in011-in010);
	      double in10x = in100 + dx*(in101-in100);
	      double in11x = in110 + dx*(in111-in110);
	      
	      double in0yx = in00x + dy*(in01x-in00x);
	      double in1yx = in10x + dy*(in11x-in10x);
	      
	      *O++ = (T) (in0yx + dz*(in1yx-in0yx));
	    }
	}
    }
} 


template <class T>
void Convolution<T>::coeffBin(T const *IN,int x1,int y1,int z1,T *volOut,int x2,int y2,int z2)
{

  float *OUTZ=new float[x1*y1*z2];
  if ( OUTZ == NULL) {
    cerr << "Failed to allocate memory!" << endl;
    exit(1);
  }

  float *OUTY=new float[x1*y2*z2];
  if ( OUTY == NULL) {
    cerr << "Failed to allocate memory!" << endl;
    exit(1);
  }

  T *OUTX=volOut;

  int x,ii,dii,ff,dff, outlenX,outlenY,outlenZ,inlenX,inlenY,inlenZ;
  int a,b,c;
  double val;
  int i,j,k,l;  

  inlenX=x1;
  inlenY=y1;
  inlenZ=z1;
  outlenX=x2;
  outlenY=y2;
  outlenZ=z2;

  // Interpolation in z-dimension;
  a=inlenZ/outlenZ;
  b=inlenZ%outlenZ;
  ii=0;  // ii indexes into bin;
  ff=outlenZ/2; // ff is fractional remainder;
  x=b*_Oversample;
  dii=x/outlenZ;  // dii is ii increment;
  dff=x%outlenZ;  // dff is ff increment;

  k=0;
  c=_kernLengthZ/2;

  float sum=0;
  for(i=0;i<inlenX;i++)
    {
      for(j=0;j<inlenY;j++)
	{
	  for(x=0;x<outlenZ;x++)
	    {
	      // compute convolution centered at current position;
	      val =0;

	      for (l=mmax(0,k-(c-1));l<=k;l++)
		{
		  int t = k-l;
		  sum+=_kernZ[t*_Oversample+ii];
		}

	      for (l=k+1; l<=mmin(inlenZ-1,k+c);l++)
		{
		  int t = l-k;
		  sum+=_kernZ[t*_Oversample-ii];
		}

	      if(ii==0 && k>=c)
		{
		  sum+=_kernZ[c*_Oversample+ii];
		}

	      for (l=mmax(0,k-(c-1));l<=k;l++)
		{
		  int t = k-l;
		  val+=IN[i+inlenX*(j+inlenY*l)]*_kernZ[t*_Oversample+ii]/sum;
		}

	      for (l=k+1; l<=mmin(inlenZ-1,k+c);l++)
		{
		  int t = l-k;
		  val+=IN[i+inlenX*(j+inlenY*l)]*_kernZ[t*_Oversample-ii]/sum;
		}

	      if(ii==0 && k>=c)
		{
		  val+=IN[i+inlenX*(j+inlenY*(k-c))]*_kernZ[c*_Oversample+ii]/sum;
		}
	      sum=0;

	      OUTZ[i+inlenX*(j+inlenY*x)]=val;

	      // Bresenham-like algorithm to recenter kernel;
	      if((ff+=dff)>=outlenZ) // check if fractional part overflows;
		{
		  ff-=outlenZ; // normalize;
		  ii++;  // increment integer part;
		}
	      if((ii+=dii)>=_Oversample) // check if integer part overflows;
		{  
		  ii = ii-_Oversample; // normalize;
		  k++; // increment input pointer;
		}
	      k+=a;
	    }
	  k=0;
	  ii=0;
	}
      k=0;
      ii=0;
    }

  // Interpolation in y-dimension;
  inlenZ=outlenZ;
  a=inlenY/outlenY;
  b=inlenY%outlenY;
  ii=0;  // ii indexes into bin;
  ff=outlenY/2; // ff is fractional remainder;
  x=b*_Oversample;
  dii=x/outlenY;  // dii is ii increment;
  dff=x%outlenY;  // dff is ff increment;

  j=0;
  c=_kernLengthY/2;

  sum=0;
  for(i=0;i<inlenX;i++)
    {
      for(k=0;k<inlenZ;k++)
	{
	  for(x=0;x<outlenY;x++)
	    {
	      val =0;
	    
	      for ( l=mmax(0,j-(c-1));l<=j;l++)
		{
		  int t = j-l;
		  sum+=_kernY[t*_Oversample+ii];
		}

	      for (l=j+1; l<=mmin(inlenY-1,j+c);l++)
		{
		  int t = l-j;
		  sum+=_kernY[t*_Oversample-ii];
		}

	      if(ii==0 && j>=c)
		{
		  sum+=_kernY[c*_Oversample+ii];
		}

	      for ( l=mmax(0,j-(c-1));l<=j;l++)
		{
		  int t = j-l;
		  val+=OUTZ[i+inlenX*(l+inlenY*k)]*_kernY[t*_Oversample+ii]/sum;
		}

	      for (l=j+1; l<=mmin(inlenY-1,j+c);l++)
		{
		  int t = l-j;
		  val+=OUTZ[i+inlenX*(l+inlenY*k)]*_kernY[t*_Oversample-ii]/sum;
		}

	      if(ii==0 && j>=c)
		{
		  val+=OUTZ[i+inlenX*((j-c)+inlenY*k)]*_kernY[c*_Oversample+ii]/sum;
		}
	      sum=0;

	      OUTY[i+inlenX*(x+outlenY*k)]=val;
	      
	      // Bresenham-like algorithm to recenter kernel;
	      if((ff+=dff)>=outlenY) // check if fractional part overflows;
		{
		  ff-=outlenY; // normalize;
		  ii++;  // increment integer part;
		}
	      if((ii+=dii)>=_Oversample) // check if integer part overflows;
		{  
		  ii-=_Oversample; // normalize;
		  j++; // increment input pointer;
		}
	      j+=a;
	    }
	  j=0;
	  ii=0;
	}
      j=0;
      ii=0;
    }
  delete [] OUTZ;

  // interpolation in x-dimension;
  inlenY=outlenY;
  a=inlenX/outlenX;
  b=inlenX%outlenX;
  ii=0;  // indexes into bin;
  ff=outlenX/2; // ff is fractional remainder;
  x=b*_Oversample;
  dii=x/outlenX;  // dii is ii increment;
  dff=x%outlenX;  // dff is ff increment;

  i=0;
  c=_kernLengthX/2;

  sum=0;
  for(k=0;k<inlenZ;k++)
    {
      for(j=0;j<inlenY;j++)
	{
	  for(x=0;x<outlenX;x++)
	    {
	      // compute convolution centered at current position;
	      val =0;

	      for (l=mmax(0,i-(c-1));l<=i;l++)
		{
		  int t = i-l;
		  sum+=_kernX[t*_Oversample+ii];
		}

	      for (l=i+1; l<=mmin(inlenX-1,i+c);l++)
		{
		  int t = l-i;
		  sum+=_kernX[t*_Oversample-ii];
		}

	      if(ii==0 && i>= c)
		{
		  sum+=_kernX[c*_Oversample+ii];
		}

	      for (l=mmax(0,i-(c-1));l<=i;l++)
		{
		  int t = i-l;
		  val+=OUTY[l+inlenX*(j+inlenY*k)]*_kernX[t*_Oversample+ii]/sum;
		}

	      for (l=i+1; l<=mmin(inlenX-1,i+c);l++)
		{
		  int t = l-i;
		  val+=OUTY[l+inlenX*(j+inlenY*k)]*_kernX[t*_Oversample-ii]/sum;
		}

	      if(ii==0 && i>= c)
		{
		  val+=OUTY[i-c+inlenX*(j+inlenY*k)]*_kernX[c*_Oversample+ii]/sum;
		}
	      sum=0;


	      if(val<0) val=0;
	      if(val>255) val=255;
	      OUTX[x+outlenX*(j+inlenY*k)]=(T) (val);
	      
	      // Bresenham-like algorithm to recenter kernel;
	      if((ff+=dff)>=outlenX) // check if fractional part overflows;
		{
		  ff-=outlenX; // normalize;
		  ii++;  // increment integer part;
		}
	      if((ii+=dii)>=_Oversample) // check if integer part overflows;
		{  
		  ii-=_Oversample; // normalize;
		  i++; // increment input pointer;
		}
	      i+=a;
	    }
	  i=0;
	  ii=0;
	}
      i=0;
      ii=0;
    }
  delete [] OUTY;

}


// initialize kernels with the function in the argument;
template <class T>
void Convolution<T>::initX(double (*func)(double))
{
  int i,size;
  size=(_kernLengthX/2)*_Oversample;

  double sum=0;

  for(i=0;i<=size;i++)
    {
      _kernX[i]=(*func)((double)i/size);
      sum+=_kernX[i];
    }

  sum=sum*2.0/(_Oversample);

  for(i=size+1;i<=_MaxSize*_Oversample;i++)
     _kernX[i]=0;
  /*
  for(i=0;i<=size;i++)
    {
      _kernX[i]/=sum;
    }
  */
}

template <class T>
void Convolution<T>::initY(double (*func)(double))
{
  int i,size;
  size=(_kernLengthY/2)*_Oversample;
  
  double sum=0;

  for(i=0;i<=size;i++)
    {
      _kernY[i]=(*func)((double)i/size);
      sum+=_kernY[i];
    }
  for(i=size+1;i<=_MaxSize*_Oversample;i++)
    _kernY[i]=0;

  sum= sum*2.0/(_Oversample);
  /*
  for(i=0;i<=size;i++)
    {
      _kernY[i]/=sum;
    }
  */
}

template <class T>
void Convolution<T>::initZ(double (*func)(double))
{
  int i,size;
  size=_kernLengthZ/2*_Oversample;
 
  double sum=0;

  for(i=0;i<=size;i++)
    {
      _kernZ[i]=(*func)((double)i/size);
      sum+=_kernZ[i];
    }

  for(i=size+1;i<=_MaxSize*_Oversample;i++)
    _kernZ[i]=0;

  sum=sum*2.0/(_Oversample);

  /*
  for(i=0;i<=size;i++)
    {
      _kernZ[i]/=sum;
    }
  */
}



// interpolating functions used in init;
template <class T>
double Convolution<T>::linear( double x)
{
  if(fabs(x)>=1)
    return 0;
  else
    return 1-fabs(x);
}

template <class T>
double Convolution<T>::cubic( double x)
{
  if(fabs(x)>=1)
    return 0;
  else if(fabs(x)<1 && fabs(x)>=0.5)
    return 0.5*fabs(2*2*2*x*x*x)-5*0.5*fabs(2*2*x*x)+8*0.5*fabs(2*x)-2;
  else
    return (0.5+2)*fabs(2*2*2*x*x*x)-(0.5+3)*fabs(2*2*x*x)+1;
}

template <class T>
double Convolution<T>::gauss( double x)
{
  if(fabs(x)>=1)
    return 0;
  else
    return exp(-3*3*x*x/2);  
}


#ifdef GNU_COMPILER
template class Convolution<unsigned char>;
template class Convolution<unsigned short>;
template class Convolution<short>;
template class Convolution<int>;
template class Convolution<float>;
#endif

#ifdef SGI_COMPILER
#pragma instantiate Convolution<unsigned char>
#pragma instantiate Convolution<unsigned short>
#pragma instantiate Convolution<short>
#pragma instantiate Convolution<int>
#pragma instantiate Convolution<float>
#endif


#ifdef RS6K

#pragma define(Convolution<unsigned char>)
#pragma define(Convolution<unsigned short>)
#pragma define(Convolution<short>)
#pragma define(Convolution<int>)
#pragma define(Convolution<float>)

#endif
