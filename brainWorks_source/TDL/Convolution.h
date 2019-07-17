///////////////////////////////////////////////////////////////////////////
//
// File: Convolution.h 
//
// Author: Guoling Tong
//
// Purpose: Class Convolution header file, providing routines for convolutions. 
//         
// 
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: Convolution.h,v 1.1 2004/11/15 04:44:07 joeh Exp $
//
// Revision History
//
// $Log: Convolution.h,v $
// Revision 1.1  2004/11/15 04:44:07  joeh
// first check in of BrainWorks
//
// Revision 1.2  1999/10/08 18:42:38  sarang
// revision
//
// Revision 1.1  1999/10/01 15:28:03  kem
// Initial revision
//
// Revision 1.2  1999/09/29 14:43:47  tong
// continue testing
//
// Revision 1.1  1999/07/20 16:24:19  RAZORDB
// Initial revision
//
// Revision 1.2  1999/07/15 22:50:15  rst
//
///////////////////////////////////////////////////////////////////////////
#include <OS.h>

#ifndef __CONVOLUTION_H__
#define __CONVOLUTION_H__

// constants used in coefficient bins resampling

#define _MaxSize 10
#define _Oversample 512  // subdivisions per pixel;


template <class T>
class Convolution 
{

 public:

  // Resampling types;
  enum ConvolutionType { TRI, FFT,LINEAR, CUBIC, GAUSS };

  // interface for the convolution;
  static void convolution(ConvolutionType flag,T const *volIn,int x1,int y1,int z1,
		       T *volOut,int x2,int y2,int z2);

  static void convolution(ConvolutionType flag,float para, T const *volIn,int x1,int y1,int z1,
			  T *volOut);
  
 protected:

  static int _kernLengthX;  // kernel length in three dimenstions respectively;
  static int _kernLengthY;
  static int _kernLengthZ;

  static float _RatioX, _RatioY, _RatioZ; // size ratio of input and output;

   // kernel used for interpolation;
  static float _kernX[_MaxSize*_Oversample+1];
  static float _kernY[_MaxSize*_Oversample+1];
  static float _kernZ[_MaxSize*_Oversample+1];

  // Trilinear resampling (not requiring kern)
  static void trilinear(T const *volIn,int x1,int y1,int z1,T *volOut,int x2,int y2,int z2);

  //  FFT resampling ( not requiring kern )
  static void fft(T const *volIn,int x1,int y1,int z1,T *volOut,int x2,int y2,int z2);


  // initialize kernel with the function in the argument;
  static void initX(double (*func)(double));
  static void initY(double (*func)(double));
  static void initZ(double (*func)(double));
  static void initialize(ConvolutionType flag);

  //  Resampling with coefficient bins;
  static void coeffBin(T const *volIn,int x1,int y1,int z1,T *volOut,int x2,int y2,int z2);

  // interpolating functions used currently, more can be added here if necessary;
  static double linear( double );
  static double cubic( double );
  static double gauss( double );


 private:

  static const char *const _ClassRevision;
  static const char *const _ClassName;
  static const int         _ClassVersion;
  
};

template <class T> const char *const Convolution<T>::_ClassRevision = "$Id: Convolution.h,v 1.1 2004/11/15 04:44:07 joeh Exp $";
template <class T> const char *const Convolution<T>::_ClassName = "Convolution";
template <class T> const int Convolution<T>::_ClassVersion = 1; 


#endif // __CONVOLUTION_H__

