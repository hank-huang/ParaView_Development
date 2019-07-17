///////////////////////////////////////////////////////////////////////////
//
// IntellX, L.L.C., Proprietary
// The contents of this file comprise confidential, proprietary and/or
// trade secret information which is the property of IntellX LLC.
// Dissemination, distribution or other disclosure of this information
// is strictly prohibited without the express permission of IntellX LLC.
// Copyright (c) 1999 IntellX, L.L.C.
// All rights reserved
//
// File: ITX_FFT.h
//
// Author: Guoling Tong
//
// Purpose: Header file for the class FFT.C which gives a uniform interface for
//          the SGI complib.sgimath fft transforms and fftw transforms.
// 
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: ITX_FFT.h,v 1.1 2004/11/15 04:44:08 joeh Exp $
//
// Revision History
//
// $Log: ITX_FFT.h,v $
// Revision 1.1  2004/11/15 04:44:08  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:27:46  kem
// Initial revision
//
// Revision 1.1  1999/07/09 18:28:27  RAZORDB
// Initial revision
//
///////////////////////////////////////////////////////////////////////////
#include <OS.h>
#include <iostream.h>

#ifndef ITX_FFT_H
#define ITX_FFT_H

// if compiled on SGI machine, use SGI_fft library.
#ifdef SGI_FFT

#include <fft.h>

typedef float FFT_real;
typedef complex FFT_complex;

//otherwise, use fftw library.
#else

#define FFTW_ENABLE_FLOAT

#include <fftw.h>
#include <rfftw.h>

//typedef float FFT_real;
//typedef complex FFT_complex;
typedef fftw_real FFT_real;
typedef fftw_complex FFT_complex;

#endif


class ITX_FFT 
{
 public:

  enum FFT_type { R2C, C2R, C2C };
  
  // initialize the size and plan of the transformation;
   ITX_FFT();	
  ~ITX_FFT();  // destroy the plan;

  void initialize(int dimension,FFT_type tType,int xSize,int ySize=1,int zSize=1);
  // complex to complex transformation of 1,2,3 dimension and complex to real transform of 1, 2,3d;
  void fftForward(FFT_complex *in);
  void fftInverse(FFT_complex *in);

  // real to complex transformation of 1,2,3d;
  void fftForward(FFT_real *in);

  void fftInverse(FFT_real *in){
    cerr << "Inverse transform of real data not implemented!" << endl; }


 private:


  int _isInitialized;
  int _dim;  // dimension of the array;
  FFT_type _flag;   // type of the transform, i.e. R2C,etc.
  int _nx,_ny,_nz; // size of the transformation;

#ifdef SGI_FFT
  FFT_real *_rCoeff; // plans for SGI_fft library;
  FFT_complex *_cCoeff;
#else
  // plans for fftw library
  fftw_plan _cPlanF1,_cPlanB1;   // plans for 1D C2C Forward and Backward;
  fftwnd_plan _cPlanFn,_cPlanBn; // plans for nD C2C Forward and Backward; 
  rfftw_plan _rPlanF1,_rPlanB1;  // plans for 1D R2C and 1D C2R;
  rfftwnd_plan _rPlanFn,_rPlanBn;// plans for nD R2C and nD C2R;
#endif  

};

#endif
