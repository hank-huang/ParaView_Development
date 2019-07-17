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
// File: ITX_FFT.C 
//
// Author: Guoling Tong
//
// Purpose: ITX_FFT class that gives a uniform interface for
//          the SGI complib.sgimath fft transforms and fftw transforms.
// 
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: ITX_FFT.C,v 1.1 2004/11/15 04:44:08 joeh Exp $
//
// Revision History
//
// $Log: ITX_FFT.C,v $
// Revision 1.1  2004/11/15 04:44:08  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:28:52  kem
// Initial revision
//
// Revision 1.1  1999/07/09 18:28:14  RAZORDB
// Initial revision
//
///////////////////////////////////////////////////////////////////////////


#include <iostream.h>
#include <stdlib.h>
#include <TDL/ITX_FFT.h>


// constructor
ITX_FFT::ITX_FFT()
{
	_isInitialized = 0;
}
	
void ITX_FFT::initialize(int dimension, FFT_type tType, int xSize,int ySize,int zSize)
{

//Set up the variables
	_dim 	= dimension;
	_flag 	= tType;
	_nx 	= xSize;
	_ny	= ySize;
	_nz	= zSize;

  // Error checking for the arguments;

  if(dimension<=0 || dimension>=4)
    {
      cerr<<"Error:dimension should be between 1,2 or 3!"<<endl;
      exit(1);
    }
  if(tType != R2C && tType!= C2R && tType != C2C)
    {
      cerr<<"Error: FFT type should be R2C, C2R or C2C!"<<endl;
      exit(1);
    }
  if(xSize<=0)
    {
      cerr<<"Error: xSize should be a positive integer!"<<endl;
      exit(1);
    }
  if(ySize<=0)
    {
      cerr<<"Error: ySize should be a positive integer!"<<endl;
      exit(1);
    }
  if(zSize<=0)
    {
      cerr<<"Error: zSize should be a positive integer!"<<endl;
      exit(1);
    }

#ifdef SGI_FFT

  _rCoeff=NULL;
  _cCoeff=NULL;

#else
  
  _cPlanF1=_cPlanB1=NULL;
  _cPlanBn=_cPlanFn=NULL;
  _rPlanF1=_rPlanB1=NULL;
  _rPlanFn=_rPlanBn=NULL;

#endif

  _isInitialized = 1;

#ifdef SGI_FFT
  if(_dim==1)
    {
      if(_flag==R2C || _flag==C2R)
	{
	  _rCoeff=new FFT_real[_nx+15];
	  _rCoeff=scfft1dui(_nx,NULL); 
	}
      
      if(_flag==C2C)
	{
	  _cCoeff=new FFT_complex[_nx+15];
	  _cCoeff=cfft1di(_nx,NULL);
	}
    }
  
  if(_dim==2)
    {
      if(_flag==R2C || _flag==C2R)
	{
	  _rCoeff=new FFT_real[_nx+15+2*(_ny+15)];
	  _rCoeff=scfft2dui(_nx,_ny,NULL); 
	}
      
      if(_flag==C2C)
	{
	  _cCoeff=new FFT_complex[_nx+15+_ny+15];
	  _cCoeff=cfft2di(_nx,_ny,NULL); 
  	}
    }
  
  if(_dim==3)
    {
      if(_flag==R2C || _flag==C2R)
	{
	  _rCoeff=new FFT_real[_nx+15+2*(_ny+15)+2*(_nz+15)];
	  _rCoeff=scfft3dui(_nx,_ny,_nz,NULL); 
	}
      
      if(_flag==C2C)
	{
	  _cCoeff=new FFT_complex[_nx+15+_ny+15+_nz+15];
	  _cCoeff=cfft3di(_nx,_ny,_nz,NULL); 
	}
    }
   
#else
  
  if(_dim==1)
    {
      if(_flag==R2C || _flag==C2R)
	{
	  _rPlanF1=rfftw_create_plan(_nx,FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE);
	  _rPlanB1=rfftw_create_plan(_nx,FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE);
	}
      
      if(_flag==C2C)
	{
	  _cPlanF1=fftw_create_plan(_nx,FFTW_FORWARD,FFTW_ESTIMATE);
	  _cPlanB1=fftw_create_plan(_nx,FFTW_BACKWARD,FFTW_ESTIMATE); 
	}
    }
  
  if(_dim==2)
    {
      if(_flag==R2C || _flag==C2R)
	{
	  _rPlanFn=rfftw2d_create_plan(_ny,_nx,FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE | FFTW_IN_PLACE);
	  _rPlanBn=rfftw2d_create_plan(_ny,_nx,FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE | FFTW_IN_PLACE);
	}
      
      if(_flag==C2C)
	{
	  _cPlanFn=fftw2d_create_plan(_ny,_nx,FFTW_FORWARD,FFTW_ESTIMATE | FFTW_IN_PLACE); 
	  _cPlanBn=fftw2d_create_plan(_ny,_nx,FFTW_BACKWARD,FFTW_ESTIMATE | FFTW_IN_PLACE); 
	}
    }
  
  if(_dim==3)
    {
      if(_flag==R2C || _flag==C2R)
	{
	  _rPlanFn=rfftw3d_create_plan(_nz,_ny,_nx,FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE | FFTW_IN_PLACE);
	  _rPlanBn=rfftw3d_create_plan(_nz,_ny,_nx,FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE | FFTW_IN_PLACE);
	}
      
      if(_flag==C2C)
	{
	  _cPlanFn=fftw3d_create_plan(_nz,_ny,_nx,FFTW_FORWARD,FFTW_ESTIMATE | FFTW_IN_PLACE);
	  _cPlanBn=fftw3d_create_plan(_nz,_ny,_nx,FFTW_BACKWARD,FFTW_ESTIMATE | FFTW_IN_PLACE);
	}
    }

#endif

}


// destructor
ITX_FFT::~ITX_FFT()
{
#ifdef SGI_FFT

  if(_flag==R2C || _flag==C2R)
    {
      delete [] _rCoeff;
      _rCoeff=NULL;
    }
  
  if(_flag==C2C)
    {
      delete [] _cCoeff;
      _cCoeff=NULL;
    }
  
#else
  if(_dim==1)
    {
      if(_flag==R2C || _flag==C2R)
	{
	  rfftw_destroy_plan(_rPlanF1);
	  rfftw_destroy_plan(_rPlanB1);
	}
      
      if(_flag==C2C)
	{
	  fftw_destroy_plan(_cPlanF1);
	  fftw_destroy_plan(_cPlanB1);
	}
    }
  
  if(_dim==2 || _dim==3)
    {
      if(_flag==R2C || _flag==C2R)
	{
	  rfftwnd_destroy_plan(_rPlanFn);
	  rfftwnd_destroy_plan(_rPlanBn);
	}
      
      if(_flag==C2C)
	{
	  fftwnd_destroy_plan(_cPlanFn);
	  fftwnd_destroy_plan(_cPlanBn);
	}
    }
  
#endif

}


// complex to complex Forward transformation of 1,2,3d;
void ITX_FFT::fftForward(FFT_complex *in)
{
  int i,j;

	if (_isInitialized == 0) {
	cerr << "fftForward called before Initialization" << endl;
	exit(1);
	}

#ifdef SGI_FFT

  int nxPlusStride=2*((_nx+2)/2);
  
  if(_dim==1)
    {
      cfft1d(-1,_nx,in,1,_cCoeff);
    }
  if(_dim==2)
    {      
      cfft2d(-1,_nx,_ny,in,_nx,_cCoeff);
    }
  if(_dim==3)
    {
      cfft3d(-1,_nx,_ny,_nz,in,_nx,_ny,_cCoeff);
    }
  
#else
 
  if(_dim==1)
    {
      FFT_complex *out;
      out=new FFT_complex[_nx];
      
      fftw_one(_cPlanF1,in,out);
      for(int i=0;i<_nx;i++)
	in[i]=out[i];
      
      delete [] out; out=NULL;
    }
  if(_dim==2 || _dim==3)
    {      
      fftwnd_one(_cPlanFn,in,NULL);
    }
  
#endif
      
}


// complex to complex backward complex to real transformation of 1,2,3d;
void ITX_FFT::fftInverse(FFT_complex *in)
{
  int i,j;

	if (_isInitialized == 0) {
	cerr << "fftForward called before Initialization" << endl;
	exit(1);
	}

  if(_flag==C2C)
    {
#ifdef SGI_FFT

      int nxPlusStride=2*((_nx+2)/2);
      
      if(_dim==1)
	{
	  cfft1d(1,_nx,in,1,_cCoeff);
	}
      if(_dim==2)
	{      
	  cfft2d(1,_nx,_ny,in,_nx,_cCoeff);
	}
      if(_dim==3)
	{
	  cfft3d(1,_nx,_ny,_nz,in,_nx,_ny,_cCoeff);
	}
      
#else
 
      if(_dim==1)
	{
	  FFT_complex *out;
	  out=new FFT_complex[_nx];

	  fftw_one(_cPlanB1,in,out);
	  for(int i=0;i<_nx;i++)
	    in[i]=out[i];

	  delete [] out; out=NULL;
	}
      if(_dim==2 || _dim==3)
	{      
	  fftwnd_one(_cPlanBn,in,NULL);
	}
  
#endif
      
    }
 
  if(_flag == C2R || _flag == R2C )
    {
      
#ifdef SGI_FFT
      
      int nxPlusStride=2*((_nx+2)/2);
      
      if(_dim==1)
	{
	  FFT_real *temp=(FFT_real *)in;	  
	  csfft1du(1,_nx,temp,1,_rCoeff);	  
	}      
	  
      if(_dim==2)
	{ 
	  FFT_real *temp=(FFT_real *)in;
	  csfft2du(1,_nx,_ny,temp,nxPlusStride,_rCoeff);
	}
      if(_dim==3)
	{
	  FFT_real *temp=(FFT_real *)in;
	  csfft3du(1,_nx,_ny,_nz,temp,nxPlusStride,_ny,_rCoeff);
	}
      
#else
 
      if(_dim==1)
	{
	  FFT_real *out;
	  out=new FFT_real[_nx];
	  FFT_real *temp;
	  temp=new FFT_real[_nx];

	  FFT_real *temp2=(FFT_real *)in;
 
	  temp[0]=in[0].re;
	  if(_nx%2==0)
	    {
	      temp[_nx/2]=in[_nx/2].re;
	    }
	  else
	    {
	      temp[_nx/2]=in[_nx/2].re;
	      temp[_nx/2+1]=in[_nx/2].im;
	    }
	  for(i=1;i<_nx/2;i++)
	    {
	      temp[i]=in[i].re;
	      temp[_nx-i]=in[i].im;
	    }
	  
	  rfftw_one(_rPlanB1,temp,temp2);
	  
	  delete [] temp; temp=NULL;
	}      

      if(_dim==2 || _dim==3)
	{      
	  rfftwnd_one_complex_to_real(_rPlanBn,in,NULL);
	}
      
#endif

    }
}

// real to complex transformation of 1,2,3d;
void ITX_FFT::fftForward(FFT_real *in)
{

	if (_isInitialized == 0) {
	cerr << "fftForward called before Initialization" << endl;
	exit(1);
	}

  if(_flag==R2C || _flag==C2R)
    {
#ifdef SGI_FFT
            
      int nxPlusStride=2*((_nx+2)/2);
      
      if(_dim==1)
	{
	  scfft1du(-1,_nx,in,1,_rCoeff);
	}
      if(_dim==2)
	{      
	  scfft2du(-1,_nx,_ny,in,nxPlusStride,_rCoeff);
	}
      if(_dim==3)
	{
	  scfft3du(-1,_nx,_ny,_nz,in,nxPlusStride,_ny,_rCoeff);
	}
      
#else
      
      if(_dim==1)
	{
	  FFT_real *out;
	  out=new FFT_real[_nx];

	  rfftw_one(_rPlanF1,in,out);

	  in[0]=out[0];
	  in[1]=0;

	  if(_nx%2==0)
	    {
	      in[_nx]=out[_nx/2];
	    }
	  else
	    {
	      in[_nx-1]=out[_nx/2];
	      in[_nx]=out[_nx-_nx/2];
	    }
	  for(int i=1;i<_nx/2;i++)
	    {
	      in[2*i]=out[i];
	      in[2*i+1]=out[_nx-i];
	    }

	  delete [] out; out=NULL;

	}

      if(_dim==2 || _dim==3)
	{      
	  rfftwnd_one_real_to_complex(_rPlanFn,in,NULL);
	}
      
#endif

    }

}


