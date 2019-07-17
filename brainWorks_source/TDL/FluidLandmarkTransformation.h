#ifndef __FluidLandmarkTransformation_h_
#define __FluidLandmarkTransformation_h_
//////////////////////////////////////////////////////////////////////////
//
// IntellX, L.L.C., Proprietary
// Do not reproduce without permission in writing.
// Copyright (c) 1998 IntellX, L.L.C.
// All rights reserved
//
// File: FluidLandmarkTransformation.h
//
// Author: Keith Doolittle
//	   11/98
//
// FluidLandmarkTransformation class 
//
//////////////////////////////////////////////////////////////////////////
#include <OS.h>
#include <string.h>
#include <ADL/Array1D.h>
#include <ADL/Array2D.h>
#include <ADL/Array3D.h>
#include <ADL/Matrix.h>
#include <TDL/AffineLandmarkTransformation.h>
#include <TDL/FieldTransformation.h>

#define NUMT	5

class FluidLandmarkTransformation : public FieldTransformation {
  friend class FluidFluidTransformation;
public:

  // Ctor and Dtor.
  FluidLandmarkTransformation(); 
  ~FluidLandmarkTransformation(); 
    
  int  save(const char * filename );  // virtual in Transformation class
  int  load(const char * filename );  // virtual in Transformation class

  int initialize(const Matrix<double> &atlasPoints,
		 const Matrix<double> &patientPoints,
		int AtlasSizeX,  int AtlasSizeY,  int AtlasSizeZ,
		double AtlasScaleX, double AtlasScaleY, double AtlasScaleZ,
		int PatientSizeX,  int PatientSizeY,  int PatientSizeZ,
		double PatientScaleX, double PatientScaleY, double PatientScaleZ
		);
    
  int calculate();             // virtual in Transformation class

  void print();

  // virtual in FieldTransformation class
  void getHField(Array3D<float>*);
  void getHField(Array3D<float> *, double, double, double,
        double, double, double);


  int setAlpha(double a)	{ _alpha = a; return 0; }
  int setSigma(double a)	{ _sigma = a; return 0; }
  int setBeta(double a)		{ _beta = a; return 0; }
  int setNumberIterations(int niter)
				{ if(niter>0) { _numiter = niter; return 0; }
				  else return 1;
				}

protected:

  // Returns 0 upon success, and non-zero otherwise
  int saveClassParameters ( ofstream & );
  int loadClassParameters ( ifstream & );
  ostream & printClassParameters ( ostream & );

private: 
  static const char * const mClassRevision;
  static const char * const mClassName;
  static const int          mClassVersion;
  
  double _alpha;
  double _sigma;
  double _beta;
  double _numiter;
  double _numt;

  // Private member variables
  Matrix<double> _atlasPoints;
  Matrix<double> _patientPoints;
  Matrix<double> K1,K2,mat;
  int _numPoints;

  int _nx,_ny,_nz;

  Array3D<float> _TotalHField;
  Array2D<double> pathx,pathy,pathz;
  AffineLandmarkTransformation _mAffine;

  Array1D<double> v1x,v1y,v1z;
  Array1D<double> v2x,v2y,v2z;
  Array2D<double> v3x,v3y,v3z;
  Array1D<double> wtx,wty,wtz;

  int findPath();
  void calculateWeights(double *x1, double *y1, double *z1,
			double *x2, double *y2, double *z2);

};

#endif
