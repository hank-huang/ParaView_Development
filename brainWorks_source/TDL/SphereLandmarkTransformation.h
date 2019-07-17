#ifndef __SphereLandmarkTransformation_h_
#define __SphereLandmarkTransformation_h_
//////////////////////////////////////////////////////////////////////////
//
// IntellX, L.L.C., Proprietary
// Do not reproduce without permission in writing.
// Copyright (c) 1998 IntellX, L.L.C.
// All rights reserved
//
// File: FluidLandmarkTransformation.h
//
// Author: Muge Bakircioglu
//	   08/99
//
// SphereLandmarkTransformation class 
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
#include <TDL/Surface.h>
#include <TDL/SurfaceUtils.h>
#include <TDL/Point.h>
#include <TDL/Timer.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <math.h>

//#define NUMT	10
#define NUMT 3
class SphereLandmarkTransformation : public FieldTransformation {
  
public:

  // Ctor and Dtor.
  SphereLandmarkTransformation(); 
  ~SphereLandmarkTransformation(); 
    
  int  save(const char * filename );  // virtual in Transformation class
  int  load(const char * filename );  // virtual in Transformation class

  int initialize(const Matrix<double> &atlasPnts,
		 const Matrix<double> &patientPnts,
		 Surface &atlasSurf, 
		 Point &patientNP, Point &atlasNP);
	  
  int calculate() { Surface out; calculate(out); return 0; }  
  int calculate(Surface &outSurf);      // virtual in Transformation class

  void print();

  // virtual in FieldTransformation class
       void getHField(Array3D<float>*) { return; } 
  void getHField(Array3D<float> *, double, double, double,
		   double, double, double) { return ; }

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
  
  double _sigma;

  double _beta;
  double _numiter;
  double _numt;

  // Private member variables
  Matrix<double> _atlasPoints;
  Matrix<double> _patientPoints;
  Surface _atlasSurf;
  Matrix<double> K1,K2,mat;
  int _numPoints;
  int _numVert;
  Matrix<double> Rot;
  double targRad;  
  
   Array2D<double> pathx, pathy, pathz;

  Array1D<double> v1x,v1y,v1z;
  Array1D<double> v2x,v2y,v2z;
  Array2D<double> v3x,v3y,v3z;
  Array1D<double> wtu,wtv;

  int findPath();
  void calculateWeights(double *x1, double *y1, double *z1,
			double *x2, double *y2, double *z2);

};

#endif
