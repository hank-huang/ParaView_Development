#ifndef __FLUIDFLUIDTRANSFORMATION_H__
#define __FLUIDFLUIDTRANSFORMATION_H__

///////////////////////////////////////////////////////////////////////////
//
// IntellX, L.L.C., Proprietary
// The contents of this file comprise confidential, proprietary and/or
// trade secret information which is the property of IntellX LLC.
// Dissemination, distribution or other disclosure of this information is
// strictly prohibited without the express permission of IntellX LLC.
// Copyright (c) 1997 IntellX, L.L.C.
// All rights reserved
//
// File:  FluidFluidTransformation.h
//
// Author:  Kevin E. Mark
//
// Purpose:  This class computes a fluid landmark transformation and a composed
//           fluid transformation on a subvolume.
//
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: FluidFluidTransformation.h,v 1.1 2004/11/15 04:44:08 joeh Exp $
//
// Revision History
//
// $Log: FluidFluidTransformation.h,v $
// Revision 1.1  2004/11/15 04:44:08  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:28:24  kem
// Initial revision
//
// Revision 1.2  1999/07/09 17:48:58  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.1  1999/01/08 17:29:11  RAZORDB
// Initial revision
//
//
///////////////////////////////////////////////////////////////////////////
#include <OS.h>
#include <string.h>
#include <ADL/Array1D.h>
#include <ADL/Array2D.h>
#include <ADL/Array3D.h>
#include <ADL/GreyVolume.h>
#include <ADL/Matrix.h>
#include <TDL/FieldTransformation.h>
#include <TDL/FluidTransformation.h>
#include <TDL/FluidLandmarkTransformation.h>

class FluidFluidTransformation : public FieldTransformation {

public:
  FluidFluidTransformation(); 
  ~FluidFluidTransformation(); 
    
  int initialize(GreyVolume<unsigned char> const &atlas, 
		 GreyVolume<unsigned char> const &patient, 
		 Matrix<double> const &atlasLandmarks,   // in voxels
		 Matrix<double> const &patientLandmarks, // in voxels
		 double atlasSubvolumeOffsetX, // in mm
		 double atlasSubvolumeOffsetY, // in mm
		 double atlasSubvolumeOffsetZ, // in mm
		 int atlasSubvolumeSizeX,      // in voxels
		 int atlasSubvolumeSizeY,      // in voxels
		 int atlasSubvolumeSizeZ,      // in voxels
		 double atlasSubvolumeResX,    // in mm/voxel
		 double atlasSubvolumeResY,    // in mm/voxel
		 double atlasSubvolumeResZ,   // in mm/voxel
		 const Array3D<unsigned char> *maskPtr = NULL);

  int calculate();                    // virtual in Transformation class

  int  save( const char * filename ); // virtual in Transformation class
  int  load( const char * filename ); // virtual in Transformation class
  void print();

  int setNumberIterations( int numOuter, int numInner = 1 );
  int setLaplacianWeight( double weight );
  int setGradWeight( double weight );
  int setAdditiveWeight( double weight );
  int setMaxPerturbation( double pmax ); 
  int setMagicStoppingNumber(double x); 
  int setCalculateTotalJacobian( int flag );
  int setInitialTransformation( const char * filename );
  int setInitNumBasis(int initNumBasis);

protected:
  // Returns 0 upon success, and non-zero otherwise
  int saveClassParameters ( ofstream & );
  int loadClassParameters ( ifstream & );
  ostream& printClassParameters ( ostream & );
  
private: 

  int findGoodFactorization(int x);
  int isDivisibleBy235(int x);

  void getHField(Array3D<float> *hf)
	{ getHField(hf,1.0,1.0,1.0,0.0,0.0,0.0); }

  void getHField(Array3D<float> *hField,  // virtual in FieldTransformation class 
		 double resolutionX, 
		 double resolutionY, 
		 double resolutionZ, 
		 double offsetX, 
		 double offsetY, 
		 double offsetZ);

  void extractSubvolume(const Array3D<unsigned char> &volume,
			const int subvolumeOffsetInVoxelsX,
			const int subvolumeOffsetInVoxelsY,
			const int subvolumeOffsetInVoxelsZ,
			const double magFactorX,
			const double magFactorY,
			const double magFactorZ,
			Array3D<unsigned char> &subvolume);
  static const char * const mClassRevision;
  static const char * const mClassName;
  static const int          mClassVersion;

  FluidLandmarkTransformation _fluidLandmarkTransform;
  FluidTransformation _fluidTransform;

  const Array3D<unsigned char> *_atlasPtr;
  const Array3D<unsigned char> *_patientPtr;
  const Array3D<unsigned char> *_maskPtr;
  const Matrix<double> *_atlasLandmarksPtr;
  const Matrix<double> *_patientLandmarksPtr;
  const Array3D<unsigned char> *_atlasSubvolumePtr;
  int _atlasSubSizeX;
  int _atlasSubSizeY;
  int _atlasSubSizeZ;
  double _atlasSubResX;       // atlas subvolume resolution in mm/voxel
  double _atlasSubResY;
  double _atlasSubResZ;
  double _atlasSubOffsetX;    // atlas subvolume offset relative to
  double _atlasSubOffsetY;    // full atlas volume in mm
  double _atlasSubOffsetZ;
  int _patientSubSizeX;
  int _patientSubSizeY;
  int _patientSubSizeZ;
  double _patientSubResX;     // patient subvolume resolution in mm/voxel
  double _patientSubResY;
  double _patientSubResZ;
  double _patientSubOffsetX;  // patient subvolume offset relative to 
  double _patientSubOffsetY;  // full patient volume in mm
  double _patientSubOffsetZ;
};

#endif // __FLUIDFLUIDTRANSFORMATION_H__
