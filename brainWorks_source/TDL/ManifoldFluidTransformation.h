#ifndef __MANIFOLDFLUIDTRANSFORMATION_H__
#define __MANIFOLDFLUIDTRANSFORMATION_H__

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
// File:  ManifoldFluidTransformation.h
//
// Author:  Kevin E. Mark
//
// Purpose:  This class computes a manifold transformation and a composed
//           fluid transformation on a subvolume.
//
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: ManifoldFluidTransformation.h,v 1.1 2004/11/15 04:44:08 joeh Exp $
//
// Revision History
//
// $Log: ManifoldFluidTransformation.h,v $
// Revision 1.1  2004/11/15 04:44:08  joeh
// first check in of BrainWorks
//
// Revision 1.2  1999/11/09 21:28:54  kem
// Add methods for computing h field values given a point
//
// Revision 1.1  1999/10/01 15:28:22  kem
// Initial revision
//
// Revision 1.12  1999/09/28 19:35:58  RAZORDB
// Add flag in initialize to specify subvolume in atlas or patient
//
// Revision 1.11  1999/09/28 15:31:42  RAZORDB
// Add matchIntensities() method
//
// Revision 1.10  1999/07/09 17:49:34  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.9  1999/06/03 16:37:10  RAZORDB
// Add adjustIntensities() method
//
// Revision 1.8  1998/12/18 18:14:35  RAZORDB
// Fix for 7.2.1 compiler Warning 1356
//
// Revision 1.7  1998/12/10 19:03:42  RAZORDB
// IDL/AW merge
//
// Revision 1.6  1998/11/03 21:29:48  RAZORDB
// Add support for mask volume
//
// Revision 1.5  1998/10/14 23:22:33  kem
// Add setInitNumBasis() for elastic mode
//
// Revision 1.4  1998/06/02 15:24:28  rst
// Changed fluid to support elastic
//
// Revision 1.3  1998/04/17 14:53:04  kem
// Move trilinearHField() to FieldTransformation
//
// Revision 1.2  1998/04/14 16:26:08  kem
// Add methods to find good factorizations for the fft's
//
// Revision 1.1  1998/04/09 19:09:20  kem
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
#include <TDL/ManifoldTransformation.h>

class ManifoldFluidTransformation : public FieldTransformation {

public:
  ManifoldFluidTransformation(); 
  ~ManifoldFluidTransformation(); 
    
  int initialize(GreyVolume<unsigned char> const &atlas, 
		 GreyVolume<unsigned char> const &patient, 
		 Array2D<double> const &atlasLandmarks,   // in voxels
		 Array2D<double> const &patientLandmarks, // in voxels
		 Array1D<double> const &landmarkVariance,
		 double atlasSubvolumeOffsetX, // in mm
		 double atlasSubvolumeOffsetY, // in mm
		 double atlasSubvolumeOffsetZ, // in mm
		 int atlasSubvolumeSizeX,      // in voxels
		 int atlasSubvolumeSizeY,      // in voxels
		 int atlasSubvolumeSizeZ,      // in voxels
		 double atlasSubvolumeResX,    // in mm/voxel
		 double atlasSubvolumeResY,    // in mm/voxel
		 double atlasSubvolumeResZ,    // in mm/voxel
		 int isSubvolumeInAtlas = 1,   // if 0 then subvolume is in patient 
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

  void getHField(const float patientInVoxelsX,
		 const float patientInVoxelsY,
		 const float patientInVoxelsZ,
		 float &atlasInVoxelsX,
		 float &atlasInVoxelsY,
		 float &atlasInVoxelsZ);

  void computePatientSubvolume(
    const ManifoldTransformation &manifoldTransform,
    const double &atlasSubOffsetX, 
    const double &atlasSubOffsetY, 
    const double &atlasSubOffsetZ,
    const int &atlasSubSizeX, 
    const int &atlasSubSizeY, 
    const int &atlasSubSizeZ,
    const double &atlasSubResX, 
    const double &atlasSubResY, 
    const double &atlasSubResZ, 
    const int &patientSizeX,
    const int &patientSizeY,
    const int &patientSizeZ,
    const double &patientResX,
    const double &patientResY,
    const double &patientResZ,
    int &patientSubOffsetInVoxelsX,
    int &patientSubOffsetInVoxelsY,
    int &patientSubOffsetInVoxelsZ,
    int &patientSubSizeX,
    int &patientSubSizeY,
    int &patientSubSizeZ);

  void extractSubvolume(const Array3D<unsigned char> &volume,
			const int subvolumeOffsetInVoxelsX,
			const int subvolumeOffsetInVoxelsY,
			const int subvolumeOffsetInVoxelsZ,
			const double magFactorX,
			const double magFactorY,
			const double magFactorZ,
			Array3D<unsigned char> &subvolume);

  void adjustIntensities(Array3D<unsigned char> &volume);
  void matchIntensities(Array3D<unsigned char> &vol1,
			Array3D<unsigned char> &vol2);

  int invApply(const Surface &inSurface, // should override the method in
	       Surface *outSurface,      // FieldTransformation
	       double Rx, double Ry, double Rz,
	       double Ox, double Oy, double Oz);



  static const char * const mClassRevision;
  static const char * const mClassName;
  static const int mClassVersion;

  ManifoldTransformation _manifoldTransform;
  FluidTransformation _fluidTransform;

  const Array3D<unsigned char> *_atlasPtr;
  const Array3D<unsigned char> *_patientPtr;
  const Array3D<unsigned char> *_maskPtr;
  const Array2D<double> *_atlasLandmarksPtr;
  const Array2D<double> *_patientLandmarksPtr;
  const Array1D<double> *_landmarkVariancePtr;
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
  int _isSubvolumeInAtlas;    // 0 indicates input subvolume parameters in
                              // patient coordinates
};

#endif // __MANIFOLDFLUIDTRANSFORMATION_H__
