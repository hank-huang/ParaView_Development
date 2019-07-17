#ifndef __FIELDTRANSFORMATION_H__
#define __FIELDTRANSFORMATION_H__
///////////////////////////////////////////////////////////////////////////
//
// IntellX, L.L.C., Proprietary
// Do not reproduce without permission in writing.
// Copyright (c) 1997 IntellX, L.L.C.
// All rights reserved
//
// File: FieldTransformation.h
//
// Author: IntellX
//
// Purpose: This is the base class ofr derived transformation classes that
//  use a deformation field as the parameter that describes the transformation
//
// Note: If you are changing this class.  The version number of this class
// should change as more or less data is saved so that we may be able to
// load in old versions.
//
// To Do:
//   1) deal with failed call to loadTransformationClassParameters
//
// Things to think about...
//   1) Maybe put in code that can read in older version of the class
// 
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: FieldTransformation.h,v 1.1 2004/11/15 04:44:07 joeh Exp $
//
// Revision History
//
// $Log: FieldTransformation.h,v $
// Revision 1.1  2004/11/15 04:44:07  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:28:14  kem
// Initial revision
//
// Revision 1.19  1999/07/09 17:48:48  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.18  1998/12/18 18:14:21  RAZORDB
// Fix for 7.2.1 compiler Warning 1356
//
// Revision 1.17  1998/12/10 19:03:32  RAZORDB
// IDL/AW merge
//
// Revision 1.16  1998/05/06 19:08:31  rst
// Change invApply prototype
//
// Revision 1.15  1998/04/17 14:52:31  kem
// Add inverse apply for a point in patient space
//
// Revision 1.14  1998/04/09 19:07:06  kem
// Change apply() for volumes
//
// Revision 1.13  1998/03/05 18:04:27  rst
// new applys
//
// Revision 1.12  1997/12/12 21:50:51  csb
// Merge 1.7.1.4 and 1.11
//
// Revision 1.11  1997/12/12 21:45:29  csb
// Change class version to 3 because added parameters (resolutions)
//
// Revision 1.10  1997/12/12 21:22:35  csb
// Reverted from version 1.8
//
// Revision 1.8  1997/12/12 19:20:30  abed
// Incorporating resolution changes
//
// Revision 1.7  1997/10/02 20:37:44  rst
// Fixed apply() routines for Surfaces and Points
//
// Revision 1.6  1997/09/16 21:10:48  kem
// Change atlas and patient dimension variables to protected
//
// Revision 1.5  1997/09/16 20:52:50  kem
// Add atlas and patient dimension variables
//
// Revision 1.4  1997/09/16 19:31:28  kem
// Remove getUField() method
//
// Revision 1.3  1997/08/25 17:39:20  kem
// Add default value for getHField()
//
// Revision 1.2  1997/08/22 20:55:20  rst
// Add apply() for Surfaces and Points
//
// Revision 1.1  1997/08/01 19:49:19  csb
// Initial revision
//
// Revision 1.5  1997/08/01 15:38:48  csb
// Incorportate Abed's, fluid, elastic and field
//
// Revision 1.4  1997/06/25 17:35:31  csb
// Split into .h and .c
//
// Revision 1.3  1997/06/24 22:24:31  csb
// Added Revision, changed order of saving and loading
//
// Revision 1.2  1997/06/02 21:58:44  csb
// Changes to protected functions
//
// Revision 1.1  1997/05/30 23:59:15  csb
// Initial revision
//
//////////////////////////////////////////////////////////////////////////
#include <OS.h>
#include <stream.h>
#include <ADL/Array3D.h>
#include <ADL/Pnt.h>
#include <TDL/Transformation.h>
#include <TDL/Point.h>
#include <TDL/Surface.h>

class FieldTransformation : public Transformation {

public:

  FieldTransformation() { 
    _haveHField = 0; 
    _atlasSizeX = 0;
    _atlasSizeY = 0;
    _atlasSizeZ = 0;
    _patientSizeX = 0;
    _patientSizeY = 0;
    _patientSizeZ = 0;
    _atlasResX = 1.0; 
    _atlasResY = 1.0; 
    _atlasResZ = 1.0; 
    _patientResX = 1.0; 
    _patientResY = 1.0; 
    _patientResZ = 1.0; 
  }

  virtual ~FieldTransformation() { /* empty */ } // virtual destructor

  //
  // Define the function apply() here as it the same for all derived classes
  //

  // apply() for volumes
/*   int apply(const Array3D<unsigned char> &inVolume, */
/*             Array3D<unsigned char> *outVolume, */
/*             const int interpolateFlag = 1); */

  // apply() for surfaces
  int apply (const Surface &inSurface,
	     Surface *outSurface);

  // apply() for surfaces with resolution
  int apply(const Surface &inSurface,
	    Surface *outSurface,
	    double Rx, double Ry, double Rz,
	    double Ox, double Oy, double Oz);

  // apply() for points
  int apply (const Point &inPoint,
	     Point *outPoint);

  int apply(const Array3D<unsigned char> &inVolume, 
		Array3D<unsigned char> *outVolume, 
		double Rx = 1.0, double Ry = 1.0, double Rz = 1.0, 
		double sOx = 0.0, double sOy = 0.0, double sOz = 0.0, 
		const int interpolateFlag = 1); 

#if 0
  int apply(const Array3D<unsigned char> &inVolume, 
	Array3D<unsigned char> *outVolume, 

	// Output resolution.
	double Rx, double Ry, double Rz, 

	// Input volume resolution.
	double Rtx, double Rty, double Rtz, 

	// input Seg. resolution.
	double Rgx, double Rgy, double Rgz, 

	// Input Seg. Offset wrt Input Volume
	double gOx, double gOy, double gOz, 

	// output cube offset w.r.t Study.
	double sOx, double sOy, double sOz, 

	const int interpolateFlag); 
#else
  // The right way to apply to volumes
  int apply(const Array3D<unsigned char> &inputVolume, 
	    double inputResX,     // inputVolume resolution in mm/voxel
	    double inputResY,
	    double inputResZ,
	    double inputOffsetX,  // inputVolume offset relative to atlas volume
	    double inputOffsetY,  // in mm
	    double inputOffsetZ,
	    Array3D<unsigned char> *outputVolume,
	    double outputResX,    // outputVolume resolution in mm/voxel
	    double outputResY,
	    double outputResZ,
	    double outputOffsetX, // outputVolume offset relative to patient volume
	    double outputOffsetY, // in mm
	    double outputOffsetZ,
	    const int interpolateFlag); 
#endif

  // function to get a copy of the h-field
  // (will generate if necessary)
  Array3D<float> gimmeHField();

  int setHField( Array3D<float> & hField);

  // apply inverse, to drag segmentation from a patient into the atlas
  int invApply(const Array3D<unsigned char> &inVolume,
	       Array3D<unsigned char> *outVolume,
	       double Rx, double Ry, double Rz, 
	       double sOx, double sOy, double sOz);
//	       double Rox, double Roy, double Roz,
//	       double Oox, double Ooy, double Ooz);

  // inverse for surfaces (from patient back into atlas)
  int invApply(const Surface &inSurface,
	       Surface *outSurface,
	       double Rx, double Ry, double Rz,
	       double Ox, double Oy, double Oz);

  // inverse apply for a point
  // takes a point in patient voxel coordinates and maps it to
  // its corresponding atlas voxel coordinate
  int invApply(const Pnt &patientVoxelPoint, Pnt &atlasVoxelPoint);

protected:

  void trilinearHField(const Array3D<float> &hField,
		       float x, 
		       float y,
		       float z,
		       float *Hx, 
		       float *Hy, 
		       float *Hz);

  // This should be overridden by the derived class that knows how to get
  // an h-field
  virtual void getHField(Array3D<float> *) = 0;

  virtual void getHField(Array3D<float> *, 
	double Rx, double Ry, double Rz, 
	double Ox, double Oy, double Oz) = 0; 

  int saveClassParameters( ofstream & );
  int loadClassParameters( ifstream & );
  ostream & printClassParameters ( ostream & );

  int _atlasSizeX;
  int _atlasSizeY;
  int _atlasSizeZ;
  int _patientSizeX;
  int _patientSizeY;
  int _patientSizeZ;

  double _atlasResX; 
  double _atlasResY; 
  double _atlasResZ; 
  double _patientResX; 
  double _patientResY; 
  double _patientResZ; 

private:
  static const char * const mClassRevision;
  static const char * const mClassName;
  static const int          mClassVersion;
  
  int _haveHField;
  Array3D<float> _hField;
};

#endif // __FIELDTRANSFORMATION_H__
