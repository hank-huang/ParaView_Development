#ifndef __REMOVESURFACETRANSFORMATION_H__
#define __REMOVESURFACETRANSFORMATION_H__
///////////////////////////////////////////////////////////////////////////
//
// IntellX, L.L.C., Proprietary
// The contents of this file comprise confidential, proprietary and/or
// trade secret information which is the property of IntellX LLC.
// Dissemination, distribution or other disclosure of this information
// is strictly prohibited without the express permission of IntellX LLC.
// Copyright (c) 1998 IntellX, L.L.C.
// All rights reserved
//
// File: RemoveSurfaceTransformation.h
//
// Author: Rob Teichman
//  algorithm from Sarang Joshi
//
// Purpose: Remove Surface Transformation class definition
//
// Note: If you are changing this class.  The version number of this class
// should change as more or less data is saved so that we may be able to
// load in old versions.
// 
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: RemoveSurfaceTransformation.h,v 1.1 2004/11/15 04:44:08 joeh Exp $
//
// Revision History
//
// $Log: RemoveSurfaceTransformation.h,v $
// Revision 1.1  2004/11/15 04:44:08  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:28:30  kem
// Initial revision
//
// Revision 1.4  1999/07/09 17:49:51  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.3  1998/12/18 18:14:48  RAZORDB
// Fix for 7.2.1 compiler Warning 1356
//
// Revision 1.2  1998/12/10 19:03:51  RAZORDB
// IDL/AW merge
//
// Revision 1.1  1998/05/18 19:25:56  rst
// Initial revision
//
//
///////////////////////////////////////////////////////////////////////////

#include <OS.h>
#include <string.h>
#include <ADL/Array3D.h>
#include <ADL/GreyVolume.h>
#include <TDL/Surface.h>
#include <TDL/Point.h>
#include <TDL/FieldTransformation.h>

class RemoveSurfaceTransformation : public FieldTransformation
{

public:

  // Ctor and Dtor.
  RemoveSurfaceTransformation(); 
  ~RemoveSurfaceTransformation(); 
  
  // save and load are virtuals in Transformation
  int save( const char * filename );
  int load( const char * filename );

  // calculate
  int calculate();

  // initialize
  // input a surface to remove and a bone mask (may be empty,
  // but must be sized properly)
  // Returns 0 upon success, and non-zero otherwise
  int initialize(Surface const &surf,
		 GreyVolume<unsigned char> const &boneMask,
		 double seedx = 0.0,
		 double seedy = 0.0,
		 double seedz = 0.0);

  // set the seed (the point the tumor will be shrunk into)
  // if not set explicitly, the algorithm will use the centriod
  // of the surface
  int setSeed(Point const &seed);
  int setSeed(double sx, double sy, double sz)
    { Point temp(sx, sy, sz); return setSeed(temp); }

  void print();


protected:

  // Returns 0 upon success, and non-zero otherwise
  int saveClassParameters ( ofstream & );
  int loadClassParameters ( ifstream & );
  ostream & printClassParameters ( ostream & );


private:
  
  static const char * const mClassRevision;
  static const char * const mClassName;
  static const int          mClassVersion;

  // virtual in FieldTransformation class
  void getHField(Array3D<float> *hf)
	{ getHField(hf,1.0,1.0,1.0,0.0,0.0,0.0); }

  void getHField(Array3D<float> *, double, double, double, 
	double, double, double);

  // h-field
  Array3D<float> mHField;

  // seed: 0,0,0 implies use centriod
  Point mSeed;

  // initial surface
  Surface mSurf;

  // bone mask volume
  GreyVolume<unsigned char> mBoneMask;

};
#endif // __REMOVESURFACETRANSFORMATION_H__
