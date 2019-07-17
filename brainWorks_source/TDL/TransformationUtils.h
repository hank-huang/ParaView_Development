#ifndef __TRANSFORMATION_UTILS_H__
#define __TRANSFORMATION_UTILS_H__

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
// File: TransformationUtils.h
//
// Author: csb
//
// Purpose: provide routines that all transformation classes share
// 
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: TransformationUtils.h,v 1.1 2004/11/15 04:44:09 joeh Exp $
//
// Revision History
//
// $Log: TransformationUtils.h,v $
// Revision 1.1  2004/11/15 04:44:09  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:27:51  kem
// Initial revision
//
// Revision 1.5  1998/10/29 00:19:30  RAZORDB
// Add flag in interpolate() for non-h fields
//
// Revision 1.4  1998/06/01 15:38:29  csb
// add compactdisplay() function
//
// Revision 1.3  1998/04/08 20:51:48  kem
// Add loadLandmarks() method
//
// Revision 1.2  1998/01/17 19:07:04  csb
// Add loadvolume() and resample for shorts
//
// Revision 1.1  1997/12/12 15:59:03  csb
// Initial revision
//
//
///////////////////////////////////////////////////////////////////////////
#include <OS.h>
#include <ADL/Array3D.h>
#include <ADL/GreyVolume.h>

class TransformationUtils {

public:
  static int
  interpolate( const Array3D<unsigned char> & input,
                     Array3D<unsigned char> & output,
               const Array3D<float>         & hfield);

  static int
  interpolate( float * field_at_xyz_Array,
	       const Array3D<float> & fieldInterlaced, // nz * ny * nx*3
	       double z,
	       double y,
	       double x,
	       int isHFieldFlag = 1);

  static int
  upsampleInterlacedField( const Array3D<float> & downField,
                                 Array3D<float> & upField);

  static int
  readTransformationType( const char * filename,
                                char * type);

  //////////////////////////////////////
  // Trilinear interpolation functions  
  //////////////////////////////////////
  static unsigned char
  trilinear( const Array3D<unsigned char> &X,
             double z,
             double y,
             double x,
             unsigned char bkgrnd );

  static float
  trilinear( const Array3D<float> &X,
             double z,
             double y,
             double x,
             float bkgrnd );

  static void
  trilinear( Array1D<double> *vx,
             const Array3D<float> &X,
             const Array3D<float> &Y, 
             const Array3D<float> &Z,
             double z, 
             double y,
             double x,
             float bkgrnd );

  //////////////////////////
  // Gradient functions
  //////////////////////////
  static void
  gradient(       Array1D<double> * GradientPtr, 
            const Array3D<float>  & Function,
                  int x, int y, int z, 
                  double dx, double dy, double dz);

  static void 
  gradient(        Array1D<double>        * GradientPtr, 
            const  Array3D<unsigned char> & Function,
                   int x, int y, int z, 
                   double dx, double dy, double dz);



  //////////////////////////
  // Resample functions
  //////////////////////////
  static void
  resample(unsigned char const *volin,
           int x1, 
           int y1,
           int z1,
           unsigned char *volout,
           int x2,
           int y2,
           int z2);

  static void
  resample(unsigned short const *volin,
           int x1, 
           int y1,
           int z1,
           unsigned short *volout,
           int x2,
           int y2,
           int z2);

  /////////////////////////////
  // Volume loading functions
  /////////////////////////////
  static int
  loadVolume(const char *filename,
             GreyVolume<unsigned short> &volume);

  static int
  loadVolume(const char *filename,
             GreyVolume<unsigned char> &volume);

  static int loadLandmarks(const char *filename,
			   Array2D<double> &landmarks);

   void
   compactdisplay(Vector<double> vec);
    

};

#endif // __TRANSFORMATION_UTILS_H__
