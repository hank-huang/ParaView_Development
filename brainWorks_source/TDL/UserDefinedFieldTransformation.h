#ifndef __USERDEFINEDFIELDTRANSFORMATION_H__
#define __USERDEFINEDFIELDTRANSFORMATION_H__
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
// File: UserDefinedFieldTransformation.h
//
// Author: Rob Teichman
//
// Purpose: User Defined Field Transformation class definition
//
// Note: If you are changing this class.  The version number of this class
// should change as more or less data is saved so that we may be able to
// load in old versions.
// 
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: UserDefinedFieldTransformation.h,v 1.1 2004/11/15 04:44:09 joeh Exp $
//
// Revision History
//
// $Log: UserDefinedFieldTransformation.h,v $
// Revision 1.1  2004/11/15 04:44:09  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:28:25  kem
// Initial revision
//
// Revision 1.4  1999/07/09 17:50:31  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.3  1998/12/18 18:15:21  RAZORDB
// Fix for 7.2.1 compiler Warning 1356
//
// Revision 1.2  1998/12/10 19:04:09  RAZORDB
// IDL/AW merge
//
// Revision 1.1  1998/02/13 22:33:04  rst
// Initial revision
//
//
///////////////////////////////////////////////////////////////////////////

#include <OS.h>
#include <string.h>
#include <ADL/Array1D.h>
#include <ADL/Array2D.h>
#include <ADL/Array3D.h>
#include <ADL/Matrix.h>
#include <TDL/FieldTransformation.h>

class UserDefinedFieldTransformation : public FieldTransformation
{

public:

  // Ctor and Dtor.
  UserDefinedFieldTransformation(); 
  ~UserDefinedFieldTransformation(); 
  
  // save and load are virtuals in Transformation
  int save( const char * filename );
  int load( const char * filename );

  // calculate
  int calculate();

  // initialize take in the h-Field that defines this transformation
  // Returns 0 upon success, and non-zero otherwise
  int initialize(Array3D<float> const &hField);

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
  void getHField(Array3D<float> *hf);
  void getHField(Array3D<float> *, double, double, double, 
	double, double, double);

  // user supplied h-field
  Array3D<float> mHField;

};

#endif // __USERDEFINEDFIELDTRANSFORMATION_H__
