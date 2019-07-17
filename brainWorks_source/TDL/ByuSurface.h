#ifndef __BYUSURFACE_H__
#define __BYUSURFACE_H__

///////////////////////////////////////////////////////////////////////////
//
// IntellX, L.L.C., Proprietary
// Do not reproduce without permission in writing.
// Copyright (c) 1997 IntellX, L.L.C.
// All rights reserved
//
// File: ByuSurface.h
//
// Author: Rob Teichman
//
// Purpose: ByuSurface Class Definition
//   This class is a derrived class from the Surface base class.
// It provides routines for input and output of BYU format surfaces.
// 
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: ByuSurface.h,v 1.1 2004/11/15 04:44:07 joeh Exp $
//
// Revision History
//
// $Log: ByuSurface.h,v $
// Revision 1.1  2004/11/15 04:44:07  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:27:41  kem
// Initial revision
//
// Revision 1.2  1999/09/21 21:00:09  RAZORDB
// TDL update
//
// Revision 1.1  1997/08/22 20:35:37  rst
// Initial revision
//
// Revision 1.3  1997/08/22 15:05:11  rst
// Cleaning up files
//
//
///////////////////////////////////////////////////////////////////////////
#include <OS.h>
#include <iostream.h>
#include <TDL/Surface.h>


class ByuSurface : public Surface {

public:

  // default constructor
  ByuSurface() :
    Surface() {}

  // copy constructor
  ByuSurface(ByuSurface &S) :
    Surface(S) {}

  // destructor
  ~ByuSurface() {}

  // instantiate functions for i/o
  bool read (istream &inFile = cin);
  bool read (const char *fileName);
  bool write (ostream &outFile = cout) const;
  bool write (const char *fileName) const;

  // standard assignment
  ByuSurface& operator = (ByuSurface const &S);

  // assignment from base class (for conversion)
  ByuSurface& operator = (Surface const &S);

};


#endif // __BYUSURFACE_H__
