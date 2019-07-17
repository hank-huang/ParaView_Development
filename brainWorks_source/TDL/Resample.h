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
// File: Resample.h 
//
// Author: Guoling Tong
//
// Purpose: Resample header file, providing a uniform interface for resampling. 
//         
// 
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: Resample.h,v 1.2 2004/11/28 03:37:52 joeh Exp $
//
// Revision History
//
// $Log: Resample.h,v $
// Revision 1.2  2004/11/28 03:37:52  joeh
// New Major Version will no longer compile on AIX
// but will compile on newer versions of gcc
// tested on gcc 3.2 gcc 3.3 and gcc 3.4
// It should run on Redhat 9.0/Fedora Core 1/Fedora Core 2/Fedora Core 3
//
// Revision 1.1  1999/10/01 15:27:59  kem
// Initial revision
//
// Revision 1.4  1999/09/29 14:45:03  tong
// checked in by kem
//
// Revision 1.3  1999/07/20 15:55:34  RAZORDB
// Unified Resampling
//
// Revision 1.2  1999/07/15 22:50:15  rst
// Unified Resampling
//
///////////////////////////////////////////////////////////////////////////

#ifndef __RESAMPLE_H__
#define __RESAMPLE_H__

#include <OS.h>
#include <TDL/Convolution.h>

template <class T>
class Resample : public Convolution<T>
{

 public:

  // a single interface for the resampling;
  static void resample(typename Convolution<T>::ConvolutionType flag,T const *volIn,
		       int x1,int y1,int z1,T *volOut,int x2,int y2,int z2);

 private:

  static const char *const _ClassRevision;
  static const char *const _ClassName;
  static const int         _ClassVersion;

};

template <class T> const char *const Resample<T>::_ClassRevision = "$Id: Resample.h,v 1.2 2004/11/28 03:37:52 joeh Exp $";
template <class T> const char *const Resample<T>::_ClassName = "Resample";
template <class T> const int Resample<T>::_ClassVersion = 1; 


#endif // __RESAMPLE_H__

