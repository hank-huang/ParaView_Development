#ifndef __LINE_H__
#define __LINE_H__

///////////////////////////////////////////////////////////////////////////
//
// IntellX, L.L.C., Proprietary
// Do not reproduce without permission in writing.
// Copyright (c) 1997 IntellX, L.L.C.
// All rights reserved
//
// File: Line.h
//
// Author: Rob Teichman
//
// Purpose: Header file for Line Class
// 
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: Line.h,v 1.1 2004/11/15 04:44:08 joeh Exp $
//
// Revision History
//
// $Log: Line.h,v $
// Revision 1.1  2004/11/15 04:44:08  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:27:36  kem
// Initial revision
//
// Revision 1.2  1999/05/12 17:25:30  RAZORDB
// update
//
// Revision 1.1  1997/08/22 20:35:34  rst
// Initial revision
//
// Revision 1.2  1997/08/22 20:24:08  rst
// Cleaning up files
//
// Revision 1.1  1997/08/05 14:27:57  rst
// Initial revision
//
// Revision 1.1  1997/07/15  19:28:18  rst
// Initial revision
//
//
///////////////////////////////////////////////////////////////////////////

#include <OS.h>
#include <iostream.h>
#include <stdio.h>
#define ADL_COMPLEX_DATATYPE
#include <ADL/Array1D.h>
#include <ADL/Array2D.h>
#include <TDL/Point.h>


class Line {

public:

  // default constructor
  Line() :
    mNumPoints(0) { sprintf(name,"<line>"); }

  // copy constructor
  Line(Line const &L);

  // destructor
  virtual ~Line() {}

  // assigment operator
  Line & operator= (Line const &L);

  // clear a line back to no segments
  void clearLine ()
  { mNumPoints = 0; mStartPoint.setDim(0); mEndPoint.setDim(0); }

  // add a segment to the line
  void addSegment (Point const start, Point const end);

  // combine two lines
  Line & operator+= (Line const &L);
  Line operator+ (Line const &P) const;

  // KWD: Array1D compatibility
  bool operator!= (Line const &L) const { return(false); }

  // number of segments in current line
  int numSegments() const
  { return mNumPoints; }

  // cumulative length of line
  float length() const;

  // print out a line
  void print(ostream &os) const;

  // get a segment
  bool getSegment(Point &start, Point &end, int num);

  // set a segments value
  bool setSegment(Point start, Point end, int num);

  // read/write
  virtual bool read(const char *fname);
  virtual bool read(istream &inFile=cin);
  virtual bool write(const char *fname) const;
  virtual bool write(ostream &outFile=cout) const;

  void setName(const char *nm);

protected:
  char name[256];

  Array1D<Point> mStartPoint;
  Array1D<Point> mEndPoint;
  int mNumPoints;

};

// Non-member print function
ostream & operator<< (ostream &os, Line const &L);

#endif // LINE_CLASS_H
