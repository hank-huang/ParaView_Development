#ifndef __POINT_H__
#define __POINT_H__

///////////////////////////////////////////////////////////////////////////
//
// File: Point.h
//
// Author: Rob Teichman
//
// Purpose: Point Class definition
//  A single point in R^3 and related operations.
// 
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: Point.h,v 1.1 2004/11/15 04:44:08 joeh Exp $
//
// Revision History
//
// $Log: Point.h,v $
// Revision 1.1  2004/11/15 04:44:08  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:27:37  kem
// Initial revision
//
// Revision 1.8  1999/07/09 17:47:07  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.7  1998/05/18 19:13:01  rst
// Changing float to double for accuracy required by eigen shapes
//
// Revision 1.5  1998/03/26 17:37:01  rst
// operator/= and operator*=
//
// Revision 1.3  1997/10/20 19:28:49  rst
// bug in call to seeded invFieldTrans
//
// Revision 1.2  1997/08/26 15:43:19  kem
// Change invFieldTrans() hField argument to reference
//
// Revision 1.1  1997/08/22 20:35:36  rst
// Initial revision
//
// Revision 1.3  1997/08/22 20:24:13  rst
// Cleaning up files
//
// Revision 1.2  1997/08/20 19:28:15  rst
// Added invFieldTransformation routine
//
// Revision 1.1  1997/08/05 14:28:08  rst
// Initial revision
//
///////////////////////////////////////////////////////////////////////////
#include <OS.h>
#include <iostream.h>
#include <math.h>
#include <ADL/Array3D.h>


class Point {

public:

  // Default constructor.
  Point() :
    mX(0), mY(0), mZ(0) {}

  // Create point init'd to another Point.
  Point(Point const &P) :
    mX(P.mX), mY(P.mY), mZ(P.mZ) {}

  // Create point init'd to coords specified
  Point(double x, double y, double z) :
    mX(x), mY(y), mZ(z) {}
  
  // Default Destructor.
  ~Point() {}

  // KWD
  Point(double A) : mX(A), mY(A), mZ(A) {}
  Point(int A) : mX((double)A), mY((double)A), mZ((double)A) {}

  // Set a point to a location.
  void set(double x, double y, double z);


  // Get value of one dimension of a Point
  double x() const
  { return (mX); }
  double y() const
  { return (mY); }
  double z() const
  { return (mZ); }
  
  // Set a point from another point
  Point & operator= (Point const &P);
  Point & operator= (double P);
  Point & operator= (int P);

  // Equality/inequality of points.  Two points are equal
  // if they are equal in all dimensions.
  bool operator == (Point const &P) const;
  bool operator != (Point const &P) const;

  // KWD : Array1D<Point> compatibility
  bool operator > (Point const &P) const { return(this->norm() > P.norm()); }
  bool operator < (Point const &P) const { return(this->norm() < P.norm()); }
  
  // Returns the length of the vector
  double norm() const
  { return (sqrt(this->innerProd(*this))); }

  // Basic math functions
  //
  // Addition
  Point & operator+= (Point const &P);
  Point operator+ (Point const &P) const;

  // Subtraction
  Point & operator-= (Point const &P);
  Point operator- (Point const &P) const;

  // inner product
  double innerProd (Point const &P) const;

  // cross product
  Point cross (Point const &P) const;

  // division by a scalar
  Point & operator/= (double const &f);
  Point operator/ (double const &f) const;
  Point operator/ (Point const &f) const;

  // multiplication by a scalar
  Point & operator*= (double const &f);
  Point operator* (double const &f) const;

  // KWD : for Array1D compatibility
  Point operator* (const Point &P) const { return(*this); }

  // inverse field transformation on a point
  bool invFieldTrans (const Array3D<float> &H);
  bool invFieldTrans (const Array3D<float> &H, const Point seed);

  // Print the value of a point.
  void print(ostream &os) const;
  

private:

  // Data Elements
  double mX,
    mY,
    mZ;

};


//
// Definition of non-member print function for cout
//
ostream & operator<< (ostream &os, Point const &P);



#endif // __POINT_H__
