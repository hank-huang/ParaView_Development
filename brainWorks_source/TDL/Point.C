///////////////////////////////////////////////////////////////////////////
//
// IntellX, L.L.C., Proprietary
// Do not reproduce without permission in writing.
// Copyright (c) 1997 IntellX, L.L.C.
// All rights reserved
//
// File: Point.C
//
// Author: Rob Teichman
//
// Purpose: Point class body
// 
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: Point.C,v 1.1 2004/11/15 04:44:08 joeh Exp $
//
// Revision History
//
// $Log: Point.C,v $
// Revision 1.1  2004/11/15 04:44:08  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:28:44  kem
// Initial revision
//
// Revision 1.14  1999/07/09 17:47:39  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.13  1998/05/18 19:13:02  rst
// Changing float to double for accuracy required by eigen shapes
//
// Revision 1.11  1998/04/22 18:53:41  rst
// fix out of array bounds access
//
// Revision 1.10  1998/03/26 17:37:02  rst
// operator/= and operator*=
//
// Revision 1.8  1998/03/05 17:57:46  rst
// Change BOXSZ to 20
//
// Revision 1.7  1997/10/23 19:34:04  rst
// fixed bug in invFieldTrans
//
// Revision 1.6  1997/10/20 19:28:51  rst
// bug in call to seeded invFieldTrans
//
// Revision 1.5  1997/10/10 19:46:41  rst
// Add neighborhood generation and optimize
//
// Revision 1.4  1997/10/02 20:21:56  rst
// fix spike problem
//
// Revision 1.3  1997/09/19 15:48:51  rst
// optimized two loops in invFieldTrans
//
// Revision 1.2  1997/08/26 15:46:29  kem
// Optimize invFieldTrans()
//
// Revision 1.1  1997/08/22 20:36:49  rst
// Initial revision
//
// Revision 1.3  1997/08/22 20:24:22  rst
// Cleaning up files
//
// Revision 1.2  1997/08/20 19:28:20  rst
// Added invFieldTransformation routine
//
// Revision 1.1  1997/08/05 14:26:01  rst
// Initial revision
//
///////////////////////////////////////////////////////////////////////////

#include <iostream.h>
#include <math.h>
#include <ADL/Lapack.h>
#include <ADL/Array1D.h>
#include <TDL/Point.h>

//
// macros for common operations
//

#define sqr(x) ((x) * (x))
#define max(x,y) (((x) > (y)) ? (x) : (y))
#define min(x,y) (((x) < (y)) ? (x) : (y))

static bool debug = false;

//
// Definitions of member functions
//

//
// Function: set
//
// Description:
//	Set the value of a point to the values specified
//
//
void Point::set(double x, double y, double z)
{
  mX = x;
  mY = y;
  mZ = z;
}

//
// Assignment operator (=)
//
Point & Point::operator= (Point const &P)
{
  mX = P.mX;
  mY = P.mY;
  mZ = P.mZ;
  return *this;
}

Point & Point::operator= (double A)
{
  mX = A;
  mY = A;
  mZ = A;
  return *this;
}

Point & Point::operator= (int A)
{
  mX = (double)A;
  mY = (double)A;
  mZ = (double)A;
  return *this;
}

//
// Equality operator (==)
//
// Returns: boolean valued result of comparison
// Equality fails if any of the coordinates is not equal
//
bool Point::operator==(Point const &P) const
{
  return (mX == P.mX && mY == P.mY && mZ == P.mZ);
}

//
// Inequality operator (!=)
//
// Returns: boolean valued result of comparison
// Fails only if all coordinates match
//
bool Point::operator!=(Point const &P) const
{
  return (mX != P.mX || mY != P.mY || mZ != P.mZ);
}


//
// Basic math functions
//
// Addition
Point & Point::operator+= (Point const &P)
{
  mX += P.mX;
  mY += P.mY;
  mZ += P.mZ;
  return *this;
}

Point Point::operator+ (Point const &P) const
{
  Point temp(*this);
  temp += P;
  return temp;
}

// Subtraction
Point & Point::operator-= (Point const &P)
{
  mX -= P.mX;
  mY -= P.mY;
  mZ -= P.mZ;
  return *this;
}

Point Point::operator- (Point const &P) const
{
  Point temp(*this);
  temp -= P;
  return temp;
}

//
// Inner Product
//
double Point::innerProd (Point const &P) const
{
  return(mX*P.mX + mY*P.mY + mZ*P.mZ);
}

//
// Cross Product
//
Point Point::cross (Point const &P) const
{
  Point temp;
  temp.mX = mY * P.mZ - mZ * P.mY;
  temp.mY = mZ * P.mX - mX * P.mZ;
  temp.mZ = mX * P.mY - mY * P.mX;
  return temp;
}

//
// division by a scaler
//
Point & Point::operator/= (double const &f)
{
  mX /= f;
  mY /= f;
  mZ /= f;
  return *this;
}

Point Point::operator/ (double const &f) const
{
  Point temp(*this);
  temp /= f;
  return temp;
}

Point Point::operator/ (Point const &f) const
{
// do nothing for Array1D compatability only
  return(f);
}



//
// multiplication by a scaler
//

Point & Point::operator*= (double const &f)
{
  mX *= f;
  mY *= f;
  mZ *= f;
  return *this;
}

Point Point::operator* (double const &f) const
{
  Point temp(*this);
  temp *= f;
  return temp;
}



///////////////////////////////////////////////////////////////////////////
//
// invFieldTrans
//
// find the inverted location of the given Point as specified by the H field
//
// two versions are supplied.  The first version just takes the H-Field as
// a parameter, and finds the inverse.  The second version takes the H-field
// and a seed, and finds the inverse assuming it lies within a BOXSZ box
// around the seed.
//
///////////////////////////////////////////////////////////////////////////
bool Point::invFieldTrans (const Array3D<float> &H)
{
  const int SKIP = 4;		// subsampling spacing for first pass

  int i,j,k,t;

  // copy point components to local variables for convenience
  double px = x(),
    py = y(),
    pz = z();

  // variables to keep track of the min point for the first pass
  int xMin, yMin, zMin;
  double minDist = HUGE;
  double dist;


  //
  // First do a rough pass over the volume and find the closest point
  //
  const float * const * const *hFieldPtrPtrPtr = H.address();

  for (k=0; k<H.getZsize(); k+=SKIP)
    for (j=0; j<H.getYsize(); j+=SKIP)
      for (i=0; i<H.getXsize()/3; i+=SKIP)
      {
	double distX = hFieldPtrPtrPtr[k][j][3*i] - px;
	double distY = hFieldPtrPtrPtr[k][j][3*i+1] - py;
	double distZ = hFieldPtrPtrPtr[k][j][3*i+2] - pz;

	// Compute the squared distace to the vertex from the current
	// lattice point
	dist = (distX*distX + distY*distY + distZ*distZ);

	// if dist == 0 we are done!!
	if (dist == 0)
	{
	  set(i, j, k);
	  cout << "Exact Match Found!!" << endl;
	  return true;
	}

	// Check if this point is closer than the previous closest
	if (dist < minDist)
	{
	  xMin = i;
	  yMin = j;
	  zMin = k;
	  minDist = dist;
	}
      }

  Point seed(xMin, yMin, zMin);
  return invFieldTrans(H, seed);
}


bool Point::invFieldTrans (const Array3D<float> &H, const Point seed)
{
  const int NBHDSZ = 36;       	// number of neighbors to keep track of
  const int BOXSZ = 20;		// size of box for final search pass
  const int HBOXSZ = BOXSZ/2;	// 1/2 the box size

  int i, j, k, t;
  double dist;
  double px = x(),
    py = y(),
    pz = z();
  // these factors were added to scale into the unit cube
//   int xs = H.getXsize()/3 - 1,
//     ys = H.getYsize() - 1,
//     zs = H.getZsize() - 1;

  // use integer valued seed voxel
  // if seed lies outside the volume, use the closest point
  int xSeed = (int) (seed.x() + 0.5),
    ySeed = (int) (seed.y() + 0.5),
    zSeed = (int) (seed.z() + 0.5);

  if (xSeed < 0) xSeed = 0;
  if (xSeed > H.getXsize()/3 - 1) xSeed = H.getXsize()/3 - 1;

  if (ySeed < 0) ySeed = 0;
  if (ySeed > H.getYsize() - 1) ySeed = H.getYsize() - 1;

  if (zSeed < 0) zSeed = 0;
  if (zSeed > H.getZsize() - 1) zSeed = H.getZsize() - 1;


  //
  // Search over a box of size BOXSZ in each direction around the seed
  // to find the N minimum distance points (or an exact match)
  //

  // tables of values used in determining nearest grid point to desired point
  Array1D<int> xMins(NBHDSZ),
    yMins(NBHDSZ),
    zMins(NBHDSZ);
  Array1D<float> minDists(NBHDSZ);
  xMins = 0;
  yMins = 0;
  zMins = 0;
  minDists = HUGE;

  const float * const * const *hFieldPtrPtrPtr = H.address();
  for (k=max(0, zSeed-HBOXSZ); k<min(H.getZsize(), zSeed+HBOXSZ); k++)
    for (j=max(0, ySeed-HBOXSZ); j<min(H.getYsize(), ySeed+HBOXSZ); j++)
      for (i=max(0, xSeed-HBOXSZ); i<min(H.getXsize()/3, xSeed+HBOXSZ); i++)
      {
	double distX = hFieldPtrPtrPtr[k][j][3*i] - px;
	double distY = hFieldPtrPtrPtr[k][j][3*i+1] - py;
	double distZ = hFieldPtrPtrPtr[k][j][3*i+2] - pz;

	// Compute the distace to the points
	dist = (distX*distX + distY*distY + distZ*distZ);

        // if dist == 0 we are done!!
	if (dist == 0)
	{
// 	  set((double)i/xs, (double)j/ys, (double)k/zs);
	  set(i, j, k);
	  cout << "Exact Match Found!!" << endl;
	  return true;
	}

	// check if dist is less than any saved distance, if so
	// insert point into table in order
	if (dist < minDists[NBHDSZ-1])
	{
	  t = NBHDSZ-2;
	  while (t >= 0 && dist < minDists[t])
	  {
	    minDists[t+1] = minDists[t];
	    xMins[t+1] = xMins[t];
	    yMins[t+1] = yMins[t];
	    zMins[t+1] = zMins[t];
	    t--;
	  }
	  t++;
	  minDists[t] = dist;
	  xMins[t] = i;
	  yMins[t] = j;
	  zMins[t] = k;
	}
      }

  //
  // check if any of the minimum distance points lie on the
  // boundary of the box, indicating that the seed was
  // probably bad
  //
  if (debug)
    for (i=0; i<NBHDSZ; i++)
    {
      if (xMins[i] <= xSeed-HBOXSZ ||
	  xMins[i] >= xSeed+HBOXSZ)
	cout << "X boundary: i=" << i << ", x=" << xMins[i] 
	     << ", xSeed=" << xSeed << endl;

      if (yMins[i] <= ySeed-HBOXSZ ||
	  yMins[i] >= ySeed+HBOXSZ)
	cout << "Y boundary: i=" << i << ", y=" << yMins[i]
	     << ", ySeed=" << ySeed << endl;

      if (zMins[i] <= zSeed-HBOXSZ ||
	  zMins[i] >= zSeed+HBOXSZ)
	cout << "Z boundary: i=" << i << ", z=" << zMins[i]
	     << ", zSeed=" << zSeed << endl;
    }


  //
  // Now do the interpolation.
  // build a representation of the desired point using the NBHDSZ lattice
  // points found that are closest to it.
  //

  double A[4][NBHDSZ];
  double B[3][NBHDSZ];
  double Q[4][4];
  double AA[3][4];

  // Fill the matricies
  for (i=0; i<NBHDSZ; i++)
  {
//     B[0][i] = (double) xMins[i]/xs;
//     B[1][i] = (double) yMins[i]/ys;
//     B[2][i] = (double) zMins[i]/zs;
    B[0][i] = xMins[i];
    B[1][i] = yMins[i];
    B[2][i] = zMins[i];

    A[0][i] = H[zMins[i]][yMins[i]][3*xMins[i]] - px;
    A[1][i] = H[zMins[i]][yMins[i]][3*xMins[i]+1] - py;
    A[2][i] = H[zMins[i]][yMins[i]][3*xMins[i]+2] - pz;
    A[3][i] = 1.0;
  }

  // Build the matrix Q = AA'
  for (i=0; i<4; i++)
    for (j=0; j<4; j++)
    {
      Q[i][j] = 0;
      for (k=0; k<NBHDSZ; k++)
	Q[i][j] += A[i][k] * A[j][k];
    }

  // Build AA = BA'
  // From the math, this is actually AB', but we build (AB')' = BA'
  // because the fortran (LAPACK) routine called assumes Column major
  // instead of row major (like C) order
  for(i=0; i<3; i++)
    for(j=0; j<4; j++)
    {
      AA[i][j] = 0;
      for(k=0; k<NBHDSZ; k++)
	AA[i][j] += B[i][k] * A[j][k];
    }

  // Print the matrix Q to test
//   cout << "Point::invFieldTrans Q Matrix is:\n  ";
//   for(i=0;i<4;i++)
//   {
//     for(j=0;j<4;j++)
//       cout << Q[i][j] << " ";
//     cout << endl << "  ";
//   }
//   cout << endl;

  // Set up variables for LAPACK routine
  // Use RCOND of -1 to let LAPACK do what it likes
  int INFO = 1,
    LDA = 4,
    LDB = 4,
    M = 4,
    N = 4,
    NRHS = 3,
    RANK;
  double RCOND = -1;
  double work[100];
  int JPVT[4];

  for(i=0; i<4; i++)
    JPVT[i] = 0;

  // call routine to compute the least squares approx
  //
  // Since this is a fortran routine, all matrices are really transposed.
  // Thus, the formulation of the problem shows that the result of this
  // minimization is actually G', however since it is G' in Fortran, it
  // is G here, so there is no need to transpose the result
  //
  Lapack::xgelsx(M, N, NRHS, &Q[0][0], LDA, &AA[0][0], LDB, JPVT,
		 RCOND, &RANK, work, &INFO);

  if (debug)
    cout << "info = " << INFO << ", RANK = " << RANK << endl;

  // check for successful completion
  if (INFO != 0)
  {
    cerr << "Point::invFieldTrans: error in xgelsx, info = " << INFO << endl;
    return false;
  }

  // output AA for test
//   cout << "Point::invFieldTrans AA (result) is:\n  ";
//   for (i=0; i<3; i++)
//   {
//     for (j=0; j<3; j++)
//       cout << AA[i][j] << " ";
//     cout << "\n  ";
//   }
//   cout << endl;

  set(AA[0][3], AA[1][3], AA[2][3]);

// only do this check when working in the unit cube
//  if (AA[0][3] > 1.0 || AA[1][3] > 1.0 || AA[2][3] > 1.0)
//    cout << "Error! value > 1.0" << endl;

  if (debug)
    cout << "Result: " << *this <<  endl;

  return true;
}


//
// Print the value of a point to an ostream.
//
void Point::print(ostream &os) const
{
  os << mX << " " << mY << " " << mZ;
}


//
// non-member print function
//
ostream & operator<< (ostream &os, Point const &P)
{
  P.print(os);
  return os;
}
