#ifndef __SURFACE_UTILS__
#define __SURFACE_UTILS__

///////////////////////////////////////////////////////////////////////////
//
// IntellX, L.L.C., Proprietary
// The contents of this file comprise confidential, proprietary and/or
// trade secret information which is the property of IntellX LLC.
// Dissemination, distribution or other disclosure of this information
// is strictly prohibited without the express permission of IntellX LLC.
// Copyright (c) 1997 IntellX, L.L.C.
// All rights reserved
//
// File: SurfaceUtils.h
//
// Author: Rob Teichman
//
// Purpose: Utility routines for operating on Surface Class objects
// 
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: SurfaceUtils.h,v 1.4 2008/08/08 20:34:00 mbowers Exp $
//
// Revision History
//
// $Log: SurfaceUtils.h,v $
// Revision 1.4  2008/08/08 20:34:00  mbowers
// Add a function to get the path length (for the script CurveGen script function).
//
// Revision 1.3  2006/05/31 19:31:15  mbowers
// Add a surface area utility function.
//
// Revision 1.2  2004/11/20 06:21:45  lei
// *** empty log message ***
//
// Revision 1.1  2004/11/15 04:44:08  joeh
// first check in of BrainWorks
//
// Revision 1.2  1999/10/11 16:25:59  sarang
// Revision
//
// Revision 1.1  1999/10/01 15:28:02  kem
// Initial revision
//
// Revision 1.10  1999/09/21 22:01:22  RAZORDB
// TDL update
//
// Revision 1.9  1999/05/17 15:41:19  RAZORDB
// Add findAffine member
//
// Revision 1.8  1999/05/12 17:25:44  RAZORDB
// update
//
// Revision 1.7  1999/01/19 17:21:55  RAZORDB
// kwd
//
// Revision 1.6  1998/11/16 19:07:47  kwd
// add curvature generation
//
// Revision 1.5  1998/05/18 19:07:47  rst
// apply for Manifold directly on Surfaces (no h-field)
//
// Revision 1.4  1998/05/12 20:10:30  sarang
// Add functions needed for Schiz.
//
// Revision 1.3  1998/03/05 17:59:50  rst
// metrics
//
// Revision 1.2  1997/11/26 22:39:34  rst
// Made all routines static, new routines
//
// Revision 1.1  1997/11/21 21:27:32  rst
// Initial revision
//
//
///////////////////////////////////////////////////////////////////////////

#include <OS.h>
#include <ADL/GreyVolume.h>
#include <ADL/Array1D.h>
#include <TDL/ManifoldTransformation.h>
#include <TDL/Surface.h>

//
// SurfaceUtils class
//

class SurfaceUtils {
public:

  static double surfaceArea(Surface &S,
               float xPixelSize,
               float yPixelSize,
               float zPixelSize);

  static ItXECode redistribute(Surface &S,int niter, double alpha);

  static void smooth(Surface &S, int windowSize);

  static void deformSphereToPoints(int numFixedPoints,Pnt *points,Surface &S);

  static float plgndr(int l,int m,float x);

  static float plgndr_der(int l, float x);
	
  static float BuildCov(float A,int N);

  static int  deformToPoints(Surface &S,int numFixedPoints, Pnt *points);

  static void unitSphere (Surface &,int maxlevel = 4);

  static int maxCoordDist (Surface const &S1,
			   Surface const &S2,
			   float *xMaxDist,
			   float *yMaxDist,
			   float *zMaxDist);

  static float maxVertexDist (Surface const &S1,
			      Surface const &S2);

  static float meanDist (Surface const &S1,
			 Surface const &S2);

  static int meanVectorDiff (Surface const &S1,
			     Surface const &S2,
			     Point *mean);

  static int findRotation(Surface const &surf1,
			  Surface const &surf2,
			  Matrix<float> *R,
			  Vector<float> *trans);

  static int findAffine (Surface const &surf1,
			 Surface const &surf2,
			 Matrix<float> *R,
			 Vector<float> *trans);

  //change of coordinate functions

  static void rectToSph(Point const &P,double &theta, double &phi);
  static void sphToRect(Point &P, double const &rad, double const &theta,
			double const &phi); 	
  static void rectToStereo(Point const &P,double &u, double &v);
  static void stereoToRect(Point &P, double const &u, double const &v);
  
  static void findStereoBasis(Point const &P, Point &basis_u, Point &basis_v);	

  //find the N-segment geodesic (arc) between two points on sphere
      
  static void rotatePtToPt(Point const &P,Point const &Q, 
			   Array1D<Point> &path, int const &N); 

  static void  alignNP(Point &NP, Surface &surf);
  static void  unalignNP(Point &P, Surface &surf);  
  static void alignNP(Point &NP, Matrix<double> &landmks);
  static void findSphereRotation(Matrix<double> &A1, Matrix<double>&B1,
Matrix<double> *R);
  // manifold transformation
  
  static int applyManiTrans(ManifoldTransformation &mani, Surface &S);
  
  static ItXECode smoothSurfaceAvg (Surface &surf, Surface &S,int numiter);
  static ItXECode normalizeSurface (Surface &surf, Surface &S,double tv);
  static ItXECode projectSurfaceToSphere(Surface &S);
  static ItXECode projectSurfaceToSphere(Surface &S, Surface &out);
  static ItXECode closedSurfaceToSphere(Surface &surf, Surface &S,int numiter);
  
  // 'delpts' array must be of size S.numPoints()
  //          delpts[i] != 0 means delete point i
  static ItXECode removeSurfacePoints(Surface &S, Array1D<u_char> &delpts);
  static ItXECode keepSurfacePoints(Surface &S, Array1D<u_char> &delpts); // Lei 01/22/2004
  // 'delfac' array must be of size S.numPolys()
  //          delfac[i] != 0 means delete poly i
  static ItXECode removeSurfacePolys(Surface &S, Array1D<u_char> &delfac);

  // find the length of the path specified by the indeces in the array
  static ItXECode findPathLength(
    Surface &S, Array1D<int> &path, double& length);

  //
  // Remove all duplicate facets and un-used points
  //
  static ItXECode cleanSurface(Surface &S);
private:
  //
  // smoothing helper fcn
  //
  static ItXECode Ill1(Matrix<float> &r, Matrix<float> &b, float thresh);

  //
  // set array to 1 if that vertice is along an edge
  //
  static void tagEdges(Surface &S, Array1D<int> &edges);

  //
  // return the midpoint on a line between two points
  //
  static Point midpoint (Point const &a,
			 Point const &b)
  { return ((a+b)/2); }

  //
  // normalize a point to unit magnitude
  //
  static Point normalize (Point const &p)
  { return (p / p.norm()); }

};

#endif // __SURFACE_UTILS__
