///////////////////////////////////////////////////////////////////////////
//
// File: TransformationUtils.C
//
// Author: csb (Guoling Tong modified the FFT routine to make it work with
//              ITX_FFT)
//         
// Purpose: provide routines that all transformation classes share
// 
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: TransformationUtils.C,v 1.1 2004/11/15 04:44:09 joeh Exp $
//
// Revision History
//
// $Log: TransformationUtils.C,v $
// Revision 1.1  2004/11/15 04:44:09  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:28:54  kem
// Initial revision
//
// Revision 1.10  1999/07/09 17:53:50  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.9  1998/10/29 00:19:55  RAZORDB
// Add method for interpolation of gradient field
//
// Revision 1.8  1998/07/15 21:19:50  kem
// Take out floor() in interpolate methods
//
// Revision 1.7  1998/07/08 22:01:09  kem
// Fix prototype for compactdisplay()
//
// Revision 1.6  1998/06/01 15:38:31  csb
// add compactdisplay() function
//
// Revision 1.5  1998/05/12 22:53:56  abed
// Change the Resample routine.
//
// Revision 1.4  1998/04/08 20:52:01  kem
// Add loadLandmarks() method
//
// Revision 1.3  1998/03/16 21:00:55  csb
// print errors to cerr instead of cout and ...
//
// Revision 1.2  1998/01/17 19:10:44  csb
// Add loadvolume() and resample for shorts
//
// Revision 1.1  1997/12/12 15:58:48  csb
// Initial revision
//
//
/////////////////////////////////////////////////////////////////////////////
#include <fstream.h>
#include <iostream.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <ADL/Array1D.h>
#include <ADL/Array2D.h>
#include <ADL/Array3D.h>
#include <ADL/GreyVolume.h>

#include <TDL/TransformationUtils.h>
#include <TDL/AnalyzeVol.h>
#include <TDL/ITX_FFT.h>

/////////////////////////////////////////////////////////////////////////////////////////
// These functions are used by the trilinear template functions
// Check if Cube indices are in-range.
//
// Abed M. Hammoud
// 09Oct1995-ESSRL
// Copyright 1993-1995
/////////////////////////////////////////////////////////////////////////////////////////
inline int inRange(int zmin, int zmax, int ymin, int ymax, int xmin, int xmax, 
	int z, int y, int x) {
   return (x >= xmin && x < xmax && y >= ymin && y < ymax && z >= zmin && z < zmax); 
}

inline double evalExp(double fx,  double v0,  double v1) {
   return v0 + fx * (v1 - v0); 
}



/*****************************************************************************************
* Tri-linear interpolation.
*
* The routine employs Tri-linear interpolation to overcome rounding errors.
* For a pixel at (z, y, x) in the final volume Y, find the nearest
* 8 pixels to the pixel (z0, y0, x0) in the original volume X and 
* interpolate. For points outside of original image support, 
* use background.
*
* Reference : Tri-linear interpolation. Graphics Gym X, pp 521.
* Steve Hill.
* 
* Abed M. Hammoud
* 09Oct1995-ESSRL
* Copyright 1993-1995
*****************************************************************************************/
template <class T>
double
trilinearTemplate( const Array3D<T> & X,
                         double       z,
                         double       y,
                         double       x,
                         T            bkgrnd )
{
  int zmin = 0;
  int ymin = 0;
  int xmin = 0;
  int zmax = X.getZsize();
  int ymax = X.getYsize();
  int xmax = X.getXsize();

   int z0 = int(floor(z)), z1 = z0 + 1; 
   int y0 = int(floor(y)), y1 = y0 + 1; 
   int x0 = int(floor(x)), x1 = x0 + 1; 
   int xsize = X.getXsize(), ysize = X.getYsize(); 
   double d000, d001, d010, d011, d100, d101, d110, d111; 

   /*****************************************************************
   * (0 <= fx, fy, fz <= 1), fractional position between data points.
   *****************************************************************/

   double fz = z - z0, fy = y - y0, fx = x - x0; 

   /**********************************************
   * Tri-linear interpolation for interior points.
   **********************************************/

   if (x0 >= xmin && x1 < xmax && y0 >= ymin && y1 < ymax && 
       		z0 >= zmin && z1 < zmax) {

      // Use pointer arithmetic.
      T const *dp = &X[z0][y0][x0]; 

      /****************************
      * Data used in interpolation.
      ****************************/

      d000 = dp[0]; d100 = dp[1]; 
      dp += xsize; 

      d010 = dp[0]; d110 = dp[1]; 
      dp += xsize * ysize; 

      d011 = dp[0]; d111 = dp[1]; 
      dp -= xsize; 

      d001 = dp[0]; d101 = dp[1]; 

   } else {

      /**********************************************
      * Tri-linear interpolation for boundary points.
      * Points outside of cube are set to bkgrnd.
      **********************************************/

      d000 = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z0, y0, x0) 
		? X[z0][y0][x0] : bkgrnd; 
      d001 = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z1, y0, x0) 
		? X[z1][y0][x0] : bkgrnd; 
      d010 = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z0, y1, x0) 
		? X[z0][y1][x0] : bkgrnd; 
      d011 = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z1, y1, x0) 
		? X[z1][y1][x0] : bkgrnd; 

      d100 = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z0, y0, x1) 
		? X[z0][y0][x1] : bkgrnd; 
      d101 = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z1, y0, x1) 
		? X[z1][y0][x1] : bkgrnd; 
      d110 = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z0, y1, x1) 
		? X[z0][y1][x1] : bkgrnd; 
      d111 = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z1, y1, x1) 
		? X[z1][y1][x1] : bkgrnd; 

   }

   /********************
   * Interpolation, (1).
   ********************/

   double dx00 = evalExp(fx, d000, d100); 
   double dx01 = evalExp(fx, d001, d101); 
   double dx10 = evalExp(fx, d010, d110); 
   double dx11 = evalExp(fx, d011, d111); 

   /********************
   * Interpolation, (2).
   ********************/

   double dxy0 = evalExp(fy, dx00, dx10); 
   double dxy1 = evalExp(fy, dx01, dx11); 

   /********************
   * Interpolation, (3).
   ********************/

   return evalExp(fz, dxy0, dxy1); 

}

/*****************************************************************************************
* Tri-linear interpolation.
*
* The routine employs Tri-linear interpolation to overcome rounding errors.
* For a pixel at (z, y, x) in the final volume Y, find the nearest
* 8 pixels to the pixel (z0, y0, x0) in the original volume I and 
* interpolate. For points outside of original image support, 
* use background.
*
* This routine returns the values of the X, Y, Z cubes at the location
* (z, y, x). X(z, y, x), Y(z, y, x) and Z(z, y, x) are stored in
* vx[0], vx[1], and vx[2] respectively.
*
* Reference : Tri-linear interpolation. Graphics Gym X, pp 521.
* Steve Hill.
* 
* Abed M. Hammoud
* 09Oct1995-ESSRL
* Copyright 1993-1995
*****************************************************************************************/
template <class T>
void trilinearTemplate(       Array1D<double> * vx,
                        const Array3D<T>      & X,
                        const Array3D<T>      & Y, 
                        const Array3D<T>      & Z,
                              double            z,         
                              double            y,
                              double            x,
                              T                 bkgrnd )
{

  int zmin = 0;
  int ymin = 0;
  int xmin = 0;
  int zmax = X.getZsize();
  int ymax = X.getYsize();
  int xmax = X.getXsize();

   int z0 = int(floor(z)), z1 = z0 + 1; 
   int y0 = int(floor(y)), y1 = y0 + 1; 
   int x0 = int(floor(x)), x1 = x0 + 1; 
   int xsize = X.getXsize(), ysize = X.getYsize(); 
   double d000x, d001x, d010x, d011x, d100x, d101x, d110x, d111x; 
   double d000y, d001y, d010y, d011y, d100y, d101y, d110y, d111y; 
   double d000z, d001z, d010z, d011z, d100z, d101z, d110z, d111z; 

   /*****************************************************************
   * (0 <= fx, fy, fz <= 1), fractional position between data points.
   *****************************************************************/

   double fz = z - z0, fy = y - y0, fx = x - x0; 

   /**********************************************
   * Tri-linear interpolation for interior points.
   **********************************************/

   if (x0 >= xmin && x1 < xmax && y0 >= ymin && y1 < ymax && 
       		z0 >= zmin && z1 < zmax) {

      // Use pointer arithmetic.
      T const *dpx = &X[z0][y0][x0]; 
      T const *dpy = &Y[z0][y0][x0]; 
      T const *dpz = &Z[z0][y0][x0]; 

      /****************************
      * Data used in interpolation.
      ****************************/

      d000x = dpx[0]; d100x = dpx[1]; 
      d000y = dpy[0]; d100y = dpy[1]; 
      d000z = dpz[0]; d100z = dpz[1]; 
      dpx += xsize; 
      dpy += xsize; 
      dpz += xsize; 

      d010x = dpx[0]; d110x = dpx[1]; 
      d010y = dpy[0]; d110y = dpy[1]; 
      d010z = dpz[0]; d110z = dpz[1]; 
      dpx += xsize * ysize; 
      dpy += xsize * ysize; 
      dpz += xsize * ysize; 

      d011x = dpx[0]; d111x = dpx[1]; 
      d011y = dpy[0]; d111y = dpy[1]; 
      d011z = dpz[0]; d111z = dpz[1]; 
      dpx -= xsize; 
      dpy -= xsize; 
      dpz -= xsize; 

      d001x = dpx[0]; d101x = dpx[1]; 
      d001y = dpy[0]; d101y = dpy[1]; 
      d001z = dpz[0]; d101z = dpz[1]; 

   } else {

      /**********************************************
      * Tri-linear interpolation for boundary points.
      * Points outside of cube are set to bkgrnd.
      **********************************************/

      d000x = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z0, y0, x0) 
		? X[z0][y0][x0] : bkgrnd; 
      d000y = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z0, y0, x0) 
		? Y[z0][y0][x0] : bkgrnd; 
      d000z = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z0, y0, x0) 
		? Z[z0][y0][x0] : bkgrnd; 

      d001x = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z1, y0, x0) 
		? X[z1][y0][x0] : bkgrnd; 
      d001y = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z1, y0, x0) 
		? Y[z1][y0][x0] : bkgrnd; 
      d001z = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z1, y0, x0) 
		? Z[z1][y0][x0] : bkgrnd; 

      d010x = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z0, y1, x0) 
		? X[z0][y1][x0] : bkgrnd; 
      d010y = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z0, y1, x0) 
		? Y[z0][y1][x0] : bkgrnd; 
      d010z = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z0, y1, x0) 
		? Z[z0][y1][x0] : bkgrnd; 

      d011x = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z1, y1, x0) 
		? X[z1][y1][x0] : bkgrnd; 
      d011y = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z1, y1, x0) 
		? Y[z1][y1][x0] : bkgrnd; 
      d011z = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z1, y1, x0) 
		? Z[z1][y1][x0] : bkgrnd; 

      d100x = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z0, y0, x1) 
		? X[z0][y0][x1] : bkgrnd; 
      d100y = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z0, y0, x1) 
		? Y[z0][y0][x1] : bkgrnd; 
      d100z = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z0, y0, x1) 
		? Z[z0][y0][x1] : bkgrnd; 

      d101x = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z1, y0, x1) 
		? X[z1][y0][x1] : bkgrnd; 
      d101y = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z1, y0, x1) 
		? Y[z1][y0][x1] : bkgrnd; 
      d101z = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z1, y0, x1) 
		? Z[z1][y0][x1] : bkgrnd; 

      d110x = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z0, y1, x1) 
		? X[z0][y1][x1] : bkgrnd; 
      d110y = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z0, y1, x1) 
		? Y[z0][y1][x1] : bkgrnd; 
      d110z = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z0, y1, x1) 
		? Z[z0][y1][x1] : bkgrnd; 

      d111x = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z1, y1, x1) 
		? X[z1][y1][x1] : bkgrnd; 
      d111y = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z1, y1, x1) 
		? Y[z1][y1][x1] : bkgrnd; 
      d111z = inRange(zmin, zmax, ymin, ymax, xmin, xmax, z1, y1, x1) 
		? Z[z1][y1][x1] : bkgrnd; 

   }

   /********************
   * Interpolation, (1).
   ********************/

   double dx00x = evalExp(fx, d000x, d100x); 
   double dx00y = evalExp(fx, d000y, d100y); 
   double dx00z = evalExp(fx, d000z, d100z); 

   double dx01x = evalExp(fx, d001x, d101x); 
   double dx01y = evalExp(fx, d001y, d101y); 
   double dx01z = evalExp(fx, d001z, d101z); 

   double dx10x = evalExp(fx, d010x, d110x); 
   double dx10y = evalExp(fx, d010y, d110y); 
   double dx10z = evalExp(fx, d010z, d110z); 

   double dx11x = evalExp(fx, d011x, d111x); 
   double dx11y = evalExp(fx, d011y, d111y); 
   double dx11z = evalExp(fx, d011z, d111z); 

   /********************
   * Interpolation, (2).
   ********************/

   double dxy0x = evalExp(fy, dx00x, dx10x); 
   double dxy0y = evalExp(fy, dx00y, dx10y); 
   double dxy0z = evalExp(fy, dx00z, dx10z); 

   double dxy1x = evalExp(fy, dx01x, dx11x); 
   double dxy1y = evalExp(fy, dx01y, dx11y); 
   double dxy1z = evalExp(fy, dx01z, dx11z); 

   /********************
   * Interpolation, (3).
   ********************/

   (*vx)[0] = evalExp(fz, dxy0x, dxy1x); 
   (*vx)[1] = evalExp(fz, dxy0y, dxy1y); 
   (*vx)[2] = evalExp(fz, dxy0z, dxy1z); 

}


/*****************************************************************************************
* G = Gradient(X, x, y, dx, dy) returns the numerical partial
* derivatives of Cube X at point (x, y, z), in array G, such that 
* G[0] = dX / dx, G[1] = dX / dy, G[2] = dX / dz.  
*
* dx, dy, and dz are scalars containing the sample spacing in the 
* X, Y, and Z directions.
*
* Its assumed that the volume of interest resides in the subvolume:
*
* 	xmin <= x < xmax
* 	ymin <= x < ymax
* 	zmin <= x < zmax
*
* Gradient(X) assumes dx = dy = dz = 1.  
* 
* Abed Hammoud
* IntellX
*****************************************************************************************/
template <class T>
void
gradientTemplate(       Array1D<double> * G, 
                  const Array3D<T>      & F,
                        int               x,
                        int               y,
                        int               z, 
                        double            dx,
                        double            dy,
                        double            dz  )
{

  int zmin = 0;
  int ymin = 0;
  int xmin = 0;
  int zmax = F.getZsize();
  int ymax = F.getYsize();
  int xmax = F.getXsize();

/***************
* House keeping.
***************/

   if (G->getNelm() != 3)
      throwError("gradient: output array have wrong size."); 

   *G = 0.0; 

/***************************************************************
* Interior points, (Gradient is continuous from right and left).
***************************************************************/

   if (x > xmin && x < (xmax - 1)) 
      (*G)[0] = ((double) F[z][y][x + 1] - (double) F[z][y][x - 1]) / (2.0 * dx); 

/*****************
* Boundary points.
*****************/

   else if (xmax == 1)		// For a degenerate cube.
      (*G)[0] = 0.0; 

   else if (x == xmin)
      (*G)[0] = ((double) F[z][y][x + 1] - (double) F[z][y][x]) / dx; 

   else if (x == (xmax - 1))
      (*G)[0] = ((double) F[z][y][x] - (double) F[z][y][x - 1]) / dx; 

/******************************************************************
* For points where the function is not defined its not a good
* idea to return zero for the gradient. This is because a gradient
* of zero might indicate a maxima/minima point for the function.
******************************************************************/

   else 
      throwError("gradient(): The gradient is not defined at this point."); 

/***************************************************************
* Interior points, (Gradient is continuous from right and left).
***************************************************************/

   if (y > ymin && y < (ymax - 1)) 
      (*G)[1] = ((double) F[z][y + 1][x] - (double) F[z][y - 1][x]) / (2 * dy); 

/*****************
* Boundary points.
*****************/

   else if (ymax == 1)		// For a degenerate cube.
      (*G)[1] = 0.0; 

   else if (y == ymin)
      (*G)[1] = ((double) F[z][y + 1][x] - (double) F[z][y][x]) / dy; 

   else if (y == (ymax - 1))
      (*G)[1] = ((double) F[z][y][x] - (double) F[z][y - 1][x]) / dy; 

/******************************************************************
* For points where the function is not defined its not a good
* idea to return zero for the gradient. This is because a gradient
* of zero might indicate a maxima/minima point for the function.
******************************************************************/

   else 
      throwError("gradient(): The gradient is not defined at this point."); 

/***************************************************************
* Interior points, (Gradient is continuous from right and left).
***************************************************************/

   if (z > zmin && z < (zmax - 1)) 
      (*G)[2] = ((double) F[z + 1][y][x] - (double) F[z - 1][y][x]) / (2 * dz); 

/*****************
* Boundary points.
*****************/

   else if (zmax == 1)		// For a degenerate cube.
      (*G)[2] = 0.0; 

   else if (z == zmin)
      (*G)[2] = ((double) F[z + 1][y][x] - (double) F[z][y][x]) / dz; 

   else if (z == (zmax - 1))
      (*G)[2] = ((double) F[z][y][x] - (double) F[z - 1][y][x]) / dz; 

/******************************************************************
* For points where the function is not defined its not a good
* idea to return zero for the gradient. This is because a gradient
* of zero might indicate a maxima/minima point for the function.
******************************************************************/

   else 
      throwError("gradient(): The gradient is not defined at this point."); 

}

#pragma instantiate void gradientTemplate(Array1D<double> *G, Array3D<float> const &F,         int x, int y, int z, double dx, double dy, double dz)
#pragma instantiate void gradientTemplate(Array1D<double> *G, Array3D<unsigned char> const &F, int x, int y, int z, double dx, double dy, double dz)

#pragma instantiate double trilinearTemplate(const Array3D<unsigned char>&, double, double, double, unsigned char)
#pragma instantiate double trilinearTemplate(const Array3D<float>&,         double, double, double, float)

#pragma instantiate void   trilinearTemplate(Array1D<double>*, const Array3D<float>&, const Array3D<float>&, const Array3D<float>&, double, double, double, float)
                                                                                                                                             

/////////////////////////////////////////////////////////////////////////////////////////
// Trilinearly interpolate into an unsigned char Array3D
/////////////////////////////////////////////////////////////////////////////////////////
unsigned char
TransformationUtils::trilinear( const Array3D<unsigned char> &X,
                                double z,
                                double y,
                                double x,
                                unsigned char bkgrnd )
{
  // Round to nearest unsigned char
  return (unsigned char) ( trilinearTemplate( X, z, y, x, bkgrnd) + 0.5);
}


/////////////////////////////////////////////////////////////////////////////////////////
// Trilinearly interpolate into an float Array3D
/////////////////////////////////////////////////////////////////////////////////////////
float
TransformationUtils::trilinear( const Array3D<float> &X,
                                double z,
                                double y,
                                double x,
                                float bkgrnd )
{
  return (float) trilinearTemplate( X, z, y, x, bkgrnd);
}


/////////////////////////////////////////////////////////////////////////////////////////
// Trilinear function that trilinearly interpolates 3 Array3d<float>'s at once.
/////////////////////////////////////////////////////////////////////////////////////////
void
TransformationUtils::trilinear( Array1D<double> *vx,
                                const Array3D<float> &X,
                                const Array3D<float> &Y, 
                                const Array3D<float> &Z,
                                double z, 
                                double y,
                                double x,
                                float bkgrnd )
{
  // Call template trilinear that takes in three floats
  trilinearTemplate( vx,
                     X,
                     Y,
                     Z,
                     z, 
                     y,
                     x,
                     bkgrnd);
  
}


int
TransformationUtils::interpolate( const Array3D<unsigned char> & input,
                                        Array3D<unsigned char> & output,
                                  const Array3D<float>         & hfield)
{
  if ( input.getXsize() != output.getXsize() ||
       input.getYsize() != output.getYsize() ||
       input.getZsize() != output.getZsize() )
  {
    cerr << "Input and output sizes must be the same for this routine." << endl;
    return 1;
  }

  if ( (input.getXsize()*3) != hfield.getXsize() ||
       input.getYsize() != hfield.getYsize() ||
       input.getZsize() != hfield.getZsize() )
  {
    cerr << "Input and h-field sizes must be the same for this routine." << endl;
    return 2;
  }

  unsigned char * outputDataPtr = output.data();
  const float   * hfieldDataPtr = hfield.data();
  int numX = input.getXsize();
  int numY = input.getYsize();
  int numZ = input.getZsize();

  int index;
  int total = numX * numY * numZ;

  for (index = 0; index < total; index++)
  {
    *outputDataPtr = TransformationUtils::trilinear( input,
                                                     (double) hfieldDataPtr[2],
                                                     (double) hfieldDataPtr[1],
                                                     (double) hfieldDataPtr[0],
                                                     (unsigned char) 0);
    outputDataPtr += 1;
    hfieldDataPtr += 3;
  }

  return 0;
}



int
TransformationUtils::interpolate(float * field_at_xyz_Array,
				 const Array3D<float> & fieldInterlaced, // nz * ny * nx*3
				 double z,
				 double y,
				 double x,
				 int isHFieldFlag)
{
#if 0
  if (NULL == field_at_xyz_Array)
  {
    cerr << "field_at_xyz_Array not allocated" << endl;
    return 1;
  }

  if ( 0 != (fieldInterlaced.getXsize() % 3) )
  {
    cerr << "fieldInterlaced must be an interlaced field" << endl;
    return 2;
  }
#endif

  /*****************************************************************
   * (0 <= fx, fy, fz <= 1), fractional position between data points.
   *****************************************************************/
  int x0 = (int) x;
  int y0 = (int) y;
  int z0 = (int) z;
  
  if (z < 0)
    z0--;
  if (y < 0)
    y0--;
  if (x < 0)
    x0--;

  //  int z0 = floor(z), y0 = floor(y), x0 = floor(x); 
  int z1 = z0 + 1, y1 = y0 + 1, x1 = x0 + 1; 
  float fz = z - z0, fy = y - y0, fx = x - x0; 

  int xsize = fieldInterlaced.getXsize() / 3;
  int ysize = fieldInterlaced.getYsize();
  int zsize = fieldInterlaced.getZsize(); 

  // These variables hold the 8 nearest neighbors to the point x, y, z in the X, Y, and
  // Z dimension
  float fieldX000, fieldX001, fieldX010, fieldX011, fieldX100, fieldX101, fieldX110, fieldX111; 
  float fieldY000, fieldY001, fieldY010, fieldY011, fieldY100, fieldY101, fieldY110, fieldY111; 
  float fieldZ000, fieldZ001, fieldZ010, fieldZ011, fieldZ100, fieldZ101, fieldZ110, fieldZ111; 

  const float *const *const *fieldPtrPtrPtr = fieldInterlaced.address();

  int compute = 0;

  /**********************************************
   * Tri-linear interpolation for interior points.
   **********************************************/
  if ( x0 >= 0 && x1 < xsize && 
       y0 >= 0 && y1 < ysize && 
       z0 >= 0 && z1 < zsize) {
    
    // Use pointer arithmetic.
    float const *fieldPtr;
    
    /****************************
     * Data used in interpolation.
     ****************************/
    fieldPtr = &fieldPtrPtrPtr[z0][y0][x0 * 3]; 
    fieldX000 = fieldPtr[0]; fieldX001 = fieldPtr[3]; 
    fieldY000 = fieldPtr[1]; fieldY001 = fieldPtr[4]; 
    fieldZ000 = fieldPtr[2]; fieldZ001 = fieldPtr[5]; 

    fieldPtr = &fieldPtrPtrPtr[z0][ y1 ][x0 * 3]; 
    fieldX010 = fieldPtr[0]; fieldX011 = fieldPtr[3]; 
    fieldY010 = fieldPtr[1]; fieldY011 = fieldPtr[4]; 
    fieldZ010 = fieldPtr[2]; fieldZ011 = fieldPtr[5]; 

    fieldPtr = &fieldPtrPtrPtr[ z1 ][ y0 ][x0 * 3]; 
    fieldX100 = fieldPtr[0]; fieldX101 = fieldPtr[3]; 
    fieldY100 = fieldPtr[1]; fieldY101 = fieldPtr[4]; 
    fieldZ100 = fieldPtr[2]; fieldZ101 = fieldPtr[5]; 

    fieldPtr = &fieldPtrPtrPtr[ z1 ][ y1 ][x0 * 3]; 
    fieldX110 = fieldPtr[0]; fieldX111 = fieldPtr[3]; 
    fieldY110 = fieldPtr[1]; fieldY111 = fieldPtr[4]; 
    fieldZ110 = fieldPtr[2]; fieldZ111 = fieldPtr[5]; 

    compute = 1;

  } else {

    /**********************************************
     * Tri-linear interpolation for boundary points.
     * Points outside of Array3D are set to bkgrnd.
     **********************************************/
    compute = 0;

    // set these to zero because we will only change the ones that fall 
    // inside the image
    if (isHFieldFlag) {

#if 1
    fieldX000 = x; fieldX001 = x; 
    fieldY000 = y; fieldY001 = y; 
    fieldZ000 = z; fieldZ001 = z; 
                                
    fieldX010 = x; fieldX011 = x; 
    fieldY010 = y; fieldY011 = y; 
    fieldZ010 = z; fieldZ011 = z; 
                                
    fieldX100 = x; fieldX101 = x; 
    fieldY100 = y; fieldY101 = y; 
    fieldZ100 = z; fieldZ101 = z; 
                                
    fieldX110 = x; fieldX111 = x; 
    fieldY110 = y; fieldY111 = y; 
    fieldZ110 = z; fieldZ111 = z; 
#else
    fieldX000 = x; fieldX001 = x+1; // is this the right way?!
    fieldY000 = y; fieldY001 = y; 
    fieldZ000 = z; fieldZ001 = z; 
                                
    fieldX010 = x;   fieldX011 = x+1; 
    fieldY010 = y+1; fieldY011 = y+1; 
    fieldZ010 = z;   fieldZ011 = z; 
                                
    fieldX100 = x;   fieldX101 = x+1; 
    fieldY100 = y;   fieldY101 = y; 
    fieldZ100 = z+1; fieldZ101 = z+1; 
                                
    fieldX110 = x;   fieldX111 = x+1; 
    fieldY110 = y+1; fieldY111 = y+1; 
    fieldZ110 = z+1; fieldZ111 = z+1; 
#endif
    } else {

      fieldX000 = fieldX001 = 0; 
      fieldY000 = fieldY001 = 0; 
      fieldZ000 = fieldZ001 = 0; 
      
      fieldX010 = fieldX011 = 0; 
      fieldY010 = fieldY011 = 0; 
      fieldZ010 = fieldZ011 = 0; 
      
      fieldX100 = fieldX101 = 0; 
      fieldY100 = fieldY101 = 0; 
      fieldZ100 = fieldZ101 = 0; 
    
      fieldX110 = fieldX111 = 0; 
      fieldY110 = fieldY111 = 0; 
      fieldZ110 = fieldZ111 = 0; 
    }


    int inBoundsX0 = (x0 >= 0) && (x0 < xsize);
    int inBoundsY0 = (y0 >= 0) && (y0 < ysize);
    int inBoundsZ0 = (z0 >= 0) && (z0 < zsize);
    int inBoundsX1 = (x1 >= 0) && (x1 < xsize);
    int inBoundsY1 = (y1 >= 0) && (y1 < ysize);
    int inBoundsZ1 = (z1 >= 0) && (z1 < zsize);

    if (inBoundsZ0 && inBoundsY0 && inBoundsX0) {
      fieldX000 = fieldPtrPtrPtr[z0][y0][x0*3  ];
      fieldY000 = fieldPtrPtrPtr[z0][y0][x0*3+1];
      fieldZ000 = fieldPtrPtrPtr[z0][y0][x0*3+2];
      compute = 1;
    }
    if (inBoundsZ0 && inBoundsY0 && inBoundsX1) {
      fieldX001 = fieldPtrPtrPtr[z0][y0][x1*3  ];
      fieldY001 = fieldPtrPtrPtr[z0][y0][x1*3+1];
      fieldZ001 = fieldPtrPtrPtr[z0][y0][x1*3+2];
      compute = 1;
    }
    if (inBoundsZ0 && inBoundsY1 && inBoundsX0) {
      fieldX010 = fieldPtrPtrPtr[z0][y1][x0*3  ];
      fieldY010 = fieldPtrPtrPtr[z0][y1][x0*3+1];
      fieldZ010 = fieldPtrPtrPtr[z0][y1][x0*3+2];
      compute = 1;
    }
    if (inBoundsZ0 && inBoundsY1 && inBoundsX1) {
      fieldX011 = fieldPtrPtrPtr[z0][y1][x1*3  ];
      fieldY011 = fieldPtrPtrPtr[z0][y1][x1*3+1];
      fieldZ011 = fieldPtrPtrPtr[z0][y1][x1*3+2];
      compute = 1;
    }
    if (inBoundsZ1 && inBoundsY0 && inBoundsX0) {
      fieldX100 = fieldPtrPtrPtr[z1][y0][x0*3  ];
      fieldY100 = fieldPtrPtrPtr[z1][y0][x0*3+1];
      fieldZ100 = fieldPtrPtrPtr[z1][y0][x0*3+2];
      compute = 1;
    }
    if (inBoundsZ1 && inBoundsY0 && inBoundsX1) {
      fieldX101 = fieldPtrPtrPtr[z1][y0][x1*3  ];
      fieldY101 = fieldPtrPtrPtr[z1][y0][x1*3+1];
      fieldZ101 = fieldPtrPtrPtr[z1][y0][x1*3+2];
      compute = 1;
    }
    if (inBoundsZ1 && inBoundsY1 && inBoundsX0) {
      fieldX110 = fieldPtrPtrPtr[z1][y1][x0*3  ];
      fieldY110 = fieldPtrPtrPtr[z1][y1][x0*3+1];
      fieldZ110 = fieldPtrPtrPtr[z1][y1][x0*3+2];
      compute = 1;
    }
    if (inBoundsZ1 && inBoundsY1 && inBoundsX1) {
      fieldX111 = fieldPtrPtrPtr[z1][y1][x1*3  ];
      fieldY111 = fieldPtrPtrPtr[z1][y1][x1*3+1];
      fieldZ111 = fieldPtrPtrPtr[z1][y1][x1*3+2];
      compute = 1;
    }
  }

  if (compute) {
    /********************
     * Interpolation, (1).
     ********************/
  
    float fieldX00x = fieldX000 + fx*(fieldX001-fieldX000);
    float fieldX01x = fieldX010 + fx*(fieldX011-fieldX010);
    float fieldX10x = fieldX100 + fx*(fieldX101-fieldX100);
    float fieldX11x = fieldX110 + fx*(fieldX111-fieldX110);

    float fieldX0yx = fieldX00x + fy*(fieldX01x-fieldX00x);
    float fieldX1yx = fieldX10x + fy*(fieldX11x-fieldX10x);
  
    field_at_xyz_Array[0] = fieldX0yx + fz*(fieldX1yx-fieldX0yx);

    float fieldY00x = fieldY000 + fx*(fieldY001-fieldY000);
    float fieldY01x = fieldY010 + fx*(fieldY011-fieldY010);
    float fieldY10x = fieldY100 + fx*(fieldY101-fieldY100);
    float fieldY11x = fieldY110 + fx*(fieldY111-fieldY110);

    float fieldY0yx = fieldY00x + fy*(fieldY01x-fieldY00x);
    float fieldY1yx = fieldY10x + fy*(fieldY11x-fieldY10x);
  
    field_at_xyz_Array[1] = fieldY0yx + fz*(fieldY1yx-fieldY0yx);

    float fieldZ00x = fieldZ000 + fx*(fieldZ001-fieldZ000);
    float fieldZ01x = fieldZ010 + fx*(fieldZ011-fieldZ010);
    float fieldZ10x = fieldZ100 + fx*(fieldZ101-fieldZ100);
    float fieldZ11x = fieldZ110 + fx*(fieldZ111-fieldZ110);

    float fieldZ0yx = fieldZ00x + fy*(fieldZ01x-fieldZ00x);
    float fieldZ1yx = fieldZ10x + fy*(fieldZ11x-fieldZ10x);
  
    field_at_xyz_Array[2] = fieldZ0yx + fz*(fieldZ1yx-fieldZ0yx);

  }
  else
  {
    if (isHFieldFlag) {
      field_at_xyz_Array[0] = x;
      field_at_xyz_Array[1] = y;
      field_at_xyz_Array[2] = z;
    } else { // assume zero when outside boundaries
      field_at_xyz_Array[0] = 0;
      field_at_xyz_Array[1] = 0;
      field_at_xyz_Array[2] = 0;
    }
  }

  return 0;
}


void
TransformationUtils::gradient(       Array1D<double> * GradientPtr, 
                               const Array3D<float>  & Function,
                                     int x, int y, int z, 
                                     double dx, double dy, double dz)
{
  gradientTemplate( GradientPtr, Function, x, y, z, dx, dy, dz);
}

void 
TransformationUtils::gradient(        Array1D<double>        * GradientPtr, 
                               const  Array3D<unsigned char> & Function,
                                      int x, int y, int z, 
                                      double dx, double dy, double dz)
{
  gradientTemplate( GradientPtr, Function, x, y, z, dx, dy, dz);
}





/*****************************************************************************************
* INTERPFT, 3-D interpolation using a FFT method. Y = INTERPFT(X, zsize, ysize, xsize) 
* returns a 3D array Y of dimensions (zsize, ysize, xsize) obtained by interpolation 
* in the Fourier transform of X.
*
* Abed Hammoud
* 24Jan1997-IntellX
* Copyright 1993-1997
*****************************************************************************************/
// void
// TransformationUtils::resample(unsigned char const *volin,
//                               int x1, 
//                               int y1,
//                               int z1,
//                               unsigned char *volout,
//                               int x2,
//                               int y2,
//                               int z2)
void TransformationUtils::resample(unsigned short const *I, 
				   int ix,
				   int iy,
				   int iz,
				   unsigned short *O,
				   int ox,
				   int oy,
				   int oz){

 
/***********************************
* Dimensions of input output arrays.
***********************************/

//    int iz = I.getZsize(), iy = I.getYsize(), ix = I.getXsize(); 
//    int oz = O.getZsize(), oy = O.getYsize(), ox = O.getXsize(); 

/***************
* Special cases.
***************/

//    if (oz == 0 || oy == 0 || ox == 0) {
//       O = Array3D<unsigned short>(); 
//       return O; 
//    }

/*************************************************************************
* If necessary, make (oz) to be a multiple of (iz), i.e, make (oz) > (iz).
* Also, make (oy) to be a multiple of (iy), i.e, make (oy) > (iy),
* and make (ox) to be a multiple of (ix), i.e make (ox) > (ix).
* This enables us to do decimation by downsampling in the spatial domain.
*************************************************************************/
 

   int zinc = 1; 
   if (oz <= iz) {
      zinc = (iz / oz) + 1; 
      oz *= zinc; 		// New (oz).
   }   

   int yinc = 1; 
   if (oy <= iy) {
      yinc = (iy / oy) + 1; 
      oy *= yinc; 		// New (oy).
   }   

   int xinc = 1; 
   if (ox <= ix) {
      xinc = (ix / ox) + 1; 
      ox *= xinc; 		// New (ox).
   }   

/***************************************************************
* Work space array. The array is made big enough to hold the +ve
* frequency half of the FFT of the real output sequence.
***************************************************************/

   int oxc = 2 * ((ox + 2) / 2);

   Array3D<float> W(oz, oy, oxc); 

/******************
* Clear work space.
******************/

   W = 0.0f; 

/****************************
* Make a copy of input array.
****************************/

   int i, j, k; 
   int xstride = oxc - ix; 
   int ystride = oxc * (oy - iy); 
   int Xstride = 2*((ix+2)/2)-ix;

   //   float *wptr = W.data(); 

   //   unsigned short const *iptr = I.data(); 
   unsigned short const *iptr = I; 

   //   for (i = iz; i; i--, wptr += ystride)
   //     for (j = iy; j; j--, wptr += xstride)
   //       for (k = ix; k; k--)
   //	      *wptr++ = *iptr++; 


   FFT_real *wptr=new FFT_real[iz*iy*(ix+Xstride)];
   FFT_real *temp=wptr;
   for (i = iz; i; i--)
      for (j = iy; j; j--, wptr += Xstride)
	 for (k = ix; k; k--)
	    *wptr++ = *iptr++; 
   wptr=temp;

/*****************************
* Forward real to complex FFT.
*****************************/

   //   float *coeff = fftpack::fft3di(ix, iy, iz, (float *) 0); 
   //   float *coeff = ::scfft3dui(ix, iy, iz, (float *) 0); 
   //   fftpack::rcfft3d(-1, ix, iy, iz, W.data(), oxc, oy, coeff); 
   //   ::scfft3du(-1, ix, iy, iz, W.data(), oxc, oy, coeff); 
   //   free(coeff);
   ITX_FFT fftF;
   fftF.initialize(3, ITX_FFT::R2C,ix,iy,iz);
   fftF.fftForward(wptr);

   temp=(FFT_real *)W.data();
   FFT_real *temp2=wptr;
   for (i = 0; i<iz; i++)
     for (j = 0; j<iy; j++)
       for (k=0; k<ix+Xstride; k++)
	 {
	   W[i][j][k] = *wptr++;
	 }
   wptr=temp2;
   delete [] wptr;

/******************************
* Pad with zeros in the middle.
******************************/

// Nyquist frequiencis.
   int znyqst = iz / 2; 		
   int ynyqst = iy / 2; 		

// Move 3 blocks of the FFT around and insert zeros.
   int ip, jp; 
   int ixc = 2 * ((ix + 2) / 2); 
   for (i = 0; i <= znyqst; i++) {
      for (j = iy - 1, jp = oy - 1; j > ynyqst; j--, jp--) {
	 for (k = 0; k < ixc; k++) {
	    W[i][jp][k] = W[i][j][k]; 
	    W[i][j][k] = 0.0f; 
	 }
      }
   }
   for (i = iz - 1, ip = oz - 1; i > znyqst; i--, ip--) {
      for (j = iy - 1, jp = oy - 1; j > ynyqst; j--, jp--) {
	 for (k = 0; k < ixc; k++) {
	    W[ip][jp][k] = W[i][j][k]; 
	    W[i][j][k] = 0.0f; 
	 }
      }
      for (j = 0; j <= ynyqst; j++) {
	 for (k = 0; k < ixc; k++) {
	    W[ip][j][k] = W[i][j][k]; 
	    W[i][j][k] = 0.0f; 
	 }
      }
   }

/*******************************************************
* For even number of elements in the z and y dimensions.
*******************************************************/

// Z = Znyqst plane.
   if (!(iz % 2)) {
      ip = oz - znyqst; 
      for (j = 0; j < iy; j++) {
	 for (k = 0; k < ixc; k++) {
	    W[znyqst][j][k] /= 2.0f; 
	    W[ip][j][k] = W[znyqst][j][k]; 
	 }
      }
   }

// Y = Ynyqst plane.
   if (!(iy % 2)) {
      jp = oy - ynyqst; 
      for (i = 0; i< iz; i++) {
	 for (k = 0; k < ixc; k++) {
	    W[i][ynyqst][k] /= 2.0f; 
	    W[i][jp][k] = W[i][ynyqst][k];
	 }
      }
   }

/*****************************
* Complex to real inverse FFT.
*****************************/

   //   coeff = fftpack::fft3di(ox, oy, oz, (float *) 0);
   //   coeff = ::scfft3dui(ox, oy, oz, (float *) 0);
   //   fftpack::crfft3d(1, ox, oy, oz, W.data(), oxc, oy, coeff);
   //   ::csfft3du(1, ox, oy, oz, W.data(), oxc, oy, coeff);
   //   free(coeff);
   ITX_FFT fftI;
   fftI.initialize(3,ITX_FFT::C2R,ox,oy,oz);
   fftI.fftInverse((FFT_complex *)W.data());


/**********************************************
* Copy into output Array. Skip extra points.
* Also, normalize by length of input sequence.
**********************************************/
//    unsigned short *optr = O.data(); 
   unsigned short *optr = O; 
   float alpha = 1.0f / (iz * iy * ix); 
   for (k = 0; k < oz; k += zinc)
      for (j = 0; j < oy; j += yinc)
	 for (i = 0; i < ox; i += xinc)
	    *optr++ = (unsigned short) ((W[k][j][i] > 0.0f) ? 
				(alpha * W[k][j][i]) : 0.0f); 

/******************
* Return new array.
******************/

//   return O; 

}


/*****************************************************************************************
* INTERPFT, 3-D interpolation using a FFT method. Y = INTERPFT(X, zsize, ysize, xsize) 
* returns a 3D array Y of dimensions (zsize, ysize, xsize) obtained by interpolation 
* in the Fourier transform of X.
*
* Abed Hammoud
* 24Jan1997-IntellX
* Copyright 1993-1997
*****************************************************************************************/
//Array3D<unsigned char>& interpft(Array3D<unsigned char> const &I, Array3D<unsigned char> &O) {

   void TransformationUtils::resample(unsigned char const *I,
				      int ix,
				      int iy,
				      int iz,
				      unsigned char *O,
				      int ox,
				      int oy,
				      int oz) {


/***********************************
* Dimensions of input output arrays.
***********************************/

//    int iz = I.getZsize(), iy = I.getYsize(), ix = I.getXsize(); 
//    int oz = O.getZsize(), oy = O.getYsize(), ox = O.getXsize(); 

/***************
* Special cases.
***************/

//    if (oz == 0 || oy == 0 || ox == 0) {
//       O = Array3D<unsigned char>(); 
//       return O; 
//    }

/*************************************************************************
* If necessary, make (oz) to be a multiple of (iz), i.e, make (oz) > (iz).
* Also, make (oy) to be a multiple of (iy), i.e, make (oy) > (iy),
* and make (ox) to be a multiple of (ix), i.e make (ox) > (ix).
* This enables us to do decimation by downsampling in the spatial domain.
*************************************************************************/

   int zinc = 1; 
   if (oz <= iz) {
      zinc = (iz / oz) + 1; 
      oz *= zinc; 		// New (oz).
   }   

   int yinc = 1; 
   if (oy <= iy) {
      yinc = (iy / oy) + 1; 
      oy *= yinc; 		// New (oy).
   }   

   int xinc = 1; 
   if (ox <= ix) {
      xinc = (ix / ox) + 1; 
      ox *= xinc; 		// New (ox).
   }   

/***************************************************************
* Work space array. The array is made big enough to hold the +ve
* frequency half of the FFT of the real output sequence.
**************************************************************ls
 */

   int oxc = 2 * ((ox + 2) / 2); 
   Array3D<float> W(oz, oy, oxc); 

/******************
* Clear work space.
******************/

   W = 0.0f; 

/****************************
* Make a copy of input array.
****************************/

   int i, j, k; 
   int xstride = oxc - ix; 
   int ystride = oxc * (oy - iy); 
   int Xstride = 2*((ix+2)/2)-ix;

   //   float *wptr = W.data(); 

   //   unsigned char const *iptr = I.data(); 
   unsigned char const *iptr = I; 

   //   for (i = iz; i; i--, wptr += ystride)
   //     for (j = iy; j; j--, wptr += xstride)
   //        for (k = ix; k; k--)
   //	       *wptr++ = *iptr++; 



   FFT_real *wptr=new FFT_real[iz*iy*Xstride];
   FFT_real *temp=wptr;
   for (i = iz; i; i--)
      for (j = iy; j; j--, wptr += Xstride)
	 for (k = ix; k; k--)
	    *wptr++ = *iptr++; 
   wptr=temp;


/*****************************
* Forward real to complex FFT.
*****************************/

   //   float *coeff = fftpack::fft3di(ix, iy, iz, (float *) 0); 
   //   float *coeff = ::scfft3dui(ix, iy, iz, (float *) 0); 
   //   fftpack::rcfft3d(-1, ix, iy, iz, W.data(), oxc, oy, coeff); 
   //   ::scfft3du(-1, ix, iy, iz, W.data(), oxc, oy, coeff); 
   //   free(coeff);
   
   ITX_FFT fftF;
   fftF.initialize(3,ITX_FFT::R2C,ix,iy,ix);
   fftF.fftForward(wptr);

   temp=(FFT_real *)W.data();
   FFT_real *temp2=wptr;
   for (i = 0 ; i<iz;i++)
     for (j = 0; j<iy; j++)
       for (k = 0;k<ix+Xstride; k++)
	 {
	   W[i][j][k] = *wptr++; 
	 }
   wptr=temp2;
   delete [] wptr;

/******************************
* Pad with zeros in the middle.
******************************/

// Nyquist frequiencis.
   int znyqst = iz / 2; 		
   int ynyqst = iy / 2; 		

// Move 3 blocks of the FFT around and insert zeros.
   int ip, jp; 
   int ixc = 2 * ((ix + 2) / 2); 
   for (i = 0; i <= znyqst; i++) {
      for (j = iy - 1, jp = oy - 1; j > ynyqst; j--, jp--) {
	 for (k = 0; k < ixc; k++) {
	    W[i][jp][k] = W[i][j][k]; 
	    W[i][j][k] = 0.0f; 
	 }
      }
   }
   for (i = iz - 1, ip = oz - 1; i > znyqst; i--, ip--) {
      for (j = iy - 1, jp = oy - 1; j > ynyqst; j--, jp--) {
	 for (k = 0; k < ixc; k++) {
	    W[ip][jp][k] = W[i][j][k]; 
	    W[i][j][k] = 0.0f; 
	 }
      }
      for (j = 0; j <= ynyqst; j++) {
	 for (k = 0; k < ixc; k++) {
	    W[ip][j][k] = W[i][j][k]; 
	    W[i][j][k] = 0.0f; 
	 }
      }
   }

/*******************************************************
* For even number of elements in the z and y dimensions.
*******************************************************/

// Z = Znyqst plane.
   if (!(iz % 2)) {
      ip = oz - znyqst; 
      for (j = 0; j < iy; j++) {
	 for (k = 0; k < ixc; k++) {
	    W[znyqst][j][k] /= 2.0f; 
	    W[ip][j][k] = W[znyqst][j][k]; 
	 }
      }
   }

// Y = Ynyqst plane.
   if (!(iy % 2)) {
      jp = oy - ynyqst; 
      for (i = 0; i < iz; i++) {
	 for (k = 0; k < ixc; k++) {
	    W[i][ynyqst][k] /= 2.0f; 
	    W[i][jp][k] = W[i][ynyqst][k];
	 }
      }
   }

/*****************************
* Complex to real inverse FFT.
*****************************/

   //   coeff = fftpack::fft3di(ox, oy, oz, (float *) 0);
   //   coeff = ::scfft3dui(ox, oy, oz, (float *) 0);
   //   fftpack::crfft3d(1, ox, oy, oz, W.data(), oxc, oy, coeff);
   //   ::csfft3du(1, ox, oy, oz, W.data(), oxc, oy, coeff);
   //   free(coeff);
   ITX_FFT fftI;
   fftI.initialize(3,ITX_FFT::C2R,ox,oy,oz);
   fftI.fftInverse((FFT_complex *)W.data());


/**********************************************
* Copy into output Array. Skip extra points.
* Also, normalize by length of input sequence.
**********************************************/

   float tmp; 
//   unsigned char *optr = O.data(); 
   unsigned char *optr = O; 
   float alpha = 1.0f / (iz * iy * ix); 
   for (k = 0; k < oz; k += zinc) {
      for (j = 0; j < oy; j += yinc) {
	 for (i = 0; i < ox; i += xinc) {
	    tmp = alpha * W[k][j][i]; 
	    if (tmp < 0.0f) tmp = 0.0f; 
	    if (tmp > 255.0f) tmp = 255.0f; 
	    *optr++ = (unsigned char) tmp; 
	 }
      }
   }

/******************
* Return new array.
******************/

//   return O; 

}

#if 0

/* Resample an input volume to the size of the size specified
   Note the input size and the output size can be of any legnth!!
   Use the FFT to upsample and downsamle the volume!!
   Assumes the volumes are unsigned char!!
   USAGE: resample involume xin yin zin outvolume xout yout zout
*/

// void
// TransformationUtils::resample(unsigned char const *volin,
//                               int x1, 
//                               int y1,
//                               int z1,
//                               unsigned char *volout,
//                               int x2,
//                               int y2,
//                               int z2)
// {
//   float tmp1;
//   int i,j,k;
//   float alpha;
  
//   // Dimensions of workvolume1
//   int z3 = z1;
//   int y3 = y1;
//   int x3 = x1+2;
//   int xy3 = x3*y3;

//   // Dimensions of workvolume2
//   int z4 = z2;
//   int y4 = y2;
//   int x4 = x2+2;
//   int xy4 = x4*y4;
  
//   // create work space for fft
//   float *workvolume1 = (float *)calloc(z3*y3*x3,sizeof(float));
//   float *coef1 = (float *)calloc(((x1+15)+2*(y1+15)+2*(z1+15)),sizeof(float));

//   // copy input volume into workvolume with padding for fft
//   unsigned char const *cptr = volin;
//   float *fptr = workvolume1;
//   for(i=0; i<z1; i++) {
//     for(j=0; j<y1; j++) {
//       for(k=0; k<x1; k++) {
// 	*fptr++ = (float) *cptr++;
//       }
//       fptr++;
//       fptr++;
//     }
//   }	
  
//   cout << z1 << " " << y1 << " " << x1 << " " << endl;
//   cout << z2 << " " << y2 << " " << x2 << " " << endl;

//   // compute the fft
//   sfft3dui(x1, y1, z1, coef1);
//   sfft3du(-1, x1, y1, z1, workvolume1, x3, y3, coef1);
  
//   // create work space for fft
//   float *workvolume2 = (float *)calloc(z4*y4*x4,sizeof(float));
//   float *coef2 = (float *)calloc(((x2+15)+2*(y2+15)+2*(z2+15)),sizeof(float));
  
//   /* Now fill the first part of the array by the transform */
//   int zm = (z1 < z2) ? z1 : z2;
//   int ym = (y1 < y2) ? y1 : y2;
//   int xm = (x1 < x2) ? x1 : x2;
  
//   alpha = 1.0/sqrt((float)(z1*y1*x1));

//   int k40, k30, j40, j30;
//   int k41, k31, j41, j31;
//   for (k=0; k<zm/2; k++) {
//     k40 = k*xy4;
//     k30 = k*xy3;
//     k41 = (z2-1-k)*xy4;
//     k31 = (z1-1-k)*xy3;
//     for (j=0; j<ym/2; j++) {
//       j40 = j*x4;
//       j30 = j*x3;
//       j41 = (y2-1-j)*x4;
//       j31 = (y1-1-j)*x3;
//       for(i=0; i<2*((xm+2)/2); i++){
//   	workvolume2[i+j40+k40] = alpha*workvolume1[i+j30+k30];
//   	workvolume2[i+j41+k40] = alpha*workvolume1[i+j31+k30];
//   	workvolume2[i+j40+k41] = alpha*workvolume1[i+j30+k31];
//   	workvolume2[i+j41+k41] = alpha*workvolume1[i+j31+k31];
//       }
//     }
//   }

//   sfft3dui(x2, y2, z2, coef2);
//   sfft3du(1, x2, y2, z2, workvolume2, x4, y4, coef2);

//   // scale and copy workvolume2 into output volume
//   fptr = workvolume2;
//   unsigned char *cptr0 = volout;
//   for (i=0; i<z2; i++) {
//     for (j=0; j<y2; j++) {
//       for (k=0; k<x2; k++) {
// 	tmp1 = alpha*(*fptr++);
// 	if (tmp1 < 0 ) tmp1 = 0;
// 	if (tmp1 > 255 ) tmp1 = 255;
// 	*cptr0++ = (unsigned char) tmp1;
//       }
//       fptr++;
//       fptr++;
//     }
//   }

//   free(workvolume1);
//   free(workvolume2);
//   free(coef1);
//   free(coef2);
// }

// void
// TransformationUtils::resample(unsigned short const *volin,
//                               int x1, 
//                               int y1,
//                               int z1,
//                               unsigned short *volout,
//                               int x2,
//                               int y2,
//                               int z2)
// {
//   float tmp1;
//   int i,j,k;
//   float alpha;
  
//   // Dimensions of workvolume1
//   int z3 = z1;
//   int y3 = y1;
//   int x3 = x1+2;
//   int xy3 = x3*y3;

//   // Dimensions of workvolume2
//   int z4 = z2;
//   int y4 = y2;
//   int x4 = x2+2;
//   int xy4 = x4*y4;
  
//   // create work space for fft
//   float *workvolume1 = (float *)calloc(z3*y3*x3,sizeof(float));
//   float *coef1 = (float *)calloc(((x1+15)+2*(y1+15)+2*(z1+15)),sizeof(float));

//   // copy input volume into workvolume with padding for fft
//   unsigned short const *cptr = volin;
//   float *fptr = workvolume1;
//   for(i=0; i<z1; i++) {
//     for(j=0; j<y1; j++) {
//       for(k=0; k<x1; k++) {
// 	*fptr++ = (float) *cptr++;
//       }
//       fptr++;
//       fptr++;
//     }
//   }	
  
//   cout << z1 << " " << y1 << " " << x1 << " " << endl;
//   cout << z2 << " " << y2 << " " << x2 << " " << endl;

//   // compute the fft
//   sfft3dui(x1, y1, z1, coef1);
//   sfft3du(-1, x1, y1, z1, workvolume1, x3, y3, coef1);
  
//   // create work space for fft
//   float *workvolume2 = (float *)calloc(z4*y4*x4,sizeof(float));
//   float *coef2 = (float *)calloc(((x2+15)+2*(y2+15)+2*(z2+15)),sizeof(float));
  
//   /* Now fill the first part of the array by the transform */
//   int zm = (z1 < z2) ? z1 : z2;
//   int ym = (y1 < y2) ? y1 : y2;
//   int xm = (x1 < x2) ? x1 : x2;
  
//   alpha = 1.0/sqrt((float)(z1*y1*x1));

//   int k40, k30, j40, j30;
//   int k41, k31, j41, j31;
//   for (k=0; k<zm/2; k++) {
//     k40 = k*xy4;
//     k30 = k*xy3;
//     k41 = (z2-1-k)*xy4;
//     k31 = (z1-1-k)*xy3;
//     for (j=0; j<ym/2; j++) {
//       j40 = j*x4;
//       j30 = j*x3;
//       j41 = (y2-1-j)*x4;
//       j31 = (y1-1-j)*x3;
//       for(i=0; i<2*((xm+2)/2); i++){
//   	workvolume2[i+j40+k40] = alpha*workvolume1[i+j30+k30];
//   	workvolume2[i+j41+k40] = alpha*workvolume1[i+j31+k30];
//   	workvolume2[i+j40+k41] = alpha*workvolume1[i+j30+k31];
//   	workvolume2[i+j41+k41] = alpha*workvolume1[i+j31+k31];
//       }
//     }
//   }

//   sfft3dui(x2, y2, z2, coef2);
//   sfft3du(1, x2, y2, z2, workvolume2, x4, y4, coef2);

//   // scale and copy workvolume2 into output volume
//   fptr = workvolume2;
//   unsigned short *cptr0 = volout;
//   for (i=0; i<z2; i++) {
//     for (j=0; j<y2; j++) {
//       for (k=0; k<x2; k++) {
// 	tmp1 = alpha*(*fptr++);
// 	if (tmp1 < 0 ) tmp1 = 0;
// 	if (tmp1 > 255 ) tmp1 = 255;
// 	*cptr0++ = (unsigned short) tmp1;
//       }
//       fptr++;
//       fptr++;
//     }
//   }

//   free(workvolume1);
//   free(workvolume2);
//   free(coef1);
//   free(coef2);
// }

/////////////////////////////////////////////////////////////////////////////
// 
// upsampleField
// 
// Takes 3 field arrays (down{X,Y,Z}Field) (h or u it doesn't matter) and upsamples
// them into up{X,Y,Z}Field using trilinear interpolation
//
// Note: this operation can just as easily "downsample", it doesn't care
// 
// Returns 2 if downXField's size is 0x0x0
// Returns 3 if upXField's size is 0x0x0
//
// Assumes that the y and z's field sizes are the same as the x's
// Assumes that fields are in pixels coordinates and not mm
//
/////////////////////////////////////////////////////////////////////////////
int
TransformationUtils::upsampleFields( const Array3D<float> & downXField,
                                     const Array3D<float> & downYField,
                                     const Array3D<float> & downZField,
                                           Array3D<float> & upXField,
                                           Array3D<float> & upYField,
                                           Array3D<float> & upZField )
{

  double downSizeX = downXField.getXsize();
  double downSizeY = downXField.getYsize();
  double downSizeZ = downXField.getZsize();
  cout << "Original: " << downSizeX << "x" << downSizeY << "x" << downSizeZ << endl;

  if (0 == downSizeX || 0 == downSizeY || 0 == downSizeZ )
  {
    cerr << "At least one dimension for down displacement field is zero" << endl;
    return 2;
  }
  
  double upsampledSizeX = upXField.getXsize();
  double upsampledSizeY = upXField.getYsize();
  double upsampledSizeZ = upXField.getZsize();
  cout << "Upsampled: " << upsampledSizeX << "x" << upsampledSizeY << "x" << upsampledSizeZ << endl;

  if (0 == upsampledSizeX || 0 == upsampledSizeY || 0 == upsampledSizeZ )
  {
    cerr << "At least one dimension for the upsampled displacement field is zero" << endl;
    cerr << "I don't believe that these arrays have been initialized correctly" << endl;
    return 3;
  }
  
  // If these are the same size then just copy them
  if ( downSizeX == upsampledSizeX &&
       downSizeY == upsampledSizeY &&
       downSizeZ == upsampledSizeZ)
  {
    cout << "Sizes are the same so just copying" << endl;
    // Not much work needs to be done here, they are the same dimensions
    upXField = downXField;
    upYField = downYField;
    upZField = downZField;
  }
  else
  {
    // Things are a little more complicated here, we need to upsample (or
    // downsample is the case may be) so that up{X,Y,Z}Field have the
    // correct dimensions.  We also need to scale them accordingly as these
    // fields are in pixel coordinates not mm
    Array1D<double> U(3); 
    double ScratchX = 0;
    double ScratchY = 0;
    double ScratchZ = 0;

    double downToUpsampledScaleX = upsampledSizeX / downSizeX; // Upsampledx = downx * downToUpsampledScaleX;
    double downToUpsampledScaleY = upsampledSizeY / downSizeY; 
    double downToUpsampledScaleZ = upsampledSizeZ / downSizeZ; 
    double upsampledToDownScaleX =  1.0 / downToUpsampledScaleX; // downx = upsampledx * upsampledToDownScaleX;
    double upsampledToDownScaleY =  1.0 / downToUpsampledScaleY;
    double upsampledToDownScaleZ =  1.0 / downToUpsampledScaleZ;

    // Loop through each point in the upsampled image, find the corresponding 
    // point in down image
    for (int z = 0; z < upsampledSizeZ; z++)
    {
      // Scale so we are in the down object's coordinate system
      ScratchZ = (double) z * upsampledToDownScaleZ;

      for (int y = 0; y < upsampledSizeY; y++)
      {
        ScratchY = (double) y * upsampledToDownScaleY;

        for (int x = 0; x < upsampledSizeX; x++)
        {
          ScratchX = (double) x * upsampledToDownScaleX;

          // Cannot use _nx, _ny, _nz as they have not been saved with the transformation file!
          TransformationUtils::trilinear( &U,
                                          downXField, downYField, downZField,
                                          ScratchZ, ScratchY, ScratchX,
                                          (float) 0.0);  // i don't like using 0 here, 0 is not a correct
                                                         // background value, but then what is?

          upXField[z][y][x] = U[0] * downToUpsampledScaleX;
          upYField[z][y][x] = U[1] * downToUpsampledScaleY;
          upZField[z][y][x] = U[2] * downToUpsampledScaleZ;
        }
      }
    }
  }

  int x, y, z;
  cout << "Calculating min, max, and sum of new fields" << endl;
  Array1D<double> minup(3), maxup(3), sumup(3);
  minup[0] = maxup[0] = upXField[0][0][0];
  minup[1] = maxup[1] = upYField[0][0][0];
  minup[2] = maxup[2] = upZField[0][0][0];
  sumup = (double) 0;
  for (z = 0; z < upsampledSizeZ; z++)
  {
    for (y = 0; y < upsampledSizeY; y++)
    {
      for (x = 0; x < upsampledSizeX; x++)
      {
        if (upXField[z][y][x] < minup[0]) minup[0] = upXField[z][y][x];
        if (upYField[z][y][x] < minup[1]) minup[1] = upYField[z][y][x];
        if (upZField[z][y][x] < minup[2]) minup[2] = upZField[z][y][x];

        if (upXField[z][y][x] > maxup[0]) maxup[0] = upXField[z][y][x];
        if (upYField[z][y][x] > maxup[1]) maxup[1] = upYField[z][y][x];
        if (upZField[z][y][x] > maxup[2]) maxup[2] = upZField[z][y][x];

        sumup[0] += upXField[z][y][x];
        sumup[1] += upYField[z][y][x];
        sumup[2] += upZField[z][y][x];
      }
    }
  }
  cout << "mins " << minup[0] << " " << minup[1] << " " << minup[2] << endl;
  cout << "maxs " << maxup[0] << " " << maxup[1] << " " << maxup[2] << endl;
  cout << "sums " << sumup[0] << " " << sumup[1] << " " << sumup[2] << endl;

  cout << "Calculating min, max, and sum of old fields" << endl;
  Array1D<double> mindown(3), maxdown(3), sumdown(3);
  mindown[0] = maxdown[0] = upXField[0][0][0];
  mindown[1] = maxdown[1] = upYField[0][0][0];
  mindown[2] = maxdown[2] = upZField[0][0][0];
  sumdown = (double) 0;
  for (z = 0; z < downSizeZ; z++)
  {
    for (y = 0; y < downSizeY; y++)
    {
      for (x = 0; x < downSizeX; x++)
      {
        if (downXField[z][y][x] < mindown[0]) mindown[0] = downXField[z][y][x];
        if (downYField[z][y][x] < mindown[1]) mindown[1] = downYField[z][y][x];
        if (downZField[z][y][x] < mindown[2]) mindown[2] = downZField[z][y][x];

        if (downXField[z][y][x] > maxdown[0]) maxdown[0] = downXField[z][y][x];
        if (downYField[z][y][x] > maxdown[1]) maxdown[1] = downYField[z][y][x];
        if (downZField[z][y][x] > maxdown[2]) maxdown[2] = downZField[z][y][x];

        sumdown[0] += downXField[z][y][x];
        sumdown[1] += downYField[z][y][x];
        sumdown[2] += downZField[z][y][x];
      }
    }
  }
  cout << "mins " << mindown[0] << " " << mindown[1] << " " << mindown[2] << endl;
  cout << "maxs " << maxdown[0] << " " << maxdown[1] << " " << maxdown[2] << endl;
  cout << "sums " << sumdown[0] << " " << sumdown[1] << " " << sumdown[2] << endl;


  return 0;
}
#endif

/////////////////////////////////////////////////////////////////////////////
//
// Reads in the transformation type from the file filename
// and returns it in the array type.
//
/////////////////////////////////////////////////////////////////////////////
int
TransformationUtils::readTransformationType( const char * filename,
                                                   char * type)
{
  const char * const functionName= "readTransformationType";
  
  ifstream inStream (filename, ios::in);
  
  if (! inStream )
  {
    cerr << "Could not open file: \"" << filename << "\"" << endl;
    return 2;
  }

  int ReadInClassNameLength;
  inStream.read((char *) &ReadInClassNameLength, sizeof (ReadInClassNameLength));
  ReadInClassNameLength = big_endian_int(ReadInClassNameLength);
  // Let's not get crazy, even _I_ don't make class names this big
  if (ReadInClassNameLength <= 0 || ReadInClassNameLength > 1000) { 
    cerr << "I don't think that this is a transformation file." << endl;
    cerr << "Apparent class length is " << ReadInClassNameLength << endl;
    // Throw error here!
    return 1;
  }
  
  // Alloc memory for the class name and read it in
  char * ReadInClassName = new char[ReadInClassNameLength + 1];
  if ( ! ReadInClassName ) {
    cerr << functionName << ": ";
    cerr << "Failed to alloc " << ReadInClassNameLength+1 << " bytes." << endl;
    // Throw error here!
    return 2;
  }
  inStream.read((char *) ReadInClassName, ReadInClassNameLength);
  ReadInClassName[ReadInClassNameLength] = '\0'; // End-of-string so I can use strcmp()
  
  // Copy into the memory that the caller points to
  strcpy(type, ReadInClassName);

  // Close file
  inStream.close();

  // Free temporary memory
  delete ReadInClassName;

  return 0;
}



int
TransformationUtils::upsampleInterlacedField( const Array3D<float> & downField,
                                                    Array3D<float> & upField)
{

  double downSizeX = downField.getXsize() / 3;
  double downSizeY = downField.getYsize();
  double downSizeZ = downField.getZsize();
  cout << "Original: " << downSizeX << "x" << downSizeY << "x" << downSizeZ << endl;

  if (0 == downSizeX || 0 == downSizeY || 0 == downSizeZ )
  {
    cerr << "At least one dimension for down displacement field is zero" << endl;
    return 2;
  }
  
  double upsampledSizeX = upField.getXsize() / 3;
  double upsampledSizeY = upField.getYsize();
  double upsampledSizeZ = upField.getZsize();
  cout << "Upsampled: " << upsampledSizeX << "x" << upsampledSizeY << "x" << upsampledSizeZ << endl;

  if (0 == upsampledSizeX || 0 == upsampledSizeY || 0 == upsampledSizeZ )
  {
    cerr << "At least one dimension for the upsampled displacement field is zero" << endl;
    cerr << "I don't believe that these arrays have been initialized correctly" << endl;
    return 3;
  }
  
  // If these are the same size then just copy them
  if ( downSizeX == upsampledSizeX && downSizeY == upsampledSizeY && downSizeZ == upsampledSizeZ)
  {
    cout << "Sizes are the same so just copying" << endl;
    // Not much work needs to be done here, they are the same dimensions
    upField = downField;
  }
  else
  {
    // Things are a little more complicated here, we need to upsample (or
    // downsample is the case may be) so that up{X,Y,Z}Field have the
    // correct dimensions.  We also need to scale them accordingly as these
    // fields are in pixel coordinates not mm
    double downToUpsampledScaleX = upsampledSizeX / downSizeX; // Upsampledx = downx * downToUpsampledScaleX;
    double downToUpsampledScaleY = upsampledSizeY / downSizeY; 
    double downToUpsampledScaleZ = upsampledSizeZ / downSizeZ; 
    double upsampledToDownScaleX =  1.0 / downToUpsampledScaleX; // downx = upsampledx * upsampledToDownScaleX;
    double upsampledToDownScaleY =  1.0 / downToUpsampledScaleY;
    double upsampledToDownScaleZ =  1.0 / downToUpsampledScaleZ;

    // Loop through each point in the upsampled image, find the corresponding 
    // point in down image
    float interpolatedValueArray[3];
    double newX = 0;
    double newY = 0;
    double newZ = 0;

    for (int z = 0; z < upsampledSizeZ; z++)
    {
      // Scale so we are in the down object's coordinate system
      newZ = (double) z * upsampledToDownScaleZ;

      for (int y = 0; y < upsampledSizeY; y++)
      {
        newY = (double) y * upsampledToDownScaleY;

        for (int x = 0; x < upsampledSizeX; x++)
        {
          newX = (double) x * upsampledToDownScaleX;

          // Cannot use _nx, _ny, _nz as they have not been saved with the transformation file!
          TransformationUtils::interpolate( interpolatedValueArray,
                                            downField,
                                            (double) newZ,
                                            (double) newY,
                                            (double) newX);

          upField[z][y][x*3    ] = interpolatedValueArray[0] * downToUpsampledScaleX;
          upField[z][y][x*3 + 1] = interpolatedValueArray[1] * downToUpsampledScaleY;
          upField[z][y][x*3 + 2] = interpolatedValueArray[2] * downToUpsampledScaleZ;
        }
      }
    }
  }

#if 1
  int x, y, z;
  cout << "Calculating min, max, and sum of new fields" << endl;
  Array1D<double> minup(3), maxup(3), sumup(3);
  minup[0] = maxup[0] = upField[0][0][0];
  minup[1] = maxup[1] = upField[0][0][1];
  minup[2] = maxup[2] = upField[0][0][2];
  sumup = (double) 0;
  for (z = 0; z < upsampledSizeZ; z++)
  {
    for (y = 0; y < upsampledSizeY; y++)
    {
      for (x = 0; x < upsampledSizeX; x++)
      {
        if (upField[z][y][x*3+0] < minup[0]) minup[0] = upField[z][y][x*3+0];
        if (upField[z][y][x*3+1] < minup[1]) minup[1] = upField[z][y][x*3+1];
        if (upField[z][y][x*3+2] < minup[2]) minup[2] = upField[z][y][x*3+2];
              
        if (upField[z][y][x*3+0] > maxup[0]) maxup[0] = upField[z][y][x*3+0];
        if (upField[z][y][x*3+1] > maxup[1]) maxup[1] = upField[z][y][x*3+1];
        if (upField[z][y][x*3+2] > maxup[2]) maxup[2] = upField[z][y][x*3+2];

        sumup[0] += upField[z][y][x*3+0];
        sumup[1] += upField[z][y][x*3+1];
        sumup[2] += upField[z][y][x*3+2];
      }
    }
  }
  cout << "mins " << minup[0] << " " << minup[1] << " " << minup[2] << endl;
  cout << "maxs " << maxup[0] << " " << maxup[1] << " " << maxup[2] << endl;
  cout << "sums " << sumup[0] << " " << sumup[1] << " " << sumup[2] << endl;

  cout << "Calculating min, max, and sum of old fields" << endl;
  Array1D<double> mindown(3), maxdown(3), sumdown(3);
  mindown[0] = maxdown[0] = upField[0][0][0];
  mindown[1] = maxdown[1] = upField[0][0][1];
  mindown[2] = maxdown[2] = upField[0][0][2];
  sumdown = (double) 0;
  for (z = 0; z < downSizeZ; z++)
  {
    for (y = 0; y < downSizeY; y++)
    {
      for (x = 0; x < downSizeX; x++)
      {
        if (downField[z][y][x*3+0] < mindown[0]) mindown[0] = downField[z][y][x*3+0];
        if (downField[z][y][x*3+1] < mindown[1]) mindown[1] = downField[z][y][x*3+1];
        if (downField[z][y][x*3+2] < mindown[2]) mindown[2] = downField[z][y][x*3+2];
                                                               
        if (downField[z][y][x*3+0] > maxdown[0]) maxdown[0] = downField[z][y][x*3+0];
        if (downField[z][y][x*3+1] > maxdown[1]) maxdown[1] = downField[z][y][x*3+1];
        if (downField[z][y][x*3+2] > maxdown[2]) maxdown[2] = downField[z][y][x*3+2];

        sumdown[0] += downField[z][y][x*3+0];
        sumdown[1] += downField[z][y][x*3+1];
        sumdown[2] += downField[z][y][x*3+2];
      }
    }
  }
  cout << "mins " << mindown[0] << " " << mindown[1] << " " << mindown[2] << endl;
  cout << "maxs " << maxdown[0] << " " << maxdown[1] << " " << maxdown[2] << endl;
  cout << "sums " << sumdown[0] << " " << sumdown[1] << " " << sumdown[2] << endl;
#endif

  return 0;
}

int
TransformationUtils::loadVolume(const char *filename,
                                GreyVolume<unsigned short> &volume)
{
  cout << "Loading " << filename << "...";

  AnalyzeVol tmpVol;

  if (!tmpVol.load(filename))
    return 0;

  int xsize = tmpVol.getXsize();
  int ysize = tmpVol.getYsize();
  int zsize = tmpVol.getZsize();

  volume.setDim(zsize, ysize, xsize);
  volume.setVoxelXScale(tmpVol.getVoxelXScale());
  volume.setVoxelYScale(tmpVol.getVoxelYScale());
  volume.setVoxelZScale(tmpVol.getVoxelZScale());

  // transpose volume
  int n = volume.getNelm();
  unsigned short *inData = tmpVol.data();
  unsigned short *outData = volume.data();

  int i;
  for (i=0; i<n; i++)
  {
    outData[n-1-i] = inData[i];
  }

  return 1;
}


int
TransformationUtils::loadVolume(const char *filename,
                                GreyVolume<unsigned char> &volume)
{
  cout << "Loading " << filename << "...";

  AnalyzeVol tmpVol;

  if (!tmpVol.load(filename))
    return 0;
  int xsize = tmpVol.getXsize();
  int ysize = tmpVol.getYsize();
  int zsize = tmpVol.getZsize();

  volume.setDim(zsize, ysize, xsize);
  volume.setVoxelXScale(tmpVol.getVoxelXScale());
  volume.setVoxelYScale(tmpVol.getVoxelYScale());
  volume.setVoxelZScale(tmpVol.getVoxelZScale());

  // transpose volume
  int n = volume.getNelm();
  unsigned short *inData = tmpVol.data();
  unsigned char *outData = volume.data();

  int i;
  // find maxVal to normalize volume
  unsigned short maxVal = 0;
  for (i=0; i<n; i++)
    if (inData[i] > maxVal)
      maxVal = inData[i];

  float intensityFactor = 1.0;
  if (maxVal > 255)
  {
    intensityFactor = 255.0/maxVal;  // scale intensity values to 0..255
  }
  cout << "maxVal: " << (int) maxVal << " intensityFactor: " << intensityFactor;

  // normalize and transpose volume
  float x;
  for (i=0; i<n; i++)
  {
    x = inData[i]*intensityFactor;
    outData[n-1-i] = (unsigned char) x;  // transpose volume
  }

  return 1;
}


int
TransformationUtils::loadLandmarks(const char *filename,
				   Array2D<double> &landmarks)
{
  cout << "Loading " << filename << "...";

  ifstream file(filename, ios::in);
  if (!file)
    return 0;

  char tmp[40];
  file >> tmp;
  if (strcmp(tmp, "Landmarks-1.0") != 0) {
    cerr << "Incorrect landmark file format" << endl;
    return 0;
  }

  int numLandmarks;
  file >> numLandmarks >> ws;

  Array2D<double> activeLandmarks;
  activeLandmarks.setDim(numLandmarks, 3);

  float tmpVariance;
  double x, y, z;
  int isActive;
  int numActiveLandmarks = 0;
  int i, j;
  for (i=0; i<numLandmarks; i++) {
    file.get(tmp, 13, '\n');
    for (j=1; tmp[j]!='\"'; j++);
    tmp[j] = 0;
    file >> x;
    file >> y;
    file >> z;
    file >> isActive;
    file >> tmpVariance >> ws; // should be real variance
    if (isActive) {
      activeLandmarks[numActiveLandmarks][0] = x;
      activeLandmarks[numActiveLandmarks][1] = y;
      activeLandmarks[numActiveLandmarks][2] = z;
      numActiveLandmarks++;
    } else {
      //      _tempPoints[i].makeInactive();
      //      XmToggleButtonSetState(_pntAtlas[i], False, False);
    }
  }

  // copy active landmarks into output landmarks array
  landmarks.setDim(numActiveLandmarks, 3);
  for (i=0; i<numActiveLandmarks; i++) {
    landmarks[i][0] = activeLandmarks[i][0];
    landmarks[i][1] = activeLandmarks[i][1];
    landmarks[i][2] = activeLandmarks[i][2];
  }

#if 0
  // read in landmark lines (currently not in use)
  int n;
  int m, j;
  Pnt p;
  file >> n >> ws;
  for (i=0; i<n; i++) {
    file.get(tmp, 13, '\n');
    for (j=1; tmp[j]!='\"'; j++);
    tmp[j] = 0;
    XmTextFieldSetString(_lineLabel[i], &tmp[1]);
    file >> m;
    for (j=0; j<m; j++) {
      file >> p.x();
      file >> p.y();
      file >> p.z();
      _tempOrigLines[i].push_back(p);
    }
    file >> ws;
    if (m > 0) {
      XmToggleButtonSetState(_lineAtlas[i], True, False);
      lineResample(_tempOrigLines[i], _tempLines[i]);
    } else
      XmToggleButtonSetState(_lineAtlas[i], False, False);
  }
#endif

  return 1;
}

////////////////////////////////////////////////////////////////////////////////
//
// Function: compactdisplay()
//
// Purpose: displays vector in one line with trailing zeros only if required
//
// Inputs: vector to display
//
// Outputs: none
//
////////////////////////////////////////////////////////////////////////////////
void
TransformationUtils::compactdisplay(Vector<double> vec)
{
  long f = cout.flags(ios::floatfield);
  cout << "[";
  for (int i=0; i<vec.getNelm() - 1; i++)
  {
    cout << vec[i] << " ";
  }
  
  cout << vec[ vec.getNelm() -1 ] << "]";
  //cout.flags(f);
}
