#ifndef __LINEARTRANSFORMATION_H__
#define __LINEARTRANSFORMATION_H__

//////////////////////////////////////////////////////////////////////////
//
// LinearTransformation class
// 
// This is the base class for derived transformation classes that
// define a transformation with a matrix multiply and a translation vector.
//
// Note: If you are changing this class.  The version number of this class
// should change is more or less data is saved so that we may be able to
// load in old versions.
//
// To Do:
//   1) deal with failed call to loadTransformationClassParameters
//
// Things to think about...
//   1) Maybe put in code that can read in older version of the class
//
//////////////////////////////////////////////////////////////////////////
//
// Revision Log:
//
// $Log: LinearTransformation.h,v $
// Revision 1.1  2004/11/15 04:44:08  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:28:11  kem
// Initial revision
//
// Revision 1.12  1999/07/09 17:49:26  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.11  1999/01/07 16:25:13  RAZORDB
// Minor fixes (rst)
//
// Revision 1.10  1998/12/18 18:14:32  RAZORDB
// Fix for 7.2.1 compiler Warning 1356
//
// Revision 1.9  1998/10/08 20:34:12  csb
// Add support for storing transformations in millimeters
//
// Revision 1.8  1998/07/23 22:41:53  csb
// Add getTransformMatrix() accessor function
//
// Revision 1.7  1998/04/17 19:48:25  rst
// add stub for Surface apply with resolution
//
// Revision 1.6  1998/04/08 21:08:06  kem
// Change apply() for volumes
//
// Revision 1.5  1997/12/23 20:39:48  csb
// return value in apply
//
// Revision 1.4  1997/12/12 23:16:06  csb
// change prototype for getHField() and apply()
//
// Revision 1.3  1997/10/02 20:37:45  rst
// Fixed apply() routines for Surfaces and Points
//
// Revision 1.2  1997/09/11 22:08:15  rst
// added surface apply
//
// Revision 1.1  1997/08/01 19:49:22  csb
// Initial revision
//
// Revision 1.4  1997/08/01 15:34:19  csb
// Incorportate Abed's, fluid, elastic and field
//
// Revision 1.3  1997/06/25 22:52:12  kem
// Add apply() method for a list of points
//
// Revision 1.2  1997/06/25 17:37:46  kem
// Split file into .h and .C
//
//
//
//////////////////////////////////////////////////////////////////////////
#include <OS.h>
#include <stream.h>
#include <TDL/Transformation.h>
#include <ADL/Matrix.h>
#include <ADL/Vector.h>
#include <TDL/Surface.h>

class LinearTransformation : public Transformation {
public:
  LinearTransformation();
  virtual ~LinearTransformation() {} // virtual destructor
    
  void setInMMFlag(int flag);
  int getInMMFlag();

  // Define the function apply() here as it the same for all derived classes

  // apply transformation to a volume
  int apply(const Array3D<unsigned char> &inVolume,
	    Array3D<unsigned char> *outVolume,
	    double = 1.0, double = 1.0, double = 1.0, 
	    double = 0.0, double = 0.0, double = 0.0, 
	    const int interpolateFlag = 1);

  int apply(const Array3D<unsigned char> &inputVolume, 
	    double inputResX,     // inputVolume resolution in mm/voxel
	    double inputResY,
	    double inputResZ,
	    double inputOffsetX,  // inputVolume offset relative to atlas
	    double inputOffsetY,  // in mm
	    double inputOffsetZ,
	    Array3D<unsigned char> *outputVolume,
	    double outputResX,    // outputVolume resolution in mm/voxel
	    double outputResY,
	    double outputResZ,
	    double outputOffsetX, // outputVolume offset relative to patient
	    double outputOffsetY, // in mm
	    double outputOffsetZ,
	    const int interpolateFlag) { return 0; }

  // apply transformation to a list of points
  int apply(const Matrix<double> &inPointList,
	    Matrix<double> &outPointList);
    
  // apply transformation to a surface
  int apply(const Surface &inSurface,
	    Surface *outSurface);

  // apply for surface with resolution (not implemented)
  int apply(const Surface &inSurface,
	    Surface *outSurface,
	    double outputResX,
	    double outputResY,
	    double outputResZ,
	    double outputOffsetX,
	    double outputOffsetY,
	    double outputOffsetZ)
  { return 0; }

  // apply transformation to a Point (not implemented)
  int apply(const Point &inPoint, Point *outPoint);

  //
  // apply inverse of transformation
  //

  // apply inverse to a volume
  virtual int invApply(const Array3D<unsigned char> &inputVolume, 
		       Array3D<unsigned char> *outputVolume,
		       double outputResX,    // outputVolume resolution in mm/voxel
		       double outputResY,
		       double outputResZ,
		       double outputOffsetX, // outputVolume offset relative to patient
		       double outputOffsetY, // in mm
		       double outputOffsetZ)
    { cerr << "invApply not implemented!!" << endl; return 1; };

  // apply inverse to surfaces, with resolution
  virtual int invApply (const Surface &inSurface,
			Surface *outSurface,
			double outputResX,    // desired resolution of outSurface
			double outputResY,    // in mm/voxel
			double outputResZ,
			double outputOffsetX, // desired offset of outSurface
			double outputOffsetY, // relative to patient volume, in mm
			double outputOffsetZ)
    { cerr << "invApply not implemented!!" << endl; return 1; };




  double meanSquaredError(const Matrix<double> &inPointList,
			  const Matrix<double> &outPointList);

  int save( const char * filename );                  // virtual in Transformation class
  int load( const char * filename );                  // virtual in Transformation class
  void print();

  int getTransformMatrix( Matrix<double> * matrixPtr );

  // initialize defined in derived class
  // calculate defined in derived class

  inline Matrix<double> &matrix()      { return(mAffineMatrix); }
  inline Vector<double> &translation() { return(mTranslationVector); }

protected:
  int       saveClassParameters( ofstream & );
  int       loadClassParameters( ifstream & );
  ostream & printClassParameters ( ostream & );

  Matrix<double> mAffineMatrix;      // A in y = Ax + b
  Vector<double> mTranslationVector; // b in y = Ax + b
  int mInMMFlag;
    

private:
  static const char * const mClassRevision;
  static const char * const mClassName;
  static const int          mClassVersion;
};
#endif // __LINEARTRANSFORMATION_H__
