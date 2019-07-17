#ifndef __TRANSFORMATION_H__
#define __TRANSFORMATION_H__

//////////////////////////////////////////////////////////////////////////
//
// Transformation class
// 
// This is the base class for all of the derived transformation classes
//
// Note: If you are changing this class.  The version number of this class
// should change is more or less data is saved so that we may be able to
// load in old versions.
//
// To Do:
//   1) Handle errors in loading incorrect versions better
//   2) Don't know what the apply parameters are yet
//
// Things to think about...
//   1) look at error codes when saving?
//   2) should the save,load and print functions be called the same
//      name so the children can call them easily and call their
//      function the same thing?
//   3) Maybe put in code that can read in older version of the class
//
//////////////////////////////////////////////////////////////////////////
//
// Revision Log:
// 
// $Log: Transformation.h,v $
// Revision 1.1  2004/11/15 04:44:09  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:28:17  kem
// Initial revision
//
// Revision 1.9  1999/07/09 17:50:23  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.8  1999/01/07 16:25:15  RAZORDB
// Minor fixes (rst)
//
// Revision 1.7  1998/12/18 18:15:12  RAZORDB
// Fix for 7.2.1 compiler Warning 1356
//
// Revision 1.6  1998/12/10 19:04:04  RAZORDB
// IDL/AW merge
//
// Revision 1.5  1998/04/17 19:48:11  rst
// Adding apply for surfaces with resolution
//
// Revision 1.4  1998/04/08 21:01:16  kem
// Change virtual apply() prototype for volumes
//
// Revision 1.3  1997/12/12 23:19:03  csb
// change prototype for 2 apply()s
//
// Revision 1.2  1997/10/02 20:37:46  rst
// Fixed apply() routines for Surfaces and Points
//
// Revision 1.1  1997/08/01 19:49:27  csb
// Initial revision
//
// Revision 1.6  1997/08/01 15:28:58  csb
// Incorportate Abed's, fluid, elastic and field
//
// Revision 1.5  1997/06/25 17:35:39  csb
// Split into .h and .c
//
// Revision 1.4  1997/06/24 22:24:13  csb
// Added Revision, changed order of saving and loading
//
// Revision 1.3  1997/06/24 20:48:58  csb
// Changes
//
// Revision 1.2  1997/06/02 21:50:20  csb
// Changes
//
// Revision 1.1  1997/05/30 23:59:16  csb
// Initial revision
//
//
//////////////////////////////////////////////////////////////////////////
#include <OS.h>
#include <stream.h>
#include <ADL/Array3D.h>
#include <TDL/Point.h>
#include <TDL/Surface.h>
#include <TDL/IDLdefines.h>
#include <TDL/ItXWorkClass.h>

class Transformation : public ItXWorkClass {
public:
  Transformation() : mVerboseModeFlag(0),
                     mOkToCalculateFlag(0),
                     mOkToApplyFlag(0)          { /* empty */ };

  virtual ~Transformation() {}; // virtual destructor

  // Return success values as ints
  virtual int calculate() = 0;
  virtual int save( const char * filename ) = 0;
  virtual int load( const char * filename ) = 0;

  // apply function for volumes
  virtual int apply(const Array3D<unsigned char> &inVolume,
                    Array3D<unsigned char> *outVolume,
		    double = 1.0, double = 1.0, double = 1.0, 
		    double = 0.0, double = 0.0, double = 0.0, 
                    const int interpolateFlag = 1) = 0;

  virtual int apply(const Array3D<unsigned char> &inputVolume, 
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
		    const int interpolateFlag) = 0; 
#if 0
  virtual int apply(const Array3D<unsigned char> &inVolume, 
		Array3D<unsigned char> *outVolume, 
		double Rx, double Ry, double Rz, 
		double Rtx, double Rty, double Rtz, 	// Template resolution.
		double Rgx, double Rgy, double Rgz, 	// Segmentation resolution.
		double gOx, double gOy, double gOz, 	// Segmentation Offset w.r.t T
		double sOx, double sOy, double sOz, 	// output cube offset w.r.t S.
		const int interpolateFlag) = 0; 
#endif

  // apply function for surfaces
  virtual int apply (const Surface &inSurface,
	             Surface *outSurface) = 0;

  // apply for surfaces with resolution
  virtual int apply (const Surface &inSurface,
	             Surface *outSurface,
		     double outputResX,    // desired resolution of outSurface
		     double outputResY,    // in mm/voxel
		     double outputResZ,
		     double outputOffsetX, // desired offset of outSurface
		     double outputOffsetY, // relative to patient volume, in mm
		     double outputOffsetZ) = 0;

  // apply function for points
  virtual int apply (const Point &inPoint,
		     Point *outPoint) = 0;

  //
  // apply inverse routines
  //

  // apply inverse to a volume
  virtual int invApply(const Array3D<unsigned char> &inputVolume, 
		       Array3D<unsigned char> *outputVolume,
		       double outputResX,    // outputVolume resolution in mm/voxel
		       double outputResY,
		       double outputResZ,
		       double outputOffsetX, // outputVolume offset relative to patient
		       double outputOffsetY, // in mm
		       double outputOffsetZ) = 0;

  // apply inverse to surfaces, with resolution
  virtual int invApply (const Surface &inSurface,
			Surface *outSurface,
			double outputResX,    // desired resolution of outSurface
			double outputResY,    // in mm/voxel
			double outputResZ,
			double outputOffsetX, // desired offset of outSurface
			double outputOffsetY, // relative to patient volume, in mm
			double outputOffsetZ) = 0;



  virtual void print() = 0;

  void setVerboseMode( int flag ) { mVerboseModeFlag = flag; }
  int  getVerboseMode( void )     { return mVerboseModeFlag; }

protected:
  // Returns 0 upon success, and non-zero otherwise
  int saveClassParameters ( ofstream & );
  int loadClassParameters ( ifstream & );
  ostream & printClassParameters ( ostream & );

  // Check and set flags
  int isOkToCalculate() const { return mOkToCalculateFlag; }
  int isOkToApply() const     { return mOkToApplyFlag; }
  void setOkToCalculate(int val) { mOkToCalculateFlag = val; }
  void setOkToApply(int val)     { mOkToApplyFlag = val; }

private:
  int mVerboseModeFlag;         // Print out algorithmic info during calculate()
                                // This should show information about the algorithm
                                // but not slow down the calculation significantly
  int mOkToCalculateFlag;       // This is set when the necessary things are
                                // done so that calculate() may be called
  int mOkToApplyFlag;           // This is set when we are ready to apply a
                                // transformation
  static const char * const mClassRevision;
  static const char * const mClassName;
  static const int          mClassVersion;
}; // class Transformation

#endif // __TRANSFORMATION_H__

