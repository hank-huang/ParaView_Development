#ifndef __RIGIDLANDMARKTRANSFORMATION_H__
#define __RIGIDLANDMARKTRANSFORMATION_H__

//////////////////////////////////////////////////////////////////////////
//
//
//////////////////////////////////////////////////////////////////////////
//
// Revision Log:
//
// $Log: RigidLandmarkTransformation.h,v $
// Revision 1.1  2004/11/15 04:44:08  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:28:23  kem
// Initial revision
//
// Revision 1.3  1999/07/09 17:49:59  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.2  1998/12/18 18:15:01  RAZORDB
// Fix for 7.2.1 compiler Warning 1356
//
// Revision 1.1  1997/08/01 19:49:26  csb
// Initial revision
//
// Revision 1.4  1997/06/25 17:37:52  kem
// Split file into .h and .C
//
//
//
//////////////////////////////////////////////////////////////////////////

#include <OS.h>
#include <stream.h>
#include <TDL/LinearTransformation.h>

class RigidLandmarkTransformation : public LinearTransformation {
public:

  RigidLandmarkTransformation() {
    _atlasPoints = NULL; 
    _patientPoints = NULL;
    _scaleMatrix.setDim(3, 3);
    _scaleMatrix.eye();
  };
    
  ~RigidLandmarkTransformation() {};
    
  int save( const char * filename );                  // virtual in Transformation class
  int load( const char * filename );                  // virtual in Transformation class
  int calculate();             // virtual in Transformation class
  void print();

  int initialize(const Matrix<double> * const atlasPoints,
		 const Matrix<double> * const patientPoints,
		 const float &atlasVoxelScaleX = 1.0,
		 const float &atlasVoxelScaleY = 1.0,
		 const float &atlasVoxelScaleZ = 1.0,
		 const float &patientVoxelScaleX = 1.0,
		 const float &patientVoxelScaleY = 1.0,
		 const float &patientVoxelScaleZ = 1.0); // non-virtual

protected:
  // Returns 0 upon success, and non-zero otherwise
  int saveClassParameters ( ofstream & );
  int loadClassParameters ( ifstream & );
  ostream & printClassParameters ( ostream & );

private:
  static const char * const mClassRevision;
  static const char * const mClassName;
  static const int          mClassVersion;

  const Matrix<double> * _atlasPoints;
  const Matrix<double> * _patientPoints;
  Matrix<double> _scaleMatrix;
};

#endif // __RIGIDLANDMARKTRANSFORMATION_H__
