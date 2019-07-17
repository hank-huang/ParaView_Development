#ifndef __AFFINELANDMARKTRANSFORMATION_H__
#define __AFFINELANDMARKTRANSFORMATION_H__

//////////////////////////////////////////////////////////////////////////
//
//
//////////////////////////////////////////////////////////////////////////
//
// Revision Log:
//
// $Log: AffineLandmarkTransformation.h,v $
// Revision 1.1  2004/11/15 04:44:07  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:28:15  kem
// Initial revision
//
// Revision 1.4  1999/09/21 17:06:49  RAZORDB
// TDL update
//
// Revision 1.3  1999/07/09 17:48:22  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.2  1998/12/18 18:14:15  RAZORDB
// Fix for 7.2.1 compiler Warning 1356
//
// Revision 1.1  1997/08/01 19:49:15  csb
// Initial revision
//
// Revision 1.2  1997/06/25 17:37:40  kem
// Split file into .h and .C
//
//
//
//////////////////////////////////////////////////////////////////////////
#include <OS.h>
#include <stream.h>
#include <TDL/LinearTransformation.h>

class AffineLandmarkTransformation : public LinearTransformation {
public:

  AffineLandmarkTransformation() {_atlasPoints = NULL; _patientPoints = NULL;};
    
  ~AffineLandmarkTransformation() {};
    
  int save( const char * filename );                  // virtual in Transformation class
  int load( const char * filename );                  // virtual in Transformation class
  int calculate();             // virtual in Transformation class
  void print();

  int initialize(const Matrix<double> * const atlasPoints,
		 const Matrix<double> * const patientPoints);            // not-virtual

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
};

#endif // __AFFINELANDMARKTRANSFORMATION_H__
