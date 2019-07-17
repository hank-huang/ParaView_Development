#ifndef __AFFINEMOMENTTRANSFORMATION_H__
#define __AFFINEMOMENTTRANSFORMATION_H__

//////////////////////////////////////////////////////////////////////////
//
//
//////////////////////////////////////////////////////////////////////////
//
// Revision Log:
//
// $Log: AffineMomentTransformation.h,v $
// Revision 1.1  2004/11/15 04:44:07  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:28:26  kem
// Initial revision
//
// Revision 1.3  1999/07/09 17:48:30  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.2  1998/12/18 18:14:17  RAZORDB
// Fix for 7.2.1 compiler Warning 1356
//
// Revision 1.1  1997/08/01 19:49:16  csb
// Initial revision
//
// Revision 1.1  1997/06/25 17:40:09  kem
// Initial revision
//
//
//
//////////////////////////////////////////////////////////////////////////
#include <OS.h>
#include <stream.h>
#include <TDL/LinearTransformation.h>

class AffineMomentTransformation : public LinearTransformation {
public:

  AffineMomentTransformation() {_atlasPtr = NULL; _patientPtr = NULL;};
    
  ~AffineMomentTransformation() {};
    
  int save( const char * filename );                  // virtual in Transformation class
  int load( const char * filename );                  // virtual in Transformation class
  int calculate();             // virtual in Transformation class
  void print();

  int initialize(const Array3D<unsigned char> * const atlasPtr,
		 const Array3D<unsigned char> * const patientPtr); // not-virtual

protected:
  // Returns 0 upon success, and non-zero otherwise
  int saveClassParameters ( ofstream & );
  int loadClassParameters ( ifstream & );
  ostream & printClassParameters ( ostream & );

private:
  static const char * const mClassName;
  static const int          mClassVersion;
  const Array3D<unsigned char> *_atlasPtr;
  const Array3D<unsigned char> *_patientPtr;
};

#endif // __AFFINEMOMENTTRANSFORMATION_H__
