#ifndef __USERDEFINEDLINEARTRANSFORMATION_H__
#define __USERDEFINEDLINEARTRANSFORMATION_H__

//////////////////////////////////////////////////////////////////////////
//
//
//////////////////////////////////////////////////////////////////////////
//
// Revision Log:
//
// $Log: UserDefinedLinearTransformation.h,v $
// Revision 1.1  2004/11/15 04:44:09  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:28:27  kem
// Initial revision
//
// Revision 1.3  1999/07/09 17:50:47  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.2  1998/12/18 18:15:24  RAZORDB
// Fix for 7.2.1 compiler Warning 1356
//
// Revision 1.1  1997/08/01 19:49:28  csb
// Initial revision
//
// Revision 1.1  1997/06/25 17:40:11  kem
// Initial revision
//
//
//
//////////////////////////////////////////////////////////////////////////

#include <OS.h>
#include <stream.h>
#include <TDL/LinearTransformation.h>

class UserDefinedLinearTransformation : public LinearTransformation {
public:

  UserDefinedLinearTransformation() {};
    
  ~UserDefinedLinearTransformation() {};
    
  int save( const char * filename );                  // virtual in Transformation class
  int load( const char * filename );                  // virtual in Transformation class
  int calculate();             // virtual in Transformation class
  void print();

  int initialize();            // not-virtual

  void setTransformation(const Matrix<double> &affineMatrix,
			 const Vector<double> &translationVector);

protected:
  // Returns 0 upon success, and non-zero otherwise
  int saveClassParameters ( ofstream & );
  int loadClassParameters ( ifstream & );
  ostream & printClassParameters ( ostream & );

private:
  static const char * const mClassName;
  static const int          mClassVersion;
};
#endif // __USERDEFINEDLINEARTRANSFORMATION_H__

