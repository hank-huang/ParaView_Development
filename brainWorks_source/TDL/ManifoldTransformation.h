#ifndef __MANIFOLDTRANSFORMATION_H__
#define __MANIFOLDTRANSFORMATION_H__

//////////////////////////////////////////////////////////////////////////
//
// ManifoldTransformation class
// 
// This is the class that performs the manifold deformation
//
// Note: If you are changing this class.  The version number of this class
// should change is more or less data is saved so that we may be able to
// load in old versions.
//
// To Do:
//
// Things to think about...
//   1) Maybe put in code that can read in older version of the class
//
//////////////////////////////////////////////////////////////////////////
//
// Revision Log:
//
// $Log: ManifoldTransformation.h,v $
// Revision 1.1  2004/11/15 04:44:08  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:28:10  kem
// Initial revision
//
// Revision 1.12  1999/07/09 17:49:43  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.11  1998/12/18 18:14:42  RAZORDB
// Fix for 7.2.1 compiler Warning 1356
//
// Revision 1.10  1998/12/10 19:03:46  RAZORDB
// IDL/AW merge
//
// Revision 1.9  1998/05/18 22:20:32  rst
// make transformPoint() public
//
// Revision 1.8  1998/04/02 17:12:48  kem
// Add ManifoldFluidTransformation as a friend
//
// Revision 1.7  1997/12/12 19:27:39  abed
// Adding resolution
//
// Revision 1.6  1997/09/17 21:58:30  kem
// Add atlas and patient dimension parameters for initialize()
//
// Revision 1.5  1997/09/16 19:47:58  kem
// Remove getUField() method
//
// Revision 1.4  1997/08/25 17:30:35  kem
// Add getHfield()
//
// Revision 1.3  1997/08/07 16:08:58  kem
// Change ElasticManifoldTransformation to ManifoldElasticTransformation
//
// Revision 1.2  1997/08/05 18:45:15  abed
// Added Accessor functions and made ElasticManifoldTransformation friend
//
// Revision 1.1  1997/08/01 19:49:23  csb
// Initial revision
//
// Revision 1.1  1997/08/01 15:51:45  csb
// Initial revision
//
//
//////////////////////////////////////////////////////////////////////////
#include <OS.h>
#include <string.h>
#include <ADL/Array1D.h>
#include <ADL/Array2D.h>
#include <ADL/Array3D.h>
#include <ADL/Matrix.h>
#include <TDL/FieldTransformation.h>

class ManifoldFluidTransformation;

class ManifoldTransformation : public FieldTransformation {

  // Friends.
  friend class ManifoldFluidTransformation; 

public:

  // Ctor and Dtor.
  ManifoldTransformation(); 
  ~ManifoldTransformation(); 
    
  int  save( const char * filename );  // virtual in Transformation class
  int  load( const char * filename );  // virtual in Transformation class

  int initialize(const Array2D<double> &atlasPoints,
		 const Array2D<double> &patientPoints,
		 const Array1D<double> &variancePoints,
		 const int atlasSizeX, const int atlasSizeY, const int atlasSizeZ,
		 const int patientSizeX, const int patientSizeY, const int patientSizeZ, 
		 const double atlasResX = 1.0, 
		 const double atlasResY = 1.0, 
		 const double atlasResZ = 1.0, 
		 const double patientResX = 1.0, 
		 const double patientResY = 1.0, 
		 const double patientResZ = 1.0);
    
  int calculate();             // virtual in Transformation class

  void print();

  int getNumPoints() const; 
  Matrix<double> getGreensCoeff() const; 
  Matrix<double> getAffineMatrix() const; 
  Array2D<double> getPatientPoints() const; 
  Vector<double> getTranslationVector() const; 

  // Calculate transformed coordinates
  void transformPoint(const double x0, const double x1, const double x2,
		      double &y0, double &y1, double &y2);

protected:

  // Returns 0 upon success, and non-zero otherwise
  int saveClassParameters ( ofstream & );
  int loadClassParameters ( ifstream & );
  ostream & printClassParameters ( ostream & );

private: 
  static const char * const mClassRevision;
  static const char * const mClassName;
  static const int          mClassVersion;
  void getHField(Array3D<float> *arr)
	{ getHField(arr,1.0,1.0,1.0,0.0,0.0,0.0); }

  void getHField(Array3D<float> *, double Rx, double Ry, 
	double Rz, double Ox, double Oy, double Oz); 

  double trilin(Array3D<unsigned char> const &volume, 
	       double z, double y, double x);

  // Private member variables
  Array2D<double> _atlasPoints;
  Array2D<double> _patientPoints;
  Array1D<double> _variancePoints;
  int _numPoints;

  double _affineMatrix00;
  double _affineMatrix01;
  double _affineMatrix02;
  double _affineMatrix10;
  double _affineMatrix11;
  double _affineMatrix12;
  double _affineMatrix20;
  double _affineMatrix21;
  double _affineMatrix22;
  Vector<double> _translationVector;
  Matrix<double> _greensCoeff;

};
#endif // __MANIFOLDTRANSFORMATION_H__
