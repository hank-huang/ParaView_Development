//////////////////////////////////////////////////////////////////////////
//
// ManifoldTransformation class
// 
// This file contains the function definitions.
//
//////////////////////////////////////////////////////////////////////////
//
// Revision Log:
//
// $Log: ManifoldTransformation.C,v $
// Revision 1.1  2004/11/15 04:44:08  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:29:11  kem
// Initial revision
//
// Revision 1.8  1999/07/09 17:52:56  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.7  1998/12/10 19:05:42  RAZORDB
// IDL/AW merge
//
// Revision 1.6  1997/12/12 19:31:54  abed
// Adding resolution
//
// Revision 1.5  1997/09/17 22:00:43  kem
// Initialize _atlasSize? and _patientSize? variables in initialize()
//
// Revision 1.4  1997/09/16 19:46:40  kem
// Remove getUField() method
//
// Revision 1.3  1997/08/25 17:27:10  kem
// Add getHfield()
//
// Revision 1.2  1997/08/05 18:44:47  abed
// Added Accessor functions and made ElasticManifoldTransformation friend
//
// Revision 1.1  1997/08/01 19:50:06  csb
// Initial revision
//
// Revision 1.1  1997/08/01 15:51:33  csb
// Initial revision
//
//
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
//
// Include Files
//
//////////////////////////////////////////////////////////////////////////////

#include <string.h>
#include <TDL/ManifoldTransformation.h>

const char * const ManifoldTransformation::mClassRevision =     "$Id: ManifoldTransformation.C,v 1.1 2004/11/15 04:44:08 joeh Exp $";
const char * const ManifoldTransformation::mClassName = "ManifoldTransformation";
const int ManifoldTransformation::mClassVersion = 1;


/***************
* Ctor and Dtor.
***************/

ManifoldTransformation::ManifoldTransformation()
{ 
}


ManifoldTransformation::~ManifoldTransformation() { /* empty */ }



//////////////////////////////////////////////////////////////////////////////
//
// Private Methods
//
//////////////////////////////////////////////////////////////////////////////
void ManifoldTransformation::getHField(Array3D<float> *hField, 
	double Rx, double Ry, double Rz, double Ox, double Oy, double Oz) {

  cout << "getHField(new) in ManifoldTransformation" << endl;

  // get affine transformation parameters
  double tx = _translationVector[0];
  double ty = _translationVector[1];
  double tz = _translationVector[2];

  // get size of displacement field
  int fieldSizeX = hField->getXsize()/3;
  int fieldSizeY = hField->getYsize();
  int fieldSizeZ = hField->getZsize();

  double xnorm; 
  double sumX, sumY, sumZ; 
  double xreal, yreal, zreal; 
  double xDiff, yDiff, zDiff; 

//  cout << getTranslationVector() << endl; 
//  cout << getAffineMatrix() << endl; 
//  cout << getGreensCoeff() << endl; 

  cout << "\nDeforming at resolution, (in mm/voxel): " << endl; 
  cout << Rx << '\t' << Ry << '\t' << Rz << endl; 

  // compute the displacement field
  float *hFieldPtr = hField->data();
  system("date");
  cerr << "slice     of " << fieldSizeZ-1 << "\010\010\010\010\010\010\010";
  for (int z=0; z<fieldSizeZ; z++) {
    zreal = z * Rz + Oz; 
    cerr << "\010\010\010" ;
    cerr.width(3);
    cerr << z ;

    for (int y=0; y<fieldSizeY; y++) {
        yreal = y * Ry + Oy; 

      for (int x=0; x<fieldSizeX; x++) {

	xreal = x * Rx + Ox; 

	sumX = _affineMatrix00*xreal + _affineMatrix01*yreal 
		+ _affineMatrix02*zreal + tx;
	sumY = _affineMatrix10*xreal + _affineMatrix11*yreal 
		+ _affineMatrix12*zreal + ty;
	sumZ = _affineMatrix20*xreal + _affineMatrix21*yreal 
		+ _affineMatrix22*zreal + tz;

	for (int i=0; i<_numPoints; i++) {
	  xDiff = xreal - _patientPoints[i][0];
	  yDiff = yreal - _patientPoints[i][1];
	  zDiff = zreal - _patientPoints[i][2];
	  xnorm = sqrt((xDiff*xDiff) + (yDiff*yDiff) + (zDiff*zDiff));
	  sumX += xnorm * _greensCoeff[i][0];
	  sumY += xnorm * _greensCoeff[i][1];
	  sumZ += xnorm * _greensCoeff[i][2];
	}

	*hFieldPtr++ = (xreal - sumX) / _atlasResX;	// Field in voxel.
	*hFieldPtr++ = (yreal - sumY) / _atlasResY;	// Field in voxel.
	*hFieldPtr++ = (zreal - sumZ) / _atlasResZ;	// Field in voxel.

      }
    }
  }

  cerr << endl;
  system("date");

}

//////////////////////////////////////////////////////////////////////////////
//
// Protected Methods
//
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
//
// Function: saveClassParameters
//
// Purpose: Saves parameters that are owned by this class to the supplied 
//          ofstream.
//          It is envisioned that it's own parameters are saves and then
//          those of it's parent.  Since at this time this is the 
//          the final derived class, save() calls this function.
//
// Inputs: ofstream that we want to save to 
//
// Outputs: 0 upon success
//          1 if parent's save returned non-zero
//
//////////////////////////////////////////////////////////////////////////////
int
ManifoldTransformation::saveClassParameters ( ofstream & outStream)
{
  // const char * functionName = "ManifoldTransformation::saveClassParameters";

  // Write out size of class name and the class name
  int swapint;
  const int nameLength = strlen(mClassName);
  swapint = big_endian_int(nameLength);
  outStream.write((const char *) & swapint, sizeof(nameLength));
  outStream.write(mClassName,                  nameLength);

  // Write out the class version
  const int temp  = mClassVersion;    // Because the compiler may not put these variables
  swapint = big_endian_int(temp);
  outStream.write((const char *) & swapint, sizeof (temp));

  // Write out size of class revision and the class revision
  const int revisionLength = strlen(mClassRevision);
  swapint = big_endian_int(revisionLength);
  outStream.write((const char *) & swapint, sizeof(revisionLength));
  outStream.write(mClassRevision, revisionLength);

  // Write out input parameters for this class
  _atlasPoints.save(outStream);
  _patientPoints.save(outStream);
  _variancePoints.save(outStream);
  swapint = big_endian_int(_numPoints);
  outStream.write((const char *) &swapint, sizeof(_numPoints));

  // Write out the transformation results
  double swapdouble;
  swapdouble = (big_endian_double(_affineMatrix00));
  outStream.write((const char *) &swapdouble, sizeof(_affineMatrix00));
  swapdouble = (big_endian_double(_affineMatrix01));
  outStream.write((const char *) &swapdouble, sizeof(_affineMatrix01));
  swapdouble = (big_endian_double(_affineMatrix02));
  outStream.write((const char *) &swapdouble, sizeof(_affineMatrix02));
  swapdouble = (big_endian_double(_affineMatrix10));
  outStream.write((const char *) &swapdouble, sizeof(_affineMatrix10));
  swapdouble = (big_endian_double(_affineMatrix11));
  outStream.write((const char *) &swapdouble, sizeof(_affineMatrix11));
  swapdouble = (big_endian_double(_affineMatrix12));
  outStream.write((const char *) &swapdouble, sizeof(_affineMatrix12));
  swapdouble = (big_endian_double(_affineMatrix20));
  outStream.write((const char *) &swapdouble, sizeof(_affineMatrix20));
  swapdouble = (big_endian_double(_affineMatrix21));
  outStream.write((const char *) &swapdouble, sizeof(_affineMatrix21));
  swapdouble = (big_endian_double(_affineMatrix22));
  outStream.write((const char *) &swapdouble, sizeof(_affineMatrix22));
  
  _translationVector.save(outStream);
  _greensCoeff.save(outStream);

  // Save parent class's parameters
  if ( FieldTransformation::saveClassParameters( outStream ) )
  {
    return 1;
  };
  
  return 0;
}
 


//////////////////////////////////////////////////////////////////////////////
//
// Function: loadClassParameters
//
// Purpose: Load parameters that are owned by this class and all of it's parent's
//
// Inputs: ifstream that we want to write to 
//
// Outputs: Returns 0 upon success and non-zero upon failure
//          Returns 1 if read-in length of class name is negative or > 1000
//          Returns 2 if read-in class name is not correct
//          Returns 3 if read-in class version is not correct
//          Returns 4 if cannot alloc memory
//          Returns 5 if fails reading in parent class
//          Returns 6 if read-in length of class revision is negative or > 1000
//
//////////////////////////////////////////////////////////////////////////////
int
ManifoldTransformation::loadClassParameters ( ifstream & inStream)
{
  const char * const functionName = "ManifoldTransformation::loadClassParameters";

  // Read in the length of the class name from the file
  int ReadInClassNameLength;
  inStream.read((char *) &ReadInClassNameLength, sizeof (ReadInClassNameLength));
  ReadInClassNameLength = big_endian_int(ReadInClassNameLength);
  // Let's not get crazy, even _I_ don't make class names this big
  if (ReadInClassNameLength < 0 || ReadInClassNameLength > 1000) { 
    cerr << "I don't think that this is a " << mClassName << " Results File." << endl;
    cerr << "Apparent class length is " << ReadInClassNameLength << endl;
    // Throw error here!
    return 1;
  }

  // Alloc memory for the class name and read it in
  char * ReadInClassName = new char[ReadInClassNameLength + 1];
  if ( ! ReadInClassName ){
    cerr << functionName << ": ";
    cerr << "Failed to alloc " << ReadInClassNameLength+1 << " bytes." << endl;
    // Throw error here!
    return 4;
  }
  inStream.read((char *) ReadInClassName, ReadInClassNameLength);
  ReadInClassName[ReadInClassNameLength] = '\0'; // End-of-string so I can use strcmp()

  // Check to see if the file's class name is mClassName and return with failure if not
  if (0 != strcmp(ReadInClassName, mClassName))
  {
    cerr << functionName << ": ";
    cerr << "Attempted to read in Class " << ReadInClassName;
    cerr << " whereas this Class is " << mClassName << endl;
    // Maybe should backup inStream before returning
    // Maybe should skip the rest of this record before returning, probally not
    // Throw an error here!
    return 2;
  }

  // Ok, we are in the correct class, let's check if we have the correct
  // version of mClassName
  int ReadInClassVersion;
  inStream.read((char *) & ReadInClassVersion, sizeof (ReadInClassVersion));
  ReadInClassVersion = big_endian_int(ReadInClassVersion);
  if (ReadInClassVersion != mClassVersion) {
    cerr << functionName << ": ";
    cerr << "Attempted to read in " << mClassName << " Class parameters from version ";
    cerr << ReadInClassVersion << " whereas this code is version ";
    cerr << mClassVersion << endl;
    // Maybe should backup before returning
    // Throw an error here!
    return 3;
  }

  // Read in the length of the class revision from the file
  int ReadInClassRevisionLength;
  inStream.read((char *) &ReadInClassRevisionLength, sizeof (ReadInClassRevisionLength));
  ReadInClassRevisionLength = big_endian_int(ReadInClassRevisionLength);
  // Let's not get crazy, even _I_ don't make class names this big
  if (ReadInClassRevisionLength < 0 || ReadInClassRevisionLength > 1000) { 
    cerr << "I don't think that this is a " << mClassName << " Results File." << endl;
    cerr << "Apparent class revision length is " << ReadInClassRevisionLength << endl;
    // Throw error here!
    return 6;
  }

  // Alloc memory for the class revision and read it in
  char * ReadInClassRevision = new char[ReadInClassRevisionLength + 1];
  if ( ! ReadInClassRevision ){
    cerr << functionName << ": ";
    cerr << "Failed to alloc " << ReadInClassRevisionLength+1 << " bytes." << endl;
    // Throw error here!
    return 4;
  }
  inStream.read((char *) ReadInClassRevision, ReadInClassRevisionLength);
  // End-of-string so I can use strcmp()
  ReadInClassRevision[ReadInClassRevisionLength] = '\0'; 
  // Do nothing with it for now, mainly for debugging outside class


  // Ok, things look good, let's read in the parameters that mClassName owns

  // Read in input parameters for this class
  _atlasPoints.load(inStream);
  _patientPoints.load(inStream);
  _variancePoints.load(inStream);
  inStream.read((char *) &_numPoints, sizeof(_numPoints));
  _numPoints = big_endian_int(_numPoints);
  // Read in the transformation results
  inStream.read((char *) &_affineMatrix00, sizeof(_affineMatrix00));
  inStream.read((char *) &_affineMatrix01, sizeof(_affineMatrix01));
  inStream.read((char *) &_affineMatrix02, sizeof(_affineMatrix02));
  inStream.read((char *) &_affineMatrix10, sizeof(_affineMatrix10));
  inStream.read((char *) &_affineMatrix11, sizeof(_affineMatrix11));
  inStream.read((char *) &_affineMatrix12, sizeof(_affineMatrix12));
  inStream.read((char *) &_affineMatrix20, sizeof(_affineMatrix20));
  inStream.read((char *) &_affineMatrix21, sizeof(_affineMatrix21));
  inStream.read((char *) &_affineMatrix22, sizeof(_affineMatrix22));

  //do byteswapping if needed
  _affineMatrix00=(big_endian_double(_affineMatrix00));
  _affineMatrix01=(big_endian_double(_affineMatrix01));
  _affineMatrix02=(big_endian_double(_affineMatrix02));
  _affineMatrix10=(big_endian_double(_affineMatrix10));
  _affineMatrix11=(big_endian_double(_affineMatrix11));
  _affineMatrix12=(big_endian_double(_affineMatrix12));
  _affineMatrix20=(big_endian_double(_affineMatrix20));
  _affineMatrix21=(big_endian_double(_affineMatrix21));
  _affineMatrix22=(big_endian_double(_affineMatrix22));
  
_translationVector.load(inStream);
  _greensCoeff.load(inStream);
 
  if ( FieldTransformation::loadClassParameters( inStream ) )
  {
    return 5;
  };

  // Return with success
  return 0;
}



//////////////////////////////////////////////////////////////////////////////
//
// Function: save
//
// Purpose: Save parameters for this class (including all base classes)
//
// Inputs: filename
//
// Outputs: Return 0 upon success
//          Returns 1 if cannot open file
//          Returns 2 if other failure
//
//////////////////////////////////////////////////////////////////////////////
int
ManifoldTransformation::save( const char * filename )
{
  ofstream outStream (filename, ios::out);

  if ( ! outStream )
  {
    cerr << "Failed to open file \"" << filename << "\"" << endl;
    return 1;
  }
  if ( saveClassParameters( outStream ) )
  {
    return 2;
  };
  outStream.close();

  // cout << "save() in ManifoldTransformation" << endl;
  return 0;
}



//////////////////////////////////////////////////////////////////////////////
//
// Function: printClassParameters
//
// Purpose: Calls parent's print function and then prints parameters that are owned
//          by this class to the supplied ostream.
//          It is envisioned that the child of this class calls this function and
//          then prints out it's own parameters.
//
// Inputs: ostream that we want to print to 
//
// Outputs: reference to the passed in ostream
//
//////////////////////////////////////////////////////////////////////////////
ostream & 
ManifoldTransformation::printClassParameters ( ostream & outStream)
{

  FieldTransformation::printClassParameters ( outStream );

  outStream << mClassName << "::mClassVersion: "  << mClassVersion << endl;
  outStream << mClassName << "::mClassRevision: " << mClassRevision << endl;
  //  outStream << mClassName << "::sigma: " << sigma << endl;

  return outStream;
}



//////////////////////////////////////////////////////////////////////////////
//
// Public Methods
//
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
//
// Function: load
//
// Purpose: Load parameters for this class (including all base classes)
//
// Inputs: filename
//
// Outputs: Return 0 upon success
//          Returns 1 if cannot open file
//          Returns 2 if other failure
//
//////////////////////////////////////////////////////////////////////////////
int
ManifoldTransformation::load( const char * filename )
{
  ifstream inStream (filename, ios::in);
  if ( ! inStream )
  {
    cerr << "Failed to open file \"" << filename << "\"" << endl;
    return 1;
  }

  if ( loadClassParameters( inStream ) )
  {
    return 2;
  }

  // cout << "load() in ManifoldTransformation" << endl;
  inStream.close();

  return 0;
}

//////////////////////////////////////////////////////////////////////////////
//
// Function: initialize
//
// Purpose: Initialize class so that calculate may run
//
// Inputs: atlas landmarks
//         patient landmarks
//         landmark variance
//
// Outputs: Sets mOkToCalculateFlag to 1 if successfull
//
//////////////////////////////////////////////////////////////////////////////
int ManifoldTransformation::initialize(const Array2D<double> &atlasPoints,
				       const Array2D<double> &patientPoints,
				       const Array1D<double> &variancePoints,
				       const int atlasSizeX,
				       const int atlasSizeY,
				       const int atlasSizeZ,
				       const int patientSizeX,
				       const int patientSizeY,
				       const int patientSizeZ, 
				       const double atlasResX, 
				       const double atlasResY, 
				       const double atlasResZ, 
				       const double patientResX, 
				       const double patientResY, 
				       const double patientResZ) {

  cout << "initialize() in ManifoldTransformation" << endl;

  // error checking
  if ((atlasPoints.getNrow() != patientPoints.getNrow()) ||
      (atlasPoints.getNcol() != patientPoints.getNcol())) {
    return 1;
  }
  if ((atlasPoints.getNcol() != 3) ||
      (patientPoints.getNcol() != 3)) {
    return 2;
  }
  if ((atlasPoints.getNrow() < 4) ||
      (patientPoints.getNrow() < 4)) {
    return 3;
  }

  _atlasSizeX = atlasSizeX;
  _atlasSizeY = atlasSizeY;
  _atlasSizeZ = atlasSizeZ;
  _patientSizeX = patientSizeX;
  _patientSizeY = patientSizeY;
  _patientSizeZ = patientSizeZ;

  _atlasResX = atlasResX; 
  _atlasResY = atlasResY; 
  _atlasResZ = atlasResZ; 
  _patientResX = patientResX; 
  _patientResY = patientResY; 
  _patientResZ = patientResZ; 

  // Take resolution into account.
  _atlasPoints = atlasPoints; 
  _patientPoints = patientPoints; 
  _variancePoints = variancePoints; 

  // copy landmarks
  _numPoints = _atlasPoints.getNrow();

  for (int i = 0; i < _numPoints; i++) {

     // Atlas points.
     _atlasPoints[i][0] *= _atlasResX;
     _atlasPoints[i][1] *= _atlasResY;
     _atlasPoints[i][2] *= _atlasResZ;

     // Patient points.
     _patientPoints[i][0] *= _patientResX; 
     _patientPoints[i][1] *= _patientResY; 
     _patientPoints[i][2] *= _patientResZ; 

  }

  setOkToCalculate( 1 );
   
  return 0;

}



//////////////////////////////////////////////////////////////////////////////
//
// Function: calculate
//
// Purpose: This method computes the manifold displacement field
//
// Inputs: none
//
// Outputs: Turns mOkToApplyFlag to 1 if successful
//          Returns 0 if sucessful
//          Returns 1 if not ready to calculate
//          Returns 2 if no solution could be found for manifoldMatrix
//          Returns 3 if solution is not good enough as measured by
//                       the norm of the error matrix
//
//////////////////////////////////////////////////////////////////////////////
int
ManifoldTransformation::calculate()
{
  cout << "calculate() in ManifoldTransformation" << endl;

  cout.precision(7);
  cerr.precision(7);

  if ( ! isOkToCalculate() )
    return 1;

  // construct A & B matrices
  Matrix<double> B(_numPoints + 4, 3);
  Matrix<double> A(_numPoints + 4, _numPoints + 4);

  int i;
  int j;

  ShowWorking("Manifold calculating");

  for (i=0; i < _numPoints; i++) {

    // Diagonal elements of A.
    A[i][i] = _variancePoints[i]; 
    
    double patienti0 = _patientPoints[i][0];
    double patienti1 = _patientPoints[i][1];
    double patienti2 = _patientPoints[i][2];

    // Off-diagonal elements of A.
    for (j = i + 1; j < _numPoints; j++) {
      double point0diff = patienti0 - _patientPoints[j][0];
      double point1diff = patienti1 - _patientPoints[j][1];
      double point2diff = patienti2 - _patientPoints[j][2];
      A[i][j] = sqrt(point0diff*point0diff +
		     point1diff*point1diff +
		     point2diff*point2diff);
      A[j][i] = A[i][j]; 
    }

    // Extra columns and rows of A.
    A[i][j] = A[j][i] = patienti0;
    A[i][j + 1] = A[j + 1][i] = patienti1;
    A[i][j + 2] = A[j + 2][i] = patienti2;
    A[i][j + 3] = A[j + 3][i] = 1; 

    // Elements of B.
    B[i][0] = (patienti0 - _atlasPoints[i][0]); 
    B[i][1] = (patienti1 - _atlasPoints[i][1]); 
    B[i][2] = (patienti2 - _atlasPoints[i][2]); 

  }

  // Extra diag(A) and rest of B.
  for (i=_numPoints; i < _numPoints + 4; i++) {
    B[i][0] = B[i][1] = B[i][2] = 0; 
    for (j = i; j < _numPoints + 4; j++)
      A[i][j] = A[j][i] = 0; 
  }

  // Print condition number of A.
//  cout << "condition number of A = " << cond(A) << endl; 

  // Solve for transformation parameters.
  Matrix<double> manifoldMatrix(A.getNcol(), B.getNcol());
  int solveResult = MatrixUtils<double>::solve(A, B, &manifoldMatrix,"nobalance");
  if (solveResult == 0)
    return 2;
  
  // Check the manifoldMatrix
  Matrix<double> errorMatrix = A * manifoldMatrix - B;
  double errorMatrixNorm = MatrixUtils<double>::norm(errorMatrix);
  cout << "errorMatrixNorm = " << errorMatrixNorm << endl;
  if (errorMatrixNorm > 1.0e-4)
    return 3;


  ///////////////////////////////////////////////////////
  //
  //  Extract parameters from the manifoldMatrix
  //
  ///////////////////////////////////////////////////////

  // Green's function coefficients (Beta)
  _greensCoeff.setDim(_numPoints, 3); 
  int k;
  for(k=0; k < _numPoints; k++) {
    _greensCoeff[k][0] = manifoldMatrix[k][0]; 
    _greensCoeff[k][1] = manifoldMatrix[k][1]; 
    _greensCoeff[k][2] = manifoldMatrix[k][2]; 
  }

  // affine transformation matrix
  _affineMatrix00 = manifoldMatrix[_numPoints][0]; 
  _affineMatrix01 = manifoldMatrix[_numPoints+1][0];
  _affineMatrix02 = manifoldMatrix[_numPoints+2][0]; 
  _affineMatrix10 = manifoldMatrix[_numPoints][1]; 
  _affineMatrix11 = manifoldMatrix[_numPoints+1][1];
  _affineMatrix12 = manifoldMatrix[_numPoints+2][1]; 
  _affineMatrix20 = manifoldMatrix[_numPoints][2]; 
  _affineMatrix21 = manifoldMatrix[_numPoints+1][2];
  _affineMatrix22 = manifoldMatrix[_numPoints+2][2];

  // extract the translation vector
  _translationVector.setDim(3); 
  _translationVector[0] = manifoldMatrix[_numPoints + 3][0]; 
  _translationVector[1] = manifoldMatrix[_numPoints + 3][1]; 
  _translationVector[2] = manifoldMatrix[_numPoints + 3][2]; 

//  cout << _translationVector << endl;

  // compute mean absolute error distance between transformed patient
  // points and atlas points
  double meanError = 0.0;
  double y0, y1, y2;
  for (k=0; k<_numPoints; k++) {
    transformPoint((double) _patientPoints[k][0], 
		   (double) _patientPoints[k][1], 
		   (double) _patientPoints[k][2],
		   y0, y1, y2);
    cout << "(" << _atlasPoints[k][0] << ", " 
		<< _atlasPoints[k][1] << ", "
		<< _atlasPoints[k][2] << ") ---> "
		<< "(" << y0 << ", " << y1 << ", " << y2 << ")" << endl; 
    y0 -= _atlasPoints[k][0];
    y1 -= _atlasPoints[k][1];
    y2 -= _atlasPoints[k][2];
    meanError += sqrt(y0*y0 + y1*y1 + y2*y2);
  }
  meanError /= _numPoints;
  cout.precision(7);
  cout << "average landmark error = " << meanError << endl;

  // Say we have a transformation and we may now apply it
  setOkToApply( 1 );

  return 0;
}



//////////////////////////////////////////////////////////////////////////////
//
// Function: transformPoint()
//
// Purpose: This procedure takes a point in patient space and finds the
//          corresponding point in atlas space
//
// Inputs: point (x0, x1, x2)  in patient space to be transformed
//
// Outputs: transformed point (y0, y1, y2) in atlas space
//
//////////////////////////////////////////////////////////////////////////////
void ManifoldTransformation::transformPoint(const double x0, 
					    const double x1, 
					    const double x2,
					    double &y0,
					    double &y1,
					    double &y2) {

  y0 = x0 - (_affineMatrix00*x0 + 
	     _affineMatrix01*x1 + 
	     _affineMatrix02*x2 + _translationVector[0]);
  y1 = x1 - (_affineMatrix10*x0 + 
	     _affineMatrix11*x1 + 
	     _affineMatrix12*x2 + _translationVector[1]);
  y2 = x2 - (_affineMatrix20*x0 + 
	     _affineMatrix21*x1 + 
	     _affineMatrix22*x2 + _translationVector[2]);

  for (int i=0; i<_numPoints; i++) {
    double x0diff = x0 - _patientPoints[i][0];
    double x1diff = x1 - _patientPoints[i][1];
    double x2diff = x2 - _patientPoints[i][2];
    double norm = sqrt((x0diff*x0diff) + (x1diff*x1diff) + (x2diff*x2diff));
    y0 -= norm * _greensCoeff[i][0];
    y1 -= norm * _greensCoeff[i][1];
    y2 -= norm * _greensCoeff[i][2];
  }
}

/***************************************************************************
* Tri-linear interpolation.
*
* The routine employs Tri-linear interpolation to overcome rounding errors.
* For a pixel at (z, y, x) in the final volume Y, find the nearest
* 8 pixels to the pixel (z0, y0, x0) in the original volume I and 
* interpolate. For points outside of original image support, 
* use background.
*
* Abed M. Hammoud
****************************************************************************/
double ManifoldTransformation::trilin(Array3D<unsigned char> const &volume, 
				     double z, double y, double x) {

  /***************
   * House Keeping.
   ***************/
  
  int z0 = (int)floor(z), y0 = (int)floor(y), x0 = (int)floor(x); 
  int z1 = z0 + 1, y1 = y0 + 1, x1 = x0 + 1; 
  unsigned char d000, d001, d010, d011, d100, d101, d110, d111; 
  int xsize = volume.getXsize(), ysize = volume.getYsize(), zsize = volume.getZsize(); 

  int inBoundsX0;
  int inBoundsY0;
  int inBoundsZ0;
  int inBoundsX1;
  int inBoundsY1;
  int inBoundsZ1;
  int compute;
  const unsigned char *const *const *volumePtrPtrPtr = volume.address();

  /*****************************************************************
   * (0 <= fx, fy, fz <= 1), fractional position between data points.
   *****************************************************************/

  double fz = z - z0, fy = y - y0, fx = x - x0; 

  /**********************************************
   * Tri-linear interpolation for interior points.
   **********************************************/
  
  if (x0 >= 0 && x1 < xsize && y0 >= 0 && y1 < ysize && 
      z0 >= 0 && z1 < zsize) {
    
    // Use pointer arithmetic.
    unsigned char const *dp = &volumePtrPtrPtr[z0][y0][x0]; 
    
    /****************************
     * Data used in interpolation.
     ****************************/

    d000 = dp[0]; d001 = dp[1]; 
    dp += xsize; 
    
    d010 = dp[0]; d011 = dp[1]; 
    dp += xsize * ysize; 
    
    d110 = dp[0]; d111 = dp[1]; 
    dp -= xsize; 
    
    d100 = dp[0]; d101 = dp[1]; 
    
    compute = 1;
  } else {

    /**********************************************
     * Tri-linear interpolation for boundary points.
     * Points outside of Array3D are set to bkgrnd.
     **********************************************/
    compute = 0;
    d000 = d001 = d010 = d011 = 0;
    d100 = d101 = d110 = d111 = 0;
    inBoundsX0 = (x0 >= 0) && (x0 < xsize);
    inBoundsY0 = (y0 >= 0) && (y0 < ysize);
    inBoundsZ0 = (z0 >= 0) && (z0 < zsize);
    inBoundsX1 = (x1 >= 0) && (x1 < xsize);
    inBoundsY1 = (y1 >= 0) && (y1 < ysize);
    inBoundsZ1 = (z1 >= 0) && (z1 < zsize);

    if (inBoundsZ0 && inBoundsY0 && inBoundsX0) {
      d000 = volumePtrPtrPtr[z0][y0][x0];
      compute = 1;
    }
    if (inBoundsZ0 && inBoundsY0 && inBoundsX1) {
      d001 = volumePtrPtrPtr[z0][y0][x1];
      compute = 1;
    }
    if (inBoundsZ0 && inBoundsY1 && inBoundsX0) {
      d010 = volumePtrPtrPtr[z0][y1][x0];
      compute = 1;
    }
    if (inBoundsZ0 && inBoundsY1 && inBoundsX1) {
      d011 = volumePtrPtrPtr[z0][y1][x1];
      compute = 1;
    }
    if (inBoundsZ1 && inBoundsY0 && inBoundsX0) {
      d100 = volumePtrPtrPtr[z1][y0][x0];
      compute = 1;
    }
    if (inBoundsZ1 && inBoundsY0 && inBoundsX1) {
      d101 = volumePtrPtrPtr[z1][y0][x1];
      compute = 1;
    }
    if (inBoundsZ1 && inBoundsY1 && inBoundsX0) {
      d110 = volumePtrPtrPtr[z1][y1][x0];
      compute = 1;
    }
    if (inBoundsZ1 && inBoundsY1 && inBoundsX1) {
      d111 = volumePtrPtrPtr[z1][y1][x1];
      compute = 1;
    }
  }

  if (compute) {
    /********************
     * Interpolation, (1).
     ********************/
  
    double d00x = d000 + fx*(d001-d000);
    double d01x = d010 + fx*(d011-d010);
    double d10x = d100 + fx*(d101-d100);
    double d11x = d110 + fx*(d111-d110);

    double d0yx = d00x + fy*(d01x-d00x);
    double d1yx = d10x + fy*(d11x-d10x);
  
    return (d0yx + fz*(d1yx-d0yx));
  } else 
    return 0;
}


////////////////////////////////////////////////////////////////////////////////
//
// Function: print
//
// Purpose: Prints parameters
//
// Inputs: To be determined
//
// Outputs: To be determined
//
////////////////////////////////////////////////////////////////////////////////
void ManifoldTransformation::print() {

  cout << "print() in ManifoldTransformation" << endl;
  ManifoldTransformation::printClassParameters(cout);

}



//////////////////////////////////////////////////////////////////////////////
// Accessor functions.
//
//
//
//////////////////////////////////////////////////////////////////////////////

Matrix<double> ManifoldTransformation::getAffineMatrix() const {

   Matrix<double> M(3, 3); 

   M[0][0] = _affineMatrix00; 
   M[0][1] = _affineMatrix01; 
   M[0][2] = _affineMatrix02; 

   M[1][0] = _affineMatrix10; 
   M[1][1] = _affineMatrix11; 
   M[1][2] = _affineMatrix12; 

   M[2][0] = _affineMatrix20; 
   M[2][1] = _affineMatrix21; 
   M[2][2] = _affineMatrix22; 

   return M; 

}

Matrix<double> ManifoldTransformation::getGreensCoeff() const {
     return _greensCoeff; 
}

Vector<double> ManifoldTransformation::getTranslationVector() const {
   return _translationVector; 
}

int ManifoldTransformation::getNumPoints() const {
   return _numPoints; 
}

Array2D<double>  ManifoldTransformation::getPatientPoints() const {
   return _patientPoints; 
}

