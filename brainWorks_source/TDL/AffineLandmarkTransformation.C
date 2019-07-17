#include <string.h>
#include <TDL/AffineLandmarkTransformation.h>

//////////////////////////////////////////////////////////////////////////
//
//
//////////////////////////////////////////////////////////////////////////
//
// Revision Log:
//
// $Log: AffineLandmarkTransformation.C,v $
// Revision 1.1  2004/11/15 04:44:07  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:29:20  kem
// Initial revision
//
// Revision 1.3  1999/07/09 17:50:56  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.2  1998/12/10 19:04:14  RAZORDB
// IDL/AW merge
//
// Revision 1.1  1997/08/01 19:49:58  csb
// Initial revision
//
// Revision 1.2  1997/08/01 14:54:38  csb
// Incorportate Abed's, fluid, elastic and field
//
// Revision 1.1  1997/06/25 17:39:42  kem
// Initial revision
//
//
//
//////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//
// Private Methods
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//
// Protected Methods
//
////////////////////////////////////////////////////////////////////////////////

const char * const AffineLandmarkTransformation::mClassRevision = "AffineLandmarkTransformation.h,v 1.1";
const char * const AffineLandmarkTransformation::mClassName = "AffineLandmarkTransformation";
const int          AffineLandmarkTransformation::mClassVersion = 1;


////////////////////////////////////////////////////////////////////////////////
//
// Function: saveClassParameters
//
// Purpose: Saves parameters that are owned by this class to the supplied ofstream.
//          It is envisioned that the child of this class calls this function and
//          then saves it's own parameters. Since at this time this is the 
//          the final derived class, save() calls this function.
//
// Inputs: ofstream that we want to save to 
//
// Outputs: 0 upon success
//          1 if parent's save returned non-zero
//
////////////////////////////////////////////////////////////////////////////////
int
AffineLandmarkTransformation::saveClassParameters ( ofstream & outStream)
{
  int swapint; //temp variable that holds swap info
  

  // const char * functionName = "AffineLandmarkTransformation::saveClassParameters";

  // Write out size of class name and the class name
  const int nameLength = strlen(mClassName);
  swapint = big_endian_int(nameLength);
  outStream.write((const char *) & swapint, sizeof(nameLength));
  outStream.write(mClassName,                  nameLength);

  // Write out the class version
  const int temp  = mClassVersion;          // Because the compiler may not put these variables
  swapint = big_endian_int(temp);
  outStream.write((const char *) & swapint, sizeof (temp));

  // Write out size of class revision and the class revision
  const int revisionLength = strlen(mClassRevision);
  swapint = big_endian_int(revisionLength);
  outStream.write((const char *) & swapint, sizeof(revisionLength));
  outStream.write(mClassRevision, revisionLength);

  // Save parent class's parameters
  if ( LinearTransformation::saveClassParameters( outStream ) )
  {
    return 1;
  };
  
  return 0;
}
 

////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////
int
AffineLandmarkTransformation::loadClassParameters ( ifstream & inStream)
{
  const char * const functionName = "AffineLandmarkTransformation::loadClassParameters";

  // Read in the length of the class name from the file
  int ReadInClassNameLength;
  inStream.read((char *) &ReadInClassNameLength, sizeof (ReadInClassNameLength));
  ReadInClassNameLength = big_endian_int(ReadInClassNameLength);
  // Let's not get crazy, even _I_ don't make class names this big
  if (ReadInClassNameLength < 0 || ReadInClassNameLength > 1000) { 
    cerr << "I don't think that this is a AffineLandmarkTransformation Results File." << endl;
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
  ReadInClassRevision[ReadInClassRevisionLength] = '\0'; // End-of-string so I can use strcmp()
  // Do nothing with it for now, mainly for debugging outside class

  if ( LinearTransformation::loadClassParameters( inStream ) )
  {
    return 5;
  };


  // Return with success
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////
int
AffineLandmarkTransformation::save( const char * filename )
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

  // cout << "save() in AffineLandmarkTransformation" << endl;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////
ostream & 
AffineLandmarkTransformation::printClassParameters ( ostream & outStream)
{
  LinearTransformation::printClassParameters ( outStream );
  outStream << mClassName << "::mClassVersion: " << mClassVersion << endl;
  outStream << mClassName << "::mClassRevision: " << mClassVersion << endl;

  return outStream;
}


////////////////////////////////////////////////////////////////////////////////
//
// Public Methods
//
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////
int
AffineLandmarkTransformation::load( const char * filename )
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

  // cout << "load() in AffineLandmarkTransformation" << endl;
  inStream.close();

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//
// Function: calculate
//
// Purpose: This method computes the affine transformation Y = X*A' + B
//          that minimizes squared error
//
// Inputs: Two sets of landmark points
//
// Outputs: Turns mOkToApplyFlag to 1 if successfull
//          Returns 0 if successful
//          Returns 1 if transformation not initialized
//          Returns 2 if M matrix is not invertible
//
////////////////////////////////////////////////////////////////////////////////
int
AffineLandmarkTransformation::calculate()
{
  if ( ! isOkToCalculate() )
    return 1;

  int pp,numPnts = _patientPoints->getNrow();

  // Check for 2-D case
  float allZ = (*_patientPoints)[0][2];
  for(pp=0;pp<numPnts;pp++)
    if(((*_patientPoints)[pp][2] != allZ)||((*_atlasPoints)[pp][2] != allZ))
	break;

  bool is2D = (pp==numPnts);
  if(is2D) cerr << "AffineLandmarkTransformation: using 2D case" << endl;

  int nsum;

  if(is2D) nsum = 2;
  else     nsum = 3;

  ShowWorking("Working...");

  // create sums to form matrices
  Vector<double> sumX(nsum);
  Vector<double> sumY(nsum);
  Matrix<double> sumXX(nsum,nsum);
  Matrix<double> sumXY(nsum,nsum);
  sumX  = 0.0;
  sumY  = 0.0;
  sumXX = 0.0;
  sumXY = 0.0;
  for (int i=0; i<numPnts; i++) {
    for (int j=0; j<nsum; j++) {
      sumX[j] += (*_patientPoints)[i][j];
      sumY[j] += (*_atlasPoints)[i][j];
      for (int k=0; k<nsum; k++) {
	sumXX[j][k] += (*_patientPoints)[i][j] * (*_patientPoints)[i][k];
	sumXY[j][k] += (*_patientPoints)[i][j] * (*_atlasPoints)[i][k];
      }
    }
  }

  // M = [ numPnts  sumX0    sumX1    sumX2 ;
  //       sumX0    sumX0X0  sumX0X1  sumX0X2 ;
  //       sumX1    sumX0X1  sumX1X1  sumX1X2 ;
  //       sumX2    sumX0X2  sumX1X2  sumX2X2 ]
  Matrix<double> M(nsum+1,nsum+1);

  M[0][0] = numPnts;
  M[0][1] = sumX[0];
  M[0][2] = sumX[1];
  M[1][0] = sumX[0];
  M[1][1] = sumXX[0][0];
  M[1][2] = sumXX[0][1];
  M[2][0] = sumX[1];
  M[2][1] = sumXX[0][1];
  M[2][2] = sumXX[1][1];

  if(!is2D) {
      M[0][3] = sumX[2];
      M[1][3] = sumXX[0][2];
      M[2][3] = sumXX[1][2];
      M[3][0] = sumX[2];
      M[3][1] = sumXX[0][2];
      M[3][2] = sumXX[1][2];
      M[3][3] = sumXX[2][2];
  }

  // Z = [ sumY0    sumY1    sumY2 ;
  //       sumX0Y0  sumX0Y1  sumX0Y2 ;
  //       sumX1Y0  sumX1Y1  sumX1Y2 ;
  //       sumX2Y0  sumX2Y1  sumX2Y2 ]
  Matrix<double> Z(nsum+1,nsum+1);
  Z = 0.;

  Z[0][0] = sumY[0];
  Z[0][1] = sumY[1];

  Z[1][0] = sumXY[0][0];
  Z[1][1] = sumXY[0][1];

  Z[2][0] = sumXY[1][0];
  Z[2][1] = sumXY[1][1];

  if(!is2D) {
    Z[0][2] = sumY[2];
    Z[1][2] = sumXY[0][2];
    Z[2][2] = sumXY[1][2];
    Z[3][0] = sumXY[2][0];
    Z[3][1] = sumXY[2][1];
    Z[3][2] = sumXY[2][2];
  } 

  ShowWorking("Working...");

  // resultMatrix = M inverse * Z

  Matrix<double> result(nsum+1,nsum+1);
  Matrix<double> Minverse(nsum+1,nsum+1);

  if (MatrixUtils<double>::inv(M, &Minverse) == 0) 
    return 2;

  MatrixUtils<double>::multiply(Minverse,Z,&result);

  // resultMatrix = [ b0   b1   b2 ;
  //                  a00  a10  a20 ;
  //                  a01  a11  a21 ;
  //                  a02  a12  a22 ]

  // compute translation vector from result
  mTranslationVector[0] = result[0][0];
  mTranslationVector[1] = result[0][1];
  if(is2D)
    mTranslationVector[2] = 0.;
  else
    mTranslationVector[2] = result[0][2];

  // compute affine matrix from result
  mAffineMatrix[0][0] = result[1][0];
  mAffineMatrix[0][1] = result[2][0];

  mAffineMatrix[1][0] = result[1][1];
  mAffineMatrix[1][1] = result[2][1];

  if(is2D) {
    mAffineMatrix[2][0] = 0.;
    mAffineMatrix[2][1] = 0.;
    mAffineMatrix[2][2] = 1.;
    mAffineMatrix[0][2] = 0.;
    mAffineMatrix[1][2] = 0.;
  } else {
    mAffineMatrix[0][2] = result[3][0];
    mAffineMatrix[1][2] = result[3][1];
    mAffineMatrix[2][0] = result[1][2];
    mAffineMatrix[2][1] = result[2][2];
    mAffineMatrix[2][2] = result[3][2];
  }

  // Say we have a transformation and we may now apply it
  setOkToApply( 1 );

  ShowWorking("DONE.");
  ShowWorking((char*)NULL);

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//
// Function: initialize
//
// Purpose: Initialize class so that calculate may run
//
// Inputs: Two sets of landmark points
//
// Outputs: Sets mOkToCalculateFlag to 1 if successfull
//          Returns 0 if successful
//          Returns 1 if matrices have different dimensions
//          Returns 2 if matrices do not have 3 columns (3-D points)
//          Returns 3 if matrices do not have at least 4 rows (4 points)
//
////////////////////////////////////////////////////////////////////////////////
int
AffineLandmarkTransformation::initialize(const Matrix<double> * const atlasPoints, 
				const Matrix<double> * const patientPoints)
{
  if ((atlasPoints->getNrow() != patientPoints->getNrow()) ||
      (atlasPoints->getNcol() != patientPoints->getNcol())) {
    return 1;
  }
  if ((atlasPoints->getNcol() != 3) ||
      (patientPoints->getNcol() != 3)) {
    return 2;
  }
  if ((atlasPoints->getNrow() < 4) ||
      (patientPoints->getNrow() < 4)) {
    return 3;
  }

  _atlasPoints = atlasPoints;
  _patientPoints = patientPoints;

  setOkToCalculate( 1 );

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
void
AffineLandmarkTransformation::print()
{
  // cout << "print() in AffineLandmarkTransformation" << endl;
  AffineLandmarkTransformation::printClassParameters(cout);
}
