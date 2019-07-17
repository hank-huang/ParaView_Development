#include <string.h>
#include <TDL/AffineMomentTransformation.h>

const char * const AffineMomentTransformation::mClassName = "AffineMomentTransformation";
const int          AffineMomentTransformation::mClassVersion = 1;


//////////////////////////////////////////////////////////////////////////
//
//
//////////////////////////////////////////////////////////////////////////
//
// Revision Log:
//
// $Log: AffineMomentTransformation.C,v $
// Revision 1.1  2004/11/15 04:44:07  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:29:31  kem
// Initial revision
//
// Revision 1.3  1999/07/09 17:51:05  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.2  1998/12/10 19:04:21  RAZORDB
// IDL/AW merge
//
// Revision 1.1  1997/08/01 19:49:59  csb
// Initial revision
//
// Revision 1.2  1997/08/01 15:00:16  csb
// Incorportate Abed's, fluid, elastic and field
//
// Revision 1.1  1997/06/25 17:39:43  kem
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
AffineMomentTransformation::saveClassParameters ( ofstream & outStream)
{
  // const char * functionName = "AffineMomentTransformation::saveClassParameters";

  int swapint;
  // Write out size of class name and the class name
  const int nameLength = strlen(mClassName);
  swapint = big_endian_int(nameLength);
  outStream.write((const char *) & swapint, sizeof(nameLength));
  outStream.write(mClassName,                  nameLength);

  // Write out the class version
  const int temp  = mClassVersion;          // Because the compiler may not put these variables
  swapint = big_endian_int(temp);
  outStream.write((const char *) & swapint, sizeof (temp));

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
//          Returns 1 if read-in class name is negative or > 1000
//          Returns 2 if read-in class name is not correct
//          Returns 3 if read-in class version is not correct
//          Returns 4 if cannot alloc memory
//          Returns 5 if fails reading in parent class
//
////////////////////////////////////////////////////////////////////////////////
int
AffineMomentTransformation::loadClassParameters ( ifstream & inStream)
{
  const char * const functionName = "AffineMomentTransformation::loadClassParameters";

  // Read in the length of the class name from the file
  int ReadInClassNameLength;
  inStream.read((char *) &ReadInClassNameLength, sizeof (ReadInClassNameLength));
  // Let's not get crazy, even _I_ don't make class names this big
  ReadInClassNameLength = big_endian_int(ReadInClassNameLength);
  if (ReadInClassNameLength < 0 || ReadInClassNameLength > 1000) { 
    cerr << "I don't think that this is a AffineMomentTransformation Results File." << endl;
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
AffineMomentTransformation::save( const char * filename )
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

  // cout << "save() in AffineMomentTransformation" << endl;
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
AffineMomentTransformation::printClassParameters ( ostream & outStream)
{
  LinearTransformation::printClassParameters ( outStream );
  outStream << mClassName << "::mClassVersion: " << mClassVersion << endl;

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
AffineMomentTransformation::load( const char * filename )
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

  // cout << "load() in AffineMomentTransformation" << endl;
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
AffineMomentTransformation::calculate()
{
  if (! isOkToCalculate() )
    return 1;

  ShowWorking(0.0);

  Vector<double> meanAtlas(3);
  Vector<double> meanPatient(3);
  Matrix<double> covAtlas(3,3);
  Matrix<double> covPatient(3,3);
  Matrix<double> A(3,3);
  Matrix<double> B(3,3);
  
  int xAtlas = _atlasPtr->getXsize();
  int yAtlas = _atlasPtr->getYsize();
  int zAtlas = _atlasPtr->getZsize();

  int xPatient = _patientPtr->getXsize();
  int yPatient = _patientPtr->getYsize();
  int zPatient = _patientPtr->getZsize();
    
  int x, y, z;

  int tmp;
  const unsigned char * const * const *atl = _atlasPtr->address();

  ShowWorking(0.1);

  // compute mean of the atlas
  double normAtlas = 0;
  meanAtlas = 0.0;
  for (z=0; z<zAtlas; z++)
    for (y=0; y<yAtlas; y++)
      for (x=0; x<xAtlas; x++) {
	tmp = (int) atl[z][y][x];
	meanAtlas[0] += x*tmp;
	meanAtlas[1] += y*tmp;
	meanAtlas[2] += z*tmp;
	normAtlas += tmp;
      }
  meanAtlas /= normAtlas;


  ShowWorking(0.4);

  // compute covariance of the atlas
  double xm, ym, zm;
  double xtmp, ytmp, ztmp;
  covAtlas = 0.0;
  for (z=0; z<zAtlas; z++) {
    zm = z-meanAtlas[2];
    for (y=0; y<yAtlas; y++) {
      ym = y-meanAtlas[1];
      for (x=0; x<xAtlas; x++) {
	xm = x-meanAtlas[0];
	tmp = (int) atl[z][y][x];
	xtmp = xm*tmp;
	ytmp = ym*tmp;
	ztmp = zm*tmp;
	covAtlas[0][0] += xm * xtmp;
	covAtlas[0][1] += xm * ytmp;
	covAtlas[0][2] += xm * ztmp;
	covAtlas[1][1] += ym * ytmp;
	covAtlas[1][2] += ym * ztmp;
	covAtlas[2][2] += zm * ztmp;
      }
    }
  }
  covAtlas[1][0] = covAtlas[0][1];
  covAtlas[2][0] = covAtlas[0][2];
  covAtlas[2][1] = covAtlas[1][2];
  covAtlas /= normAtlas;

  ShowWorking(0.6);

  const unsigned char * const * const *pat = _patientPtr->address();
  // compute mean of the patient
  double normPatient = 0.0;
  meanPatient = 0.0;
  for (z=0; z<zPatient; z++)
    for (y=0; y<yPatient; y++)
      for (x=0; x<xPatient; x++) {
	tmp = (int) pat[z][y][x];
	meanPatient[0] += x*tmp;
	meanPatient[1] += y*tmp;
	meanPatient[2] += z*tmp;
	normPatient += tmp;
      }
  meanPatient /= normPatient;

  ShowWorking(0.8);

  // compute covariance of the patient
  covPatient = 0.0;
  for (z=0; z<zPatient; z++) {
    zm = z-meanPatient[2];
    for (y=0; y<yPatient; y++) {
      ym = y-meanPatient[1];
      for (x=0; x<xPatient; x++) {
	xm = x-meanPatient[0];
	tmp = (int) pat[z][y][x];
	xtmp = xm*tmp;
	ytmp = ym*tmp;
	ztmp = zm*tmp;
	covPatient[0][0] += xm * xtmp;
	covPatient[0][1] += xm * ytmp;
	covPatient[0][2] += xm * ztmp;
	covPatient[1][1] += ym * ytmp;
	covPatient[1][2] += ym * ztmp;
	covPatient[2][2] += zm * ztmp;
      }
    }
  }
  covPatient[1][0] = covPatient[0][1];
  covPatient[2][0] = covPatient[0][2];
  covPatient[2][1] = covPatient[1][2];
  covPatient /= normPatient;

  Matrix<double> sqrtAtlas(3,3);
  Matrix<double> sqrtPatient(3,3);

  Matrix<double> U(3,3);
  Matrix<double> S(3,3);
  Matrix<double> V(3,3);

  MatrixUtils<double>::svd(covAtlas, &U, &S, &V);
  for (x=0; x<3; x++) 
    S[x][x] = sqrt(S[x][x]);
  sqrtAtlas = U * S * V.transpose();

  MatrixUtils<double>::svd(covPatient, &U, &S, &V);
  for (x=0; x<3; x++) S[x][x] = sqrt(S[x][x]);
  sqrtPatient = U * S * V.transpose();
 
  B = sqrtPatient * MatrixUtils<double>::inv(sqrtAtlas);
  
  A = MatrixUtils<double>::inv(B);

  ShowWorking(0.9);

  Vector<double> b(3);
  b = meanAtlas - A*meanPatient;

//  cout << "Check the answer (A * K2 - K1 * Ainvtranspose) should be 0" << endl;
//  cout << A * covPatient - covAtlas * inv(A).transpose() << endl;

// #define ROTATION
// #ifdef ROTATION
//   B = sqrtPatient * sqrtAtlas.transpose();
//   svd(B, &U, &S, &V);
//   A = V*U.transpose();
//   cout << A << endl;
// #endif

  mAffineMatrix = A;
  mTranslationVector = b;

  // Say we have a transformation and we may now apply it
  setOkToApply( 1 );

  ShowWorking(-1.0);

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//
// Function: initialize
//
// Purpose: Initialize class so that calculate may run
//
// Inputs: Two volumes
//
// Outputs: Sets mOkToCalculateFlag to 1 if successfull
//          Returns 0 if successful
//          Returns 1 if volumes are empty
//
////////////////////////////////////////////////////////////////////////////////
int
AffineMomentTransformation::initialize(const Array3D<unsigned char> * const atlasPtr, 
				const Array3D<unsigned char> * const patientPtr)
{
  if (atlasPtr->isEmpty() || patientPtr->isEmpty() ) {
    return 1;
  }

  _atlasPtr = atlasPtr;
  _patientPtr = patientPtr;

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
AffineMomentTransformation::print()
{
  cout << "print() in AffineMomentTransformation" << endl;
  AffineMomentTransformation::printClassParameters(cout);
}
