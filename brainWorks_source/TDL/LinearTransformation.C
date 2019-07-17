#include <string.h>
#include <TDL/LinearTransformation.h>

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
// $Log: LinearTransformation.C,v $
// Revision 1.1  2004/11/15 04:44:08  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:29:13  kem
// Initial revision
//
// Revision 1.8  1999/07/09 17:52:29  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.7  1998/12/10 19:05:11  RAZORDB
// IDL/AW merge
//
// Revision 1.7  1998/11/16 20:34:41  kwd
// Add apply to point
//
// Revision 1.6  1998/10/08 20:34:41  csb
// Add support for storing transformations in millimeters
//
// Revision 1.5  1998/07/23 21:04:40  csb
// Add getTransformMatrix() accessor function
//
// Revision 1.4  1997/12/12 23:22:50  csb
// change prototype for 2 apply()s
//
// Revision 1.3  1997/10/02 20:37:48  rst
// Fixed apply() routines for Surfaces and Points
//
// Revision 1.2  1997/09/11 22:08:16  rst
// added surface apply
//
// Revision 1.1  1997/08/01 19:50:04  csb
// Initial revision
//
// Revision 1.3  1997/08/01 16:01:55  csb
// Incorportate Abed's, fluid, elastic and field
//
// Revision 1.2  1997/06/25 22:52:05  kem
// Add apply() method for a list of points
//
// Revision 1.1  1997/06/25 17:39:45  kem
// Initial revision
//
//
//
//////////////////////////////////////////////////////////////////////////

const char * const LinearTransformation::mClassRevision = "$Id: LinearTransformation.C,v 1.1 2004/11/15 04:44:08 joeh Exp $";
const char * const LinearTransformation::mClassName = "LinearTransformation";
const int          LinearTransformation::mClassVersion = 2;



LinearTransformation::LinearTransformation()
{ 
  mAffineMatrix.setDim(3, 3);
  mAffineMatrix.eye();
  
  mTranslationVector.setDim(3);
  mTranslationVector = 0.0;
  
  mInMMFlag = 0;
}

void
LinearTransformation::setInMMFlag(int flag)
{
  if ( 0 == flag ) 
  {
    mInMMFlag = 0;
  }
  else
  {
    mInMMFlag = 1;
  }
}

int
LinearTransformation::getInMMFlag()
{
  return mInMMFlag;
}


////////////////////////////////////////////////////////////////////////////////
//
// Protected Methods
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//
// Function: getTransformMatrix()
//
// Purpose: Returns 4x4 matrix containing the affine matrix and translation
//          matrix as calculated by this class
//
// Inputs: Takes pointer to resulting matrix, matrix must be allocated
//
// Outputs: Returns 1 if resulting matrix is not 4x4
//          Returns 2 if there is no answer yet
//
////////////////////////////////////////////////////////////////////////////////
int
LinearTransformation::getTransformMatrix( Matrix<double> * matrixPtr )
{

  if ((matrixPtr->getNrow() != 4) || (matrixPtr->getNcol() != 4))
  {
    cerr << "LinearTransformation::getTransformMatrix(): Matrix pointer must point to a 4x4" << endl;
    return 1; // return early with error
  }

  if ( ! isOkToApply() ) // i.e., calculate() is done
  {
    return 2;
  }

  matrixPtr->eye();
  
  for (int i=0; i<3; i++)
  {
    (*matrixPtr)[i][3] = mTranslationVector[i];
    for (int j=0; j<3; j++)
    {
      (*matrixPtr)[i][j] = mAffineMatrix[i][j];
    }
  }

  return 0;
}



////////////////////////////////////////////////////////////////////////////////
//
// Function: saveClassParameters
//
// Purpose: Saves parameters that are owned by this class to the supplied ofstream.
//          It is envisioned that the child of this class calls this function and
//          then saves it's own parameters.
//
// Inputs: ofstream that we want to save to 
//
// Outputs: Returns 0 if successfull
//          Returns 1 if call to parent fails
//          Returns 2 if cannot save matrix or vector
//
////////////////////////////////////////////////////////////////////////////////
int
LinearTransformation::saveClassParameters ( ofstream & outStream)
{
  // Write out size of class name and the class name
  int swapint; //temporary integer swapping variable
  
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

  // Write out parameters for this class
  if ( mAffineMatrix.save(outStream) == 0 )
  {
    return 2;
  }

  if ( mTranslationVector.save(outStream) == 0)
  {
    return 2;
  }

  outStream.write((const char *) & mInMMFlag, sizeof (mInMMFlag));

  // Save parent class's parameters
  if ( Transformation::saveClassParameters( outStream ) )
  {
    return 1;
  };


  // KWD for Faisal 11/02
  cerr << "--------------------------------------" << endl;
  cerr << "Linear Transformation Parameters:" << endl;
  cerr << "Rotation Matrix:" << endl;
  for(int j=0;j<3;j++) {
    for(int i=0;i<3;i++)
      cerr << mAffineMatrix[j][i] << " ";
    cerr << endl;
  }
  cerr << endl;
  cerr << "Translation Vector:" << endl;
  cerr << mTranslationVector[0] << " " <<
          mTranslationVector[1] << " " <<
          mTranslationVector[2] << endl;
  cerr << "--------------------------------------" << endl;
 
  return 0;
}


////////////////////////////////////////////////////////////////////////////////
//
// Function: loadClassParameters
//
// Purpose: Load parameters that are owned by this class and it's parents
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
LinearTransformation::loadClassParameters ( ifstream & inStream)
{
  const char * const functionName = "LinearTransformation::loadClassParameters";

  // Read in the length of the class name from the file
  int ReadInClassNameLength;
  inStream.read((char *) &ReadInClassNameLength, sizeof (ReadInClassNameLength));
  ReadInClassNameLength = big_endian_int(ReadInClassNameLength);
  // Let's not get crazy, even _I_ don't make class names this big
  if (ReadInClassNameLength < 0 || ReadInClassNameLength > 1000) { 
    cerr << "I don't think that this is a LinearTransformation Results File." << endl;
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
  // Difference between version 1 and 2 is mInMMFlag, so if version 1
  // we can still read it in as long as we set mInMMFlag to 0 later
  int ReadInClassVersion;
  inStream.read((char *) & ReadInClassVersion, sizeof (ReadInClassVersion));
  ReadInClassVersion = big_endian_int(ReadInClassVersion);
  if ( ! (   ReadInClassVersion == mClassVersion
          || 1 == ReadInClassVersion            ) )
  {
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

  if ( 0 == mAffineMatrix.load(inStream) )
  {
    return 2;
  }
  if ( 0 == mTranslationVector.load(inStream) )
  {
    return 2;
  }
  
  if ( 1 == ReadInClassVersion )
  {
    // Ok to read in a class version of 1. Set mInMMFlag to 0
    if ( getVerboseMode() ) 
    {
      cout << functionName << ": ";
      cout << "Setting mInMMFlag to 0" << endl;
    }
    mInMMFlag = 0;
  }
  else
  {
    inStream.read((char *) &mInMMFlag, sizeof(mInMMFlag));
    mInMMFlag = big_endian_int(mInMMFlag);
  }
  
  // Ok, things look good, let's read in the parameters that mClassName owns
  // None at this time

  if ( Transformation::loadClassParameters( inStream ) )
  {
    return 5;
  };

  // KWD for Faisal 11/02
  cerr << "--------------------------------------" << endl;
  cerr << "Linear Transformation Parameters:" << endl;
  cerr << "Rotation Matrix:" << endl;
  for(int j=0;j<3;j++) {
    for(int i=0;i<3;i++)
      cerr << mAffineMatrix[j][i] << " ";
    cerr << endl;
  }
  cerr << endl;
  cerr << "Translation Vector:" << endl;
  cerr << mTranslationVector[0] << " " <<
          mTranslationVector[1] << " " <<
          mTranslationVector[2] << endl;
  cerr << "--------------------------------------" << endl;
 

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
LinearTransformation::printClassParameters ( ostream & outStream)
{
  Transformation::printClassParameters ( outStream );
  outStream << mClassName << "::mClassVersion: " << mClassVersion << endl;
  outStream << mClassName << "::mClassRevision: " << mClassRevision << endl;

//  outStream << mClassName << "::mAffineMatrix: " << mAffineMatrix << endl;
//  outStream << mClassName << "::TranslationVector: " << mTranslationVector << endl;
  outStream << mClassName << "::InMMFlag: " << mInMMFlag << endl;

  return outStream;
}


////////////////////////////////////////////////////////////////////////////////
//
// Public Methods
//
////////////////////////////////////////////////////////////////////////////////

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
LinearTransformation::save( const char * filename )
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

  if (getVerboseMode())
    cout << "save() in LinearTransformation" << endl;

  return 0;
}

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
LinearTransformation::load( const char * filename )
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

  if (getVerboseMode())
    cout << "load() in LinearTransformation" << endl;

  inStream.close();

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
LinearTransformation::print()
{
  if (getVerboseMode())
    cout << "print() in LinearTransformation" << endl;

  LinearTransformation::printClassParameters(cout);
}


////////////////////////////////////////////////////////////////////////////////
//
// Function: apply()
//
// Purpose: Apply mAffineMatrix and mTranslationVector to a volume
//
// Inputs: 1) volume to be transformed
//         2) volume to put result
//         3) flag to trilinearly interpolate or not
//
// Outputs: transformed volume
//
// Returns: Returns 0 upon success
//          Returns 1 if volumes are the same
//          Returns 2 if input volume is empty
//          Returns 3 if output volume is empty
//          Returns 4 if not ready to apply
// 
////////////////////////////////////////////////////////////////////////////////
int
LinearTransformation::apply( const Array3D<unsigned char> &inVolume,
                             Array3D<unsigned char> * outVolumePtr,
			     double Rx, double Ry, double Rz, 
			     double sOx, double sOy, double sOz, 
                             const int interpolateFlag)
{
  if (&inVolume == outVolumePtr)
  {
    return 1;
  }
  if ( inVolume.isEmpty() )
  {
    return 2;
  }
  if ( outVolumePtr->isEmpty() )
  {
    return 3;
  }
  if ( ! isOkToApply() )
  {
    return 4;
  }

  Matrix<double> A(mAffineMatrix);
  Vector<double> b(mTranslationVector);

  // copy matrix elements
  double A00 = A[0][0], A01 = A[0][1], A02 = A[0][2];
  double A10 = A[1][0], A11 = A[1][1], A12 = A[1][2];
  double A20 = A[2][0], A21 = A[2][1], A22 = A[2][2];

  // for voxel elements we really want y = A (x+0.5) + b
  // so add A*0.5 to b
#if 0
  double b0 = b[0] + 0.5 * (A00 + A01 + A02) - 0.5;
  double b1 = b[1] + 0.5 * (A10 + A11 + A12) - 0.5;
  double b2 = b[2] + 0.5 * (A20 + A21 + A22) - 0.5;
#endif
  double b0 = b[0];
  double b1 = b[1];
  double b2 = b[2];

  // get volume sizes
  int y2Max = inVolume.getZsize();
  int y1Max = inVolume.getYsize();
  int y0Max = inVolume.getXsize();
  int x2Max = outVolumePtr->getZsize();
  int x1Max = outVolumePtr->getYsize();
  int x0Max = outVolumePtr->getXsize();

  // variables for trilinear interpolation
  double C1, C2, C3;
  unsigned char d000, d001, d010, d011;
  unsigned char d100, d101, d110, d111;
  double d00x, d01x, d10x, d11x;
  double d0xx, d1xx;

  // index variables
  int x0, x1, x2;

  // volume pointers
  const unsigned char *const *const * inVolPtrPtrPtr = inVolume.address();
  unsigned char *outVolumeDataPtr = outVolumePtr->data();

  const unsigned char *dptr;  // temporary pointer to element of inVolume

  int computeFlag;  // flag to indicate data exists for trilinear interpolation

  // y = A * x + b
  // y0 = A00 * x0 + A01 * x1 + A02 * x2 + b0
  // y1 = A10 * x0 + A11 * x1 + A12 * x2 + b1
  // y2 = A20 * x0 + A21 * x1 + A22 * x2 + b2
  
  double sum0;  // sum0 = A01 * x1 + A02 * x2 + b0
  double sum1;  // sum1 = A11 * x1 + A12 * x2 + b1
  double sum2;  // sum2 = A21 * x1 + A22 * x2 + b2
  
  double y0f;    // y0f = A00 * x0 + sum0
  double y1f;    // y1f = A10 * x0 + sum1
  double y2f;    // y2f = A20 * x0 + sum2

  int y0;        // y0 = floor(y0f);
  int y1;        // y1 = floor(y1f);
  int y2;        // y2 = floor(y2f);

  int copyInterpolate = !interpolateFlag;

  float dwork = 1.0/((float)x2Max);
  float work_perc = 0.0;

  for (x2=0; x2<x2Max; x2++,work_perc += dwork) {
    ShowWorking(work_perc);

    sum0 = A02*(x2) + b0;
    sum1 = A12*(x2) + b1;
    sum2 = A22*(x2) + b2;

    for (x1=0; x1<x1Max; x1++) {
      y0f = sum0;
      y1f = sum1;
      y2f = sum2;

      for (x0=0; x0<x0Max; x0++) {
	if (copyInterpolate) {

	  // trunc indices
	  y0 = (int)floor(y0f);
	  y1 = (int)floor(y1f);
	  y2 = (int)floor(y2f);
	  
	  // check for indices in bounds
	  if (y0>=0 && y0<y0Max &&
	      y1>=0 && y1<y1Max &&
	      y2>=0 && y2<y2Max) {
	    *outVolumeDataPtr++ = inVolPtrPtrPtr[y2][y1][y0]; // copy nearest voxel
	  } else
	    *outVolumeDataPtr++ = 0;

	} else { // use trilinear interpolation
#if 0
#else
#endif
	  // truncate indices
	  y0 = (int)floor(y0f);
	  y1 = (int)floor(y1f);
	  y2 = (int)floor(y2f);
	
	  // get data from volume 1
	  if (y0>=0 && y0<(y0Max-1) &&
	      y1>=0 && y1<(y1Max-1) &&
	      y2>=0 && y2<(y2Max-1)) {  // indices lie within the volume
	    computeFlag = 1;
	    dptr = &inVolPtrPtrPtr[y2][y1][y0];
	    d000 = *dptr;
	    d001 = dptr[1]; dptr+=y0Max;
	    d010 = *dptr;
	    d011 = dptr[1]; dptr+=y0Max*y1Max;
	    d110 = *dptr;
	    d111 = dptr[1]; dptr-=y0Max;
	    d100 = *dptr;
	    d101 = dptr[1];
	  } else { // at least one index is outside the volume
	    computeFlag = 0;
	    d000 = d001 = d010 = d011 = d100 = d101 = d110 = d111 = 0;
	    int y01 = y0+1;
	    int y11 = y1+1;
	    int y21 = y2+1;

	    int y0InBounds = (y0 >= 0 && y0 < y0Max);
	    int y1InBounds = (y1 >= 0 && y1 < y1Max);
	    int y2InBounds = (y2 >= 0 && y2 < y2Max);
	    int y01InBounds = (y01 >= 0 && y01 < y0Max);
	    int y11InBounds = (y11 >= 0 && y11 < y1Max);
	    int y21InBounds = (y21 >= 0 && y21 < y2Max);

	    // use indices that are in bounds
	    if (y2InBounds && y1InBounds && y0InBounds) {
	      d000 = inVolPtrPtrPtr[y2][y1][y0];
	      computeFlag = 1;
	    }
	    if (y2InBounds && y1InBounds && y01InBounds) {
	      d001 = inVolPtrPtrPtr[y2][y1][y01];
	      computeFlag = 1;
	    }
	    if (y2InBounds && y11InBounds && y0InBounds) {
	      d010 = inVolPtrPtrPtr[y2][y11][y0];
	      computeFlag = 1;
	    }
	    if (y2InBounds && y11InBounds && y01InBounds) {
	      d011 = inVolPtrPtrPtr[y2][y11][y01];
	      computeFlag = 1;
	    }
	    if (y21InBounds && y1InBounds && y0InBounds) {
	      d100 = inVolPtrPtrPtr[y21][y1][y0];
	      computeFlag = 1;
	    }
	    if (y21InBounds && y1InBounds && y01InBounds) {
	      d101 = inVolPtrPtrPtr[y21][y1][y01];
	      computeFlag = 1;
	    }
	    if (y21InBounds && y11InBounds && y0InBounds) {
	      d110 = inVolPtrPtrPtr[y21][y11][y0];
	      computeFlag = 1;
	    }
	    if (y21InBounds && y11InBounds && y01InBounds) {
	      d111 = inVolPtrPtrPtr[y21][y11][y01];
	      computeFlag = 1;
	    }
	  }

	  if (computeFlag) {
	    // do the trilinear interpolation
	    C1 = y0f - y0;
	    C2 = y1f - y1;
	    C3 = y2f - y2;

	    d00x = d000 + C1*(d001 - d000);
	    d01x = d010 + C1*(d011 - d010);
	    d10x = d100 + C1*(d101 - d100);
	    d11x = d110 + C1*(d111 - d110);
	    
	    d0xx = d00x + C2*(d01x - d00x);
	    d1xx = d10x + C2*(d11x - d10x);
	    
	    // set interpolated voxel 
	    *outVolumeDataPtr++ = (unsigned char)(d0xx + C3*(d1xx - d0xx));
	  } else 
	    // set voxel to zero
	    *outVolumeDataPtr++ = 0U;

	}  // end interpolation

	y0f += A00;
	y1f += A10;
	y2f += A20;
      }  // end x0 loop
      sum0 += A01;
      sum1 += A11;
      sum2 += A21;
    }  // end x1 loop
  }  // end x2 loop

  return 0;
}



////////////////////////////////////////////////////////////////////////////////
//
// Function: apply()
//
// Purpose: Apply mAffineMatrix and mTranslationVector to a list of points
//
// Inputs: 1) list of points to be transformed
//         2) list of points to put result
//
// Outputs: transformed list of points
//
// Returns: Returns 0 upon success
//          Returns 1 if either input list or output list does not have
//                       three columns
//          Returns 2 if number of points in lists are not equal
//          Returns 3 if not ready to apply
// 
////////////////////////////////////////////////////////////////////////////////
int
LinearTransformation::apply(const Matrix<double> &inPointList,
			    Matrix<double> &outPointList)
{
  if ((inPointList.getNcol() != 3) ||
      (outPointList.getNcol() != 3)) 
  {
    return 1;
  }
  if (inPointList.getNrow() != outPointList.getNrow())
  {
    return 2;
  }
  if ( ! isOkToApply() )
  {
    return 3;
  }

  // y = Ax + b ->  x = Ainverse * y - Ainverse * b
  Matrix<double> Ainv(3, 3);
  Vector<double> invTranslationVector(3);

  Ainv = MatrixUtils<double>::inv(mAffineMatrix);

  invTranslationVector = -1.0 * (Ainv * mTranslationVector);

//  outPointList = Ainv.transpose() * inPointList;
// KWD : bug (?) Matrix::multiply always fails 
//       since it should be inPointList * Ainv.transpose()

  Matrix<double> AinvT(3, 3);
  Ainv.transpose(&AinvT);

  MatrixUtils<double>::multiply(inPointList,AinvT,&outPointList);

  int numPnts = inPointList.getNrow();
  for (int i=0; i<numPnts; i++)
    for (int j=0; j<3; j++)
      outPointList[i][j] += invTranslationVector[j];

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//
// Function: apply()
//
// Purpose: Apply mAffineMatrix and mTranslationVector to a surface
//
// Inputs: 1) surface to be transformed
//         2) surface to put result
//
// Outputs: transformed surface
//
// Returns: Returns true (1) upon success, flase (0) on error
// 
////////////////////////////////////////////////////////////////////////////////
int LinearTransformation::apply(const Surface &inSurface,
				Surface *outSurface)
{

  if (getVerboseMode())
    cout << "apply() Surface in LinearTransformation" << endl;

  // make sure it is ok to apply transformation
  if (!isOkToApply())
    throwError("LinearTransformation::apply(Surface): cannot apply transformation");

  // copy surface to prepare for transformation
  *outSurface = inSurface;

  // compute inverse of affine matrix for apply to surface
  // first convert to floats
  Matrix<float> A(3,3);
  Vector<float> B(3);

  for (int i=0; i<3; i++)
  {
    B[i] = mTranslationVector[i];
    A[i][0] = mAffineMatrix[i][0];
    A[i][1] = mAffineMatrix[i][1];
    A[i][2] = mAffineMatrix[i][2];
  }

  Matrix<float> A_inv = MatrixUtils<float>::inv(A);
  Vector<float> B_inv = B * -1.0f;
  B_inv = A_inv*B_inv;

  // generate deformed surface
  outSurface->affine(A_inv, B_inv);

  return 1;
}

////////////////////////////////////////////////////////////////////////////////
//
// Function: apply()
//
// Purpose: Apply mAffineMatrix and mTranslationVector to a Point
//
// Inputs: 1) Point to be transformed
//         2) ptr to Point to put result
//
// Outputs: transformed Point
//
// Returns: Returns true (1) upon success, flase (0) on error
// 
////////////////////////////////////////////////////////////////////////////////
int LinearTransformation::apply(const Point &inPoint,
				Point *outPoint)
{

  if (getVerboseMode())
    cout << "apply() Point in LinearTransformation" << endl;

  // make sure it is ok to apply transformation
  if (!isOkToApply())
    throwError("LinearTransformation::apply(Point): cannot apply transformation");

  // kwd
  // y = Ax + b ->  x = Ainverse * y - Ainverse * b
  Matrix<double> Ainv(3, 3);
  Vector<double> invTranslationVector(3);

  Ainv = MatrixUtils<double>::inv(mAffineMatrix);
  invTranslationVector = -1.0 * (Ainv * mTranslationVector);

  Matrix<double> pMat(1,3);
  Matrix<double> nMat(1,3);

  pMat[0][0] = inPoint.x();
  pMat[0][1] = inPoint.y();
  pMat[0][2] = inPoint.z();

  nMat = pMat * Ainv.transpose();
  
  outPoint->set( nMat[0][0] + invTranslationVector[0],
		nMat[0][1] + invTranslationVector[1],
		nMat[0][2] + invTranslationVector[2]);
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//
// Function: meanSquaredError()
//
// Purpose: compute mean squared error between a transformed set of points
//          and another set of points
//
// Inputs: 1) list of points to be transformed
//         2) list of points for comparison
//
// Outputs: mean squared error
//
// Returns: Returns mean squared error
// 
////////////////////////////////////////////////////////////////////////////////
double
LinearTransformation::meanSquaredError(const Matrix<double> &atlasPoints,
				       const Matrix<double> &patientPoints)
{
#if 0
  if ((inPointList.getNcol() != 3) ||
      (outPointList.getNcol() != 3)) 
  {
    return 1;
  }
  if (inPointList.getNrow() != outPointList.getNrow())
  {
    return 2;
  }
  if ( ! isOkToApply() )
  {
    return 3;
  }
#endif

  Matrix<double> atlasTransformed(atlasPoints.getDim());
  apply(atlasPoints, atlasTransformed);

  atlasTransformed -= patientPoints;

  int numPnts = patientPoints.getNrow();
  double meanSquaredError = 0.0;
  for (int i=0; i<numPnts; i++)
    for (int j=0; j<3; j++)
      meanSquaredError += atlasTransformed[i][j]*atlasTransformed[i][j];
  meanSquaredError /= numPnts;

  return meanSquaredError;
}
