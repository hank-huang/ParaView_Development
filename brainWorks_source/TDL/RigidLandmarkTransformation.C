#include <string.h>
#include <TDL/RigidLandmarkTransformation.h>

//////////////////////////////////////////////////////////////////////////
//
//
//////////////////////////////////////////////////////////////////////////
//
// Revision Log:
//
// $Log: RigidLandmarkTransformation.C,v $
// Revision 1.1  2004/11/15 04:44:08  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:29:24  kem
// Initial revision
//
// Revision 1.4  1999/07/09 17:53:27  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.3  1998/12/10 19:05:49  RAZORDB
// IDL/AW merge
//
// Revision 1.2  1997/09/12 20:52:42  kem
// Add check for zero patient voxel scales in initialize()
//
// Revision 1.1  1997/08/01 19:50:08  csb
// Initial revision
//
// Revision 1.2  1997/08/01 15:06:53  csb
// Incorportate Abed's, fluid, elastic and field
//
// Revision 1.1  1997/06/25 17:39:47  kem
// Initial revision
//
//
//
//////////////////////////////////////////////////////////////////////////

const char * const RigidLandmarkTransformation::mClassRevision = "$Id: RigidLandmarkTransformation.C,v 1.1 2004/11/15 04:44:08 joeh Exp $";
const char * const RigidLandmarkTransformation::mClassName = "RigidLandmarkTransformation";
const int          RigidLandmarkTransformation::mClassVersion = 1;


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
RigidLandmarkTransformation::saveClassParameters ( ofstream & outStream)
{
  // const char * functionName = "RigidLandmarkTransformation::saveClassParameters";
  int swapint;
  // Write out size of class name and the class name
  const int nameLength = strlen(mClassName);
  swapint = big_endian_int(nameLength);
  outStream.write((const char *) & swapint, sizeof(nameLength));
  outStream.write(mClassName,                  nameLength);

   // Write out the class version
  const int temp  = mClassVersion;          // Because the compiler may not put these variables in memory
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
RigidLandmarkTransformation::loadClassParameters ( ifstream & inStream)
{
  const char * const functionName = "RigidLandmarkTransformation::loadClassParameters";

  // Read in the length of the class name from the file
  int ReadInClassNameLength;
  inStream.read((char *) &ReadInClassNameLength, sizeof (ReadInClassNameLength));
  ReadInClassNameLength = big_endian_int(ReadInClassNameLength);
  // Let's not get crazy, even _I_ don't make class names this big
  if (ReadInClassNameLength < 0 || ReadInClassNameLength > 1000) { 
    cerr << "I don't think that this is a RigidLandmarkTransformation Results File." << endl;
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
RigidLandmarkTransformation::save( const char * filename )
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

  // cout << "save() in RigidLandmarkTransformation" << endl;
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
RigidLandmarkTransformation::printClassParameters ( ostream & outStream)
{
  LinearTransformation::printClassParameters ( outStream );
  outStream << mClassName << "::mClassVersion: " << mClassVersion << endl;
  outStream << mClassName << "::mClassRevision: " << mClassRevision << endl;

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
RigidLandmarkTransformation::load( const char * filename )
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

  // cout << "load() in RigidLandmarkTransformation" << endl;
  inStream.close();

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//
// Function: calculate
//
// Purpose: This method computes a rigid transformation that 'minimizes
//          squared error' (supposedly)
//
// Inputs:  Two sets of landmark points
//
// Outputs: Turns mOkToApplyFlag to 1 if successfull
//          Returns 0 if successful
//          Returns 1 if transformation not initialized
//
////////////////////////////////////////////////////////////////////////////////
int
RigidLandmarkTransformation::calculate()
{
  if (! isOkToCalculate() )
    return 1;

  Vector<double> patientMean(3);
  Vector<double> atlasMean(3);


  // should do some error checking!!!


  // use ADL routine to find transformation
  ShowWorking("Working...");
  MatrixUtils<double>::findTransform(*_patientPoints, *_atlasPoints, 
		&mAffineMatrix, &patientMean, &atlasMean);

  ShowWorking("Working...");
  // scale transformation parameters
  mAffineMatrix = _scaleMatrix * mAffineMatrix;

  // compute translation vector from ADL result
  mTranslationVector = atlasMean - mAffineMatrix * patientMean;

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
//         The voxel scales for the atlas and the patient
//
// Outputs: Sets mOkToCalculateFlag to 1 if successfull
//          Returns 0 if successful
//          Returns 1 if matrices have different dimensions
//          Returns 2 if matrices do not have 3 columns (3-D points)
//          Returns 3 if matrices do not have at least 4 rows (4 points)
//
////////////////////////////////////////////////////////////////////////////////
int
RigidLandmarkTransformation::initialize(const Matrix<double> * const atlasPoints, 
					const Matrix<double> * const patientPoints,
					const float &atlasVoxelScaleX,
					const float &atlasVoxelScaleY,
					const float &atlasVoxelScaleZ,
					const float &patientVoxelScaleX,
					const float &patientVoxelScaleY,
					const float &patientVoxelScaleZ)
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

  if ((atlasVoxelScaleX != 0.0) &&
      (atlasVoxelScaleY != 0.0) &&
      (atlasVoxelScaleZ != 0.0) &&
      (patientVoxelScaleX != 0.0) &&
      (patientVoxelScaleY != 0.0) &&
      (patientVoxelScaleZ != 0.0)) {
    _scaleMatrix[0][0] = patientVoxelScaleX/atlasVoxelScaleX;
    _scaleMatrix[1][1] = patientVoxelScaleY/atlasVoxelScaleY;
    _scaleMatrix[2][2] = patientVoxelScaleZ/atlasVoxelScaleZ;
  } else {
    _scaleMatrix.eye();
  }

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
RigidLandmarkTransformation::print()
{
  // cout << "print() in RigidLandmarkTransformation" << endl;
  RigidLandmarkTransformation::printClassParameters(cout);
}
