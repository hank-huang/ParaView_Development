#include <string.h>
#include <TDL/UserDefinedLinearTransformation.h>

//////////////////////////////////////////////////////////////////////////
//
//
//////////////////////////////////////////////////////////////////////////
//
// Revision Log:
//
// $Log: UserDefinedLinearTransformation.C,v $
// Revision 1.1  2004/11/15 04:44:09  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:29:32  kem
// Initial revision
//
// Revision 1.3  1999/07/09 17:54:31  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.2  1997/10/20 16:34:27  csb
// reverse order of save and load so same as others
//
// Revision 1.1  1997/08/01 19:50:11  csb
// Initial revision
//
// Revision 1.2  1997/08/01 15:31:05  csb
// Incorportate Abed's, fluid, elastic and field
//
// Revision 1.1  1997/06/25 17:39:48  kem
// Initial revision
//
//
//
//////////////////////////////////////////////////////////////////////////

const char * const UserDefinedLinearTransformation::mClassName = "UserDefinedLinearTransformation";
const int          UserDefinedLinearTransformation::mClassVersion = 1;


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
UserDefinedLinearTransformation::saveClassParameters ( ofstream & outStream)
{
  // const char * functionName = "UserDefinedLinearTransformation::saveClassParameters";

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
UserDefinedLinearTransformation::loadClassParameters ( ifstream & inStream)
{
  const char * const functionName = "UserDefinedLinearTransformation::loadClassParameters";

  // Read in the length of the class name from the file
  int ReadInClassNameLength;
  inStream.read((char *) &ReadInClassNameLength, sizeof (ReadInClassNameLength));
  ReadInClassNameLength = big_endian_int(ReadInClassNameLength);
  // Let's not get crazy, even _I_ don't make class names this big
  if (ReadInClassNameLength < 0 || ReadInClassNameLength > 1000) { 
    cerr << "I don't think that this is a UserDefinedLinearTransformation Results File." << endl;
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
  ReadInClassVersion=big_endian_int(ReadInClassVersion);
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
UserDefinedLinearTransformation::save( const char * filename )
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

  // cout << "save() in UserDefinedLinearTransformation" << endl;
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
UserDefinedLinearTransformation::printClassParameters ( ostream & outStream)
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
UserDefinedLinearTransformation::load( const char * filename )
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

  // cout << "load() in UserDefinedLinearTransformation" << endl;
  inStream.close();

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//
// Function: calculate
//
// Purpose: This method does nothing since the linear transformation must
//          be defined by the user
//
// Inputs: none
//
// Outputs: Returns 0 if successful
//
////////////////////////////////////////////////////////////////////////////////
int
UserDefinedLinearTransformation::calculate()
{
  return 0;
}


////////////////////////////////////////////////////////////////////////////////
//
// Function: initialize
//
// Purpose: This method does nothing since the linear transformation must
//          be defined by the user
//
// Inputs: none
//
// Outputs: Returns 0 if successful
//
////////////////////////////////////////////////////////////////////////////////
int
UserDefinedLinearTransformation::initialize() 
{
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
UserDefinedLinearTransformation::print()
{
  cout << "print() in UserDefinedLinearTransformation" << endl;
  UserDefinedLinearTransformation::printClassParameters(cout);
}


////////////////////////////////////////////////////////////////////////////////
//
// Function: setTransformation
//
// Purpose: This method sets the rotation and translation matrices of
//          the linear transformation
//
// Inputs: rotation matrix
//         translation vector
//
// Outputs: none
//
////////////////////////////////////////////////////////////////////////////////
void
UserDefinedLinearTransformation::setTransformation(const Matrix<double> &A,
						   const Vector<double> &b)
{
  mAffineMatrix = A;
  mTranslationVector = b;
  setOkToApply( 1 );
}
