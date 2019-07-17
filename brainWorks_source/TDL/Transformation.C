//////////////////////////////////////////////////////////////////////////
//
// Transformation class
// 
// This file contains function definitions.
//
//////////////////////////////////////////////////////////////////////////
//
// Revision Log:
// 
// $Log: Transformation.C,v $
// Revision 1.1  2004/11/15 04:44:09  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:29:21  kem
// Initial revision
//
// Revision 1.2  1999/07/09 17:53:36  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.1  1997/08/01 19:50:10  csb
// Initial revision
//
// Revision 1.1  1997/06/25 17:36:18  csb
// Initial revision
//
//
//
//////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//
// Include Files
//
////////////////////////////////////////////////////////////////////////////////

#include <string.h>
#include <TDL/Transformation.h>

const char * const Transformation::mClassRevision = "Transformation.h,v 1.1";
const char * const Transformation::mClassName     = "Transformation";
const int          Transformation::mClassVersion  = 2;

////////////////////////////////////////////////////////////////////////////////
//
// Function: printClassParameters
//
// Purpose: Prints parameters that are owned by this class to the supplied ostream
//          It is envisioned that the child of this class calls this function and
//          then prints out it's own parameters.
//
// Inputs: ostream that we want to print to 
//
// Outputs: reference to the passed in ostream
//
////////////////////////////////////////////////////////////////////////////////
ostream & 
Transformation::printClassParameters ( ostream & outStream)
{
  outStream << mClassName << "::mClassVersion: "           << mClassVersion           << endl;
  outStream << mClassName << "::mClassRevision: "          << mClassRevision           << endl;
  outStream << mClassName << "::mVerboseModeFlag: "        << mVerboseModeFlag        << endl;
  outStream << mClassName << "::mOkToCalculateFlag: "      << mOkToCalculateFlag      << endl;
  outStream << mClassName << "::mOkToApplyFlag: "          << mOkToApplyFlag          << endl;

  return outStream;
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
// Outputs: Returns 0 upon success
//
////////////////////////////////////////////////////////////////////////////////
int
Transformation::saveClassParameters ( ofstream & outStream)
{
  // const char * functionName = "Transformation::saveClassParameters";
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

  // Write out parameters for this class
  swapint = big_endian_int(mVerboseModeFlag);
  outStream.write((char *) &swapint ,             sizeof (mVerboseModeFlag));
  swapint = big_endian_int(mOkToCalculateFlag);
  outStream.write((char *) & swapint,      sizeof (mOkToCalculateFlag));
  swapint = big_endian_int(mOkToApplyFlag);
  outStream.write((char *) & swapint,      sizeof (mOkToApplyFlag));

  return 0;
}
 

////////////////////////////////////////////////////////////////////////////////
//
// Function: loadClassParameters
//
// Purpose: Load parameters that are owned by this class
//
// Inputs: ifstream that we want to write to 
//
// Outputs: Returns 0 upon success and non-zero upon failure
//          Returns 1 if read-in class name is negative or > 1000
//          Returns 2 if read-in class name is not correct
//          Returns 3 if read-in class version is not correct
//          Returns 4 if cannot alloc memory
//          Returns 5 if read-in class revision is negative or > 1000
//
////////////////////////////////////////////////////////////////////////////////
int
Transformation::loadClassParameters ( ifstream & inStream)
{
  const char * const functionName = "Transformation::loadClassParameters";

  // Read in the length of the class name from the file
  int ReadInClassNameLength;
  inStream.read((char *) &ReadInClassNameLength, sizeof (ReadInClassNameLength));
  ReadInClassNameLength = big_endian_int(ReadInClassNameLength);
  // Let's not get crazy, even _I_ don't make class names this big
  if (ReadInClassNameLength < 0 || ReadInClassNameLength > 1000) { 
    cerr << "I don't think that this is a " << mClassName << " Results File." << endl;
    cerr << "Apparent class name length is " << ReadInClassNameLength << endl;
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
  // 
  // We can read version 1 safely for now [csb:19970624.1506MST]
  int ReadInClassVersion;
  inStream.read((char *) & ReadInClassVersion, sizeof (ReadInClassVersion));
  ReadInClassVersion = big_endian_int(ReadInClassVersion);
  if ( ! (    ReadInClassVersion == mClassVersion
           || ReadInClassVersion == 1             ) ) {
    cerr << functionName << ": ";
    cerr << "Attempted to read in " << mClassName << " Class parameters from version ";
    cerr << ReadInClassVersion << " whereas this code is version ";
    cerr << mClassVersion << endl;
    // Maybe should backup before returning
    // Throw an error here!
    return 3;
  }

  if (1 == ReadInClassVersion)
  {
    int TempInt;
    // Ok, things look good, let's read in the parameters that mClassName owns
    inStream.read((char *) & TempInt,            sizeof (TempInt)); // To simulate reading in mHaveTransformationFlag
    inStream.read((char *) & mVerboseModeFlag,   sizeof (mVerboseModeFlag));
    mVerboseModeFlag = big_endian_int(mVerboseModeFlag);
    inStream.read((char *) & mOkToCalculateFlag, sizeof (mOkToCalculateFlag));
    mOkToCalculateFlag = big_endian_int(mOkToCalculateFlag);
    inStream.read((char *) & mOkToApplyFlag,     sizeof (mOkToApplyFlag));
    mOkToApplyFlag = big_endian_int(mOkToApplyFlag);
  }
  else
  {
    // Read in the length of the class revision from the file
    int ReadInClassRevisionLength;
    inStream.read((char *) &ReadInClassRevisionLength, sizeof (ReadInClassRevisionLength));
    ReadInClassRevisionLength = big_endian_int(ReadInClassRevisionLength);
    // Let's not get crazy, even _I_ don't make class names this big
    if (ReadInClassRevisionLength < 0 || ReadInClassRevisionLength > 1000) { 
      cerr << "I don't think that this is a " << mClassName << " Results File." << endl;
      cerr << "Apparent class revision length is " << ReadInClassRevisionLength << endl;
      // Throw error here!
      return 5;
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


    // Ok, things look good, let's read in the parameters that mClassName owns
    inStream.read((char *) & mVerboseModeFlag,   sizeof (mVerboseModeFlag));
    mVerboseModeFlag = big_endian_int(mVerboseModeFlag);
    inStream.read((char *) & mOkToCalculateFlag, sizeof (mOkToCalculateFlag));
    mOkToCalculateFlag = big_endian_int(mOkToCalculateFlag);
    inStream.read((char *) & mOkToApplyFlag,     sizeof (mOkToApplyFlag));
    mOkToApplyFlag = big_endian_int(mOkToApplyFlag);
  } 
  //}
// Return with success
return 0;
}


///////////////////////////////////////////////////////////////////////
// Old code that I may want someday
///////////////////////////////////////////////////////////////////////

  // This is the size of all the parameters that we want to save
  // const int mClassSizeOfParameters = ( sizeof(size_t) // size returned from strlen(mClassName)
  //                                      + 14          // size "Transformation", there must be a better way
  //                                      + sizeof(int) // This variable
  //                                      + sizeof(mClassVersion) 
  //                                      + sizeof(int) // mHaveTransformationFlag
  //                                      + sizeof(int) // mVerboseModeFlag
  //                                      + sizeof(int) // mOkToCalculateFlag
  //                                      + sizeof(int) // mOkToApplyFlag
  //                                    );

//   outStream << mClassName << "::mClassSizeOfParameters: " << mClassSizeOfParameters << endl;

//  const int temp2 = mClassSizeOfParameters; // in memory and we need a pointer for write
//  outStream.write((const char *) & temp2,             sizeof (temp2));

  // // Ok, we are in the correct class version, double check that the size hasn't changed.
  // int ReadInSize;
  // inStream.read((char *) & ReadInSize, sizeof (ReadInSize));
  // if (ReadInSize != mClassSizeOfParameters) {
  //   cerr << functionName << ": ";
  //   cerr << "Size of file is " << ReadInSize;
  //   cerr << " whereas size of this class is " << mClassSizeOfParameters << endl;
  //   // Maybe should backup before returning
  //   // Throw an error here!
  //   return 4;
  // }

  // cout << "Class version from file: " << ReadInClassVersion << " == " << mClassVersion<< endl;
  // cout << "Class size from file: " << ReadInSize << " == " << mClassSizeOfParameters << endl;












