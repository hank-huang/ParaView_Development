///////////////////////////////////////////////////////////////////////////
//
// IntellX, L.L.C., Proprietary
// The contents of this file comprise confidential, proprietary and/or
// trade secret information which is the property of IntellX LLC.
// Dissemination, distribution or other disclosure of this information
// is strictly prohibited without the express permission of IntellX LLC.
// Copyright (c) 1998 IntellX, L.L.C.
// All rights reserved
//
// File: UserDefinedFieldTransformation.C
//
// Author: Rob Teichman
//
// Purpose: User Defined Field Transformation class body
//
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: UserDefinedFieldTransformation.C,v 1.1 2004/11/15 04:44:09 joeh Exp $
//
// Revision History
//
// $Log: UserDefinedFieldTransformation.C,v $
// Revision 1.1  2004/11/15 04:44:09  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:29:30  kem
// Initial revision
//
// Revision 1.3  1999/07/09 17:54:19  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.2  1998/12/10 19:06:05  RAZORDB
// IDL/AW merge
//
// Revision 1.1  1998/02/13 22:32:33  rst
// Initial revision
//
//
///////////////////////////////////////////////////////////////////////////

#include <TDL/UserDefinedFieldTransformation.h>

const char * const UserDefinedFieldTransformation::mClassRevision =     "$Id: UserDefinedFieldTransformation.C,v 1.1 2004/11/15 04:44:09 joeh Exp $";
const char * const UserDefinedFieldTransformation::mClassName = "UserDefinedFieldTransformation";
const int UserDefinedFieldTransformation::mClassVersion = 0;


//
// Constructor and Destructor
//
UserDefinedFieldTransformation::UserDefinedFieldTransformation()
{
  if (getVerboseMode())
    cout << "Calling UserDefinedFieldTransformation constructor" << endl;
}


UserDefinedFieldTransformation::~UserDefinedFieldTransformation()
{
  if (getVerboseMode())
    cout << "Calling UserDefinedFieldTransformation destructor" << endl;
}


///////////////////////////////////////////////////////////////////////////
//
// Private Methods
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//
// Function: getHField
//
// Purpose: To be determined
//
// Inputs: To be determined
//
// Outputs: To be determined 
//
///////////////////////////////////////////////////////////////////////////
void UserDefinedFieldTransformation::getHField(Array3D<float> *hField)
{

  if (getVerboseMode())
    cout << "getHField() in UserDefinedFieldTransformation" << endl;

  *hField = mHField;

}

void UserDefinedFieldTransformation::getHField(Array3D<float> *hField, 
	double Rx, double Ry, double Rz, double Ox, double Oy, double Oz)
{
if((hField->getXsize() == mHField.getXsize()) &&
   (hField->getYsize() == mHField.getYsize()) &&
   (hField->getZsize() == mHField.getZsize()) &&
   (Rx == _patientResX) && (Ry == _patientResY) && (Rz == _patientResZ) &&
   (Ox == 0.0) && (Oy ==0.0) && (Oz == 0.0)) {
	*hField = mHField;
	} else {

 // get fluid H field dimensions
  int hFieldSizeX = hField->getXsize()/3;
  int hFieldSizeY = hField->getYsize();
  int hFieldSizeZ = hField->getZsize();

  int THfilednX = mHField.getXsize()/3;
  int THfilednY = mHField.getYsize();
  int THfilednZ = mHField.getZsize();

  int x, y, z, k;
  float *hFieldPtr = hField->data();
  cerr << "slice     of " << hFieldSizeZ-1 << "\010\010\010\010\010\010\010";
  for (z = 0; z < hFieldSizeZ; z++) {
    cerr << "\010\010\010" ;
    cerr.width(3);
    cerr << z;
    for (y = 0; y < hFieldSizeY; y++) {
      for (x = 0; x < hFieldSizeX; x++) {
	float patientX = (Rx * x + Ox)/_patientResX;
        float patientY = (Ry * y + Oy)/_patientResY;
	float patientZ = (Rz * z + Oz)/_patientResZ;

	float Hx,Hy,Hz;

	if ((patientX >= 0) && (patientX < THfilednX-1) &&
            (patientY >= 0) && (patientY < THfilednY-1) &&
            (patientZ >= 0) && (patientZ < THfilednZ-1)) {
	trilinearHField(mHField, patientX, patientY, patientZ,
                          &Hx, &Hy, &Hz);
	} else {
	  Hx = patientX*_patientResX/_atlasResX;
          Hy = patientY*_patientResY/_atlasResY;
          Hz = patientZ*_patientResZ/_atlasResZ;
        }

	*hFieldPtr++ = Hx;
        *hFieldPtr++ = Hy;
        *hFieldPtr++ = Hz;
        }
      }
    }
  }
}



///////////////////////////////////////////////////////////////////////////
//
// Protected Methods
//
///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////
//
// Function: saveClassParameters
//
// Purpose: Saves parameters that are owned by this class to the
//          supplied ofstream.
//
// Inputs: ofstream that we want to save to 
//
// Outputs: 0 upon success
//          1 if parent's save returned non-zero
//
///////////////////////////////////////////////////////////////////////////
int UserDefinedFieldTransformation::saveClassParameters ( ofstream & outStream)
{

/************************************************
* Write out size of class name and the class name
************************************************/
  int swapint;
  const int nameLength = strlen(mClassName);
  swapint = big_endian_int(nameLength);
  outStream.write((char const *) &swapint, sizeof(nameLength));
  outStream.write((char const *) mClassName, nameLength);

/****************************
* Write out the class version
****************************/

  const int temp  = mClassVersion;       
  swapint = big_endian_int(temp);
  outStream.write((char const *) &swapint, sizeof(temp));

/********************************************************
* Write out size of class revision and the class revision
********************************************************/

  const int revisionLength = strlen(mClassRevision);
  swapint = big_endian_int(revisionLength);
  outStream.write((char const *) &swapint, sizeof(revisionLength));
  outStream.write((char const *) mClassRevision, revisionLength);

/********************************************************
* Write out class parameters
********************************************************/
  // none

/****************************
* Write out the vector field.
****************************/
   mHField.save(outStream); 

/*******************************
* Save parent class's parameters
*******************************/

  if ( FieldTransformation::saveClassParameters( outStream ) ) {
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
int UserDefinedFieldTransformation::loadClassParameters (ifstream & inStream)
{

  // Function name.
  const char * const functionName = "UserDefinedFieldTransformation::loadClassParameters";

/***************************************************
* Read in the length of the class name from the file
***************************************************/

  int ReadInClassNameLength;
  inStream.read((char *) &ReadInClassNameLength, sizeof(ReadInClassNameLength));

  ReadInClassNameLength = big_endian_int(ReadInClassNameLength);
  // Let's not get crazy, even _I_ don't make class names this big
  if (ReadInClassNameLength < 0 || ReadInClassNameLength > 1000) { 
    cerr << "I don't think that this is a " << mClassName << " Results File." << endl;
    cerr << "Apparent class length is " << ReadInClassNameLength << endl;
    // Throw error here!
    return 1;
  }

/***********************************************
* Alloc memory for the class name and read it in
***********************************************/

  char * ReadInClassName = new char[ReadInClassNameLength + 1];
  if ( ! ReadInClassName ) {
    cerr << functionName << ": ";
    cerr << "Failed to alloc " << ReadInClassNameLength+1 << " bytes." << endl;
    // Throw error here!
    // No use deleting ReadInClassName here, it didn't allocate correctly
    return 4;
  }
  inStream.read((char *) ReadInClassName, ReadInClassNameLength);
  ReadInClassName[ReadInClassNameLength] = '\0'; // End-of-string so I can use strcmp()

/*********************************************************
* Check to see if the file's class name is mClassName and 
* return with failure if not
*********************************************************/

  if (0 != strcmp(ReadInClassName, mClassName)) {
    cerr << functionName << ": ";
    cerr << "Attempted to read in Class " << ReadInClassName;
    cerr << " whereas this Class is " << mClassName << endl;
    // Maybe should backup inStream before returning
    // Maybe should skip the rest of this record before returning, probally not
    // Throw an error here!
    delete ReadInClassName;
    return 2;
  }

/********************************************************************
* Ok, we are in the correct class, let's check if we have the correct
* version of mClassName
********************************************************************/

  int ReadInClassVersion;
  inStream.read((char *) &ReadInClassVersion, sizeof(ReadInClassVersion));
  ReadInClassVersion = big_endian_int(ReadInClassVersion);
  if (ReadInClassVersion != mClassVersion) {
    cerr << functionName << ": ";
    cerr << "Attempted to read in " << mClassName << " Class parameters from version ";
    cerr << ReadInClassVersion << " whereas this code is version ";
    cerr << mClassVersion << endl;
    // Maybe should backup before returning
    // Throw an error here!
    delete ReadInClassName;
    return 3;
  }

/*******************************************************
* Read in the length of the class revision from the file
*******************************************************/

  int ReadInClassRevisionLength;
  inStream.read((char *) &ReadInClassRevisionLength, sizeof(ReadInClassRevisionLength));
  ReadInClassRevisionLength = big_endian_int(ReadInClassRevisionLength);
// Let's not get crazy, even _I_ don't make class names this big
  if (ReadInClassRevisionLength < 0 || ReadInClassRevisionLength > 1000) { 
    cerr << "I don't think that this is a " << mClassName << " Results File." << endl;
    cerr << "Apparent class revision length is " << ReadInClassRevisionLength << endl;
    // Throw error here!
    delete ReadInClassName;
    return 6;
  }

/***************************************************
* Alloc memory for the class revision and read it in
***************************************************/

  char * ReadInClassRevision = new char[ReadInClassRevisionLength + 1];
  if ( ! ReadInClassRevision ){
    cerr << functionName << ": ";
    cerr << "Failed to alloc " << ReadInClassRevisionLength+1 << " bytes." << endl;
    // Throw error here!
    // No use deleting ReadInClassRevision here because it didn't malloc
    delete ReadInClassName;
    return 4;
  }
  inStream.read((char *) ReadInClassRevision, ReadInClassRevisionLength);

// End-of-string so I can use strcmp()
// Do nothing with it for now, mainly for debugging outside class

  ReadInClassRevision[ReadInClassRevisionLength] = '\0'; 

/**********************************
* Read in parameters for this class
**********************************/

// none
 
/**************************
* Read in the vector field.
**************************/

  mHField.load(inStream); 
 
/*******************************
* Load in base class parameters.
*******************************/

  if ( FieldTransformation::loadClassParameters( inStream ) ) {
    delete ReadInClassName;
    delete ReadInClassRevision;
    return 5;
  };

  // Return with success
  delete ReadInClassName;
  delete ReadInClassRevision;
  return 0;

}


///////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////
int UserDefinedFieldTransformation::save(const char *filename)
{
  if (getVerboseMode())
    cout << "save() in UserDefinedFieldTransformation" << endl;

  ofstream outStream (filename, ios::out);

  if ( ! outStream ) {
    cerr << "Failed to open file \"" << filename << "\"" << endl;
    return 1;
  }
  if ( saveClassParameters( outStream ) ) {
    return 2;
  };
  outStream.close();

  return 0;

}



///////////////////////////////////////////////////////////////////////////
//
// Function: printClassParameters
//
// Purpose: Calls parent's print function and then prints parameters that are
//          owned by this class to the supplied ostream.
//          It is envisioned that the child of this class calls this function
//          and then prints out it's own parameters.
//
// Inputs: ostream that we want to print to 
//
// Outputs: reference to the passed in ostream
//
///////////////////////////////////////////////////////////////////////////
ostream & UserDefinedFieldTransformation::printClassParameters (ostream &outStream)
{

  FieldTransformation::printClassParameters (outStream);
  outStream << mClassName << "::mClassVersion:              " << mClassVersion << endl;
  outStream << mClassName << "::mClassRevision:             " << mClassRevision << endl;

  return outStream;

}



///////////////////////////////////////////////////////////////////////////
//
// Public Methods
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////
int UserDefinedFieldTransformation::load(const char *filename)
{

  if (getVerboseMode())
    cout << "load() in UserDefinedFieldTransformation" << endl;

  ifstream inStream (filename, ios::in);
  if ( ! inStream ) {
    cerr << "Failed to open file \"" << filename << "\"" << endl;
    return 1;
  }

  if ( loadClassParameters( inStream ) ) {
    return 2;
  }

  inStream.close();

  return 0;

}



///////////////////////////////////////////////////////////////////////////
//
// Function: calculate
//
// Purpose: 
//
// Inputs: 
//
// Outputs: Turns mOkToApplyFlag to 1 if successfull
//
///////////////////////////////////////////////////////////////////////////
int UserDefinedFieldTransformation::calculate()
{

  if (getVerboseMode())
    cout << "calculate() in UserDefinedFieldTransformation" << endl;

  // calculate does nothing, so flags are set in initialize
  cout << "calculate() for UserDefinedFieldTransformation is empty." << endl;
  // Check if its ok to calculate.
//  if (!isOkToCalculate())
//    throwError("UserDefinedFieldTransformation::calculate: Unable to calculate trasformation.");

//  setOkToApply( 1 );

  return 0;

}



///////////////////////////////////////////////////////////////////////////
//
// Function: initialize
//
// Purpose: Initialize class so that calculate may run
//
// Inputs: To be determined
//
// Outputs: Sets mOkToCalculateFlag to 1 if successfull
//
///////////////////////////////////////////////////////////////////////////
int
UserDefinedFieldTransformation::initialize(Array3D<float> const &hField)
{
  if (getVerboseMode())
    cout << "initialize() in UserDefinedFieldTransformation" << endl;

  // copy hField into member variable
  mHField = hField;

  // set FieldTransformation class parameters of interest
  _atlasSizeX = _patientSizeX = hField.getXsize()/3;
  _atlasSizeY = _patientSizeY = hField.getYsize();
  _atlasSizeZ = _patientSizeZ = hField.getZsize();


  // Done with initialize.
  setOkToCalculate( 1 );

  // since calculate does nothing, set okToApply here also.
  setOkToApply( 1 );

  return 0;

}


///////////////////////////////////////////////////////////////////////////
//
// Function: print
//
// Purpose: Prints parameters
//
// Inputs: To be determined
//
// Outputs: To be determined
//
///////////////////////////////////////////////////////////////////////////
void UserDefinedFieldTransformation::print()
{
  if (getVerboseMode())
    cout << "print() in UserDefinedFieldTransformation" << endl;

  UserDefinedFieldTransformation::printClassParameters(cout);

}
