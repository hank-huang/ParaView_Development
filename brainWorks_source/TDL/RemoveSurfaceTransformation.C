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
// File: RemoveSurfaceTransformation.C
//
// Author: Rob Teichman
//
// Purpose: Remove Surface Transformation class body
//
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: RemoveSurfaceTransformation.C,v 1.1 2004/11/15 04:44:08 joeh Exp $
//
// Revision History
//
// $Log: RemoveSurfaceTransformation.C,v $
// Revision 1.1  2004/11/15 04:44:08  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:29:33  kem
// Initial revision
//
// Revision 1.2  1999/07/09 17:53:09  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.1  1998/05/18 19:25:16  rst
// Initial revision
//
//
///////////////////////////////////////////////////////////////////////////

#include <TDL/RemoveSurfaceTransformation.h>
#include <TDL/Embed.h>

const char * const RemoveSurfaceTransformation::mClassRevision = "RemoveSurfaceTransformation,v1.1";
const char * const RemoveSurfaceTransformation::mClassName = "RemoveSurfaceTransformation";
const int RemoveSurfaceTransformation::mClassVersion = 0;

//
// Constructor and Destructor
//
RemoveSurfaceTransformation::RemoveSurfaceTransformation()
{
  if (getVerboseMode())
    cout << "Calling RemoveSurfaceTransformation constructor" << endl;
}


RemoveSurfaceTransformation::~RemoveSurfaceTransformation()
{
  if (getVerboseMode())
    cout << "Calling RemoveSurfaceTransformation destructor" << endl;
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
void RemoveSurfaceTransformation::getHField(Array3D<float> *hField, 
	double Rx, double Ry, double Rz, double Ox, double Oy, double Oz)
{

  if (getVerboseMode())
    cout << "getHField() in RemoveSurfaceTransformation" << endl;

  *hField = mHField;
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
int RemoveSurfaceTransformation::saveClassParameters ( ofstream & outStream)
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
  Point swappoint = mSeed;
  double tempz, tempy, tempx;
  tempx = big_endian_double(mSeed.x());
  tempy = big_endian_double(mSeed.y());
  tempz = big_endian_double(mSeed.z());
  swappoint.set(tempx,tempy,tempz);
  outStream.write((char const *) &swappoint, sizeof(mSeed));

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
int RemoveSurfaceTransformation::loadClassParameters (ifstream & inStream)
{

  // Function name.
  const char * const functionName = "RemoveSurfaceTransformation::loadClassParameters";

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

  inStream.read((char *) &mSeed, sizeof(mSeed));
  double tempx,tempy,tempz;
  tempx = big_endian_double(mSeed.x());
  tempy = big_endian_double(mSeed.y());
  tempz = big_endian_double(mSeed.z());
  mSeed.set(tempx,tempy,tempz);

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
int RemoveSurfaceTransformation::save(const char *filename)
{
  if (getVerboseMode())
    cout << "save() in RemoveSurfaceTransformation" << endl;

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
ostream & RemoveSurfaceTransformation::printClassParameters (ostream &outStream)
{

  FieldTransformation::printClassParameters (outStream);
  outStream << mClassName << "::mClassVersion:              " << mClassVersion << endl;
  outStream << mClassName << "::mClassRevision:             " << mClassRevision << endl;

  outStream << mClassName << "::mSeed:                      " << mSeed << endl;

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
int RemoveSurfaceTransformation::load(const char *filename)
{

  if (getVerboseMode())
    cout << "load() in RemoveSurfaceTransformation" << endl;

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
int RemoveSurfaceTransformation::calculate()
{
  const int MAXITS = 150;
  const float EPS = 1.e-10;
  const float rjac = 0.95;
  const float alpha = 0.0001;

  if (getVerboseMode())
    cout << "calculate() in RemoveSurfaceTransformation" << endl;

  // Check if its ok to calculate.
  if (!isOkToCalculate())
    throwError("RemoveSurfaceTransformation::calculate: Unable to calculate trasformation.");

  //
  // build mask volume from supplied surface
  //
  int xsize = mBoneMask.getXsize(),
    ysize = mBoneMask.getYsize(),
    zsize = mBoneMask.getZsize();

  // compute seed if necessary
  if (mSeed.x() == 0 && mSeed.y() == 0 && mSeed.z() == 0)
    if (mSurf.getCentroid(&mSeed) != 0)
    {
      cout << "getCentriod failed" << endl;
      exit(-1);
    }

  if (getVerboseMode())
  {
    cout << "Building mask of surface in volume size "
	 << xsize << "x" << ysize << "x" << zsize << endl;
    cout << " using seed " << mSeed << endl;
  }

  Embed maskVol(zsize, ysize, xsize, mSurf, mSeed);

  if (getVerboseMode())
    cout << "done building mask" << endl;

  //
  // compute h-field
  //
  int n, i, j, k;
  unsigned char *mask = maskVol.data();
  unsigned char *bone = mBoneMask.data();

  if (getVerboseMode())
    cout << "allocating space for h-field" << endl;

  mHField.setDim(zsize, ysize, xsize*3);
  if (mHField.isEmpty())
  {
    cerr << "Unable to allocate memory for hField." << endl;
    exit(-1);
  }
  mHField = 0;
  float *uh = mHField.data();

  // initialize the deformation field by setting interior to seed
  if (getVerboseMode())
    cout << "initializing h-field" << endl;

  for (i=1; i<zsize-1; i++)
    for (j=1; j<ysize-1; j++)
      for (k=1; k<xsize-1; k++)
	if (mask[i*xsize*ysize+j*xsize+k] >= 30 )
	{
	  uh[3*(i*xsize*ysize+j*xsize+k)] = (float)k - mSeed.x();
	  uh[3*(i*xsize*ysize+j*xsize+k)+1] = (float)j - mSeed.y();
	  uh[3*(i*xsize*ysize+j*xsize+k)+2] = (float)i - mSeed.z();
	}

  // Now start the Laplacian Relaxation
  // This uses successive over-relaxation with zero boundary values
  if (getVerboseMode())
    cout << "Starting SOR" << endl;

  float omega = 1.0;
  float resx, resy, resz, anorm;

  for (n=0; n<MAXITS; n++)
  {
    anorm = 0.0;
    for (i=1; i<zsize-1; i++)
      for (j=1; j<ysize-1; j++)
	for (k=1; k<xsize-1; k++) 
	  /* Odd-Even ordering with the Tumor boundry condition */
	  // use 127 below instead of 255 to include boundary after fill
	  if ( ((i+j+k)%2 == n%2) && mask[i*xsize*ysize+j*xsize+k] < 127  &&
	       bone[i*xsize*ysize+j*xsize+k] == 0)
	  {
	    resx = uh[3*((i+1)*xsize*ysize+j*xsize+k)] +
	      uh[3*((i-1)*xsize*ysize+j*xsize+k)] +
	      uh[3*(i*xsize*ysize+(j+1)*xsize+k)] +
	      uh[3*(i*xsize*ysize+(j-1)*xsize+k)] +
	      uh[3*(i*xsize*ysize+j*xsize+(k+1))] +
	      uh[3*(i*xsize*ysize+j*xsize+(k-1))] -
	      (6+alpha)*uh[3*(i*xsize*ysize+j*xsize+k)]; 

	    resy = uh[3*((i+1)*xsize*ysize+j*xsize+k)+1] +
	      uh[3*((i-1)*xsize*ysize+j*xsize+k)+1] +
	      uh[3*(i*xsize*ysize+(j+1)*xsize+k)+1] +
	      uh[3*(i*xsize*ysize+(j-1)*xsize+k)+1] +
	      uh[3*(i*xsize*ysize+j*xsize+(k+1))+1] +
	      uh[3*(i*xsize*ysize+j*xsize+(k-1))+1] -
	      (6+alpha)*uh[3*(i*xsize*ysize+j*xsize+k)+1]; 

	    resz = uh[3*((i+1)*xsize*ysize+j*xsize+k)+2] +
	      uh[3*((i-1)*xsize*ysize+j*xsize+k)+2] +
	      uh[3*(i*xsize*ysize+(j+1)*xsize+k)+2] +
	      uh[3*(i*xsize*ysize+(j-1)*xsize+k)+2] +
	      uh[3*(i*xsize*ysize+j*xsize+(k+1))+2] +
	      uh[3*(i*xsize*ysize+j*xsize+(k-1))+2] -
	      (6+alpha)*uh[3*(i*xsize*ysize+j*xsize+k)+2]; 
	
	    anorm += fabs(resx) + fabs(resy) + fabs(resz);

	    uh[3*(i*xsize*ysize+j*xsize+k)] -= omega*resx/(-(6+alpha));	
	    uh[3*(i*xsize*ysize+j*xsize+k)+1] -= omega*resy/(-(6+alpha));	
	    uh[3*(i*xsize*ysize+j*xsize+k)+2] -= omega*resz/(-(6+alpha));	
	  }

    omega = (n == 0 ? 1.0/(1.0-0.5*rjac*rjac) :
	     1.0/(1.0-0.25*omega*rjac*rjac) );
    cout << "(n, omega, anorm) = (" << n << ", " << omega << ", "
	 << anorm << ")" << endl;
  }

  cout << "Final anorm = " << anorm << endl;

  // convert the u-field to h-field
  for (i=0; i<zsize; i++)
    for (j=0; j<ysize; j++)
      for (k=0; k<xsize; k++)
      {
	uh[3*(i*xsize*ysize+j*xsize+k)] = k-uh[3*(i*xsize*ysize+j*xsize+k)];
	uh[3*(i*xsize*ysize+j*xsize+k)+1] = j-uh[3*(i*xsize*ysize+j*xsize+k)+1];
	uh[3*(i*xsize*ysize+j*xsize+k)+2] = i-uh[3*(i*xsize*ysize+j*xsize+k)+2];
      }

  cout << "Deformation field generated" << endl;

  setOkToApply( 1 );

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
int RemoveSurfaceTransformation::initialize(Surface const &surf,
				GreyVolume<unsigned char> const &boneMask,
				double seedx,
			      	double seedy,
			       	double seedz)
{
  if (getVerboseMode())
    cout << "initialize() in RemoveSurfaceTransformation" << endl;

  // verify that the boneMask supplied is valid

  // take the size of the boneMask to be the size of the data volume
  // and set FieldTransformation class parameters of interest
  _atlasSizeX = _patientSizeX = boneMask.getXsize();
  _atlasSizeY = _patientSizeY = boneMask.getYsize();
  _atlasSizeZ = _patientSizeZ = boneMask.getZsize();

  // set mSeed
  mSeed.set(seedx, seedy, seedz);

  // store surface and boneMask
  mSurf = surf;
  mBoneMask = boneMask;

  // Done with initialize.
  setOkToCalculate( 1 );

  return 0;
}



///////////////////////////////////////////////////////////////////////////
//
// Function: setSeed
//
// Purpose: set seed value
//
// Inputs: To be determined
//
// Outputs: returns 0 on success
//
///////////////////////////////////////////////////////////////////////////
int
RemoveSurfaceTransformation::setSeed(Point const &seed)
{
  if (getVerboseMode())
    cout << "setSeed() in RemoveSurfaceTransformation" << endl;

  mSeed = seed;

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
void RemoveSurfaceTransformation::print()
{
  if (getVerboseMode())
    cout << "print() in RemoveSurfaceTransformation" << endl;

  RemoveSurfaceTransformation::printClassParameters(cout);

}
