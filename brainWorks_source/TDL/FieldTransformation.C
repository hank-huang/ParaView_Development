///////////////////////////////////////////////////////////////////////////
//
// IntellX, L.L.C., Proprietary
// Do not reproduce without permission in writing.
// Copyright (c) 1997 IntellX, L.L.C.
// All rights reserved
//
// File: FieldTransformation.C
//
// Author: IntellX
//
// Purpose: This file contains the function definitions (bodies)
//  for the FieldTransformation class
// 
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: FieldTransformation.C,v 1.1 2004/11/15 04:44:07 joeh Exp $
//
// Revision History
//
// $Log: FieldTransformation.C,v $
// Revision 1.1  2004/11/15 04:44:07  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:29:16  kem
// Initial revision
//
// Revision 1.18  1999/07/09 17:51:27  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.17  1999/01/08 17:30:08  RAZORDB
// Remove _haveHfield in apply() methods
//
// Revision 1.16  1999/01/07 16:25:16  RAZORDB
// Minor fixes (rst)
//
// Revision 1.15  1998/12/18 18:15:29  RAZORDB
// Change to use TransformationUtils::trilinear()
//
// Revision 1.14  1998/12/10 19:04:42  RAZORDB
// IDL/AW merge
//
// Revision 1.13  1998/09/14 20:18:13  csb
// Fix nearest neighbor so it doesn't just '(int)'
//
// Revision 1.12  1998/05/06 19:10:07  rst
// change invApply to use interpolation method
//
// Revision 1.11  1998/04/17 14:53:17  kem
// Add inverse apply for a point in patient space
//
// Revision 1.10  1998/04/09 19:10:41  kem
// Change apply() for volumes
//
// Revision 1.9  1998/03/05 18:04:02  rst
// new applys
//
// Revision 1.8  1997/12/12 21:55:19  csb
// merge 1.6.1.2 and 1.7
//
// Revision 1.7  1997/12/12 19:32:57  abed
// Incorporating resolution changes
//
// Revision 1.6  1997/10/15 16:29:19  kem
// Use interpolateFlag in apply()
//
// Revision 1.5  1997/10/02 20:37:47  rst
// Fixed apply() routines for Surfaces and Points
//
// Revision 1.4  1997/09/16 20:52:43  kem
// Add atlas and patient dimension variables
//
// Revision 1.3  1997/08/25 17:36:43  kem
// Change apply() method to use getHField()
//
// Revision 1.2  1997/08/22 20:56:57  rst
// Add apply() for Surfaces and Points
//
// Revision 1.1  1997/08/01 19:50:02  csb
// Initial revision
//
// Revision 1.2  1997/08/01 15:37:39  csb
// Incorportate Abed's, fluid, elastic and field
//
// Revision 1.1  1997/06/25 17:36:16  csb
// Initial revision
//
//
//////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//
// Include Files
//
////////////////////////////////////////////////////////////////////////////////

#include <string.h>
#include <TDL/FieldTransformation.h>
#include <TDL/TransformationUtils.h>

const char * const FieldTransformation::mClassRevision = "$Id: FieldTransformation.C,v 1.1 2004/11/15 04:44:07 joeh Exp $";
const char * const FieldTransformation::mClassName = "FieldTransformation";
const int          FieldTransformation::mClassVersion = 3;

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
//          then saves it's own parameters.
//
// Inputs: ofstream that we want to save to 
//
// Outputs: Returns 0 if successfull
//          Returns 1 if call to parent fails
//
////////////////////////////////////////////////////////////////////////////////
int
FieldTransformation::saveClassParameters ( ofstream & outStream)
{
  int swapint;
  double swapdouble;
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

  // Write out parameters for this class
  swapint = big_endian_int(_atlasSizeX);
  outStream.write((const char *) &swapint, sizeof(_atlasSizeX));
  
  swapint = big_endian_int(_atlasSizeY);
  outStream.write((const char *) &swapint, sizeof(_atlasSizeY));
  
  swapint = big_endian_int(_atlasSizeZ);
  outStream.write((const char *) &swapint, sizeof(_atlasSizeZ));
  
  swapint = big_endian_int(_patientSizeX);
  outStream.write((const char *) &swapint, sizeof(_patientSizeX));
  
  swapint = big_endian_int(_patientSizeY);
  outStream.write((const char *) &swapint, sizeof(_patientSizeY));
  
  swapint = big_endian_int(_patientSizeZ);
  outStream.write((const char *) &swapint, sizeof(_patientSizeZ));
  
  swapdouble = big_endian_double(_atlasResX);
  outStream.write((const char *) &swapdouble, sizeof(_atlasResX));
  
  swapdouble = big_endian_double(_atlasResY);
  outStream.write((const char *) &swapdouble, sizeof(_atlasResY));
  
  swapdouble = big_endian_double(_atlasResZ);
  outStream.write((const char *) &swapdouble, sizeof(_atlasResZ));
  
  swapdouble = big_endian_double(_patientResX);
  outStream.write((const char *) &swapdouble, sizeof(_patientResX));
  
  swapdouble = big_endian_double(_patientResY);
  outStream.write((const char *) &swapdouble, sizeof(_patientResY));
  
  swapdouble = big_endian_double(_patientResZ);
  outStream.write((const char *) &swapdouble, sizeof(_patientResZ));
 

  // Save parent class's parameters
  if ( Transformation::saveClassParameters( outStream ) )
  {
    return 1;
  };
 
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
FieldTransformation::loadClassParameters ( ifstream & inStream)
{
  const char * const functionName = "FieldTransformation::loadClassParameters";

  // Read in the length of the class name from the file
  int ReadInClassNameLength;
  inStream.read((char *) &ReadInClassNameLength, sizeof (ReadInClassNameLength));
  // Let's not get crazy, even _I_ don't make class names this big
  ReadInClassNameLength = big_endian_int(ReadInClassNameLength);
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
  // Let's not get crazy, even _I_ don't make class names this big
  ReadInClassRevisionLength = big_endian_int(ReadInClassRevisionLength);
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


  // Ok, things look good, let's read in the parameters that mClassName owns
  inStream.read((char *) &_atlasSizeX, sizeof(_atlasSizeX));
  _atlasSizeX = big_endian_int(_atlasSizeX);

  inStream.read((char *) &_atlasSizeY, sizeof(_atlasSizeY));
  _atlasSizeY =  big_endian_int(_atlasSizeY);

  inStream.read((char *) &_atlasSizeZ, sizeof(_atlasSizeZ));
   _atlasSizeZ =  big_endian_int(_atlasSizeZ);

  inStream.read((char *) &_patientSizeX, sizeof(_patientSizeX));
  _patientSizeX =  big_endian_int(_patientSizeX);

  inStream.read((char *) &_patientSizeY, sizeof(_patientSizeY));
  _patientSizeY =  big_endian_int(_patientSizeY);

  inStream.read((char *) &_patientSizeZ, sizeof(_patientSizeZ));
  _patientSizeZ =  big_endian_int(_patientSizeZ);

  inStream.read((char *) &_atlasResX, sizeof(_atlasResX));
  _atlasResX = big_endian_double(_atlasResX);
  
  inStream.read((char *) &_atlasResY, sizeof(_atlasResY));
  _atlasResY = big_endian_double(_atlasResY);

  inStream.read((char *) &_atlasResZ, sizeof(_atlasResZ));
  _atlasResZ = big_endian_double(_atlasResZ);

  inStream.read((char *) &_patientResX, sizeof(_patientResX));
  _patientResX = big_endian_double(_patientResX);

  inStream.read((char *) &_patientResY, sizeof(_patientResY));
  _patientResY = big_endian_double(_patientResY);

  inStream.read((char *) &_patientResZ, sizeof(_patientResZ));
  _patientResZ = big_endian_double(_patientResZ);
  
  // Read in parent's parameters
  if ( Transformation::loadClassParameters( inStream ) )
  {
    return 5;
  };

  _haveHField = 0;

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
FieldTransformation::printClassParameters ( ostream & outStream)
{
  Transformation::printClassParameters ( outStream );
  outStream << mClassName << "::mClassVersion: " << mClassVersion << endl;
  outStream << mClassName << "::mClassRevision: " << mClassRevision << endl;
  
  outStream << mClassName << "::_atlasSizeX: " << _atlasSizeX << endl;
  outStream << mClassName << "::_atlasSizeY: " << _atlasSizeY << endl;
  outStream << mClassName << "::_atlasSizeZ: " << _atlasSizeZ << endl;
  outStream << mClassName << "::_patientSizeX: " << _patientSizeX << endl;
  outStream << mClassName << "::_patientSizeY: " << _patientSizeY << endl;
  outStream << mClassName << "::_patientSizeZ: " << _patientSizeZ << endl;
  outStream << mClassName << "::_atlasResX: " << _atlasResX << endl;
  outStream << mClassName << "::_atlasResY: " << _atlasResY << endl;
  outStream << mClassName << "::_atlasResZ: " << _atlasResZ << endl;
  outStream << mClassName << "::_patientResX: " << _patientResX << endl;
  outStream << mClassName << "::_patientResY: " << _patientResY << endl;
  outStream << mClassName << "::_patientResZ: " << _patientResZ << endl;

  return outStream;

}


////////////////////////////////////////////////////////////////////////////////
//
// Public Methods
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//
// Function: apply() for volumes.
//
// Purpose: Apply H-field
//
// Inputs: To be determined
//
// Outputs: To be determined
//
////////////////////////////////////////////////////////////////////////////////
int FieldTransformation::apply(const Array3D<unsigned char> &inVolume, 
		Array3D<unsigned char> *outVolume, 
		double Rx, double Ry, double Rz, 
		double sOx, double sOy, double sOz, 
		const int interpolateFlag) {


  if (getVerboseMode())
  {
    cout << "apply() in FieldTransformation" << endl;
    system("date");
  }

  // Check if its ok to apply transformation.
  if (!isOkToApply())
    throwError("FieldTransformation::apply: cannot apply transformation"); 
      
  // Size of output volume.
  int nx = outVolume->getXsize(); 
  int ny = outVolume->getYsize(); 
  int nz = outVolume->getZsize(); 

  // Allocate space for transformation vector field
  if (!_haveHField) {
    _hField.setDim(nz, ny, nx*3);
    getHField(&_hField, Rx, Ry, Rz, sOx, sOy, sOz); 
    _haveHField = 1;
  }

  /********************************************
   * Deform input volume and store in outVolume.
   ********************************************/

  double tt; 
  int x, y, z; 

  const float *hFieldPtr = _hField.data();
  unsigned char *dfrm = outVolume->data(); 
  unsigned char bkgrnd = inVolume[0][0][0]; 

  int inx = inVolume.getXsize(); 
  int iny = inVolume.getYsize(); 
  int inz = inVolume.getZsize(); 

  float hFieldX;
  float hFieldY;
  float hFieldZ;

  if (interpolateFlag) {
    cerr << "slice     of " << nz-1 << "\010\010\010\010\010\010\010";
    for (z = 0; z < nz; z++) {
      cerr << "\010\010\010" ;
      cerr.width(3);
      cerr << z;
      
      for (y = 0; y < ny; y++) {
	for (x = 0; x < nx; x++) {
	  hFieldX = *hFieldPtr++;
	  hFieldY = *hFieldPtr++;
	  hFieldZ = *hFieldPtr++;
	  tt = TransformationUtils::trilinear(inVolume, 
				       hFieldZ, hFieldY, hFieldX, bkgrnd); 
	  *dfrm++ = (unsigned char) rint(tt); 
	  
	}
      }
    }
    cerr << endl;
  } else {
    cerr << "slice     of " << nz-1 << "\010\010\010\010\010\010\010";
    for (z = 0; z < nz; z++) {
      cerr << "\010\010\010" ;
      cerr.width(3);
      cerr << z;
      
      for (y = 0; y < ny; y++) {
	for (x = 0; x < nx; x++) {
	  hFieldX = *hFieldPtr++;
	  hFieldY = *hFieldPtr++;
	  hFieldZ = *hFieldPtr++;
	  int newX = int( hFieldX + 0.5);
	  int newY = int( hFieldY + 0.5);
	  int newZ = int( hFieldZ + 0.5);
	  if ((newX >= 0) && (newX < inx) &&
	      (newY >= 0) && (newY < iny) &&
	      (newZ >= 0) && (newZ < inz))
	    *dfrm++ = inVolume[newZ][newY][newX];
	  else
	    *dfrm++ = 0;
	}
      }
    }
    cerr << endl;
  }

  system("date");
  return 1;

}



#ifdef wait_for_mm
////////////////////////////////////////////////////////////////////////////////
//
// Function: apply() for volumes.
//
// Purpose: Apply H-field
//
// Inputs: offsets are in mm
//
// Outputs: To be determined
//
////////////////////////////////////////////////////////////////////////////////
int FieldTransformation::apply(const Array3D<unsigned char> &inVolume,
                double inResX, double inResY, double inResZ,
                double inOffX, double inOffY, double inOffZ,
                Array3D<unsigned char> *outVolume,
                double outResX, double outResY, double outResZ,
                double outOffX, double outOffY, double outOffZ,
		const int interpolatFlag)
{


  if (getVerboseMode())
  {
    cout << "apply() in FieldTransformation" << endl;
    system("date");
  }

  // Check if its ok to apply transformation.
  if (!isOkToApply())
    throwError("FieldTransformation::apply: cannot apply transformation"); 
      
  // Size of output volume.
  int nx = outVolume->getXsize(); 
  int ny = outVolume->getYsize(); 
  int nz = outVolume->getZsize(); 

  // Allocate space for transformation vector field
  _hField.setDim(nz, ny, nx*3);
  getHField(&_hField, outResX, outResY, outResZ,
		outOffX,outOffY,outOffZ):


  /********************************************
   * Deform input volume and store in outVolume.
   ********************************************/

  double tt; 
  int x, y, z; 

  const float *hFieldPtr = _hField.data();
  unsigned char *dfrm = outVolume->data(); 
  unsigned char bkgrnd = inVolume[0][0][0]; 

  int inx = inVolume.getXsize(); 
  int iny = inVolume.getYsize(); 
  int inz = inVolume.getZsize(); 

  float hFieldX;
  float hFieldY;
  float hFieldZ;

  // kwd
  float dwork = 1.0/((float)nz);
  float work_perc = 0.0;

  if (interpolateFlag) {
    for (z = 0; z < nz; z++, work_perc += dwork) {
      if(z % 5 == 0) ShowWorking(work_perc);
      for (y = 0; y < ny; y++) {
	for (x = 0; x < nx; x++) {
	  hFieldX = (*hFieldPtr++ + inOffX) / inResX;
	  hFieldY = (*hFieldPtr++ + inOffY) / inResY;
	  hFieldZ = (*hFieldPtr++ + inOffZ) / inResZ;
	  tt = TransformationUtils::trilinear(inVolume,
			 hFieldZ, hFieldY, hFieldX, bkgrnd); 
	  *dfrm++ = (unsigned char) rint(tt); 
	  
	}
      }
    }
  } else {
    for (z = 0; z < nz; z++,work_perc += dwork) {
      if(z % 5 == 0) ShowWorking(work_perc);
      for (y = 0; y < ny; y++) {
	for (x = 0; x < nx; x++) {
	  hFieldX = (*hFieldPtr++ + inOffX) / inResX;
	  hFieldY = (*hFieldPtr++ + inOffY) / inResY;
	  hFieldZ = (*hFieldPtr++ + inOffZ) / inResZ;
	  int newX = int( hFieldX + 0.5);
	  int newY = int( hFieldY + 0.5);
	  int newZ = int( hFieldZ + 0.5);
	  if ((newX >= 0) && (newX < inx) &&
	      (newY >= 0) && (newY < iny) &&
	      (newZ >= 0) && (newZ < inz))
	    *dfrm++ = inVolume[newZ][newY][newX];
	  else
	    *dfrm++ = 0;
	}
      }
    }
  }

  system("date");
  return 1;

}
#endif



////////////////////////////////////////////////////////////////////////////////
//
// Function: new apply() for volume + segmentations.
//
// Purpose: Apply H-field
//
// Inputs: To be determined
//
// Outputs: To be determined
//
////////////////////////////////////////////////////////////////////////////////
int FieldTransformation::apply(const Array3D<unsigned char> &inputVolume, 
			       double inputResX,     // in mm/voxel
			       double inputResY,
			       double inputResZ,
			       double inputOffsetX,  // offset relative to atlas
			       double inputOffsetY,  // in mm
			       double inputOffsetZ,
			       Array3D<unsigned char> *outputVolume,
			       double outputResX,    // in mm/voxel
			       double outputResY,
			       double outputResZ,
			       double outputOffsetX, // offset relative to patient
			       double outputOffsetY, // in mm
			       double outputOffsetZ,
			       const int interpolateFlag) { 

  if (getVerboseMode())
  {
    cout << "apply(New) in FieldTransformation" << endl;
    system("date");
  }

  // Check if its ok to apply transformation.
  if (!isOkToApply())
    throwError("FieldTransformation::apply: cannot apply transformation"); 
      
  // Size of output volume.
  int nx = outputVolume->getXsize(); 
  int ny = outputVolume->getYsize(); 
  int nz = outputVolume->getZsize(); 

  // Allocate space for transformation vector field
  _hField.setDim(nz, ny, nx*3);
  getHField(&_hField, outputResX, outputResY, outputResZ,
	    outputOffsetX, outputOffsetY, outputOffsetZ); 

  /********************************************
   * Deform input volume and store in outputVolume.
   ********************************************/

  double tt; 
  int x, y, z; 

  const float *hFieldPtr = _hField.data();

  unsigned char *dfrm = outputVolume->data(); 
  unsigned char bkgrnd = inputVolume[0][0][0]; 

  int inx = inputVolume.getXsize(); 
  int iny = inputVolume.getYsize(); 
  int inz = inputVolume.getZsize(); 

  float hFieldX;
  float hFieldY;
  float hFieldZ;

  double inputOffsetInVoxelsX = inputOffsetX / inputResX;
  double inputOffsetInVoxelsY = inputOffsetY / inputResY;
  double inputOffsetInVoxelsZ = inputOffsetZ / inputResZ;

  // kwd
  float dwork = 1.0/((float)nz);
  float work_perc = 0.0;

  if (interpolateFlag) {
    for (z = 0; z < nz; z++,work_perc += dwork) {
      if(z % 5 == 0) ShowWorking(work_perc);
      for (y = 0; y < ny; y++) {
	for (x = 0; x < nx; x++) {
	  hFieldX = *hFieldPtr++ * _atlasResX / inputResX - inputOffsetInVoxelsX; 
	  hFieldY = *hFieldPtr++ * _atlasResY / inputResY - inputOffsetInVoxelsY; 
	  hFieldZ = *hFieldPtr++ * _atlasResZ / inputResZ - inputOffsetInVoxelsZ; 
	  tt = TransformationUtils::trilinear(inputVolume,
			 hFieldZ, hFieldY, hFieldX, bkgrnd); 
	  *dfrm++ = (unsigned char) rint(tt); 
	}
      }
    }
  } else {
    for (z = 0; z < nz; z++,work_perc += dwork) {
      if(z % 5 == 0) ShowWorking(work_perc);
      for (y = 0; y < ny; y++) {
	for (x = 0; x < nx; x++) {
	  hFieldX = *hFieldPtr++ * _atlasResX / inputResX - inputOffsetInVoxelsX;
	  hFieldY = *hFieldPtr++ * _atlasResY / inputResY - inputOffsetInVoxelsY;
	  hFieldZ = *hFieldPtr++ * _atlasResZ / inputResZ - inputOffsetInVoxelsZ;
	  int newX = (int)floor(hFieldX);
	  int newY = (int)floor(hFieldY);
	  int newZ = (int)floor(hFieldZ);
	  if ((newX >= 0) && (newX < inx) &&
	      (newY >= 0) && (newY < iny) &&
	      (newZ >= 0) && (newZ < inz))
	    *dfrm++ = inputVolume[newZ][newY][newX];
	  else
	    *dfrm++ = 0;
	}
      }
    }
  }

  system("date");
  return 1;
}


////////////////////////////////////////////////////////////////////////////////
//
// Function: apply()
//
// Purpose: Apply h-field to a Surface
//
// Inputs: Surface to be deformed (in atlas coordinates)
//
// Outputs: deformed Surface (patient coordinates)
//
////////////////////////////////////////////////////////////////////////////////
int
FieldTransformation::apply(const Surface &inSurface,
 			   Surface *outSurface)
{

  if (getVerboseMode())
    cout << "apply() Surface in FieldTransformation" << endl;

  // Check if it is ok to apply transformation
  if (!isOkToApply())
    throwError("FieldTransformation::apply(Surface): cannot apply transformation");

  // allocate space and get H-field for transformation
  if (!_haveHField)
  {
    _hField.setDim(_patientSizeZ, _patientSizeY, _patientSizeX*3);
    getHField(&_hField);
    _haveHField = 1;
  }
    
  // copy input surface to output surface
  *outSurface = inSurface;

  // generate deformed surface
  if (!outSurface->invFieldTrans(_hField))
    throwError("FieldTransformation::apply(Surface): invFieldTrans failed");

  return 1;
}

int FieldTransformation::apply(const Surface &inSurface,
			       Surface *outSurface,
			       double Rx, double Ry, double Rz,
			       double Ox, double Oy, double Oz)
{
  if (getVerboseMode())
    cout << "apply() Surface with resolution in FieldTransformation" << endl;

  // make sure it is ok to apply
  if (!isOkToApply())
    throwError("FieldTransformation::apply(Surface/res): cannot apply Fluid transformation");

  // allocate space and get H-field for transformation
  if (!_haveHField)
  {
    _hField.setDim(_patientSizeZ, _patientSizeY, _patientSizeX*3);
    getHField(&_hField, Rx, Ry, Rz, Ox, Oy, Oz);
    _haveHField = 1;
  }
    
  // copy input surface to output surface
  *outSurface = inSurface;

  // generate deformed surface
  if (!outSurface->invFieldTrans(_hField))
    throwError("FieldTransformation::apply(Surface/res): invFieldTrans failed");

  return 1;
}


////////////////////////////////////////////////////////////////////////////////
//
// Function: apply()
//
// Purpose: Apply h-field to a Point
//
// Inputs: Point to be transformed (Point in the atlas)
//
// Outputs: transformed Point (location in the patient)
//
////////////////////////////////////////////////////////////////////////////////
int
FieldTransformation::apply(const Point &inPoint,
			   Point *outPoint)
{

  if (getVerboseMode())
    cout << "apply(Point) in FieldTransformation" << endl;

  // Check if it is ok to apply transformation
  if (!isOkToApply())
    throwError("FieldTransformation::apply(Point): cannot apply transformation");

  // allocate space and get H-field for transformation
  //  if (!_haveHField) {
    _hField.setDim(_patientSizeZ, _patientSizeY, _patientSizeX*3);
    getHField(&_hField);
    _haveHField = 1;
    //  }

  // copy input point to output point
  *outPoint = inPoint;

  // generate transformed point
  if (!outPoint->invFieldTrans(_hField))
    throwError("FieldTransformation::apply(Point): invFieldTrans failed");

  return 1;
}



////////////////////////////////////////////////////////////////////////////////
//
// Function: invApply()
//
// Purpose: Apply inverse of H-field to a volume, ie. take a segmentation in
//   the patient volume and transform it into the atlas
//
// Inputs: inVolume      Volume (segmentation) to compute the inverse of
//         outVolume     Output Volume (at requested size)
//         Rx, Ry, Rz    Resolution (mm/voxel) of the hField
//         sOx, sOy, sOz Offset of output cube in Study (voxels)
//         Rox, Roy, Roz Desired resolution of output
//         Oox, Ooy, Ooz Offset to output cube in template (voxels)
//
// Outputs:
//         Returns 0 on success and the deformed volume in outVolume
//
////////////////////////////////////////////////////////////////////////////////
int
FieldTransformation::invApply(const Array3D<unsigned char> &inVolume,
			      Array3D<unsigned char> *outVolume,
			      double  Rx, double  Ry, double  Rz, 
			      double sOx, double sOy, double sOz)
// currently the output size parameters (Rox ... Oox ...) are not used
//			      double Rox, double Roy, double Roz,
//			      double Oox, double Ooy, double Ooz)
{
  if (getVerboseMode())
  {
    cout << "invApply() in FieldTransformation" << endl;
    system("date");
  }

  // Check if its ok to apply transformation.
  if (!isOkToApply())
    throwError("FieldTransformation::invApply: cannot apply transformation"); 
      
  // Size of output volume.
  int nx = outVolume->getXsize(); 
  int ny = outVolume->getYsize(); 
  int nz = outVolume->getZsize(); 

  // Allocate space for transformation vector field
  if (!_haveHField) {
    _hField.setDim(_patientSizeZ, _patientSizeY, _patientSizeX*3);
    getHField(&_hField, Rx, Ry, Rz, sOx, sOy, sOz);
    _haveHField = 1;
  }

  /********************************************
   * Deform input volume and store in outVolume.
   ********************************************/

  double tt; 
  int i, x, y, z,
    px, py, pz; 

  unsigned char bkgrnd = inVolume[0][0][0]; 

  int inx = inVolume.getXsize(); 
  int iny = inVolume.getYsize(); 
  int inz = inVolume.getZsize(); 

  // accumulate distributed weights and counts
  Array3D<double> accum(nz, ny, nx);
  Array3D<double> count(nz, ny, nx);
  accum = 0.0;
  count = 0.0;

  int indx, indy, indz;		// indices (integer positions)
  double rx, ry, rz;		// residuals (fractional positions)
  // trilinear weights
  double dx00, dx01, dx10, dx11;
  double d000, d001, d010, d011, d100, d101, d110, d111;
  double chk, val;

//
// reverse interpolation method
//
// for each point in the image compute the trilinear weights
// and distribute the intensity as if it was generated by the
// trilinear neighbors.  Keep track of the number of contributors
// to each lattice point and normalize at the end.
//

  cerr << "slice     of " << inz-1 << "\010\010\010\010\010\010\010";
  for (z = 0; z < inz; z++)
  {
    cerr << "\010\010\010" ;
    cerr.width(3);
    cerr << z;

    for (y = 0; y < iny; y++)
    {
      for (x = 0; x < inx; x++)
      {
	val = inVolume[z][y][x];
	indx = int(_hField[z][y][3*x]);
	indy = int(_hField[z][y][3*x+1]);
	indz = int(_hField[z][y][3*x+2]);
	rx = _hField[z][y][3*x] - indx;
	ry = _hField[z][y][3*x+1] - indy;
	rz = _hField[z][y][3*x+2] - indz;

	dx00 = (1.0-ry) * (1.0-rz);
	dx01 = (1.0-ry) *    rz;
	dx10 =    ry  * (1.0-rz);
	dx11 =    ry  *    rz;

	d000 = (1.0-rx) * dx00;
	d001 = (1.0-rx) * dx01;
	d010 = (1.0-rx) * dx10;
	d011 = (1.0-rx) * dx11;
	d100 =    rx  * dx00;
	d101 =    rx  * dx01;
	d110 =    rx  * dx10;
	d111 =    rx  * dx11;


	if (indx >= 0 && indx < nx-1 &&
	    indy >= 0 && indy < ny-1 &&
	    indz >= 0 && indz < nz-1)
	{
	  // distribute intensity to 8 neighbors and
	  // increment count for later normalization
	  accum[indz  ][indy  ][indx  ] += d000 * val;
	  count[indz  ][indy  ][indx  ] += d000;
	  accum[indz+1][indy  ][indx  ] += d001 * val;
	  count[indz+1][indy  ][indx  ] += d001;
	  accum[indz  ][indy+1][indx  ] += d010 * val;
	  count[indz  ][indy+1][indx  ] += d010;
	  accum[indz+1][indy+1][indx  ] += d011 * val;
	  count[indz+1][indy+1][indx  ] += d011;
	  accum[indz  ][indy  ][indx+1] += d100 * val;
	  count[indz  ][indy  ][indx+1] += d100;
	  accum[indz+1][indy  ][indx+1] += d101 * val;
	  count[indz+1][indy  ][indx+1] += d101;
	  accum[indz  ][indy+1][indx+1] += d110 * val;
	  count[indz  ][indy+1][indx+1] += d110;
	  accum[indz+1][indy+1][indx+1] += d111 * val;
	  count[indz+1][indy+1][indx+1] += d111;

	}
	else			// handle boundary voxels
	{
	  // For simplicity, just distribute the intensity to
	  // all neighbors equally (as above), since we don't
	  // really know what value to use for the background
	  if (indx >= 0 && indx < nx)
	  {
	    if (indy >= 0 && indy < ny)
	    {
	      if (indz >= 0 && indz < nz)
	      {
		accum[indz  ][indy  ][indx  ] += d000 * val;
		count[indz  ][indy  ][indx  ] += d000;
	      }
	      if (indz+1 >= 0 && indz+1 < nz)
	      {
		accum[indz+1][indy  ][indx  ] += d001 * val;
		count[indz+1][indy  ][indx  ] += d001;
	      }
	    }
	    if (indy+1 >= 0 && indy+1 < ny)
	    {
	      if (indz >= 0 && indz < nz)
	      {
		accum[indz  ][indy+1][indx  ] += d010 * val;
		count[indz  ][indy+1][indx  ] += d010;
	      }
	      if (indz+1 >= 0 && indz+1 < nz)
	      {
		accum[indz+1][indy+1][indx  ] += d011 * val;
		count[indz+1][indy+1][indx  ] += d011;
	      }
	    }
	  }
	  if (indx+1 >= 0 && indx+1 < nx)
	  {
	    if (indy >= 0 && indy < ny)
	    {
	      if (indz >= 0 && indz < nz)
	      {
		accum[indz  ][indy  ][indx+1] += d100 * val;
		count[indz  ][indy  ][indx+1] += d100;
	      }
	      if (indz+1 >= 0 && indz+1 < nz)
	      {
		accum[indz+1][indy  ][indx+1] += d101 * val;
		count[indz+1][indy  ][indx+1] += d101;
	      }
	    }
	    if (indy+1 >= 0 && indy+1 < ny)
	    {
	      if (indz >= 0 && indz < nz)
	      {
		accum[indz  ][indy+1][indx+1] += d110 * val;
		count[indz  ][indy+1][indx+1] += d110;
	      }
	      if (indz+1 >= 0 && indz+1 < nz)
	      {
		accum[indz+1][indy+1][indx+1] += d111 * val;
		count[indz+1][indy+1][indx+1] += d111;
	      }
	    }
	  }
	}
      }
    }
  }
  cerr << endl;

  // normalize
  if (getVerboseMode())
    cout << "Normalizing..." << endl;

  int numElem = nz*ny*nx;
  double *pCount = count.data();
  double *pAccum = accum.data();
  unsigned char *dfrm = outVolume->data(); 

  for (i=0; i<numElem; i++, dfrm++, pCount++, pAccum++)
    if (*pCount != 0)
      *dfrm = (unsigned char) (*pAccum / *pCount + 0.5);
    else
      *dfrm = 0;

  if (getVerboseMode())
    system("date");

  return 1;
}


////////////////////////////////////////////////////////////////////////////////
//
// Function: invApply()
//
// Purpose: Apply inverse h-field to a Surface
//
// Inputs: Surface to be deformed (in patient space)
//
// Outputs: deformed Surface (atlas space)
//
////////////////////////////////////////////////////////////////////////////////
int FieldTransformation::invApply(const Surface &inSurface,
				  Surface *outSurface,
				  double Rx, double Ry, double Rz,
				  double Ox, double Oy, double Oz)
{
  if (getVerboseMode())
    cout << "invApply() Surface with resolution in FieldTransformation" << endl;

  // make sure it is ok to apply
  if (!isOkToApply())
    throwError("FieldTransformation::invApply(Surface/res): cannot apply transformation");

  // allocate space and get H-field for transformation
  if (!_haveHField)
  {
    _hField.setDim(_patientSizeZ, _patientSizeY, _patientSizeX*3);
    getHField(&_hField, Rx, Ry, Rz, Ox, Oy, Oz);
    _haveHField = 1;
  }
    
  // copy input surface to output surface
  *outSurface = inSurface;

  // generate deformed surface
  if (outSurface->applyFieldTrans(_hField))
    throwError("FieldTransformation::invApply(Surface/res): applyFieldTrans failed");

  return 1;
}


// this method assumes that 0 <= x < sizeX-1
// 0 <= y < sizeY-1  and  0 <= z < sizeZ-1
void FieldTransformation::trilinearHField(const Array3D<float> &hField,
						  float x, float y, float z,
						  float *Hx, 
						  float *Hy, 
						  float *Hz) {
#if 0
  const float *const *const *hFieldPtrPtrPtr = hField.address();
  int intX = (int) (x + 0.5);
  int intY = (int) (y + 0.5);
  int intZ = (int) (z + 0.5);
  
  *Hx = hFieldPtrPtrPtr[intZ][intY][intX*3];
  *Hy = hFieldPtrPtrPtr[intZ][intY][intX*3+1];
  *Hz = hFieldPtrPtrPtr[intZ][intY][intX*3+2];
#else

  int sizeX = hField.getXsize(); // equals 3 * length of x row since interlaced
  int sizeY = hField.getYsize();
  int sizeZ = hField.getZsize();
  int sizeXY = sizeX*sizeY;

  const float *_voxels = hField.data();

  int x0 = (int) x;
  int y0 = (int) y;
  int z0 = (int) z;
  int x1 = x0 + 1;
  int y1 = y0 + 1;
  int z1 = z0 + 1;
  float dx = x - x0;
  float dy = y - y0;
  float dz = z - z0;

  int x0InRange = ((x0 >= 0) && (x0 < sizeX/3));
  int y0InRange = ((y0 >= 0) && (y0 < sizeY));
  int z0InRange = ((z0 >= 0) && (z0 < sizeZ));
  int x1InRange = ((x1 >= 0) && (x1 < sizeX/3));
  int y1InRange = ((y1 >= 0) && (y1 < sizeY));
  int z1InRange = ((z1 >= 0) && (z1 < sizeZ));

  float hx000, hx001, hx010, hx011;
  float hx100, hx101, hx110, hx111;
  float hy000, hy001, hy010, hy011;
  float hy100, hy101, hy110, hy111;
  float hz000, hz001, hz010, hz011;
  float hz100, hz101, hz110, hz111;
  
  hx000 = hx001 = hx010 = hx011 = 0;
  hx100 = hx101 = hx110 = hx111 = 0;
  hy000 = hy001 = hy010 = hy011 = 0;
  hy100 = hy101 = hy110 = hy111 = 0;
  hz000 = hz001 = hz010 = hz011 = 0;
  hz100 = hz101 = hz110 = hz111 = 0;
  
  int baseIndex0 = z0*sizeXY + y0*sizeX + x0*3;
  int baseIndex1 = baseIndex0 + sizeXY;

  if (z0InRange && y0InRange && x0InRange &&
      z1InRange && y1InRange && x1InRange) {
    hx000 = _voxels[baseIndex0];
    hy000 = _voxels[baseIndex0 + 1];
    hz000 = _voxels[baseIndex0 + 2];
    hx001 = _voxels[baseIndex0 + 3];
    hy001 = _voxels[baseIndex0 + 4];
    hz001 = _voxels[baseIndex0 + 5];

    hx010 = _voxels[baseIndex0 + sizeX];
    hy010 = _voxels[baseIndex0 + sizeX + 1];
    hz010 = _voxels[baseIndex0 + sizeX + 2];
    hx011 = _voxels[baseIndex0 + sizeX + 3];
    hy011 = _voxels[baseIndex0 + sizeX + 4];
    hz011 = _voxels[baseIndex0 + sizeX + 5];

    hx100 = _voxels[baseIndex1];
    hy100 = _voxels[baseIndex1 + 1];
    hz100 = _voxels[baseIndex1 + 2];
    hx101 = _voxels[baseIndex1 + 3];
    hy101 = _voxels[baseIndex1 + 4];
    hz101 = _voxels[baseIndex1 + 5];

    hx110 = _voxels[baseIndex1 + sizeX];
    hy110 = _voxels[baseIndex1 + sizeX + 1];
    hz110 = _voxels[baseIndex1 + sizeX + 2];
    hx111 = _voxels[baseIndex1 + sizeX + 3];
    hy111 = _voxels[baseIndex1 + sizeX + 4];
    hz111 = _voxels[baseIndex1 + sizeX + 5];
  } else {
    if (z0InRange && y0InRange && x0InRange) {
      hx000 = _voxels[baseIndex0];
      hy000 = _voxels[baseIndex0 + 1];
      hz000 = _voxels[baseIndex0 + 2];
    }
    if (z0InRange && y0InRange && x1InRange) {
      hx001 = _voxels[baseIndex0 + 3];
      hy001 = _voxels[baseIndex0 + 4];
      hz001 = _voxels[baseIndex0 + 5];
    }
    if (z0InRange && y1InRange && x0InRange) {
      hx010 = _voxels[baseIndex0 + sizeX];
      hy010 = _voxels[baseIndex0 + sizeX + 1];
      hz010 = _voxels[baseIndex0 + sizeX + 2];
    }
    if (z0InRange && y1InRange && x1InRange) {
      hx011 = _voxels[baseIndex0 + sizeX + 3];
      hy011 = _voxels[baseIndex0 + sizeX + 4];
      hz011 = _voxels[baseIndex0 + sizeX + 5];
    }
    if (z1InRange && y0InRange && x0InRange) {
      hx100 = _voxels[baseIndex1];
      hy100 = _voxels[baseIndex1 + 1];
      hz100 = _voxels[baseIndex1 + 2];
    }
    if (z1InRange && y0InRange && x1InRange) {
      hx101 = _voxels[baseIndex1 + 3];
      hy101 = _voxels[baseIndex1 + 4];
      hz101 = _voxels[baseIndex1 + 5];
    }
    if (z1InRange && y1InRange && x0InRange) {
      hx110 = _voxels[baseIndex1 + sizeX];
      hy110 = _voxels[baseIndex1 + sizeX + 1];
      hz110 = _voxels[baseIndex1 + sizeX + 2];
    }
    if (z1InRange && y1InRange && x1InRange) {
      hx111 = _voxels[baseIndex1 + sizeX + 3];
      hy111 = _voxels[baseIndex1 + sizeX + 4];
      hz111 = _voxels[baseIndex1 + sizeX + 5];
    }
  }

  // interpolate x component of H field
  float d00x = hx000 + dx*(hx001 - hx000);
  float d01x = hx010 + dx*(hx011 - hx010);
  float d10x = hx100 + dx*(hx101 - hx100);
  float d11x = hx110 + dx*(hx111 - hx110);
  
  float d0yx = d00x + dy*(d01x - d00x);
  float d1yx = d10x + dy*(d11x - d10x);
  
  *Hx = (d0yx + dz*(d1yx - d0yx));


  // interpolate y component of H field
  d00x = hy000 + dx*(hy001 - hy000);
  d01x = hy010 + dx*(hy011 - hy010);
  d10x = hy100 + dx*(hy101 - hy100);
  d11x = hy110 + dx*(hy111 - hy110);
  
  d0yx = d00x + dy*(d01x - d00x);
  d1yx = d10x + dy*(d11x - d10x);
  
  *Hy = (d0yx + dz*(d1yx - d0yx));


  // interpolate z component of H field
  d00x = hz000 + dx*(hz001 - hz000);
  d01x = hz010 + dx*(hz011 - hz010);
  d10x = hz100 + dx*(hz101 - hz100);
  d11x = hz110 + dx*(hz111 - hz110);
  
  d0yx = d00x + dy*(d01x - d00x);
  d1yx = d10x + dy*(d11x - d10x);
  
  *Hz = (d0yx + dz*(d1yx - d0yx));


#endif
}


////////////////////////////////////////////////////////////////////////////////
//
// Function: invApply()
//
// Purpose: map a point in patient voxel coordinates back to its corresponding
//          atlas voxel coordinates
//
// Inputs: point to be mapped in patient space
//
// Outputs: corresponding point in atlas space
//
////////////////////////////////////////////////////////////////////////////////
int FieldTransformation::invApply(const Pnt &patientVoxelPoint,
				  Pnt &atlasVoxelPoint) {
  if (getVerboseMode())
  {
    cout << "FieldTransformation::invApply(Pnt, Pnt)" << endl;
    system("date");
  }

  // Check if its ok to apply transformation.
  if (!isOkToApply())
    throwError("FieldTransformation::invApply(Pnt, Pnt): cannot apply inverse transformation"); 
      
  double patientX = patientVoxelPoint.x();
  double patientY = patientVoxelPoint.y();
  double patientZ = patientVoxelPoint.z();

  // Allocate space for transformation vector field
  double offsetX = patientX - 1;
  double offsetY = patientY - 1;
  double offsetZ = patientZ - 1;

  double offsetXmm = offsetX * _patientResX;
  double offsetYmm = offsetY * _patientResY;
  double offsetZmm = offsetZ * _patientResZ;

  Array3D<float> hField(5, 5, 5*3);  // small sub cube around patientVoxelPoint
  getHField(&hField, _patientResX, _patientResY, _patientResZ,
	    offsetXmm, offsetYmm, offsetZmm);
  
  trilinearHField(hField,
		  1, 1, 1,
		  &atlasVoxelPoint.x(),
		  &atlasVoxelPoint.y(),
		  &atlasVoxelPoint.z());

#if 0
  // Size of output volume.
  int nx = _patientSizeX;
  int ny = _patientSizeY;
  int nz = _patientSizeZ;

  if (!_haveHField) {
    _hField.setDim(nz, ny, nx*3);
    getHField(&_hField, _patientResX, _patientResY, _patientResZ,
	      0, 0, 0); 
    _haveHField = 1;
  }

  trilinearHField(_hField,
		  patientVoxelPoint.x(),
		  patientVoxelPoint.y(),
		  patientVoxelPoint.z(),
		  &atlasVoxelPoint.x(),
		  &atlasVoxelPoint.y(),
		  &atlasVoxelPoint.z());
#endif
  

  
  system("date");
  return 1;
}




///////////////////////////////////////////////////////////////////////////
//
// Function: gimmeHField()
//
// Purpose: get the array containing the h-field
//
// Inputs: none
//
// Outputs: h-field
//
///////////////////////////////////////////////////////////////////////////
Array3D<float>
FieldTransformation::gimmeHField()
{ 
  // if h-field has not yet been computed, then compute it
  if (!_haveHField)
  {
    _hField.setDim(_patientSizeZ, _patientSizeY, _patientSizeX*3);
    getHField(&_hField);
    _haveHField = 1;
  }

  // return h-field
  return _hField;
}

int
FieldTransformation::setHField( Array3D<float> & hField)
{
  
  if ( ! isOkToApply())
  {
    cerr << "We do not have a displacment field yet. So we can not set another." << endl;
    return 1;
  }

  // Get h field from this object's sub-class
  if (!_haveHField) {
    _hField.setDim(_patientSizeZ, _patientSizeY, _patientSizeX * 3);
    getHField( &_hField );
    _haveHField = 1;
  }

  int uprc = TransformationUtils::upsampleInterlacedField( _hField,
                                                           hField);
  if (0 != uprc)
  {
    cerr << "TransformationUtils::upsampleInterlacedField failed with return code " << uprc << endl;
    return (10 + uprc);
  }

  return 0;
}
