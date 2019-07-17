///////////////////////////////////////////////////////////////////////////
//
// IntellX, L.L.C., Proprietary
// The contents of this file comprise confidential, proprietary and/or
// trade secret information which is the property of IntellX LLC.
// Dissemination, distribution or other disclosure of this information is
// strictly prohibited without the express permission of IntellX LLC.
// Copyright (c) 1997 IntellX, L.L.C.
// All rights reserved
//
// File:  FluidFluidTransformation.C
//
// Author:  Kevin E. Mark
//
// Purpose:  This class computes a fluid landmark transformation and a composed
//           fluid transformation on a subvolume.
//
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: FluidFluidTransformation.C,v 1.1 2004/11/15 04:44:07 joeh Exp $
//
// Revision History
//
// $Log: FluidFluidTransformation.C,v $
// Revision 1.1  2004/11/15 04:44:07  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:29:29  kem
// Initial revision
//
// Revision 1.2  1999/07/09 17:51:39  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.1  1999/01/08 17:29:49  RAZORDB
// Initial revision
//
//
///////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//
// Include Files
//
////////////////////////////////////////////////////////////////////////////////

#include <string.h>
#include <TDL/FluidFluidTransformation.h>
#include <TDL/UserDefinedLinearTransformation.h>
#include <ADL/GreyVolume.h>

const char * const FluidFluidTransformation::mClassRevision = "$Id: FluidFluidTransformation.C,v 1.1 2004/11/15 04:44:07 joeh Exp $";
const char * const FluidFluidTransformation::mClassName = "FluidFluidTransformation";
const int FluidFluidTransformation::mClassVersion = 1;


//////////////////////////////
// constructor and destructor
//////////////////////////////

FluidFluidTransformation::FluidFluidTransformation() {
  cout << "calling constructor" << endl;
}

FluidFluidTransformation::~FluidFluidTransformation() {
}


///////////////////
// Public Methods
///////////////////

////////////////////////////////////////////////////////////////////////////////
//
// Function: initialize
//
// Purpose: Initialize class so that calculate may run
//
// Inputs: To be determined
//
// Outputs: Sets mOkToCalculateFlag to 1 if successfull
//
////////////////////////////////////////////////////////////////////////////////
int FluidFluidTransformation::initialize(GreyVolume<unsigned char> const &atlas,
					 GreyVolume<unsigned char> const &patient,
					 Matrix<double> const &atlasLandmarks, 
					 Matrix<double> const &patientLandmarks,
					 double atlasSubOffsetX,
					 double atlasSubOffsetY,
					 double atlasSubOffsetZ,
					 int atlasSubSizeX,
					 int atlasSubSizeY,
					 int atlasSubSizeZ,
					 double atlasSubResX,
					 double atlasSubResY,
					 double atlasSubResZ,
					 const Array3D<unsigned char> *maskPtr) {

  cout << "Executing FluidFluidTransformation::initialize()" << endl;

  _atlasPtr = &atlas;
  _patientPtr = &patient;
  _maskPtr = maskPtr;
  _atlasLandmarksPtr = &atlasLandmarks;
  _patientLandmarksPtr = &patientLandmarks;

  _atlasSizeX = atlas.getXsize();
  _atlasSizeY = atlas.getYsize();
  _atlasSizeZ = atlas.getZsize();
  _atlasResX = atlas.getVoxelXScale();
  _atlasResY = atlas.getVoxelYScale();
  _atlasResZ = atlas.getVoxelZScale();

  _patientSizeX = patient.getXsize();
  _patientSizeY = patient.getYsize();
  _patientSizeZ = patient.getZsize();
  _patientResX = patient.getVoxelXScale();
  _patientResY = patient.getVoxelYScale();
  _patientResZ = patient.getVoxelZScale();

  _atlasSubOffsetX = atlasSubOffsetX;
  _atlasSubOffsetY = atlasSubOffsetY;
  _atlasSubOffsetZ = atlasSubOffsetZ;
  _atlasSubSizeX = atlasSubSizeX;
  _atlasSubSizeY = atlasSubSizeY;
  _atlasSubSizeZ = atlasSubSizeZ;
  _atlasSubResX = atlasSubResX;
  _atlasSubResY = atlasSubResY;
  _atlasSubResZ = atlasSubResZ;

  _patientSubSizeX = _atlasSubSizeX; // should be recomputed in calculate
  _patientSubSizeY = _atlasSubSizeY;
  _patientSubSizeZ = _atlasSubSizeZ;

  // ?????????????????????????????
  _patientSubResX = _atlasSubResX; // this is always the same as atlasSubRes
  _patientSubResY = _atlasSubResY; // so maybe we can remove it or _atlasSubRes
  _patientSubResZ = _atlasSubResZ;

  setOkToCalculate(1);

  return 0;
}



int FluidFluidTransformation::findGoodFactorization(int x) {
  for (int i=0; i<(x/5); i++)
    if (isDivisibleBy235(x+i))
      return (x+i);
  cout << "FluidFluidTransformation::findGoodFactorization(): "
       << "good factorization not found for " << x << endl;
  return x;
}

int FluidFluidTransformation::isDivisibleBy235(int x) {
  if (x==1)
    return 1;
  if (!(x%2)) // divisible by 2
    return isDivisibleBy235(x/2);
  else if (!(x%3)) // divisible by 3
    return isDivisibleBy235(x/3);
  else if (!(x%5)) // divisible by 5
    return isDivisibleBy235(x/5);
  else
    return 0;
}


////////////////////////////////////////////////////////////////////////////////
//
// Function: calculate
//
// Purpose: 
//
// Inputs: 
//
// Outputs: Turns mOkToApplyFlag to 1 if successfull
//          Returns 0 if sucessfull
//          Returns 1 if not ready to calculate
//
////////////////////////////////////////////////////////////////////////////////
int FluidFluidTransformation::calculate() {

  cout << "Executing FluidFluidTransformation::calculate()" << endl;

  if (!isOkToCalculate())
    return 1;

  // kwd

  ShowWorking("Initializing Fluid Landmark Transformation\n");

  // Initialize the fluid landmark transformation
  _fluidLandmarkTransform.initialize(*_atlasLandmarksPtr, *_patientLandmarksPtr, 
				     _atlasSizeX, _atlasSizeY, _atlasSizeZ,
				     _atlasResX, _atlasResY, _atlasResZ,
				     _patientSizeX, _patientSizeY, _patientSizeZ,
				     _patientResX, _patientResY, _patientResZ);

  // Calculate the fluid landmark transformation

  // kwd
  ShowWorking("Calculating Fluid Landmark Transformation\n");
  _fluidLandmarkTransform.assignWorkProc(this);
  _fluidLandmarkTransform.calculate(); 

#if 0
  // transform atlas (for testing only, can be removed)
  GreyVolume<unsigned char> transformedAtlas(_atlasSizeZ,
					     _atlasSizeY,
					     _atlasSizeX);
  _fluidLandmarkTransform.apply(*_atlasPtr,
				_atlasResX, _atlasResY, _atlasResZ,
				0, 0, 0,
				&transformedAtlas,
				_atlasResX, _atlasResY, _atlasResZ,
				0, 0, 0,
				1);
  transformedAtlas.saveAnalyze("transformedAtlas.img");
#endif

  // Compute patient subvolume offset
  float atlasSubOffsetInVoxelsX = _atlasSubOffsetX / _atlasResX;
  float atlasSubOffsetInVoxelsY = _atlasSubOffsetY / _atlasResY;
  float atlasSubOffsetInVoxelsZ = _atlasSubOffsetZ / _atlasResZ;
  cerr << "atlasSubOffset = " 
       << atlasSubOffsetInVoxelsX << " "
       << atlasSubOffsetInVoxelsY << " "
       << atlasSubOffsetInVoxelsZ << endl;
  cerr << "atlasSubSize = " 
       << _atlasSubSizeX << " "
       << _atlasSubSizeY << " "
       << _atlasSubSizeZ << endl;

#if 0
  Point patientSubOffsetInVoxels;
  _fluidLandmarkTransform.apply(atlasSubOffsetInVoxels,
				&patientSubOffsetInVoxels);
  int patientSubOffsetInVoxelsX = (int) patientSubOffsetInVoxels.x();
  int patientSubOffsetInVoxelsY = (int) patientSubOffsetInVoxels.y();
  int patientSubOffsetInVoxelsZ = (int) patientSubOffsetInVoxels.z();
#endif

  // 1. apply fluid landmark transformation to eight corners of atlas subvolume
  // 2. find min and max in three orthogonal directions to get bounding box
  // 3. center atlas subvolume-size box in bounding box

  // compute the corners of the atlas subvolume in mm
  Matrix<double> atlasCorners(8, 3); // atlas corners in mm

  atlasCorners[0][0] = _atlasSubOffsetX;
  atlasCorners[0][1] = _atlasSubOffsetY;
  atlasCorners[0][2] = _atlasSubOffsetZ;
  atlasCorners[1][0] = _atlasSubOffsetX + _atlasSubSizeX*_atlasSubResX;
  atlasCorners[1][1] = _atlasSubOffsetY;
  atlasCorners[1][2] = _atlasSubOffsetZ;
  atlasCorners[2][0] = _atlasSubOffsetX;
  atlasCorners[2][1] = _atlasSubOffsetY + _atlasSubSizeY*_atlasSubResY;
  atlasCorners[2][2] = _atlasSubOffsetZ;
  atlasCorners[3][0] = _atlasSubOffsetX + _atlasSubSizeX*_atlasSubResX; 
  atlasCorners[3][1] = _atlasSubOffsetY + _atlasSubSizeY*_atlasSubResY;
  atlasCorners[3][2] = _atlasSubOffsetZ;
  int k;
  for(k=4; k<8; k++) {
    atlasCorners[k][0] = atlasCorners[k-4][0];
    atlasCorners[k][1] = atlasCorners[k-4][1];
    atlasCorners[k][2] = atlasCorners[k-4][2] + _atlasSubSizeZ*_atlasSubResZ;
  }

  // transform atlas subvolume corners using the fluid landmark transform

  // find the min and max of each transformed coordinate
  // the min's determine the patient subvolume offset and the max's
  // determine the size of the patient subvolume
  int minX = _patientSizeX;  int maxX = 0;
  int minY = _patientSizeY;  int maxY = 0;
  int minZ = _patientSizeZ;  int maxZ = 0;

  Point atlasPoint;
  Point patientPoint;
  for (k=0; k<8; k++) {
    // convert point to voxels
    atlasPoint.set( atlasCorners[k][0] / _atlasResX,
		    atlasCorners[k][1] / _atlasResY,
		    atlasCorners[k][2] / _atlasResZ );

    // fluid landmark transform atlas point
    _fluidLandmarkTransform.apply(atlasPoint, &patientPoint);

    // update min and max
    int patientX = (int) patientPoint.x();
    int patientY = (int) patientPoint.y();
    int patientZ = (int) patientPoint.z();
    if (patientX < minX) minX = patientX;
    if (patientY < minY) minY = patientY;
    if (patientZ < minZ) minZ = patientZ;
    if (patientX > maxX) maxX = patientX;
    if (patientY > maxY) maxY = patientY;
    if (patientZ > maxZ) maxZ = patientZ;
  }

  if (minX < 0) minX = 0;
  if (minY < 0) minY = 0;
  if (minZ < 0) minZ = 0;
  if (maxX >= _patientSizeX) maxX = _patientSizeX-1;
  if (maxY >= _patientSizeY) maxY = _patientSizeY-1;
  if (maxZ >= _patientSizeZ) maxZ = _patientSizeZ-1;
  
  // compute size's such that they are a multiple of 8 to make the most of
  // the fft's
#if 0
  _patientSubSizeX = 8 * ceil((maxX-minX)*(_patientResX/_patientSubResX)/8.0);
  _patientSubSizeY = 8 * ceil((maxY-minY)*(_patientResY/_patientSubResY)/8.0);
  _patientSubSizeZ = 8 * ceil((maxZ-minZ)*(_patientResZ/_patientSubResZ)/8.0);
#else
  _patientSubSizeX = findGoodFactorization((maxX-minX)*
					   (_patientResX/_patientSubResX));
  _patientSubSizeY = findGoodFactorization((maxY-minY)*
					   (_patientResY/_patientSubResY));
  _patientSubSizeZ = findGoodFactorization((maxZ-minZ)*
					   (_patientResZ/_patientSubResZ));
#endif

  // make sure we start at an integral voxel
  int patientSubOffsetInVoxelsX = minX;
  int patientSubOffsetInVoxelsY = minY;
  int patientSubOffsetInVoxelsZ = minZ;

  // need to do this so that patientSubOffsetInVoxels is consistent
  // with the offset in mm
  _patientSubOffsetX = patientSubOffsetInVoxelsX*_patientResX;
  _patientSubOffsetY = patientSubOffsetInVoxelsY*_patientResY;
  _patientSubOffsetZ = patientSubOffsetInVoxelsZ*_patientResZ;

  cerr << "patientSubOffset = " 
       << patientSubOffsetInVoxelsX << " "
       << patientSubOffsetInVoxelsY << " "
       << patientSubOffsetInVoxelsZ << endl;
  cerr << "patientSubSize = "
       << _patientSubSizeX << " "
       << _patientSubSizeY << " "
       << _patientSubSizeZ << endl;

  // Apply fluid landmark transformation to subvolume
  GreyVolume<unsigned char> transformedAtlasSubvolume(_patientSubSizeZ,
						      _patientSubSizeY,
						      _patientSubSizeX);

   // kwd
   ShowWorking("Applying Fluid Landmark Transformation\n");

  _fluidLandmarkTransform.apply(*_atlasPtr,
				_atlasResX, _atlasResY, _atlasResZ,
				0, 0, 0,
				&transformedAtlasSubvolume,
				_atlasSubResX, _atlasSubResY, _atlasSubResZ, 
				_patientSubOffsetX, _patientSubOffsetY, 
				_patientSubOffsetZ,
				1);

  // Create patient subvolume
  GreyVolume<unsigned char> patientSubvolume(_patientSubSizeZ,
					     _patientSubSizeY,
					     _patientSubSizeX);

  extractSubvolume(*_patientPtr,
		   patientSubOffsetInVoxelsX,
		   patientSubOffsetInVoxelsY,
		   patientSubOffsetInVoxelsZ,
		   _patientResX/_atlasSubResX,
		   _patientResY/_atlasSubResY,
		   _patientResZ/_atlasSubResZ,
		   patientSubvolume);

  // Create mask subvolume
  GreyVolume<unsigned char> maskSubvolume;
  if (_maskPtr != NULL) {
    maskSubvolume.setDim(_patientSubSizeZ,
			 _patientSubSizeY,
			 _patientSubSizeX);
    extractSubvolume(*_maskPtr,
		     patientSubOffsetInVoxelsX,
		     patientSubOffsetInVoxelsY,
		     patientSubOffsetInVoxelsZ,
		     _patientResX/_atlasSubResX,
		     _patientResY/_atlasSubResY,
		     _patientResZ/_atlasSubResZ,
		     maskSubvolume);
  }

  // Save the subvolumes for verification purposes (may be removed in future)
  cerr << "Saving subvolumes...";

  char *transformedAtlasFilename = "TASV.img";
  cerr << transformedAtlasFilename << "...";
  transformedAtlasSubvolume.saveAnalyze(transformedAtlasFilename);

  char *patientFilename = "PSV.img";
  cerr << patientFilename << "...";
  patientSubvolume.saveAnalyze(patientFilename);

  if (_maskPtr != NULL) {
    char *maskFilename = "MSV.img";
    cerr << maskFilename << "...";
    maskSubvolume.saveAnalyze(maskFilename);
    cerr << "done" << endl;
  } else {
    cerr << "done" << endl;
  }

  ShowWorking("Initializing Fluid Transformation\n");

  // Initialize the fluid transformation with the subvolumes
  if (_maskPtr == NULL) {
    _fluidTransform.initialize(&transformedAtlasSubvolume, 
			       &patientSubvolume);
  } else {
    _fluidTransform.initialize(&transformedAtlasSubvolume, 
			       &patientSubvolume,
			       &maskSubvolume);
  }

  ShowWorking("Calculating Fluid Transformation\n");
  _fluidTransform.assignWorkProc(this);

  // Calculate fluid transformation on the subvolumes
  _fluidTransform.calculate(); 

  // Say we have a transformation and we may now apply it

  ShowWorking("DONE\n");
  ShowWorking((char*)NULL);

  setOkToApply(1);

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
int FluidFluidTransformation::save( const char * filename ) {

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

  // cout << "save() in FluidFluidTransformation" << endl;
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
int FluidFluidTransformation::load( const char * filename ) {
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

  // cout << "load() in FluidFluidTransformation" << endl;
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
void FluidFluidTransformation::print() {
  cout << "print() in FluidFluidTransformation" << endl;
  printClassParameters(cout) << endl;  // may crash if << endl is removed?!
}


int FluidFluidTransformation::setNumberIterations( int numOuter, int numInner ) {
  return _fluidTransform.setNumberIterations(numOuter, numInner);
}

int FluidFluidTransformation::setLaplacianWeight( double weight ) { 
   return _fluidTransform.setLaplacianWeight(weight); 
}

int FluidFluidTransformation::setGradWeight( double weight ) { 
   return _fluidTransform.setGradWeight(weight);
}

int FluidFluidTransformation::setAdditiveWeight( double weight ) { 
   return _fluidTransform.setAdditiveWeight(weight);
}

int FluidFluidTransformation::setMaxPerturbation( double pmax ) { 
   return _fluidTransform.setMaxPerturbation(pmax);
}

int FluidFluidTransformation::setMagicStoppingNumber( double x ) { 
   return _fluidTransform.setMagicStoppingNumber(x);
}

int FluidFluidTransformation::setCalculateTotalJacobian( int flag ) { 
   return _fluidTransform.setCalculateTotalJacobian(flag); 
}

int FluidFluidTransformation::setInitNumBasis( int initNumBasis ) {
  return _fluidTransform.setNumBasis(initNumBasis);
}

int FluidFluidTransformation::setInitialTransformation(const char *filename) {  
   return _fluidTransform.setInitialTransformation(filename); 
}



//////////////////////
// Protected Methods
//////////////////////

////////////////////////////////////////////////////////////////////////////////
//
// Function: saveClassParameters
//
// Purpose: Saves parameters that are owned by this class to the supplied ofstream.
//          It is envisioned that it's own parameters are saves and then
//          those of it's parent.  Since at this time this is the 
//          the final derived class, save() calls this function.
//
// Inputs: ofstream that we want to save to 
//
// Outputs: 0 upon success
//          1 if parent's save returned non-zero
//
////////////////////////////////////////////////////////////////////////////////
int FluidFluidTransformation::saveClassParameters ( ofstream & outStream) {

  // const char * functionName = "FluidFluidTransformation::saveClassParameters";
  int swapint;
  double swapdouble;
  // Write out size of class name and the class name
  const int nameLength = strlen(mClassName);
  swapint = big_endian_int(nameLength);
  outStream.write((const char *) & swapint, sizeof(nameLength));
  outStream.write(mClassName,                  nameLength);

  // Write out the class version
  const int temp  = mClassVersion;    // Because the compiler may not put these variables
  swapint =big_endian_int(temp); 
  outStream.write((const char *) & swapint, sizeof (temp));

  // Write out size of class revision and the class revision
  const int revisionLength = strlen(mClassRevision);
  swapint = big_endian_int(revisionLength);
  outStream.write((const char *) & swapint, sizeof(revisionLength));
  outStream.write(mClassRevision, revisionLength);

  // Write out parameters for this class
  _fluidLandmarkTransform.saveClassParameters( outStream );
  _fluidTransform.saveClassParameters( outStream );

  swapint = (big_endian_int(_atlasSubSizeX));
  outStream.write((const char *) &swapint, sizeof(_atlasSubSizeX));
  
  swapint = (big_endian_int(_atlasSubSizeY));
  outStream.write((const char *) &swapint, sizeof(_atlasSubSizeY));
  
  swapint = (big_endian_int(_atlasSubSizeZ));
  outStream.write((const char *) &swapint, sizeof(_atlasSubSizeZ));

  swapdouble = (big_endian_double(_atlasSubResX));
  outStream.write((const char *) &swapdouble, sizeof(_atlasSubResX));
  
  swapdouble = (big_endian_double(_atlasSubResY));
  outStream.write((const char *) &swapdouble, sizeof(_atlasSubResY));
  
  swapdouble = (big_endian_double(_atlasSubResZ));
  outStream.write((const char *) &swapdouble, sizeof(_atlasSubResZ));
  
  swapdouble = (big_endian_double(_atlasSubOffsetX));
  outStream.write((const char *) &swapdouble, sizeof(_atlasSubOffsetX));
  
  swapdouble = (big_endian_double(_atlasSubOffsetY));
  outStream.write((const char *) &swapdouble, sizeof(_atlasSubOffsetY));
  
  swapdouble = (big_endian_double(_atlasSubOffsetZ));
  outStream.write((const char *) &swapdouble, sizeof(_atlasSubOffsetZ));
  
  swapint = (big_endian_int(_patientSubSizeX));
  outStream.write((const char *) &swapint, sizeof(_patientSubSizeX));
  
  swapint = (big_endian_int(_patientSubSizeY));
  outStream.write((const char *) &swapint, sizeof(_patientSubSizeY));
  
  swapint = (big_endian_int(_patientSubSizeZ));
  outStream.write((const char *) &swapint, sizeof(_patientSubSizeZ));
  
  swapdouble = (big_endian_double(_patientSubResX));
  outStream.write((const char *) &swapdouble, sizeof(_patientSubResX));
  
  swapdouble = (big_endian_double(_patientSubResY));
  outStream.write((const char *) &swapdouble, sizeof(_patientSubResY));
  
  swapdouble = (big_endian_double(_patientSubResZ));
  outStream.write((const char *) &swapdouble, sizeof(_patientSubResZ));
  
  swapdouble = (big_endian_double(_patientSubOffsetX));
  outStream.write((const char *) &swapdouble, sizeof(_patientSubOffsetX));
  
  swapdouble = (big_endian_double(_patientSubOffsetY));
  outStream.write((const char *) &swapdouble, sizeof(_patientSubOffsetY));
  
  swapdouble = (big_endian_double(_patientSubOffsetZ));
  outStream.write((const char *) &swapdouble, sizeof(_patientSubOffsetZ));

  // Save parent class's parameters
  if ( FieldTransformation::saveClassParameters( outStream ) )
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
FluidFluidTransformation::loadClassParameters ( ifstream & inStream)
{
  const char * const functionName = 
	"FluidFluidTransformation::loadClassParameters";

  // Read in the length of the class name from the file
  int ReadInClassNameLength;
  inStream.read((char *) &ReadInClassNameLength, sizeof (ReadInClassNameLength));
   ReadInClassNameLength = big_endian_int( ReadInClassNameLength);
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
  _fluidLandmarkTransform.loadClassParameters( inStream );
  _fluidTransform.loadClassParameters( inStream );

  inStream.read((char *) &_atlasSubSizeX, sizeof(_atlasSubSizeX));
  inStream.read((char *) &_atlasSubSizeY, sizeof(_atlasSubSizeY));
  inStream.read((char *) &_atlasSubSizeZ, sizeof(_atlasSubSizeZ));
  inStream.read((char *) &_atlasSubResX, sizeof(_atlasSubResX));
  inStream.read((char *) &_atlasSubResY, sizeof(_atlasSubResY));
  inStream.read((char *) &_atlasSubResZ, sizeof(_atlasSubResZ));
  inStream.read((char *) &_atlasSubOffsetX, sizeof(_atlasSubOffsetX));
  inStream.read((char *) &_atlasSubOffsetY, sizeof(_atlasSubOffsetY));
  inStream.read((char *) &_atlasSubOffsetZ, sizeof(_atlasSubOffsetZ));
  inStream.read((char *) &_patientSubSizeX, sizeof(_patientSubSizeX));
  inStream.read((char *) &_patientSubSizeY, sizeof(_patientSubSizeY));
  inStream.read((char *) &_patientSubSizeZ, sizeof(_patientSubSizeZ));
  inStream.read((char *) &_patientSubResX, sizeof(_patientSubResX));
  inStream.read((char *) &_patientSubResY, sizeof(_patientSubResY));
  inStream.read((char *) &_patientSubResZ, sizeof(_patientSubResZ));
  inStream.read((char *) &_patientSubOffsetX, sizeof(_patientSubOffsetX));
  inStream.read((char *) &_patientSubOffsetY, sizeof(_patientSubOffsetY));
  inStream.read((char *) &_patientSubOffsetZ, sizeof(_patientSubOffsetZ));
  //do byte swapping if neccesary
  _atlasSubSizeX = (big_endian_int(_atlasSubSizeX));
  _atlasSubSizeY = (big_endian_int(_atlasSubSizeY));
  _atlasSubSizeZ = (big_endian_int(_atlasSubSizeZ));
  _atlasSubResX  = (big_endian_double(_atlasSubResX));
  _atlasSubResY  = (big_endian_double(_atlasSubResY));
  _atlasSubResZ  = (big_endian_double(_atlasSubResZ));
  _atlasSubOffsetX = (big_endian_double(_atlasSubOffsetX));
  _atlasSubOffsetY = (big_endian_double(_atlasSubOffsetY));
  _atlasSubOffsetZ = (big_endian_double(_atlasSubOffsetZ));
  _patientSubSizeX = (big_endian_int(_patientSubSizeX));
  _patientSubSizeY = (big_endian_int(_patientSubSizeY));
  _patientSubSizeZ = (big_endian_int(_patientSubSizeZ));
  _patientSubResX = (big_endian_double(_patientSubResX));
  _patientSubResY = (big_endian_double(_patientSubResY));
  _patientSubResZ = (big_endian_double(_patientSubResZ));
  _patientSubOffsetX = (big_endian_double(_patientSubOffsetX));
  _patientSubOffsetY = (big_endian_double(_patientSubOffsetY));
  _patientSubOffsetZ = (big_endian_double(_patientSubOffsetZ));

  if ( FieldTransformation::loadClassParameters( inStream ) )
  {
    return 5;
  };

  // Return with success
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
FluidFluidTransformation::printClassParameters ( ostream & outStream)
{

  FieldTransformation::printClassParameters ( outStream );

  outStream << mClassName << "::mClassVersion: "  << mClassVersion << endl;
  outStream << mClassName << "::mClassRevision: " << mClassRevision << endl;
  outStream << mClassName << "::_atlasSubSizeX: " << _atlasSubSizeX << endl;
  outStream << mClassName << "::_atlasSubSizeY: " << _atlasSubSizeY << endl;
  outStream << mClassName << "::_atlasSubSizeZ: " << _atlasSubSizeZ << endl;
  outStream << mClassName << "::_atlasSubResX: " << _atlasSubResX << endl;
  outStream << mClassName << "::_atlasSubResY: " << _atlasSubResY << endl;
  outStream << mClassName << "::_atlasSubResZ: " << _atlasSubResZ << endl;
  outStream << mClassName << "::_atlasSubOffsetX: " << _atlasSubOffsetX << endl;
  outStream << mClassName << "::_atlasSubOffsetY: " << _atlasSubOffsetY << endl;
  outStream << mClassName << "::_atlasSubOffsetZ: " << _atlasSubOffsetZ << endl;
  outStream << mClassName << "::_patientSubSizeX: " << _patientSubSizeX << endl;
  outStream << mClassName << "::_patientSubSizeY: " << _patientSubSizeY << endl;
  outStream << mClassName << "::_patientSubSizeZ: " << _patientSubSizeZ << endl;
  outStream << mClassName << "::_patientSubResX: " << _patientSubResX << endl;
  outStream << mClassName << "::_patientSubResY: " << _patientSubResY << endl;
  outStream << mClassName << "::_patientSubResZ: " << _patientSubResZ << endl;
  outStream << mClassName << "::_patientSubOffsetX: " << _patientSubOffsetX << endl;
  outStream << mClassName << "::_patientSubOffsetY: " << _patientSubOffsetY << endl;
  outStream << mClassName << "::_patientSubOffsetZ: " << _patientSubOffsetZ << endl;

  _fluidLandmarkTransform.printClassParameters( outStream );
  _fluidTransform.printClassParameters( outStream );

  return outStream;

}

////////////////////
// Private Methods
////////////////////

////////////////////////////////////////////////////////////////////////////////
//
// Function: getHField
//
// Purpose: To be determined
//
// Inputs: hField            - must be preallocated before calling this method
//         outputRes[XYZ]    - resolution of output volume in mm/voxel
//         outputOffset[XYZ] - offset of output volume relative to
//                             full volume patient coordinates in mm
//
// Outputs: H(x) field that maps from output subvolume coordinates to 
//          full volume atlas coordinates
//
////////////////////////////////////////////////////////////////////////////////

void FluidFluidTransformation::getHField(Array3D<float> *hField, 
					    double outputResX, 
					    double outputResY,
					    double outputResZ,
					    double outputOffsetX, // in mm
					    double outputOffsetY,
					    double outputOffsetZ) {
  cout << "FluidFluidTransformation::getHField() started at: " << endl;
  system("date");

  // get fluid H field
  Array3D<float> fluidHField(_patientSubSizeZ, _patientSubSizeY, 
			     _patientSubSizeX*3);
  _fluidTransform.getHField(&fluidHField);
  Array3D<float> fluidLandmarkHField(_patientSizeZ, _patientSizeY, 
				     _patientSizeX*3);
  _fluidLandmarkTransform.getHField(&fluidLandmarkHField);

  cout << "_fluidTransform.getHField() finished at: " << endl;
  system("date");
  
  // get fluid H field dimensions
  int fluidHFieldSizeX = fluidHField.getXsize()/3;
  int fluidHFieldSizeY = fluidHField.getYsize();
  int fluidHFieldSizeZ = fluidHField.getZsize();
  cout << "fluidHFieldSize = "
       << fluidHFieldSizeX << ", "
       << fluidHFieldSizeY << ", "
       << fluidHFieldSizeZ << endl;
  
  // get dimensions
  int hFieldSizeX = hField->getXsize()/3;
  int hFieldSizeY = hField->getYsize();
  int hFieldSizeZ = hField->getZsize();

#if 0   
  // get landmarks
  int numLandmarks = _manifoldTransform.getNumPoints(); 
  Array2D<double> patientLandmarks = _manifoldTransform.getPatientPoints(); 
  
  // get manifold parameters
  Matrix<double> Beta = _manifoldTransform.getGreensCoeff();
  Matrix<double> GL3 = _manifoldTransform.getAffineMatrix(); 
  double G00 = GL3[0][0], G01 = GL3[0][1], G02 = GL3[0][2]; 
  double G10 = GL3[1][0], G11 = GL3[1][1], G12 = GL3[1][2]; 
  double G20 = GL3[2][0], G21 = GL3[2][1], G22 = GL3[2][2]; 
  Vector<double> TR = _manifoldTransform.getTranslationVector(); 
  double tx = TR[0], ty = TR[1], tz = TR[2]; 
#endif
  
  // compute the transformation field H(x) = H_M(H_F(x))
  
  float fluidHx;
  float fluidHy;
  float fluidHz;
  float Hx;
  float Hy;
  float Hz;

#if 0
  double manifoldUX;
  double manifoldUY;
  double manifoldUZ;
  
  double **betaPtrPtr = Beta.address();
  double **landmarkPtrPtr = patientLandmarks.address();
#endif
  double xx, yy, zz; 
  double nrm;
  
  int x, y, z, k;
  float *hFieldPtr = hField->data();
  cerr << "slice     of " << hFieldSizeZ-1 << "\010\010\010\010\010\010\010";
  for (z = 0; z < hFieldSizeZ; z++) {
    cerr << "\010\010\010" ;
    cerr.width(3);
    cerr << z;
    for (y = 0; y < hFieldSizeY; y++) {
      for (x = 0; x < hFieldSizeX; x++) {
	
	// coordinate systems:
	// 1. output 
	// 2. patient
	// 3. patient subvolume which includes:
	//    a. post fluid subvolume wrt patient
	//    b. pre fluid subvolume wrt patient
	// 5. atlas subvolume
	// 6. atlas

	// convert x, y, z to patient coordinates in mm
	float patientX = outputResX * x + outputOffsetX;
	float patientY = outputResY * y + outputOffsetY;
	float patientZ = outputResZ * z + outputOffsetZ;
	
	// convert to fluid subvolume coordinates in patient sub voxels
	// patientSubVolX = (patientVolX - patientSubVolOffsetX)/patientSubResX
	// should be floats
	float patientSubX = (patientX - _patientSubOffsetX)/_patientSubResX;
	float patientSubY = (patientY - _patientSubOffsetY)/_patientSubResY;
	float patientSubZ = (patientZ - _patientSubOffsetZ)/_patientSubResZ;
	
	// check to see if it is inside the subvolume
	if ((patientSubX >= 0) && (patientSubX < fluidHFieldSizeX-1) &&
	    (patientSubY >= 0) && (patientSubY < fluidHFieldSizeY-1) &&
	    (patientSubZ >= 0) && (patientSubZ < fluidHFieldSizeZ-1)) {
	  
	  trilinearHField(fluidHField, patientSubX, patientSubY, patientSubZ,
			  &fluidHx, &fluidHy, &fluidHz);
	  
	  // convert to patient coordinates in mm
	  fluidHx = (fluidHx * _patientSubResX) + _patientSubOffsetX;
	  fluidHy = (fluidHy * _patientSubResY) + _patientSubOffsetY;
	  fluidHz = (fluidHz * _patientSubResZ) + _patientSubOffsetZ;
	  
	} else {
	  fluidHx = patientX;
	  fluidHy = patientY;
	  fluidHz = patientZ;
	}

	// convert to patient coordinates in voxels

	fluidHx /= _patientResX;
	fluidHy /= _patientResY;
	fluidHz /= _patientResZ;

	// do the manifold thang
	// H(x) = Hm(Hf(x))
	trilinearHField(fluidLandmarkHField, fluidHx, fluidHy, fluidHz,
			&Hx, &Hy, &Hz);
	*hFieldPtr++ = Hx;
	*hFieldPtr++ = Hy;
	*hFieldPtr++ = Hz;
	
	
#if 0
	// do the manifold thang
	// H(x) = Hm(Hf(x)) = Hf(x) - Um(Hf(x))
	manifoldUX = (G00 * Hx + G01 * Hy + G02 * Hz + tx); 
	manifoldUY = (G10 * Hx + G11 * Hy + G12 * Hz + ty); 
	manifoldUZ = (G20 * Hx + G21 * Hy + G22 * Hz + tz); 
	
	for (k=0; k<numLandmarks; k++) {
	  xx = Hx - landmarkPtrPtr[k][0];
	  yy = Hy - landmarkPtrPtr[k][1];
	  zz = Hz - landmarkPtrPtr[k][2]; 
	  nrm = sqrt(xx*xx + yy*yy + zz*zz);
	  manifoldUX += betaPtrPtr[k][0] * nrm;
	  manifoldUY += betaPtrPtr[k][1] * nrm;
	  manifoldUZ += betaPtrPtr[k][2] * nrm;
	}
	*hFieldPtr++ = (Hx - manifoldUX)/_atlasResX;
	*hFieldPtr++ = (Hy - manifoldUY)/_atlasResY;
	*hFieldPtr++ = (Hz - manifoldUZ)/_atlasResZ;
#endif
      }
    }
  }
  system("date");
  
  cout << "FluidFluidTransformation::getHField() finished at: " << endl;
  system("date"); 

}


void FluidFluidTransformation::extractSubvolume(const Array3D<unsigned char> &volume,
						   const int subOffsetInVoxelsX,
						   const int subOffsetInVoxelsY,
						   const int subOffsetInVoxelsZ,
						   const double magFactorX,
						   const double magFactorY,
						   const double magFactorZ,
						   Array3D<unsigned char> &subvolume) {

  Matrix<double> A(3,3);
  Vector<double> b(3);

  A = 0.0;
  A[0][0] = 1.0/magFactorX;
  A[1][1] = 1.0/magFactorY;
  A[2][2] = 1.0/magFactorZ;

  b[0] = subOffsetInVoxelsX;
  b[1] = subOffsetInVoxelsY;
  b[2] = subOffsetInVoxelsZ;

  UserDefinedLinearTransformation volumeToSubvolume;
  volumeToSubvolume.setTransformation(A, b);
  volumeToSubvolume.apply(volume, &subvolume);
}


