///////////////////////////////////////////////////////////////////////////
//
// FluidTransformation class
// 
// This is the class that performs the Fluid deformation
//
// Note: If you are changing this class.  The version number of this class
// should change is more or less data is saved so that we may be able to
// load in old versions.
//
// To Do:
//
// Things to think about...
//   1) Maybe put in code that can read in older version of the class
//
//////////////////////////////////////////////////////////////////////////
//
// Revision Log:
//
// $Log: FluidTransformation.C,v $
// Revision 1.1  2004/11/15 04:44:08  joeh
// first check in of BrainWorks
//
// Revision 1.2  1999/11/09 21:28:26  kem
// Add methods for computing h field values given a point
//
// Revision 1.1  1999/10/01 15:29:15  kem
// Initial revision
//
// Revision 1.22  1999/07/09 17:52:14  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.21  1999/03/16 23:37:22  kem
// Fix resolution passing
//
// Revision 1.20  1999/01/20 22:12:19  RAZORDB
// change class "Timer" to "ITXTimer"
//
// Revision 1.19  1998/12/10 19:04:54  RAZORDB
// IDL/AW merge
//
// Revision 1.18  1998/11/03 18:08:29  RAZORDB
// Fix masking of datasets
//
// Revision 1.17  1998/10/30 22:42:35  RAZORDB
// Add Masking of infiltarting Tumors and prostrate
//
// Revision 1.16  1998/10/29 00:19:43  RAZORDB
// Fix elastic inner iterations
//
// Revision 1.15  1998/10/15 21:45:14  kem
// Change tolerance for basis incrementation
//
// Revision 1.14  1998/07/15 21:19:39  kem
// Timing and optimization
//
// Revision 1.13  1998/07/10 19:18:49  kem
// Clean up unused methods and variables
//
// Revision 1.12  1998/07/10 15:57:38  sarang
// To combine Elastic and Fluid Transformation
//
// Revision 1.11  1998/06/05 20:16:02  rst
// Handle old transformation file versions
//
// Revision 1.10  1998/06/02 15:24:29  rst
// Changed fluid to support elastic
//
// Revision 1.9  1998/04/09 19:12:00  kem
// Balance Lv(x) = b(v(x))
//
// Revision 1.8  1998/02/13 22:11:17  csb
// Implemented different eigenvectors
//
// Revision 1.7  1997/12/12 22:14:22  csb
// merge 1.5.1.2 with 1.6
//
// Revision 1.6  1997/12/12 19:33:59  abed
// Change getHfield Prototype
//
// Revision 1.5  1997/11/18 21:58:33  csb
// fix negative delta problem
//
// Revision 1.4  1997/09/17 19:34:58  kem
// Initialize _atlasSize? and _patientSize? variables in initialize()
//
// Revision 1.3  1997/09/16 19:45:04  kem
// Remove getUField() method
//
// Revision 1.2  1997/08/25 17:43:00  kem
// Add getHfield()
//
// Revision 1.1  1997/08/01 19:50:03  csb
// Initial revision
//
// Revision 1.1  1997/08/01 15:49:10  csb
// Initial revision
//
// Revision 1.1  1997/08/01 15:46:55  csb
// Initial revision
//
// Revision 1.3  1997/06/24 22:24:47  csb
// Added Revision, changed order of saving and loading
//
// Revision 1.2  1997/06/02 21:58:09  csb
// Changes to function calls
//
// Revision 1.1  1997/05/30 23:59:14  csb
// Initial revision
//
//
//////////////////////////////////////////////////////////////////////////

#include <string.h>
#include <math.h>
#include <TDL/ITX_FFT.h>
#include <TDL/TDLfft.h>
#include <TDL/FluidTransformation.h>
#include <TDL/TransformationUtils.h>
#include <TDL/ManifoldTransformation.h>
#include <TDL/Timer.h>


#define USE_FFTS 1
#define COMPBODY 0
#define LLADJ 0

const char * const FluidTransformation::mClassRevision = "$Id: FluidTransformation.C,v 1.1 2004/11/15 04:44:08 joeh Exp $";
const char * const FluidTransformation::mClassName = "FluidTransformation";
const int FluidTransformation::mClassVersion = 4;

/***************
* Ctor and Dtor.
***************/

FluidTransformation::FluidTransformation():
  _niter(25),
  _nInnerIter(1),
  _numBasis(2000),
  _alpha(0.01), 
  _beta(0.01), 
  _gamma(0.0001), 
  _delta(0.01), 
  _mxper(0.2), 
  _mnjack( 0.5 ), 
  _TimMagicStoppingNumber(1e-4),
  _CalculateTotalJacobianFlag(0)
{ /* empty */ }


FluidTransformation::~FluidTransformation() {
  if (getVerboseMode())
    cout << "Calling FluidTransformation destructor" << endl;
}




////////////////////////////////////////////////////////////////////////////////
//
// Function: getHField
//
// Purpose: To be determined
//
// Inputs: To be determined
//
// Outputs: To be determined 
//
////////////////////////////////////////////////////////////////////////////////
void FluidTransformation::getHField(Array3D<float> *hField)
{
*hField = _TotalHField;
}

void FluidTransformation::getHField(Array3D<float> *hField,
        double Rx, double Ry, double Rz, double Ox, double Oy, double Oz) 
{
if((hField->getXsize() == _TotalHField.getXsize()) &&
   (hField->getYsize() == _TotalHField.getYsize()) &&
   (hField->getZsize() == _TotalHField.getZsize()) &&
   (Rx == _patientResX) && (Ry == _patientResY) && (Rz == _patientResZ) &&
   (Ox == 0.0) && (Oy ==0.0) && (Oz == 0.0)) {
	*hField = _TotalHField;
	} else {

 // get fluid H field dimensions
  int hFieldSizeX = hField->getXsize()/3;
  int hFieldSizeY = hField->getYsize();
  int hFieldSizeZ = hField->getZsize();

  int THfilednX = _TotalHField.getXsize()/3;
  int THfilednY = _TotalHField.getYsize();
  int THfilednZ = _TotalHField.getZsize();

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
	trilinearHField(_TotalHField, patientX, patientY, patientZ,
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

////////////////////////////////////////////////////////////////////////////////
//
// Function: getHField
//
// Purpose: To be determined
//
// Inputs: patient[XYZ] - point in absolute patient space in voxels
//
// Outputs: H(x) field that maps from output subvolume coordinates to 
//          full volume atlas coordinates in voxels
//
////////////////////////////////////////////////////////////////////////////////

void FluidTransformation::getHField(const float patientInVoxelsX,
				    const float patientInVoxelsY,
				    const float patientInVoxelsZ,
				    float &atlasInVoxelsX,
				    float &atlasInVoxelsY,
				    float &atlasInVoxelsZ) {
  
  // get fluid H field dimensions
  int hFieldSizeX = _TotalHField.getXsize()/3;
  int hFieldSizeY = _TotalHField.getYsize();
  int hFieldSizeZ = _TotalHField.getZsize();
  
  float Hx;
  float Hy;
  float Hz;
  
  // check to see if it is inside the subvolume
  if ((patientInVoxelsX >= 0) && (patientInVoxelsX < hFieldSizeX-1) &&
      (patientInVoxelsY >= 0) && (patientInVoxelsY < hFieldSizeY-1) &&
      (patientInVoxelsZ >= 0) && (patientInVoxelsZ < hFieldSizeZ-1)) {
	  
    trilinearHField(_TotalHField, 
		    patientInVoxelsX, patientInVoxelsY, patientInVoxelsZ,
		    &Hx, &Hy, &Hz);
  } else {
    Hx = patientInVoxelsX;
    Hy = patientInVoxelsY;
    Hz = patientInVoxelsZ;
  }

  atlasInVoxelsX = Hx;
  atlasInVoxelsY = Hy;
  atlasInVoxelsZ = Hz;
}



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
int
FluidTransformation::initialize( const Array3D<unsigned char> * atlasPtr, 
                                 const Array3D<unsigned char> * patientPtr,
				 const Array3D<unsigned char> * maskPtr,
				 float aresx, float aresy, float aresz,
				 float presx, float presy, float presz)
{
  if (getVerboseMode())
    cout << "Executing FluidTransformation::initialize()" << endl;

  _atlasSizeX = atlasPtr->getXsize();
  _atlasSizeY = atlasPtr->getYsize();
  _atlasSizeZ = atlasPtr->getZsize();

  _patientSizeX = patientPtr->getXsize();
  _patientSizeY = patientPtr->getYsize();
  _patientSizeZ = patientPtr->getZsize();

  _atlasResX = aresx;
  _atlasResY = aresy;
  _atlasResZ = aresz;

  _patientResX = presx;
  _patientResY = presy;
  _patientResZ = presz;


  // Dimensions.
  _nz = patientPtr->getZsize(); 
  _ny = patientPtr->getYsize(); 
  _nx = patientPtr->getXsize(); 
  _nxPlusStride = 2 * ((_nx + 2) / 2);
  _xstride = _nxPlusStride - _nx; 

  // Max index of bases functions.
  _xmax = (_nx + 2) / 2; 
  _ymax = (_ny + 2) / 2; 
  _zmax = (_nz + 2) / 2; 

  // Set atlas and patient pointers to point to the right place
  _AtlasPtr   = atlasPtr; 
  _PatientPtr = patientPtr; 

  // Set the Mask Pointer to the right place
  _maskPtr = maskPtr;


  // Set _deformed to be same size as patient and set to zero
  _deformed.setDim(_nz, _ny, _nx);
  _deformed = (unsigned char) 0;
  
  // set the difference operator deltas
  cout << "Defaulting deltas to 1" << endl;
      
  _Delta1 = 1.0;
  _Delta2 = 1.0;
  _Delta3 = 1.0;
  
  // difference operator vectors
  _cosW1.setDim(_nx);
  _cosW2.setDim(_ny);
  _cosW3.setDim(_nz);
  _sinW1.setDim(_nx);
  _sinW2.setDim(_ny);
  _sinW3.setDim(_nz);
  

  // Velocity vector fields.
  _Vx.setDim(_nz, _ny, _nxPlusStride); 
  _Vy.setDim(_nz, _ny, _nxPlusStride);
  _Vz.setDim(_nz, _ny, _nxPlusStride);
  _Vx = 0.f;
  _Vy = 0.f;
  _Vz = 0.f;

  // Total displacement vector fields.
  _TotalHField.setDim(_nz, _ny, _nx * 3);
  // Initialize to no displacment
  for (int z = 0; z < _nz; z++) {
    for (int y = 0; y < _ny; y++) {
      for (int x = 0; x < _nx; x++) {
	_TotalHField[z][y][3*x + 0] = x;
	_TotalHField[z][y][3*x + 1] = y;
	_TotalHField[z][y][3*x + 2] = z;
      }
    }
  }

  if ( _CalculateTotalJacobianFlag ) {
    mTotalJacobian.setDim(_nz, _ny, _nx);
    mTotalJacobian = (float) 1.0;
    
    mTotalJacobianScratch.setDim(_nz, _ny, _nx);
    mTotalJacobianScratch  = (float) 1.0;
    
    mTotalJacobianMin = 1;
    mTotalJacobianMinIndexX = 0;
    mTotalJacobianMinIndexY = 0;
    mTotalJacobianMinIndexZ = 0;
    mIncrementalJacobianMin = 1;
    mIncrementalJacobianMinIndexX = 0;
    mIncrementalJacobianMinIndexY = 0;
    mIncrementalJacobianMinIndexZ = 0;
  } else {
    // Just make sure that these are not allocated
    mTotalJacobian.setDim(0, 0, 0);
    mTotalJacobianScratch.setDim(0, 0, 0);
  }


//  mfft.initialize(3,ITX_FFT::R2C,_nx,_ny,_nz);
  
  // Generate lookup tables 
  this->genTables(); 
  
  // Done with initialize.
  setOkToCalculate( 1 );
  
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
//
////////////////////////////////////////////////////////////////////////////////
int FluidTransformation::calculate() {

  if (getVerboseMode())
    cout << "Executing FluidTransformation::calculate()" << endl;
  
  cout.precision(5);
  cout.setf(ios::scientific, ios::floatfield);
  cerr.precision(5);
  cerr.setf(ios::scientific, ios::floatfield);
  

  // Check if its ok to calculate.
  if (!isOkToCalculate())
    throwError("FluidTransformation::calculate: Unable to calculate trasformation.");
  

  // Run the engine
  
#if COMPBODY
  Array3D<float> lastBodyForceX(_nz, _ny, _nxPlusStride);
  Array3D<float> lastBodyForceY(_nz, _ny, _nxPlusStride);
  Array3D<float> lastBodyForceZ(_nz, _ny, _nxPlusStride);
  float bodyMeanSquaredErrorX;
  float bodyMeanSquaredErrorY;
  float bodyMeanSquaredErrorZ;
#endif

  // kwd
  ShowWorking("Fluid Trans::calculating()\n");

  // storage for mu's needed for elastic algorithm
  Array3D<float> muX(_nz, _ny, _nxPlusStride);
  Array3D<float> muY(_nz, _ny, _nxPlusStride);
  Array3D<float> muZ(_nz, _ny, _nxPlusStride);

  Array3D<float> uX(_nz, _ny, _nxPlusStride);
  Array3D<float> uY(_nz, _ny, _nxPlusStride);
  Array3D<float> uZ(_nz, _ny, _nxPlusStride);

  Array3D<unsigned char> defAtlas(_AtlasPtr->getDim());
  Array3D<float> newHField;
  Array3D<float> AtlasGrad;
  AtlasGrad.setDim(_TotalHField.getDim());
  newHField.setDim(_TotalHField.getDim());
  //   newHField = _TotalHField;
	
  int isDeltaSet = 0;
  double lastSquaredError=1.0;
  double squaredError;
  
  int StopYetFlag = 0;
  
  int doElasticFlag = (_nInnerIter!=1);
  int i,x,y,z;

  // Outerloop does timestep iterations, _niter (Fluid)
  // The inner loop does Elastic iterations, _nInnerIter
  for (i = 0; (i < _niter) && (!StopYetFlag); i++) {
    
#if COMPBODY
    lastBodyForceX = 0.0;
    lastBodyForceY = 0.0;
    lastBodyForceZ = 0.0;
#endif

    // Initialize to zero displacement
    for( z = 0; z < _nz; z++) {
      for(y = 0; y < _ny; y++) {
	for(x = 0; x < _nx; x++) {
	   newHField[z][y][3*x+0] = x;
	   newHField[z][y][3*x+1] = y;
	   newHField[z][y][3*x+2] = z;
	}
      }
    }

    TransformationUtils::interpolate(*_AtlasPtr, defAtlas, _TotalHField);

    // Calculate the Gradiant of the curent Template
    
    const unsigned char * const * const *defAtlasPtrPtrPtr = defAtlas.address();
    float dx0, dx1, dx2;

    for(z = 0; z < _nz; z++) {
      for(y = 0; y < _ny; y++) {
	for(x = 0; x < _nx; x++) {

	  if (x > 0 && x < (_nx-1)) {
	    dx0 = (defAtlas[z][y][x+1] - defAtlas[z][y][x-1])/2.0;
	  } else if (x == 0) {
	    dx0 = (defAtlas[z][y][x+1] - defAtlas[z][y][x]);
	  } else {
	    dx0 = (defAtlas[z][y][x]   - defAtlas[z][y][x-1]);
	  }
	  if (y > 0 && y < (_ny-1)) {
	    dx1 = (defAtlas[z][y+1][x] - defAtlas[z][y-1][x])/2.0;
	  } else if (y == 0) {
	    dx1 = (defAtlas[z][y+1][x] - defAtlas[z][y][x]);
	  } else {
	    dx1 = (defAtlas[z][y][x]   - defAtlas[z][y-1][x]);
	  }
	  if (z > 0 && z < (_nz-1)) {
	    dx2 = (defAtlas[z+1][y][x] - defAtlas[z-1][y][x])/2.0;
	  } else if (z == 0) {
	    dx2 = (defAtlas[z+1][y][x] - defAtlas[z][y][x]);
	  } else {
	    dx2 = (defAtlas[z][y][x]   - defAtlas[z-1][y][x]);
	  }	

	   AtlasGrad[z][y][3*x+0] = dx0;
	   AtlasGrad[z][y][3*x+1] = dx1;
	   AtlasGrad[z][y][3*x+2] = dx2;
	}
      }
    }

    float elasticDelta = 1.0e-4;
    double lastPosterior = 1.0e50;

    muX = 0.0;
    muY = 0.0;
    muZ = 0.0;

    for (int innerIter=0; innerIter<_nInnerIter; innerIter++) {
      
      // Calculate the body force
      genBodyForce(newHField, defAtlas,AtlasGrad, *_PatientPtr, 
		   _deformed, _Vx, _Vy, _Vz, squaredError);

#if COMPBODY
      // compare and copy the new body force to the last body force 
      compareBodyForce(_Vx, _Vy, _Vz, 
		       lastBodyForceX, lastBodyForceY, lastBodyForceZ,
		       bodyMeanSquaredErrorX,
		       bodyMeanSquaredErrorY,
		       bodyMeanSquaredErrorZ);
#endif
      
      // Increment the number of basis if not changing
      if (( fabs(lastSquaredError - squaredError )/lastSquaredError < 0.0001) &&
	  ((_numBasis < _nz-1) || (_numBasis < _ny-1) || (_numBasis < _xmax))) {
	elasticDelta = 1.0e-4;
	_numBasis+=4;
      }

      char iter_str[50];
      sprintf(iter_str,"iteration (%d %d)\n",i,innerIter);
      ShowWorking(iter_str);

      cout << "Iteration: " << i << " " << innerIter << " "
#if COMPBODY
	   << bodyMeanSquaredErrorX << " "
	   << bodyMeanSquaredErrorY << " " 
	   << bodyMeanSquaredErrorZ << " "
#endif
	   << squaredError << " NB =" << _numBasis-1 << endl;


      lastSquaredError = squaredError;
      
      if (!doElasticFlag) {

	// Calculate the velocity field.
	genVelocityField(_Vx, _Vy, _Vz); 

	if (!isDeltaSet) {
	  // Calculate which delta to use
	  _delta = FindBestDelta(&StopYetFlag); 
	  isDeltaSet = 1;
	}
      
	// Update displacement field.
	updateGlobalDisp(_delta, _Vx, _Vy, _Vz, newHField); 

      } else { // use elastic method for solving balance equation

	float priorWeight = 1.0e5;

	// Compute new mu's and displacement field given last mu,
	// body force and delta
	// body force is passed in via _Vx, _Vy, _Vz
	genDirection(muX, muY, muZ, priorWeight,
		     _Vx, _Vy, _Vz);

	elasticDelta = findOptimalDelta(muX, muY, muZ, _Vx, _Vy, _Vz,
					uX, uY, uZ, newHField,
					defAtlas, _deformed, *_PatientPtr,
					priorWeight, elasticDelta, lastPosterior);

	// lastPosterior corresponds to the posterior evaluated at elasticDelta

	// update the mu's, the u-field and the h-field
	updateMus(muX, muY, muZ, elasticDelta, _Vx, _Vy, _Vz,
		  uX, uY, uZ, newHField);

      } // end doElasticOnly branch

    } // end innerIter loop
    
    composeTransField(_TotalHField, newHField);
    _TotalHField = newHField;
    isDeltaSet = 0;
  }

  ShowWorking((char*)NULL);
  
  /******************************************************
   * Say we have a transformation and we may now apply it
   ******************************************************/
  
  setOkToApply( 1 );
  
  return 0;
}





/*****************************************************************************************
 *
 * Calculate the basis normalizing constants and the inverse variances for 
 * each of the Karhunen-Loeve coefficients.
 *
 *****************************************************************************************/
void FluidTransformation::genTables() {
  
  double scale3 = 2.0 * M_PI / _nz; 
  double scale2 = 2.0 * M_PI / _ny; 
  double scale1 = 2.0 * M_PI / _nx; 
  
  int i;
  
  for (i=0; i<_xmax; i++) 
    {
      _cosW1[i] = (2.0 * cos(_Delta1 * scale1 * (double)i) - 2.0) / (_Delta1*_Delta1);
      _sinW1[i] = sin(_Delta1 * scale1 * (double)i) / _Delta1;
    }
  for (i=0; i<_ny; i++) 
    {
      _cosW2[i] = (2.0 * cos(_Delta2 * scale2 * (double)i) - 2.0) / (_Delta2*_Delta2);
      _sinW2[i] = sin(_Delta2 * scale2 * (double)i) / _Delta2;
    }
  for (i=0; i<_nz; i++) 
    {
      _cosW3[i] = (2.0 * cos(_Delta3 * scale3 * (double)i) - 2.0) / (_Delta3*_Delta3);
      _sinW3[i] = sin(_Delta3 * scale3 * (double)i) / _Delta3;
    }
  
}



/********************************************************************************
 *
 * Generate the body force. When this routine is called, it should fill
 * the 3D volumes Vx, Vy, Vz with the value of the body force function.
 *
 ********************************************************************************/
void
FluidTransformation::genBodyForce(const Array3D<float> &hField,
				  const Array3D<unsigned char> &atlas,
				  const Array3D<float> &AtlasGrad,
				  const Array3D<unsigned char> &patient,
				  Array3D<unsigned char> &deformedAtlas,
				  Array3D<float> &bodyForceX,
				  Array3D<float> &bodyForceY,
				  Array3D<float> &bodyForceZ,
				  double &squaredError)
{
  // deform the atlas using the hField
  TransformationUtils::interpolate(atlas, deformedAtlas, hField);

  Array1D<double> G(3); 
   
  squaredError = 0.0; 

  float F[3];

  for (int z = 0; z < _nz; z++) {
    for (int y = 0; y < _ny; y++) {
      for (int x = 0; x < _nx; x++) {
	
	/******************************
	 * Calculate (T(x - ut) - S(x)).
	 ******************************/
	double deformedMinusPatient =(deformedAtlas[z][y][x]-patient[z][y][x]); 

        // mask out if needed
	if (_maskPtr && ((*_maskPtr)[z][y][x] != 0)) 
	    deformedMinusPatient = 0;

	// compute the squared error
	squaredError += deformedMinusPatient * deformedMinusPatient; 
	
	/********************************************
	 * Calculate the gradient of the ith template.
	 *********************************************/
        double hFieldX =  hField[z][y][3*x+0];
        double hFieldY =  hField[z][y][3*x+1];
        double hFieldZ =  hField[z][y][3*x+2];

	TransformationUtils::interpolate(F, AtlasGrad, 
					 hFieldZ, hFieldY, hFieldX, 0); 

	/*************************************
	 * Multiply by gradient of T(i) at (x).
	 *************************************/
	bodyForceX[z][y][x] = -deformedMinusPatient * F[0]; 
	bodyForceY[z][y][x] = -deformedMinusPatient * F[1]; 
	bodyForceZ[z][y][x] = -deformedMinusPatient * F[2]; 
      }
    }
  }
}

void
FluidTransformation::compareBodyForce(const Array3D<float> &bodyForceX,
				      const Array3D<float> &bodyForceY,
				      const Array3D<float> &bodyForceZ,
				      Array3D<float> &lastBodyForceX,
				      Array3D<float> &lastBodyForceY,
				      Array3D<float> &lastBodyForceZ,
				      float &meanSquaredErrorX,
				      float &meanSquaredErrorY,
				      float &meanSquaredErrorZ)
{
  const float * bodyForceXPtr = bodyForceX.data();
  const float * bodyForceYPtr = bodyForceY.data();
  const float * bodyForceZPtr = bodyForceZ.data();
  float * lastBodyForceXPtr = lastBodyForceX.data();
  float * lastBodyForceYPtr = lastBodyForceY.data();
  float * lastBodyForceZPtr = lastBodyForceZ.data();

  meanSquaredErrorX = 0.0;
  meanSquaredErrorY = 0.0;
  meanSquaredErrorZ = 0.0;
  for (int z = 0; z < _nz; z++) {
    for (int y = 0; y < _ny; y++) {
      for (int x = 0; x < _nx; x++) {
	float bodyDiffX = *bodyForceXPtr - *lastBodyForceXPtr;
	float bodyDiffY = *bodyForceYPtr - *lastBodyForceYPtr;
	float bodyDiffZ = *bodyForceZPtr - *lastBodyForceZPtr;
	meanSquaredErrorX += bodyDiffX * bodyDiffX;
	meanSquaredErrorY += bodyDiffY * bodyDiffY;
	meanSquaredErrorZ += bodyDiffZ * bodyDiffZ;
	*lastBodyForceXPtr++ = *bodyForceXPtr++;
	*lastBodyForceYPtr++ = *bodyForceYPtr++;
	*lastBodyForceZPtr++ = *bodyForceZPtr++;
      }
      bodyForceXPtr += _xstride; 
      bodyForceYPtr += _xstride; 
      bodyForceZPtr += _xstride; 
      lastBodyForceXPtr += _xstride; 
      lastBodyForceYPtr += _xstride; 
      lastBodyForceZPtr += _xstride; 
    }
  }
  meanSquaredErrorX /= (_nx*_ny*_nz);
  meanSquaredErrorY /= (_nx*_ny*_nz);
  meanSquaredErrorZ /= (_nx*_ny*_nz);
}


/**************************************************************************************
*
* Given the body force (driving function), calculate the 
* velocity field. The velocity field is stored as the real
* arrays Vx, Vy, Vz.
*
**************************************************************************************/
void
FluidTransformation::genVelocityField(Array3D<float> &bodyToVelocityX,
				      Array3D<float> &bodyToVelocityY,
				      Array3D<float> &bodyToVelocityZ)
{
#if USE_FFTS

  // Perform inverse fft on body force
  // scale by _nx*_ny*_nz because the fft routines do not
  
/*
KWD : dont need this now because TDLfft scales properly
  int z, y, x;
  float nxnynz = (_nx*_ny*_nz);
  for (z = 0; z < _nz; z++) {
    for (y = 0; y < _ny; y++) {
      for (x = 0; x < _nx; x++) {
	// added by kem to scale fft result correctly
	bodyToVelocityX[z][y][x] /= nxnynz;
	bodyToVelocityY[z][y][x] /= nxnynz;
	bodyToVelocityZ[z][y][x] /= nxnynz;
      }
    }
  }
*/

   // Originally we use to take inverse transform but
   // fftw does not allow you to take inverse transform of real data
   // Hence we take forward transform. This should not make in differeance
   // Sarang Joshi 6/17/99

//   mfft.fftForward(bodyToVelocityX.data());
//   mfft.fftForward(bodyToVelocityY.data());
//   mfft.fftForward(bodyToVelocityZ.data());

   mfft.FFT(bodyToVelocityX,TDLfft::Forward);
   mfft.FFT(bodyToVelocityY,TDLfft::Forward);
   mfft.FFT(bodyToVelocityZ,TDLfft::Forward);

#endif

   int ii, i, j, k; 

   float Lmatrix00, Lmatrix01, Lmatrix02;
   float Lmatrix10, Lmatrix11, Lmatrix12;
   float Lmatrix20, Lmatrix21, Lmatrix22;
   float G00, G01, G02;
   float G10, G11, G12;
   float G20, G21, G22;
   float bReal0, bReal1, bReal2;
   float bImag0, bImag1, bImag2;
   float vReal0, vReal1, vReal2;
   float vImag0, vImag1, vImag2;
   float yVector0, yVector1, yVector2;

cerr << "NUM BASIS = " << _numBasis << endl;

   for (k=0; k<_nz; k++) {
     for (int j=0; j<_ny; j++) {
       for (int i=0; i<_xmax; i++) {
	 // To use small number of Basis functions defined by _numBasis
	 // Date 29/06/98 Sarang Kem.
	 if ((i < _numBasis ) && 
	     (((j < _numBasis) || (j > _ny - _numBasis)) && 
	      ((k < _numBasis) || (k > _nz-_numBasis)))) {
	   double lambda = - _alpha * (_cosW1[i] + _cosW2[j] + _cosW3[k]) + _gamma;

	   // changed 5/27, use L Ladj instead of just L
	   // rst and sarang
#if LLADJ
	   float TLmatrix00, TLmatrix01, TLmatrix02;
	   float TLmatrix10, TLmatrix11, TLmatrix12;
	   float TLmatrix20, TLmatrix21, TLmatrix22;

	   // compute L
	   TLmatrix00 = lambda - _beta * _cosW1[i];
	   TLmatrix11 = lambda - _beta * _cosW2[j];
	   TLmatrix22 = lambda - _beta * _cosW3[k];
	   TLmatrix01 = _beta * _sinW1[i] * _sinW2[j];
	   TLmatrix02 = _beta * _sinW1[i] * _sinW3[k];
	   TLmatrix12 = _beta * _sinW2[j] * _sinW3[k];
	   TLmatrix10 = TLmatrix01; // This matrix is symmetric
	   TLmatrix20 = TLmatrix02;
	   TLmatrix21 = TLmatrix12;

	   // now compute L Ladj
	   Lmatrix00 = TLmatrix00 * TLmatrix00 + TLmatrix01*TLmatrix01 +
	     TLmatrix02 * TLmatrix02;
	   Lmatrix01 = TLmatrix00 * TLmatrix10 + TLmatrix01*TLmatrix11 +
	     TLmatrix02 * TLmatrix12;
	   Lmatrix02 = TLmatrix00 * TLmatrix20 + TLmatrix01*TLmatrix21 +
	     TLmatrix02 * TLmatrix22;
	   Lmatrix11 = TLmatrix10 * TLmatrix10 + TLmatrix11*TLmatrix11 +
	     TLmatrix12 * TLmatrix12;
	   Lmatrix12 = TLmatrix10 * TLmatrix20 + TLmatrix11*TLmatrix21 +
	     TLmatrix12 * TLmatrix22;
	   Lmatrix22 = TLmatrix20 * TLmatrix20 + TLmatrix21*TLmatrix21 +
	     TLmatrix22 * TLmatrix22;
	   Lmatrix10 = Lmatrix01; // This matrix is symmetric
	   Lmatrix20 = Lmatrix02;
	   Lmatrix21 = Lmatrix12;
#else
	   // old stuff, just compute L
	   Lmatrix00 = lambda - _beta * _cosW1[i];
	   Lmatrix11 = lambda - _beta * _cosW2[j];
	   Lmatrix22 = lambda - _beta * _cosW3[k];
	   Lmatrix01 = _beta * _sinW1[i] * _sinW2[j];
	   Lmatrix02 = _beta * _sinW1[i] * _sinW3[k];
	   Lmatrix12 = _beta * _sinW2[j] * _sinW3[k];
	   Lmatrix10 = Lmatrix01; // This matrix is symmetric
	   Lmatrix20 = Lmatrix02;
	   Lmatrix21 = Lmatrix12;
#endif

	   bReal0 = bodyToVelocityX[k][j][2*i]; 
	   bImag0 = bodyToVelocityX[k][j][2*i + 1]; 
	   
	   bReal1 = bodyToVelocityY[k][j][2*i]; 
	   bImag1 = bodyToVelocityY[k][j][2*i + 1]; 
	   
	   bReal2 = bodyToVelocityZ[k][j][2*i]; 
	   bImag2 = bodyToVelocityZ[k][j][2*i + 1];
	   
//	   LmatrixInverse = inv(Lmatrix);
	   
//	   vReal = LmatrixInverse * bReal;
//	   vImag = LmatrixInverse * bImag;

	   // Calling this function seems to take massive amounts of system time
	   // so making it inline
	   // backSolve(Lmatrix, bReal, bImag, vReal, vImag);
	   
	   // Giving that A is pos-def symetric matrix, solve Ax=b by finding
	   // cholesky decomposition GG'=A
	   // and then performing 2 back-solves, Gy=b and then G'x=y to get x.
	   // 
	   // Note: lets take the hit for now in having G be a 3x3 array
	   // instead of a 6x1 array.
	   //
	   // Solves for 2 vectors at once
	   
	   // 1. find cholesky decomposition by finding G such that GG'=A.
	   //    A must be positive definite symetric (we assume that here)
	   //    G is then lower triangular, see algorithm 4.2.1 p142-3
	   //    in Golub and VanLoan
	   // Note: these are in matlab notation 1:3
	   // [ G(1,1)   0      0    ]   [ G(1,1) G(2,1) G(3,1) ]   
	   // [ G(2,1) G(2,2)   0    ] * [   0    G(2,2) G(3,2) ] = Amatrix
	   // [ G(3,1) G(3,2) G(3,3) ]   [   0      0    G(3,3) ]

	   G00 = sqrt( Lmatrix00 );
	   G10 = Lmatrix10 / G00;
	   G20 = Lmatrix20 / G00;
	   
	   G11 = Lmatrix11 - G10 * G10;
	   G21 = Lmatrix21 - G20 * G10;
	   G11 = sqrt( G11 );
	   G21 = G21 / G11;
  
	   G22 = Lmatrix22 - ( G20*G20 + G21*G21 );
	   G22 = sqrt(G22 );
	   
	   // back-solve Gy=b to get a temporary vector y
	   // back-solve G'x=y to get answer in x
	   //
	   // Note: these are in matlab notation 1:3
	   // [ G(1,1)   0      0    ]   [ y(1) ] = b(1)
	   // [ G(2,1) G(2,2)   0    ] * [ y(2) ] = b(2)
	   // [ G(3,1) G(3,2) G(3,3) ]   [ y(3) ] = b(3)
	   //
	   // [ G(1,1) G(2,1) G(3,1) ]   [ x(1) ] = y(1)
	   // [   0    G(2,2) G(3,2) ] * [ x(2) ] = y(2)
	   // [   0      0    G(3,3) ]   [ x(3) ] = y(3)

	   yVector0 = bReal0 / G00;
	   yVector1 = ( bReal1 - G10*yVector0) / G11;
	   yVector2 = ( bReal2 - G20*yVector0 - G21*yVector1) / G22;
	   
	   vReal2 = yVector2 / G22;
	   vReal1 = ( yVector1 - G21*vReal2 ) / G11;
	   vReal0 = ( yVector0 - G10*vReal1 - G20*vReal2 ) / G00;
	   
	   yVector0 = bImag0 / G00;
	   yVector1 = ( bImag1 - G10*yVector0) / G11;
	   yVector2 = ( bImag2 - G20*yVector0 - G21*yVector1) / G22;

	   vImag2 = yVector2 / G22;
	   vImag1 = ( yVector1 - G21*vImag2 ) / G11;
	   vImag0 = ( yVector0 - G10*vImag1 - G20*vImag2 ) / G00;

	   bodyToVelocityX[k][j][2*i] = vReal0;
	   bodyToVelocityX[k][j][2*i + 1] = vImag0;
	   
	   bodyToVelocityY[k][j][2*i] = vReal1;
	   bodyToVelocityY[k][j][2*i + 1] = vImag1;
	   
	   bodyToVelocityZ[k][j][2*i] = vReal2;
	   bodyToVelocityZ[k][j][2*i + 1] = vImag2;
	 } else {
	   bodyToVelocityX[k][j][2*i] = 0;
	   bodyToVelocityX[k][j][2*i + 1] = 0;
	   
	   bodyToVelocityY[k][j][2*i] = 0;
	   bodyToVelocityY[k][j][2*i + 1] = 0;
	   
	   bodyToVelocityZ[k][j][2*i] = 0;
	   bodyToVelocityZ[k][j][2*i + 1] = 0;
	 }
       }
     }
   }


#if USE_FFTS

   // Perform inverse fft to get velocity field
//   mfft.fftInverse((FFT_complex *)bodyToVelocityX.data());
//   mfft.fftInverse((FFT_complex *)bodyToVelocityY.data());
//   mfft.fftInverse((FFT_complex *)bodyToVelocityZ.data());
 
   mfft.FFT(bodyToVelocityX,TDLfft::Reverse);
   mfft.FFT(bodyToVelocityY,TDLfft::Reverse);
   mfft.FFT(bodyToVelocityZ,TDLfft::Reverse);
#endif

}

/********************************************************************************
*
* Given the body force (driving function), calculate the 
* velocity field. The velocity field is stored as the real
* arrays Vx, Vy, Vz.
*
*********************************************************************************/
void
FluidTransformation::genDirection(Array3D<float> &muX,
				  Array3D<float> &muY,
				  Array3D<float> &muZ,
				  float priorWeight,
				  Array3D<float> &directionX,
				  Array3D<float> &directionY,
				  Array3D<float> &directionZ)
{
  // Perform inverse fft on body force

  int z, y, x;
  for (z = 0; z < _nz; z++) 
    for (y = 0; y < _ny; y++) 
      for (x = 0; x < _nx; x++) {
	directionX[z][y][x] *= -1.f;
	directionY[z][y][x] *= -1.f;
	directionZ[z][y][x] *= -1.f;
      }
  
  // scale by _nx*_ny*_nz because the fft routines do not

#ifdef notdef
KWD : dont need now with TDLfft
  int z, y, x;
  float nxnynz = (_nx*_ny*_nz);
  for (z = 0; z < _nz; z++)
    for (y = 0; y < _ny; y++)
      for (x = 0; x < _nx; x++) {
	// added by kem to scale fft result correctly
	directionX[z][y][x] /= -nxnynz;
	directionY[z][y][x] /= -nxnynz;
	directionZ[z][y][x] /= -nxnynz;
      }
#endif

        // Originally we use to take inverse transform but
        // fftw dose not allow you to take inverse transform of real data
        // Hence we take forward transform. This should not make in differeance
        // Sarang Joshi 6/17/99

//   mfft.fftForward(directionX.data());
//   mfft.fftForward(directionY.data());
//   mfft.fftForward(directionZ.data());

   mfft.FFT(directionX,TDLfft::Forward);
   mfft.FFT(directionY,TDLfft::Forward);
   mfft.FFT(directionZ,TDLfft::Forward);

   int ii, i, j, k;

   float Lmatrix00, Lmatrix01, Lmatrix02;
   float Lmatrix10, Lmatrix11, Lmatrix12;
   float Lmatrix20, Lmatrix21, Lmatrix22;
   float LL00, LL01, LL02;
   float LL10, LL11, LL12;
   float LL20, LL21, LL22;
   float muReal0, muReal1, muReal2;
   float muImag0, muImag1, muImag2;
   float LLmuReal0, LLmuReal1, LLmuReal2;
   float LLmuImag0, LLmuImag1, LLmuImag2;

   for (k=0; k<_nz; k++) {
     for (int j=0; j<_ny; j++) {
       for (int i=0; i<_xmax; i++) {
	 // To use small number of Basis functions defined by _numBasis
	 // Date 29/06/98 Sarang Kem.
	 if ((i < _numBasis ) && 
	     (((j < _numBasis) || (j > _ny - _numBasis)) && 
	      ((k < _numBasis) || (k > _nz-_numBasis)))) {
	   double lambda = - _alpha * (_cosW1[i] + _cosW2[j] + _cosW3[k]) + _gamma;

	   // old stuff, just compute L
	   Lmatrix00 = lambda - _beta * _cosW1[i];
	   Lmatrix11 = lambda - _beta * _cosW2[j];
	   Lmatrix22 = lambda - _beta * _cosW3[k];
	   Lmatrix01 = _beta * _sinW1[i] * _sinW2[j];
	   Lmatrix02 = _beta * _sinW1[i] * _sinW3[k];
	   Lmatrix12 = _beta * _sinW2[j] * _sinW3[k];
	   Lmatrix10 = Lmatrix01; // This matrix is symmetric
	   Lmatrix20 = Lmatrix02;
	   Lmatrix21 = Lmatrix12;

	   // compute L'L
	   LL00 = Lmatrix00*Lmatrix00 + Lmatrix10*Lmatrix10 + Lmatrix20*Lmatrix20;
	   LL01 = Lmatrix00*Lmatrix01 + Lmatrix10*Lmatrix11 + Lmatrix20*Lmatrix21;
	   LL02 = Lmatrix00*Lmatrix02 + Lmatrix10*Lmatrix12 + Lmatrix20*Lmatrix22;
	   LL10 = LL01;
	   LL11 = Lmatrix01*Lmatrix01 + Lmatrix11*Lmatrix11 + Lmatrix21*Lmatrix21;
	   LL12 = Lmatrix01*Lmatrix02 + Lmatrix11*Lmatrix12 + Lmatrix21*Lmatrix22;
	   LL20 = LL02;
	   LL21 = LL12;
	   LL22 = Lmatrix02*Lmatrix02 + Lmatrix12*Lmatrix12 + Lmatrix22*Lmatrix22;

	   // compute b + 2LLmu
	   // have direction = deltaMu
	   muReal0 = muX[k][j][2*i];
	   muImag0 = muX[k][j][2*i+1];
	   muReal1 = muY[k][j][2*i];
	   muImag1 = muY[k][j][2*i+1];
	   muReal2 = muZ[k][j][2*i];
	   muImag2 = muZ[k][j][2*i+1];

	   LLmuReal0 = LL00*muReal0 + LL01*muReal1 + LL02*muReal2;
	   LLmuReal1 = LL10*muReal0 + LL11*muReal1 + LL12*muReal2;
	   LLmuReal2 = LL20*muReal0 + LL21*muReal1 + LL22*muReal2;
	   LLmuImag0 = LL00*muImag0 + LL01*muImag1 + LL02*muImag2;
	   LLmuImag1 = LL10*muImag0 + LL11*muImag1 + LL12*muImag2;
	   LLmuImag2 = LL20*muImag0 + LL21*muImag1 + LL22*muImag2;

	   directionX[k][j][2*i] += 2.0*priorWeight*LLmuReal0;
	   directionX[k][j][2*i + 1] += 2.0*priorWeight*LLmuImag0;

	   directionY[k][j][2*i] += 2.0*priorWeight*LLmuReal1;
	   directionY[k][j][2*i + 1] += 2.0*priorWeight*LLmuImag1;
	   
	   directionZ[k][j][2*i] += 2.0*priorWeight*LLmuReal2;
	   directionZ[k][j][2*i + 1] += 2.0*priorWeight*LLmuImag2;

	 } else {
	   directionX[k][j][2*i] = 0;
	   directionX[k][j][2*i + 1] = 0;
	   
	   directionY[k][j][2*i] = 0;
	   directionY[k][j][2*i + 1] = 0;
	   
	   directionZ[k][j][2*i] = 0;
	   directionZ[k][j][2*i + 1] = 0;
	 }
       }
     }
   }
}



float FluidTransformation::findOptimalDelta(Array3D<float> &muX, 
					    Array3D<float> &muY, 
					    Array3D<float> &muZ, 
					    Array3D<float> &directX, 
					    Array3D<float> &directY, 
					    Array3D<float> &directZ,
					    Array3D<float> &uX, 
					    Array3D<float> &uY, 
					    Array3D<float> &uZ, 
					    Array3D<float> &hField,
					    const Array3D<unsigned char> &atlas, 
					    Array3D<unsigned char> &defAtlas, 
					    const Array3D<unsigned char> &patient,
					    float priorWeight, 
					    float initialDelta,
					    double &lastPosterior) {

  // First, we must bracket the minimum  
  float deltaA, deltaB, deltaC;
  double posteriorA, posteriorB, posteriorC;
  double priorTerm, squaredErrorTerm;
  
  const double GOLDEN = 0.38197;
  
  deltaA = 0.0;
  posteriorA = lastPosterior;

  posteriorA = computePosterior(muX, muY, muZ, deltaA, directX, directY, directZ,
				uX, uY, uZ, hField,
				atlas, defAtlas, patient,
				priorWeight, priorTerm, squaredErrorTerm);
	
  deltaB = initialDelta;
  // This procedure changes uX, uY, uZ, hField, and defAtlas.
  // The posterior equals priorTerm + likelihoodTerm.
  posteriorB = computePosterior(muX, muY, muZ, deltaB, directX, directY, directZ,
				uX, uY, uZ, hField,
				atlas, defAtlas, patient,
				priorWeight, priorTerm, squaredErrorTerm);
  
  if (posteriorB < posteriorA) { // step outwards until f(c) > f(b)
    deltaC = deltaA + (deltaB-deltaA)/GOLDEN;
    posteriorC = computePosterior(muX, muY, muZ, deltaC, directX, directY, directZ,
				  uX, uY, uZ, hField,
				  atlas, defAtlas, patient,
				  priorWeight, priorTerm, squaredErrorTerm);
    while (posteriorC <= posteriorB) {
#if 0
      // shift a and b over
      deltaA = deltaB; deltaB = deltaC;
      posteriorA = posteriorB; posteriorB = posteriorC;
      deltaC = deltaA + (deltaB-deltaA)/GOLDEN;
      posteriorC = computePosterior(muX, muY, muZ, deltaC, directX, directY, directZ,
				    uX, uY, uZ, hField,
				    atlas, defAtlas, patient,
				    priorWeight, priorTerm, squaredErrorTerm);
#else
      // shift b over
      deltaB = deltaC;
      posteriorB = posteriorC;
      deltaC = deltaA + (deltaB-deltaA)/GOLDEN;
      posteriorC = computePosterior(muX, muY, muZ, deltaC, directX, directY, directZ,
				    uX, uY, uZ, hField,
				    atlas, defAtlas, patient,
				    priorWeight, priorTerm, squaredErrorTerm);
#endif
    }
    
  } else { // f(b) >= f(a) so min is in between
    deltaC = deltaB;
    posteriorC = posteriorB;
    deltaB = deltaA + GOLDEN*(deltaC-deltaA);
    posteriorB = computePosterior(muX, muY, muZ, deltaB, directX, directY, directZ,
				  uX, uY, uZ, hField,
				  atlas, defAtlas, patient,
				  priorWeight, priorTerm, squaredErrorTerm);
    // while (posteriorB > posteriorA) {
    while ((posteriorB-posteriorA) > 1.0e-8*posteriorA) {
      deltaC = deltaB;
      posteriorC = posteriorB;
      deltaB = deltaA + GOLDEN*(deltaC-deltaA);
      posteriorB = computePosterior(muX, muY, muZ, deltaB, directX, directY, directZ,
				    uX, uY, uZ, hField,
				    atlas, defAtlas, patient,
				    priorWeight, priorTerm, squaredErrorTerm);
    }
    
  }


  // The code below does line minimization.

#if 0
  // we now have a bracketing of a < b < c such that f(b) < f(a) and f(b) < f(c)
  // now do a simple golden section search
  float sigma;
  double fsigma;
  //  while (fabs((deltaC-deltaA)/deltaA) > 1.0e-1) {
  while ((fabs((posteriorC-posteriorB)/posteriorB) > 1.0e-3) ||
	 (fabs((posteriorA-posteriorB)/posteriorB) > 1.0e-3)) {
    if (fabs(deltaC-deltaB) > fabs(deltaB-deltaA)) {
      sigma = deltaB + GOLDEN*(deltaC-deltaB);
      fsigma = computePosterior(muX, muY, muZ, sigma, directX, directY, directZ,
				uX, uY, uZ, hField,
				atlas, defAtlas, patient,
				priorWeight, priorTerm, squaredErrorTerm);
      if (fsigma < posteriorB) {
	deltaA = deltaB; deltaB = sigma;
	posteriorA = posteriorB; posteriorB = fsigma;
      } else {
	deltaC = sigma;
	posteriorC = fsigma;
      }
    } else {
      sigma = deltaB - GOLDEN*(deltaB-deltaA);
      fsigma = computePosterior(muX, muY, muZ, sigma, directX, directY, directZ,
				uX, uY, uZ, hField,
				atlas, defAtlas, patient,
				priorWeight, priorTerm, squaredErrorTerm);
      if (fsigma < posteriorB) {
	deltaC = deltaB; deltaB = sigma;
	posteriorC = posteriorB; posteriorB = fsigma;
      } else {
	deltaA = sigma;
	posteriorA = fsigma;
      }
    }

  }
#endif

  lastPosterior = posteriorB;
  return deltaB;
}






	  // This procedure changes uX, uY, uZ, newHField, and _deformed.
	  // The posterior equals priorTerm + likelihoodTerm.
double FluidTransformation::computePosterior(Array3D<float> &muX, 
					     Array3D<float> &muY, 
					     Array3D<float> &muZ, 
					     float delta, 
					     Array3D<float> &directX, 
					     Array3D<float> &directY, 
					     Array3D<float> &directZ,
					     Array3D<float> &uX, 
					     Array3D<float> &uY, 
					     Array3D<float> &uZ, 
					     Array3D<float> &newHField,
					     const Array3D<unsigned char> &atlas, 
					     Array3D<unsigned char> &defAtlas, 
					     const Array3D<unsigned char> &patient,
					     float priorWeight, 
					     double &priorTerm, 
					     double &squaredErrorTerm)
{
  const unsigned char * maskPtr;
  if ( _maskPtr != NULL ) {
    maskPtr = _maskPtr->data();
  }

  // _Vx, _Vy, and _Vz now hold the gradient direction vector
  genDisplacement(muX, muY, muZ, priorWeight,
		  delta, directX, directY, directZ, uX, uY, uZ, priorTerm);
  
  // uX, uY, and uZ now hold the displacement field
  updateGlobalDisp(1.0, uX, uY, uZ, newHField);
	  
  // deform the atlas using the hField
  TransformationUtils::interpolate(atlas, defAtlas, newHField);
  
  const unsigned char * deformedAtlasPtr = defAtlas.data(); 
  const unsigned char * patientPtr = patient.data();
  squaredErrorTerm = 0.0;
  for (int z = 0; z < _nz; z++) {
    for (int y = 0; y < _ny; y++) {
      for (int x = 0; x < _nx; x++) {
	// Calculate (T(x - ut) - S(x)).
	double deformedMinusPatient = (*deformedAtlasPtr++ - *patientPtr++); 
	// compute the squared error
	if ( _maskPtr != NULL ) { // If Mask is set.
	  if ((*maskPtr++) != 0  ) {
	    deformedMinusPatient = 0.0;
	  }
        }
	squaredErrorTerm += deformedMinusPatient * deformedMinusPatient; 

      }
    }
  }
  cout << "Posterior: " << delta << " " << squaredErrorTerm+priorTerm 
       << " = " << squaredErrorTerm << " + " << priorTerm << endl;
  
  return (priorTerm + squaredErrorTerm);
}


/*********************************************************************************
*
* Given the body force (driving function), calculate the 
* velocity field. The velocity field is stored as the real
* arrays Vx, Vy, Vz.
*
**********************************************************************************/
void
FluidTransformation::genDisplacement(Array3D<float> &muX,
				     Array3D<float> &muY,
				     Array3D<float> &muZ,
				     float priorWeight,
				     float delta,
				     Array3D<float> &directionX,
				     Array3D<float> &directionY,
				     Array3D<float> &directionZ,
				     Array3D<float> &uX,
				     Array3D<float> &uY,
				     Array3D<float> &uZ,
				     double &priorTerm) {
  float *** directionXPtr = directionX.address();
  float *** directionYPtr = directionY.address();
  float *** directionZPtr = directionZ.address();
  float *** uXPtr = uX.address();
  float *** uYPtr = uY.address();
  float *** uZPtr = uZ.address();
  float *** muXPtr = muX.address();
  float *** muYPtr = muY.address();
  float *** muZPtr = muZ.address();
  float Lmatrix00, Lmatrix01, Lmatrix02;
  float Lmatrix10, Lmatrix11, Lmatrix12;
  float Lmatrix20, Lmatrix21, Lmatrix22;
  float LuX, LuY, LuZ;
  int z, y, x;
  float uXReal, uYReal, uZReal;
  float uXImag, uYImag, uZImag;

  priorTerm = 0.0;
  for (z = 0; z < _nz; z++) {
    for (y = 0; y < _ny; y++) {
      for (x = 0; x < _xmax; x++) {
	double lambda = - _alpha * (_cosW1[x] + _cosW2[y] + _cosW3[z]) + _gamma;

	// old stuff, just compute L
	Lmatrix00 = lambda - _beta * _cosW1[x];
	Lmatrix11 = lambda - _beta * _cosW2[y];
	Lmatrix22 = lambda - _beta * _cosW3[z];
	Lmatrix01 = _beta * _sinW1[x] * _sinW2[y];
	Lmatrix02 = _beta * _sinW1[x] * _sinW3[z];
	Lmatrix12 = _beta * _sinW2[y] * _sinW3[z];
	Lmatrix10 = Lmatrix01; // This matrix is symmetric
	Lmatrix20 = Lmatrix02;
	Lmatrix21 = Lmatrix12;

	// Real part of mu
	uXReal = muXPtr[z][y][2*x] - delta * directionXPtr[z][y][2*x];
	uYReal = muYPtr[z][y][2*x] - delta * directionYPtr[z][y][2*x];
	uZReal = muZPtr[z][y][2*x] - delta * directionZPtr[z][y][2*x];

	uXPtr[z][y][2*x] = uXReal;
	uYPtr[z][y][2*x] = uYReal;
	uZPtr[z][y][2*x] = uZReal;

        LuX = Lmatrix00 * uXReal + Lmatrix01 * uYReal + Lmatrix02 * uZReal;
        LuY = Lmatrix10 * uXReal + Lmatrix11 * uYReal + Lmatrix12 * uZReal;
        LuZ = Lmatrix20 * uXReal + Lmatrix21 * uYReal + Lmatrix22 * uZReal;

	priorTerm += (double) (LuX*LuX + LuY*LuY + LuZ*LuZ);

	// Imag part of mu
	uXImag = muXPtr[z][y][2*x+1] - delta * directionXPtr[z][y][2*x+1];
	uYImag = muYPtr[z][y][2*x+1] - delta * directionYPtr[z][y][2*x+1];
	uZImag = muZPtr[z][y][2*x+1] - delta * directionZPtr[z][y][2*x+1];

	uXPtr[z][y][2*x+1] = uXImag;
	uYPtr[z][y][2*x+1] = uYImag;
	uZPtr[z][y][2*x+1] = uZImag;

        LuX = Lmatrix00 * uXImag + Lmatrix01 * uYImag + Lmatrix02 * uZImag;
        LuY = Lmatrix10 * uXImag + Lmatrix11 * uYImag + Lmatrix12 * uZImag;
        LuZ = Lmatrix20 * uXImag + Lmatrix21 * uYImag + Lmatrix22 * uZImag;

	priorTerm += (double) (LuX*LuX + LuY*LuY + LuZ*LuZ);

      }
    }
  }
  priorTerm *= (double) priorWeight;


   // Perform inverse fft to get velocity field

//   mfft.fftInverse((FFT_complex *)uX.data());
//   mfft.fftInverse((FFT_complex *)uY.data());
//   mfft.fftInverse((FFT_complex *)uZ.data());

   mfft.FFT(uX,TDLfft::Reverse);
   mfft.FFT(uY,TDLfft::Reverse);
   mfft.FFT(uZ,TDLfft::Reverse);
}

/*********************************************************************************
*
* Given the body force (driving function), calculate the 
* velocity field. The velocity field is stored as the real
* arrays Vx, Vy, Vz.
*
**********************************************************************************/
void
FluidTransformation::updateMus(Array3D<float> &muX,
			       Array3D<float> &muY,
			       Array3D<float> &muZ,
			       float delta,
			       Array3D<float> &directionX,
			       Array3D<float> &directionY,
			       Array3D<float> &directionZ,
			       Array3D<float> &uX,
			       Array3D<float> &uY,
			       Array3D<float> &uZ,
			       Array3D<float> &hField) {

  float * directionXPtr = directionX.data();
  float * directionYPtr = directionY.data();
  float * directionZPtr = directionZ.data();
  float * muXPtr = muX.data();
  float * muYPtr = muY.data();
  float * muZPtr = muZ.data();
  float * uXPtr = uX.data();
  float * uYPtr = uY.data();
  float * uZPtr = uZ.data();
  int z, y, x;

  for (z = 0; z < _nz; z++)
    for (y = 0; y < _ny; y++)
      for (x = 0; x < _nx; x++) {
	muX[z][y][x] -= (delta * directionX[z][y][x]);
	muY[z][y][x] -= (delta * directionY[z][y][x]);
	muZ[z][y][x] -= (delta * directionZ[z][y][x]);

	uX[z][y][x] = muX[z][y][x];
	uY[z][y][x] = muY[z][y][x];
	uZ[z][y][x] = muZ[z][y][x];
      }

  // Perform inverse fft to get velocity field
//   mfft.fftInverse((FFT_complex *)uX.data());
//   mfft.fftInverse((FFT_complex *)uY.data());
//   mfft.fftInverse((FFT_complex *)uZ.data());

  mfft.FFT(uX,TDLfft::Reverse);
  mfft.FFT(uY,TDLfft::Reverse);
  mfft.FFT(uZ,TDLfft::Reverse);

  // uX, uY, and uZ now hold the displacement field
  updateGlobalDisp(1.0, uX, uY, uZ, hField);
}

double
FluidTransformation::MaxNormOfVectorField(const float *a,
					  const float *b,
					  const float *c,
					  int numX,
					  int numY,
					  int numZ,
					  int strideX,
					  double *minimumPtr,
					  double *averagePtr )
{
   const float *XPtr = a;
   const float *YPtr = b;
   const float *ZPtr = c;

   int x, y, z; 
   double maximum = sqrt( *XPtr * *XPtr + *YPtr * *YPtr + *ZPtr * *ZPtr );
   int vectorFieldMaxIndexX = 0;
   int vectorFieldMaxIndexY = 0;
   int vectorFieldMaxIndexZ = 0;
   double minimum = sqrt( *XPtr * *XPtr + *YPtr * *YPtr + *ZPtr * *ZPtr );
   double average = 0.0;
   double tmp;
   for (z = 0; z < numZ; z++)
   {
     for (y = 0; y < numY; y++)
     {
       for (x = 0; x < numX; x++)
       {
         tmp = sqrt( *XPtr * *XPtr + *YPtr * *YPtr + *ZPtr * *ZPtr );

         if (tmp > maximum)
         {
           maximum = tmp; 
           vectorFieldMaxIndexX = x;
           vectorFieldMaxIndexY = y;
           vectorFieldMaxIndexZ = z;
         }
         if (tmp < minimum)
         {
           minimum = tmp; 
         }
         average += tmp;

         XPtr++;
         YPtr++;
         ZPtr++;
       }
       XPtr += strideX; 
       YPtr += strideX; 
       ZPtr += strideX; 
     }
   }

   *minimumPtr = minimum;

   *averagePtr = average / ((double)numZ*numY*numX);

  cout << "Max velocity of " << maximum 
       << " occured at (x, y, z) " << "(" << vectorFieldMaxIndexX;
  cout << ", " << vectorFieldMaxIndexY;
  cout << ", " << vectorFieldMaxIndexZ << ")" << endl;

   return maximum;
}




// Takes velocity field, displacement field, atlas, patient and delta
// and returns squared error if this delta and velocity field were
// used.

// Calculates best delta
double
FluidTransformation::FindBestDelta(int *StopFlagPtr)
{
  double MaxPertInPixels = _mxper;

  double MaxOfVelocity     = 0;
  double MinOfVelocity     = 0;
  double AverageOfVelocity = 0;

  static double OptimalDelta      = 0;


  static int FirstTimeThroughFlag = 1;


  MaxOfVelocity = MaxNormOfVectorField(_Vx.data(),
                                       _Vy.data(),
                                       _Vz.data(),
                                       _nx,
                                       _ny,
                                       _nz,
                                       _xstride,
                                       &MinOfVelocity,
                                       &AverageOfVelocity );

  if (0 == MaxOfVelocity)
  {
    cout << "Maximum of velocity is zero." << endl;
    cout << "Setting delta to zero and stopping" << endl;
    OptimalDelta = 0.0;
    *StopFlagPtr = 1;
    return OptimalDelta;
  }


  if (FirstTimeThroughFlag) 
  {
    FirstTimeThroughFlag = 0;
    OptimalDelta = MaxPertInPixels / MaxOfVelocity;
  }
  
  double maxFrobeniusNorm = calcMaxFrobeniusNorm(_Vx, _Vy, _Vz);
  cout << "Max Frobenius Norm: " << maxFrobeniusNorm << endl;

  // choose optimal delta using the maximum Frobenius norm
  OptimalDelta = _mxper / maxFrobeniusNorm;

  calcJacobian(-OptimalDelta, _Vx, _Vy, _Vz);


  if ( _CalculateTotalJacobianFlag ) 
  {
    cout << "Updating Jacobian starting" << endl;
    updateJacobian(-OptimalDelta, _Vx, _Vy, _Vz);
    cout << "Min of incremental jacobian is " << mIncrementalJacobianMin << " which occured at (";
    cout << mIncrementalJacobianMinIndexX << ", ";
    cout << mIncrementalJacobianMinIndexY << ", ";
    cout << mIncrementalJacobianMinIndexZ << ")" << endl;
    cout << "Min of total jacobian is " << mTotalJacobianMin << " which occured at (";
    cout << mTotalJacobianMinIndexX << ", ";
    cout << mTotalJacobianMinIndexY << ", ";
    cout << mTotalJacobianMinIndexZ << ")" << endl;
    
    cout << "Updating Jacobian done" << endl;
  }


  if ( MaxOfVelocity*OptimalDelta < _TimMagicStoppingNumber )
  {
    cout << "Maximum displacement (" << MaxOfVelocity*OptimalDelta;
    cout << ") is less than " << _TimMagicStoppingNumber << " so stopping" << endl;
    *StopFlagPtr = 1;
  } else {
    *StopFlagPtr = 0;	// Need to take it out.
  }

  cout << "OptimalDelta: " << OptimalDelta << endl;
  cout << "Displacement (min, max, avg): " 
       << MinOfVelocity*OptimalDelta << ", "
       << MaxOfVelocity*OptimalDelta << ", "
       << AverageOfVelocity*OptimalDelta << endl;

  return OptimalDelta;  
}


/*****************************************************************************************
*
* Update total displacement field.
*
*****************************************************************************************/
void FluidTransformation::updateGlobalDisp(const double &delta,
					   const Array3D<float> &velocityX,
					   const Array3D<float> &velocityY,
					   const Array3D<float> &velocityZ,
					   Array3D<float> &newHField)
{
//  const float *vxptr = velocityX.data(); 
//  const float *vyptr = velocityY.data(); 
//  const float *vzptr = velocityZ.data(); 
  
//  float *newHFieldPtr = newHField.data();
  
  float interpolatedArray[3];
  int x, y, z; 
  for (z = 0; z < _nz; z++) {
    for (y = 0; y < _ny; y++) {
      for (x = 0; x < _nx; x++) {
	newHField[z][y][3*x+0] = x - (-delta * velocityX[z][y][x]);
	newHField[z][y][3*x+1] = y - (-delta * velocityY[z][y][x]);
	newHField[z][y][3*x+2] = z - (-delta * velocityZ[z][y][x]);

//	*newHFieldPtr = x - (-delta * *vxptr++);
//	newHFieldPtr++;
//	*newHFieldPtr = y - (-delta * *vyptr++);
//	newHFieldPtr++;
//	*newHFieldPtr = z - (-delta * *vzptr++);
//	newHFieldPtr++;
      }

 //     vxptr += _xstride; 
 //     vyptr += _xstride; 
 //     vzptr += _xstride; 
      
    }
  }

}

/*****************************************************************************************
*
* Compose transformation fields
*
* takes in an old h-field and an update to the field, returns
* the composed field over the new field.
*
*****************************************************************************************/
void FluidTransformation::composeTransField(const Array3D<float> &oldHField,
					    Array3D<float> &newHField)
{
  float interpolatedArray[3];
  int x, y, z; 
  for (z = 0; z < _nz; z++) {
    for (y = 0; y < _ny; y++) {
      for (x = 0; x < _nx; x++) {
	
	// Deform oldHField by the incremental h-field and 
	// put the result in interpolatedArray, Then store it in the
	// newHFieldPtr so we don't have to allocate more memory
        // (yes its supposed to be z,y,x as last 3 params)
	TransformationUtils::interpolate(interpolatedArray,
					 oldHField,
					 newHField[z][y][3*x+2],
					 newHField[z][y][3*x+1],
					 newHField[z][y][3*x+0]);

	 newHField[z][y][3*x+0] = interpolatedArray[0];
	 newHField[z][y][3*x+1] = interpolatedArray[1];
	 newHField[z][y][3*x+2] = interpolatedArray[2];
      }
    }
  }

}


double FluidTransformation::calcMaxFrobeniusNorm( Array3D<float> const &fieldx, 
						  Array3D<float> const &fieldy,
						  Array3D<float> const &fieldz)
{

  Array1D<double> dfieldx(3), dfieldy(3), dfieldz(3); 
  float tempValue = 0;

  int nz = fieldx.getZsize(); 
  int ny = fieldx.getYsize(); 
  int nx = fieldx.getXsize(); 

  const float * const * const *fieldXPtrPtrPtr = fieldx.address();
  const float * const * const *fieldYPtrPtrPtr = fieldy.address();
  const float * const * const *fieldZPtrPtrPtr = fieldz.address();

  float dx0, dx1, dx2;
  float dy0, dy1, dy2;
  float dz0, dz1, dz2;
  float maxNorm = 0.0;
  float minNorm = 1.0e30;
  float avgNorm = 0.0;
#if 0
  int x, y, z; 
  for (z = 0; z < nz; z++) {
    for (y = 0; y < ny; y++) {
      for (x = 0; x < nx; x++) {

#if 1

#if 0
#else
	if (x > 0 && x < (nx-1)) {
	  dx0 = (fieldXPtrPtrPtr[z][y][x+1] - fieldXPtrPtrPtr[z][y][x-1])/2.0;
	  dy0 = (fieldYPtrPtrPtr[z][y][x+1] - fieldYPtrPtrPtr[z][y][x-1])/2.0;
	  dz0 = (fieldZPtrPtrPtr[z][y][x+1] - fieldZPtrPtrPtr[z][y][x-1])/2.0;
	} else if (x == 0) {
	  dx0 = (fieldXPtrPtrPtr[z][y][x+1] - fieldXPtrPtrPtr[z][y][x]);
	  dy0 = (fieldYPtrPtrPtr[z][y][x+1] - fieldYPtrPtrPtr[z][y][x]);
	  dz0 = (fieldZPtrPtrPtr[z][y][x+1] - fieldZPtrPtrPtr[z][y][x]);
	} else {
	  dx0 = (fieldXPtrPtrPtr[z][y][x] - fieldXPtrPtrPtr[z][y][x-1]);
	  dy0 = (fieldYPtrPtrPtr[z][y][x] - fieldYPtrPtrPtr[z][y][x-1]);
	  dz0 = (fieldZPtrPtrPtr[z][y][x] - fieldZPtrPtrPtr[z][y][x-1]);
	}
	if (y > 0 && y < (ny-1)) {
	  dx1 = (fieldXPtrPtrPtr[z][y+1][x] - fieldXPtrPtrPtr[z][y-1][x])/2.0;
	  dy1 = (fieldYPtrPtrPtr[z][y+1][x] - fieldYPtrPtrPtr[z][y-1][x])/2.0;
	  dz1 = (fieldZPtrPtrPtr[z][y+1][x] - fieldZPtrPtrPtr[z][y-1][x])/2.0;
	} else if (y == 0) {
	  dx1 = (fieldXPtrPtrPtr[z][y+1][x] - fieldXPtrPtrPtr[z][y][x]);
	  dy1 = (fieldYPtrPtrPtr[z][y+1][x] - fieldYPtrPtrPtr[z][y][x]);
	  dz1 = (fieldZPtrPtrPtr[z][y+1][x] - fieldZPtrPtrPtr[z][y][x]);
	} else {
	  dx1 = (fieldXPtrPtrPtr[z][y][x] - fieldXPtrPtrPtr[z][y-1][x]);
	  dy1 = (fieldYPtrPtrPtr[z][y][x] - fieldYPtrPtrPtr[z][y-1][x]);
	  dz1 = (fieldZPtrPtrPtr[z][y][x] - fieldZPtrPtrPtr[z][y-1][x]);
	}
	if (z > 0 && z < (nz-1)) {
	  dx2 = (fieldXPtrPtrPtr[z+1][y][x] - fieldXPtrPtrPtr[z-1][y][x])/2.0;
	  dy2 = (fieldYPtrPtrPtr[z+1][y][x] - fieldYPtrPtrPtr[z-1][y][x])/2.0;
	  dz2 = (fieldZPtrPtrPtr[z+1][y][x] - fieldZPtrPtrPtr[z-1][y][x])/2.0;
	} else if (z == 0) {
	  dx2 = (fieldXPtrPtrPtr[z+1][y][x] - fieldXPtrPtrPtr[z][y][x]);
	  dy2 = (fieldYPtrPtrPtr[z+1][y][x] - fieldYPtrPtrPtr[z][y][x]);
	  dz2 = (fieldZPtrPtrPtr[z+1][y][x] - fieldZPtrPtrPtr[z][y][x]);
	} else {
	  dx2 = (fieldXPtrPtrPtr[z][y][x] - fieldXPtrPtrPtr[z-1][y][x]);
	  dy2 = (fieldYPtrPtrPtr[z][y][x] - fieldYPtrPtrPtr[z-1][y][x]);
	  dz2 = (fieldZPtrPtrPtr[z][y][x] - fieldZPtrPtrPtr[z-1][y][x]);
	}

#endif
	tempValue = (dx0*dx0 + dx1*dx1 + dx2*dx2 +
		     dy0*dy0 + dy1*dy1 + dy2*dy2 +
		     dz0*dz0 + dz1*dz1 + dz2*dz2);
#else
	TransformationUtils::gradient(&dfieldx, fieldx, x, y, z, 1.0, 1.0, 1.0); 
	TransformationUtils::gradient(&dfieldy, fieldy, x, y, z, 1.0, 1.0, 1.0); 
	TransformationUtils::gradient(&dfieldz, fieldz, x, y, z, 1.0, 1.0, 1.0); 
	tempValue = (dfieldx[0]*dfieldx[0] + dfieldx[1]*dfieldx[1] + dfieldx[2]*dfieldx[2] +
		     dfieldy[0]*dfieldy[0] + dfieldy[1]*dfieldy[1] + dfieldy[2]*dfieldy[2] + 
		     dfieldz[0]*dfieldz[0] + dfieldz[1]*dfieldz[1] + dfieldz[2]*dfieldz[2] );
#endif	
	if (tempValue > maxNorm)
	  maxNorm = tempValue;
      }
    }
  }
#else
  for (int z=1; z<(nz-1); z++) {
    for (int y=1; y<(ny-1); y++) {
      for (int x=1; x<(nx-1); x++) {
	dx0 = (fieldXPtrPtrPtr[z][y][x+1] - fieldXPtrPtrPtr[z][y][x-1]);
	dy0 = (fieldYPtrPtrPtr[z][y][x+1] - fieldYPtrPtrPtr[z][y][x-1]);
	dz0 = (fieldZPtrPtrPtr[z][y][x+1] - fieldZPtrPtrPtr[z][y][x-1]);
	dx1 = (fieldXPtrPtrPtr[z][y+1][x] - fieldXPtrPtrPtr[z][y-1][x]);
	dy1 = (fieldYPtrPtrPtr[z][y+1][x] - fieldYPtrPtrPtr[z][y-1][x]);
	dz1 = (fieldZPtrPtrPtr[z][y+1][x] - fieldZPtrPtrPtr[z][y-1][x]);
	dx2 = (fieldXPtrPtrPtr[z+1][y][x] - fieldXPtrPtrPtr[z-1][y][x]);
	dy2 = (fieldYPtrPtrPtr[z+1][y][x] - fieldYPtrPtrPtr[z-1][y][x]);
	dz2 = (fieldZPtrPtrPtr[z+1][y][x] - fieldZPtrPtrPtr[z-1][y][x]);
	
	tempValue = sqrt(dx0*dx0 + dx1*dx1 + dx2*dx2 +
			 dy0*dy0 + dy1*dy1 + dy2*dy2 +
			 dz0*dz0 + dz1*dz1 + dz2*dz2);
	
	if (tempValue > maxNorm)
	  maxNorm = tempValue;
	if (tempValue < minNorm)
	  minNorm = tempValue;
	avgNorm += tempValue;

      }
    }
  }
  maxNorm /= 2.0;
  minNorm /= 2.0;
  avgNorm /= (2.0*(nx-2)*(ny-2)*(nz-2));
  
  cout << "Frobenius norm (min, max, avg): "
       << minNorm << ", "
       << maxNorm << ", "
       << avgNorm << endl;
#endif
  
  return maxNorm; // divided by 4 to match old results
}


/**********************************************************************************
*
* Minimum Jacobian in JAC[0], Maximum Jacobian in JAC[1].
* Average Jacobian in JAC[2].
*
***********************************************************************************/
double FluidTransformation::calcJacobian(double delta, Array3D<float> const &Ux, 
	Array3D<float> const &Uy, Array3D<float> const &Uz) {

  float j; 
  
  int nz = Ux.getZsize(); 
  int ny = Ux.getYsize(); 
  int nx = Ux.getXsize(); 
  
  const float * const * const *fieldXPtrPtrPtr = Ux.address();
  const float * const * const *fieldYPtrPtrPtr = Uy.address();
  const float * const * const *fieldZPtrPtrPtr = Uz.address();

  // changing these to floats seem to slow it down!?
  double dx0, dx1, dx2;
  double dy0, dy1, dy2;
  double dz0, dz1, dz2;

  double grad00, grad01, grad02;
  double grad10, grad11, grad12;
  double grad20, grad21, grad22;

  double minJacobian = 1.0;
  double maxJacobian = 1.0;
  double avgJacobian = 0.0;

  // Note only the internal elements of the field are being used (i.e. going
  // from 1..(n-2) instead of 0..(n-1)) since the boundary shouldn't matter 
  // much anyway.  This speeds up the loop since we don't have to test for 
  // boundary conditions.

  for (int z=1; z < (nz-1); z++) {
    for (int y=1; y < (ny-1); y++) {
      for (int x=1; x < (nx-1); x++) {

	dx0 = (fieldXPtrPtrPtr[z][y][x+1] - fieldXPtrPtrPtr[z][y][x-1])/2;
	dy0 = (fieldYPtrPtrPtr[z][y][x+1] - fieldYPtrPtrPtr[z][y][x-1])/2;
	dz0 = (fieldZPtrPtrPtr[z][y][x+1] - fieldZPtrPtrPtr[z][y][x-1])/2;
	dx1 = (fieldXPtrPtrPtr[z][y+1][x] - fieldXPtrPtrPtr[z][y-1][x])/2;
	dy1 = (fieldYPtrPtrPtr[z][y+1][x] - fieldYPtrPtrPtr[z][y-1][x])/2;
	dz1 = (fieldZPtrPtrPtr[z][y+1][x] - fieldZPtrPtrPtr[z][y-1][x])/2;
	dx2 = (fieldXPtrPtrPtr[z+1][y][x] - fieldXPtrPtrPtr[z-1][y][x])/2;
	dy2 = (fieldYPtrPtrPtr[z+1][y][x] - fieldYPtrPtrPtr[z-1][y][x])/2;
	dz2 = (fieldZPtrPtrPtr[z+1][y][x] - fieldZPtrPtrPtr[z-1][y][x])/2;

	// grad is equal to I-delta*dU
	grad00 = 1.0-delta*dx0;
	grad01 = -delta*dx1;
	grad02 = -delta*dx2;
	grad10 = -delta*dy0;
	grad11 = 1.0-delta*dy1;
	grad12 = -delta*dy2;
	grad20 = -delta*dz0;
	grad21 = -delta*dz1;
	grad22 = 1.0-delta*dz2;

	j = (grad00*(grad11*grad22 - grad12*grad21)
	     - grad01*(grad10*grad22 - grad12*grad20)
	     + grad02*(grad10*grad21 - grad11*grad20));
	
	// THE OLD WRONG WAY?!?!?
	// j = ((1 - dUx[0]) * (1 - dUy[1]) * (1 - dUz[2])) 
	// + (dUx[1] * dUy[2] * dUz[0]) + (dUx[2] * dUy[0] * dUz[1]) 
	// - (dUx[2] * (1 - dUy[1]) * dUz[0]) 
	// - (dUx[1] * dUy[0] * (1 - dUz[2])) 
	// - ((1 - dUx[0]) * dUy[2] * dUz[1]); 
	
	avgJacobian += j;
	if (j < minJacobian)
	  minJacobian = j;
	if (j > maxJacobian)
	  maxJacobian = j;
      }
    }
  }

  // Average Jacobian.
  avgJacobian /= ((nx-2)*(ny-2)*(nz-2));
  
  cout.setf(ios::fixed, ios::floatfield);
  cout << "Jacobian (min, max, avg): " 
       << minJacobian << ", "
       << maxJacobian << ", "
       << avgJacobian << " ~= 1.0" <<endl;
  cout.setf(ios::scientific, ios::floatfield);

  return minJacobian; 
}


int
FluidTransformation::setNumBasis(int numbasis)
{ 
	if (numbasis < 1 )
  {
    cerr << "FluidTransformation::setNumBasis: Number of Basis must be positive ";
    cerr << " (" << numbasis << ")" << endl;
    return 1;
  }
  else
  {
	_numBasis = numbasis;
  }
  return 0;
}

int
FluidTransformation::setNumberIterations( int numOuter, int numInner )
{ 
  if (numOuter < 1 || numInner < 1)
  {
    cerr << "FluidTransformation::setNumberIterations: Iterations must be positive ";
    cerr << " (" << numOuter << ", " << numInner << ")" << endl;
    return 1;
  }
  else
  {
    _niter = numOuter;
    _nInnerIter = numInner;
  }
  return 0;
}

int
FluidTransformation::setLaplacianWeight( double weight )
{ 
  if ( weight < 0.0)
  {
    cerr << "FluidTransformation::setLaplacianWeight: weight must be non-negative ";
    cerr << " (" << weight << ")" << endl;
    return 1;
  }
  else
  {
    _alpha = weight;
    return 0;
  }
}

int
FluidTransformation::setGradWeight( double weight )
{ // _beta = weight; 
  if ( weight < 0.0)
  {
    cerr << "FluidTransformation::setGradWeight: weight must be non-negative ";
    cerr << " (" << weight << ")" << endl;
    return 1;
  }
  else
  {
    _beta = weight;
    return 0;
  }
}

int
FluidTransformation::setAdditiveWeight( double weight )
{ // _gamma = weight; 
  if ( weight <= 0.0)
  {
    cerr << "FluidTransformation::setAdditiveWeight: weight must be positive and non-zero ";
    cerr << " (" << weight << ")" << endl;
    return 1;
  }
  else
  {
    _gamma = weight;
    return 0;
  }
}

int
FluidTransformation::setMaxPerturbation( double pmax )
{ // _mxper = pmax; 
  if ( pmax <= 0.0)
  {
    cerr << "FluidTransformation::setMaxPertubation: max pert must be > zero ";
    cerr << " (" << pmax << ")" << endl;
    return 1;
  }
  else
  {
    _mxper = pmax; 
    return 0;
  }
}

int
FluidTransformation::setMagicStoppingNumber( double x ) 
{ // _TimMagicStoppingNumber = x; 
  if ( x <= 0.0)
  {
    cerr << "FluidTransformation::setMagicStoppingNumber: stopping criteria  must be > zero ";
    cerr << " (" << x << ")" << endl;
    return 1;
  }
  else
  {
    _TimMagicStoppingNumber = x; 
    return 0;
  }
}

double
FluidTransformation::updateJacobian( double delta,
                                     Array3D<float> const &Vx, 
                                     Array3D<float> const &Vy,
                                     Array3D<float> const &Vz)
{
  // Jacobian.
  Array1D<double> Jacobian; 
  Jacobian.setDim(3); 
  
   // Array to hold min. and max. & average Jacobian.
   Jacobian[0] = Jacobian[1] = 1.0; Jacobian[2] = 0.0; 

// House keeping.
   double j; 
   Array1D<double> dVx(3), dVy(3), dVz(3); 

   int nz = mTotalJacobianScratch.getZsize(); 
   int ny = mTotalJacobianScratch.getYsize(); 
   int nx = mTotalJacobianScratch.getXsize(); 

   int x, y, z; 

   // Interpolate previous total jacobian and put it into scratch variable
   for (z = 0; z < nz; z++) {
      for (y = 0; y < ny; y++) {
	 for (x = 0; x < nx; x++) {
            // Now interpolate into the Total Jacobian image to pull out
            // the total jacobian so far for point delta*[Vx Vy Vz]
            mTotalJacobianScratch[z][y][x] = 
              TransformationUtils::trilinear( mTotalJacobian,
                                              (double) (z - delta * Vz[z][y][x]),
                                              (double) (y - delta * Vy[z][y][x]),
                                              (double) (x - delta * Vx[z][y][x]),
                                              (float) 1.0); 
	 }
      }
   }

   mIncrementalJacobianMinIndexX = 0;
   mIncrementalJacobianMinIndexY = 0;
   mIncrementalJacobianMinIndexZ = 0;
   mIncrementalJacobianMin = 1; // max value is one by definition
   mTotalJacobianMinIndexX = 0;
   mTotalJacobianMinIndexY = 0;
   mTotalJacobianMinIndexZ = 0;
   mTotalJacobianMin = 1; // max value is one by definition

   for (z = 0; z < nz; z++) {
      for (y = 0; y < ny; y++) {
	 for (x = 0; x < nx; x++) {

	    TransformationUtils::gradient(&dVx, Vx, x, y, z, 1.0, 1.0, 1.0); 
	    TransformationUtils::gradient(&dVy, Vy, x, y, z, 1.0, 1.0, 1.0); 
	    TransformationUtils::gradient(&dVz, Vz, x, y, z, 1.0, 1.0, 1.0); 

            dVx *= delta; 
            dVy *= delta;
            dVz *= delta;

	    j = ((1 - dVx[0]) * (1 - dVy[1]) * (1 - dVz[2])) 
	      + (dVx[1] * dVy[2] * dVz[0]) + (dVx[2] * dVy[0] * dVz[1]) 
	      - (dVx[2] * (1 - dVy[1]) * dVz[0]) - (dVx[1] * dVy[0] * (1 - dVz[2])) 
	      - ((1 - dVx[0]) * dVy[2] * dVz[1]); 

            // Accumulate Jacobian.
	    Jacobian[2] += j;  

            // Check against current min and max.
	    if (j < Jacobian[0]) Jacobian[0] = j; 			// minimum.
	    if (j > Jacobian[1]) Jacobian[1] = j; 			// maximum.


            if ( j < mIncrementalJacobianMin)
            {
              mIncrementalJacobianMin = j;
              mIncrementalJacobianMinIndexX = x;
              mIncrementalJacobianMinIndexY = y;
              mIncrementalJacobianMinIndexZ = z;
            }

            // Now update total Jacobian from interpolated value and multiply by j
            mTotalJacobian[z][y][x] = mTotalJacobianScratch[z][y][x] * j;

            if ( mTotalJacobian[z][y][x] < mTotalJacobianMin )
            {
              mTotalJacobianMin = mTotalJacobian[z][y][x];
              mTotalJacobianMinIndexX = x;
              mTotalJacobianMinIndexY = y;
              mTotalJacobianMinIndexZ = z;
            }
	 }
      }
   }

// Average Jacobian.
   Jacobian[2] /= (_nx * _ny * _nz); 

   return Jacobian[0]; 

}

//////////////////////////////////////////////////////////////////////////
//
// Method: FluidTransformation::setCalculateTotalJacobian( int flag )
//
// Purpose: if flag is non-zero then set up the arrays to calculate
//          the incremental jacobian.  if flag is zero, then delete
//          the memory the arrays.  Must be run after initialize()
//
// Date: [csb:19971210.0857MST]
// 
// Author: csb
// 
// Returns 0 upon success
// Returns 1 if initialize hasn't been run yet
// Returns 2 if can't allocate memory
//
//////////////////////////////////////////////////////////////////////////
int
FluidTransformation::setCalculateTotalJacobian( int flag )
{
  _CalculateTotalJacobianFlag = flag;
  
  if ( ! isOkToCalculate() )
  {
    cerr << "You must call FluidTransformation::initialize() before calling FluidTransformation::setCalculateTotalJacobian( int flag )" << endl;
    return 1;
  }
  
  if ( _CalculateTotalJacobianFlag )
  {
    mTotalJacobian.setDim(_nz, _ny, _nx);
    if (NULL == mTotalJacobian.data())
    {
      cerr << "FluidTransformation::setCalculateTotalJacobian: Failed to allocate memory" << endl;
      return 2;
    }
    mTotalJacobian = (float) 1.0;
    
    mTotalJacobianScratch.setDim(_nz, _ny, _nx);
    if (NULL == mTotalJacobianScratch.data())
    {
      cerr << "FluidTransformation::setCalculateTotalJacobian: Failed to allocate memory" << endl;
      return 2;
    }
    mTotalJacobianScratch  = (float) 1.0;
    
    mTotalJacobianMin = 1;
    mTotalJacobianMinIndexX = 0;
    mTotalJacobianMinIndexY = 0;
    mTotalJacobianMinIndexZ = 0;
    mIncrementalJacobianMin = 1;
    mIncrementalJacobianMinIndexX = 0;
    mIncrementalJacobianMinIndexY = 0;
    mIncrementalJacobianMinIndexZ = 0;
  }
  else
  {
    // Deallocate memory
    mTotalJacobian.setDim(0, 0, 0);
    mTotalJacobianScratch.setDim(0, 0, 0);
  }
  
  return 0;
}


////////////////////////////////////////////////////////////////////////////////////
//
// Returns 0 upon success
// Returns 1 if can't load in file from "filename"
// Returns 2 if can't set f-field
//
////////////////////////////////////////////////////////////////////////////////////
int
FluidTransformation::setInitialTransformation( const char * filename)
{

  if (0 == strcmp("", filename))
  {
    cout << "initial transformation filename is \"" << filename << "\".";
    cout << "Assuming this means no initial transformation" << endl;
    return 0;
  }
  
  char className[1024];
 
  int readrc = TransformationUtils::readTransformationType(filename, className);
  if (0 != readrc)
  {
    cerr << "FluidTransformation::setInitialTransformation failed. Returned " << readrc << endl;
    return 2;
  };

  // Switch here on the type of transformation
  FieldTransformation * fieldTransPtr;
  
  if (0 == strcmp(className, "FluidTransformation"))
  {
    FluidTransformation * specificfieldTransPtr = new FluidTransformation;
    fieldTransPtr = specificfieldTransPtr;
  }
  else if (0 == strcmp(className, "ManifoldTransformation"))
  {
    ManifoldTransformation * specificfieldTransPtr = new ManifoldTransformation;
    fieldTransPtr = specificfieldTransPtr;
  }
  else 
  {
    cerr << "Transformation type \"" << className << "\" not yet recognized" << endl;
    return 1;
  }
  
  cout << "Loading file \"" << filename << "\" of type \"" << className << "\"" << endl;
  int loadRC = fieldTransPtr->load(filename);
  if (0 != loadRC)
  {
    cerr << "Failed to load in transformation file \"" << filename << "\"" << endl;
    cerr << "FluidTransformation::load() returned " << loadRC << endl;
    delete fieldTransPtr;
    return (10 + loadRC);
  }

  // setUField is a FieldTransformation method
  cout << "Setting h-field ...";
  int setRC = fieldTransPtr->setHField(_TotalHField);
  if (0 != setRC)
  {
    cerr << "Failed to set the h-field to be that of file \"" << filename << "\"" << endl;
    cerr << "FluidTransformation::setHField() returned " << setRC << endl;
    delete fieldTransPtr;
    return (20 + setRC);
  }
  cout << " done" << endl;

  delete fieldTransPtr;
  cout << "Freeing fieldTransPtr" << endl;
  return 0;
}


int
FluidTransformation::dumpthem( Array3D<float> & xfield,
			       Array3D<float> & yfield,
			       Array3D<float> & zfield,
			       const char * xfile,
			       const char * yfile,
			       const char * zfile)
{
  
  FILE *fpx, *fpy, *fpz;
  fpx = fopen(xfile, "w");
  fpy = fopen(yfile, "w");
  fpz = fopen(zfile, "w");
  
  if (fpx == NULL || fpy == NULL || fpz == NULL)
    {
      cerr << "Cannot open file " << xfile << " or " << yfile << " or " << zfile << endl;
      return 1;
    }
  
  float * vxptr = xfield.data(); 
  float * vyptr = yfield.data(); 
  float * vzptr = zfield.data(); 
  float swapfloat;
  for (int z = 0; z < _nz; z++)
    {
      for (int y = 0; y < _ny; y++)
	{
	  for (int x = 0; x < _nx; x++)
	  {
	    swapfloat = big_endian_float(*vxptr);
	    fwrite(&swapfloat, 1, sizeof(float), fpx);
	    vxptr++;
	    swapfloat = big_endian_float(*vyptr);
	    fwrite(&swapfloat, 1, sizeof(float), fpy);
	    vyptr++;
	    swapfloat = big_endian_float(*vzptr);
	    fwrite(&swapfloat, 1, sizeof(float), fpz);
	    vzptr++;
	  }
	  
	  vxptr += _xstride; 
	  vyptr += _xstride; 
	  vzptr += _xstride; 
	}
    }
  fclose(fpx);
  fclose(fpy);
  fclose(fpz);
  
  return 0;
}
      


////////////////////////////////////////////////////////////////////////////////
//
// FILE HANDLING METHODS
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
int FluidTransformation::save(const char *filename) {

  ofstream outStream (filename, ios::out);
  
  if ( ! outStream ) {
    cerr << "Failed to open file \"" << filename << "\"" << endl;
    return 1;
  }
  if ( saveClassParameters( outStream ) ) {
    return 2;
  };
  outStream.close();

  if (getVerboseMode())
    cout << "save() in FluidTransformation" << endl;

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
int FluidTransformation::load(const char *filename) {

  ifstream inStream (filename, ios::in);
  if ( ! inStream ) {
    cerr << "Failed to open file \"" << filename << "\"" << endl;
    return 1;
  }

  if ( loadClassParameters( inStream ) ) {
    return 2;
  }

  if (getVerboseMode())
      cout << "load() in FluidTransformation" << endl;
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
void FluidTransformation::print() {

  if (getVerboseMode())
    cout << "print() in FluidTransformation" << endl;
  FluidTransformation::printClassParameters(cout);
  
}


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
int FluidTransformation::saveClassParameters ( ofstream & outStream) {

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

  /************************************
   * Write out parameters for this class.. Rev. Swap if neccesary
   ************************************/
  double swapdouble;
  swapint = big_endian_int(_niter);
  outStream.write((char const *) &swapint, sizeof(_niter));
  swapint = big_endian_int(_nInnerIter);
  outStream.write((char const *) &swapint, sizeof(_nInnerIter));
  swapdouble = big_endian_double(_alpha);
  outStream.write((char const *) &swapdouble, sizeof(_alpha));
  swapdouble = big_endian_double(_beta);
  outStream.write((char const *) &swapdouble,  sizeof(_beta));
  swapdouble = big_endian_double(_gamma);
  outStream.write((char const *) &swapdouble, sizeof(_gamma));
  swapdouble = big_endian_double(_delta);
  outStream.write((char const *) &swapdouble, sizeof(_delta));
  swapdouble = big_endian_double(_mxper);
  outStream.write((char const *) &swapdouble, sizeof(_mxper));
  swapdouble = big_endian_double(_mnjack);
  outStream.write((char const *) &swapdouble, sizeof(_mnjack));
  swapdouble = big_endian_double(_TimMagicStoppingNumber);
  outStream.write((char const *) &swapdouble, sizeof(_TimMagicStoppingNumber));
  swapdouble = big_endian_double(_CalculateTotalJacobianFlag);
  outStream.write((char const *) &swapdouble, sizeof(_CalculateTotalJacobianFlag));

  /****************************
   * Write out the vector field.
   ****************************/
  _TotalHField.save(outStream); 

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
int FluidTransformation::loadClassParameters (ifstream & inStream) {

  // Function name.
  const char * const functionName = "FluidTransformation::loadClassParameters";

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
  if (ReadInClassVersion == 3)
  {
    cout << "Reading old transformation: " << mClassName << " version "
	 << ReadInClassVersion << endl;
  }
  else if (ReadInClassVersion != mClassVersion)
  {
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

  inStream.read((char *) &_niter, sizeof(_niter));
  _niter = big_endian_int(_niter); //swapped it!
  if (ReadInClassVersion < 4)
    {
      //      cout << "fudging _nInnerIter for old version" << endl;
      _nInnerIter = 1;
    }
  else
    {
    inStream.read((char *) &_nInnerIter, sizeof(_nInnerIter));
    _nInnerIter = big_endian_int(_nInnerIter); //swapped it!
    }
  inStream.read((char *) &_alpha, sizeof(_alpha));
  inStream.read((char *) &_beta,  sizeof(_beta));
  inStream.read((char *) &_gamma, sizeof(_gamma));
  inStream.read((char *) &_delta, sizeof(_delta));
  inStream.read((char *) &_mxper, sizeof(_mxper));
  inStream.read((char *) &_mnjack, sizeof(_mnjack));
  inStream.read((char *) &_TimMagicStoppingNumber, sizeof(_TimMagicStoppingNumber));
  inStream.read((char *) &_CalculateTotalJacobianFlag, sizeof(_CalculateTotalJacobianFlag));
  
  //do byte swapping if neccasary
  _alpha = (big_endian_double(_alpha));
  _beta  = (big_endian_double(_beta));
  _gamma = (big_endian_double(_gamma));
  _delta = (big_endian_double(_delta));
  _mxper = (big_endian_double(_mxper));
  _mnjack= (big_endian_double (_mnjack));
  _TimMagicStoppingNumber = (int)(big_endian_double(_TimMagicStoppingNumber));
  _CalculateTotalJacobianFlag=(int)(big_endian_double(_CalculateTotalJacobianFlag));
  
/**************************
   * Read in the vector field.
   **************************/

  _TotalHField.load(inStream); 
 
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
ostream & FluidTransformation::printClassParameters (ostream &outStream) {

  FieldTransformation::printClassParameters (outStream);
  outStream << mClassName << "::mClassVersion:              " << mClassVersion << endl;
  outStream << mClassName << "::mClassRevision:             " << mClassRevision << endl;
  outStream << mClassName << "::niter:                      " << _niter << endl;
  outStream << mClassName << "::nInnerIter:                 " << _nInnerIter << endl;
  outStream << mClassName << "::alpha:                      " << _alpha << endl;
  outStream << mClassName << "::beta:                       " << _beta << endl;
  outStream << mClassName << "::gamma:                      " << _gamma << endl;
  outStream << mClassName << "::delta:                      " << _delta << endl;
  outStream << mClassName << "::mxper:                      " << _mxper << endl;
  outStream << mClassName << "::mnjack:                     " << _mnjack << endl;
  outStream << mClassName << "::MagicStopingNumber:         " << _TimMagicStoppingNumber << endl;
  outStream << mClassName << "::CalculateTotalJacobianFlag: " << _CalculateTotalJacobianFlag << endl;

  return outStream;

}
