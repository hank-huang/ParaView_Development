#ifndef __FluidTRANSFORMATION_H__
#define __FluidTRANSFORMATION_H__

//////////////////////////////////////////////////////////////////////////
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
// $Log: FluidTransformation.h,v $
// Revision 1.1  2004/11/15 04:44:08  joeh
// first check in of BrainWorks
//
// Revision 1.2  1999/11/09 21:28:48  kem
// Add methods for computing h field values given a point
//
// Revision 1.1  1999/10/01 15:28:12  kem
// Initial revision
//
// Revision 1.23  1999/07/09 17:49:17  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.22  1999/05/14 18:44:32  RAZORDB
// Reverted from version 1.20
//
// Revision 1.20  1999/03/16 23:37:06  kem
// Fix resolution passing
//
// Revision 1.19  1999/01/06 16:46:16  RAZORDB
// Add FluidFluidTransformation as a friend
//
// Revision 1.18  1998/12/18 18:14:28  RAZORDB
// Fix for 7.2.1 compiler Warning 1356
//
// Revision 1.17  1998/12/10 19:03:36  RAZORDB
// IDL/AW merge
//
// Revision 1.16  1998/11/03 18:08:41  RAZORDB
// Change mask to const array
//
// Revision 1.15  1998/10/30 22:42:15  RAZORDB
// Add Masking of infiltarting Tumors and prostrate
//
// Revision 1.14  1998/10/29 00:19:23  RAZORDB
// Add genDeltaMu() method for elastic algorithm
//
// Revision 1.13  1998/07/15 21:19:30  kem
// Timing and optimization
//
// Revision 1.12  1998/07/10 19:18:24  kem
// Clean up unused methods and variables
//
// Revision 1.11  1998/07/10 15:56:49  sarang
// To combine Elastic and Fluid Transformation
//
// Revision 1.10  1998/06/02 15:24:27  rst
// Changed fluid to support elastic
//
// Revision 1.9  1998/04/09 19:08:42  kem
// Balance Lv(x) = b(v(x))
//
// Revision 1.8  1998/02/13 22:11:10  csb
// Implemented different eigenvectors
//
// Revision 1.7  1997/12/12 22:03:35  csb
// merge 1.5.1.2 and 1.6
//
// Revision 1.6  1997/12/12 19:24:29  abed
// Change getHfield Prototype
//
// Revision 1.5  1997/09/16 19:47:39  kem
// Remove getUField() method
//
// Revision 1.4  1997/08/25 17:40:47  kem
// Add getHfield()
//
// Revision 1.3  1997/08/05 18:39:51  rst
// removed ElasticManifoldTransformation friend
//
// Revision 1.2  1997/08/05 18:35:54  rst
// make ElasticManifoldTransformation friend
//
// Revision 1.1  1997/08/01 19:49:21  csb
// Initial revision
//
// Revision 1.1  1997/08/01 15:49:23  csb
// Initial revision
//
//
//////////////////////////////////////////////////////////////////////////
#include <OS.h>
#include <stream.h>
#include <ADL/Array1D.h>
#include <TDL/ITX_FFT.h>
#include <TDL/TDLfft.h>
#include <TDL/FieldTransformation.h>

class ManifoldFluidTransformation;

class FluidTransformation : public FieldTransformation {

  friend class ManifoldFluidTransformation;
  friend class FluidFluidTransformation;

public:

  // Ctor and Dtor.
  FluidTransformation(); 
  ~FluidTransformation(); 		      // virtual in Transformation class
    
  int initialize( const Array3D<unsigned char> * atlasPtr, 
                  const Array3D<unsigned char> * patientPtr,
                  const Array3D<unsigned char> * maskPtr = NULL,
		  float aresx=1.0, float aresy=1.0, float aresz=1.0,
		  float presx=1.0, float presy=1.0, float presz=1.0); // not-virtual
  int calculate();             			// virtual in Transformation class


  int setNumberIterations( int numOuter, int numInner = 1 );
  int setNumBasis(int numbasis = 1 );
  int setLaplacianWeight( double weight );
  int setGradWeight( double weight );
  int setAdditiveWeight( double weight );
  int setMaxPerturbation( double pmax ); 
  int setMagicStoppingNumber(double x); 
  int setCalculateTotalJacobian( int flag );
  int setInitialTransformation( const char * filename );

private:

  // Generate vector field.
  void getHField(Array3D<float> *);
  void getHField(Array3D<float> *, 
		 double Rx, double Ry, double Rz, 
		 double Ox, double Oy, double Oz); 
  void getHField(const float patientInVoxelsX,
		 const float patientInVoxelsY,
		 const float patientInVoxelsZ,
		 float &atlasInVoxelsX,
		 float &atlasInVoxelsY,
		 float &atlasInVoxelsZ);
  

  // Initialize basis.
  void genTables(); 
  
  // Generate body force .
  void genBodyForce(const Array3D<float> &hField,
		    const Array3D<unsigned char> &atlas,
		    const Array3D<float> &AtlasGrad,
		    const Array3D<unsigned char> &patient,
		    Array3D<unsigned char> &deformedAtlas,
		    Array3D<float> &bodyForceX,
		    Array3D<float> &bodyForceY,
		    Array3D<float> &bodyForceZ,
		    double &squaredError); 

  void compareBodyForce(const Array3D<float> &bodyForceX,
			const Array3D<float> &bodyForceY,
			const Array3D<float> &bodyForceZ,
			Array3D<float> &lastBodyForceX,
			Array3D<float> &lastBodyForceY,
			Array3D<float> &lastBodyForceZ,
			float &meanSquaredErrorX,
			float &meanSquaredErrorY,
			float &meanSquaredErrorZ);

  // Velocity Field.
  void genVelocityField(Array3D<float> &bodyToVelocityX,
			Array3D<float> &bodyToVelocityY,
			Array3D<float> &bodyToVelocityZ); 
  
  // generate u field for elastic algorithm
  void genDirection(Array3D<float> &muX,
		    Array3D<float> &muY,
		    Array3D<float> &muZ,
		    float priorWeight,
		    Array3D<float> &directionX,
		    Array3D<float> &directionY,
		    Array3D<float> &directionZ);

  float findOptimalDelta(Array3D<float> &muX, 
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
			 double &lastPosterior);

  double computePosterior(Array3D<float> &muX, 
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
			  double &squaredErrorTerm);

  // generate displacement field given last mu, direction vector and delta
  void genDisplacement(Array3D<float> &muX,
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
		       double &priorTerm);

  // update mu's given delta and direction vector
  void updateMus(Array3D<float> &muX,
		 Array3D<float> &muY,
		 Array3D<float> &muZ,
		 float delta,
		 Array3D<float> &directionX,
		 Array3D<float> &directionY,
		 Array3D<float> &directionZ,
		 Array3D<float> &uX,
		 Array3D<float> &uY,
		 Array3D<float> &uZ,
		 Array3D<float> &hField);
  

  // Perturbations.
  double FindBestDelta(int *); 
  
  // Update total displacement field.
  void updateGlobalDisp(const double &delta,
			const Array3D<float> &velocityX,
			const Array3D<float> &velocityY,
			const Array3D<float> &velocityZ,
			Array3D<float> &newHField);

  // Compose transformation fields
  void composeTransField(const Array3D<float> &oldHField,
			 Array3D<float> &newHField);

  double MaxNormOfVectorField(const float *a,
			      const float *b,
			      const float *c,
			      int numX,
			      int numY,
			      int numZ,
			      int strideX,
			      double *minimumPtr,
			      double *averagePtr );
  
  double calcMaxFrobeniusNorm(Array3D<float> const &,
			      Array3D<float> const &, 
			      Array3D<float> const &); 
  // Check Jacobian.
  double calcJacobian(double, Array3D<float> const &, Array3D<float> const &, 
		      Array3D<float> const &); 
  
  double updateJacobian(double,
			Array3D<float> const &,
			Array3D<float> const &, 
			Array3D<float> const &); 
  
  // sine transform functions
  void sine1DTransformInitCoeff(float sineCoeff[], int n);
  void sine1DTransform(float data[], int n, float fftCoeff[], float sineCoeff[]);
  void sine3DTransform(float data[], int nx, int ny, int nz, int stride);

  int dumpthem( Array3D<float> & xfield,
		Array3D<float> & yfield,
		Array3D<float> & zfield,
		const char * xfile,
		const char * yfile,
		const char * zfile);
  
  // FILE HANDLING METHODS
public:
  int save( const char * filename );          // virtual in Transformation class
  int load( const char * filename );          // virtual in Transformation class
  void print();

protected:
  // Returns 0 upon success, and non-zero otherwise
  int saveClassParameters ( ofstream & );
  int loadClassParameters ( ifstream & );
  ostream & printClassParameters ( ostream & );



private:
  static const char * const mClassRevision;
  static const char * const mClassName;
  static const int          mClassVersion;
  
  int _niter;         // number of iterations
  int _nInnerIter;    // number of inner iterations (elastic iterations)

  int _numBasis;      // number of basis

  // Operator coefficients
  double _alpha;      // Laplacian coefficient
  double _beta;       // grad(div u) coefficient
  double _gamma;      // additive constant
  
  double _delta;      // integration delta
  
  double _mxper;      // max perturbation
  
  double _mnjack;     // min Jacobian
  
  // Dimensions of input volumes.
  int _nz, _ny, _nx, _nxPlusStride, _xstride; 
  
  // Index of basis functions.
  int _xmax, _ymax, _zmax; 
  
  // Velocity vector field.
  Array3D<float> _Vz, _Vy, _Vx; 
  
  Array3D<float> _TotalHField;


  // Variables for Inverse Method of Velocity coefficient calculation

  // deltas for the difference operator
  double _Delta1;
  double _Delta2;
  double _Delta3;
  
  // finite difference vectors
  Array1D<float> _cosW1;
  Array1D<float> _cosW2;
  Array1D<float> _cosW3;
  Array1D<float> _sinW1;
  Array1D<float> _sinW2;
  Array1D<float> _sinW3;
    
  // FFT coefficients.
   Array1D<float> _coeff; 
  
  // Image Volumes.
  Array3D<unsigned char> _deformed; 
  
  // Pointers to atlas and patient data
  const Array3D<unsigned char> * _PatientPtr; 
  const Array3D<unsigned char> * _AtlasPtr; 
  // Pointers to Mask volume
  const Array3D<unsigned char> * _maskPtr;

  
  // Tim's Magic stoping number :)
  double _TimMagicStoppingNumber; 
  
  // flag that says to calculate the total jacobian or not
  int _CalculateTotalJacobianFlag;
  
  Array3D<float> mTotalJacobian;
  Array3D<float> mTotalJacobianScratch;
  float mTotalJacobianMin;
  int mTotalJacobianMinIndexX;
  int mTotalJacobianMinIndexY;
  int mTotalJacobianMinIndexZ;
  float mIncrementalJacobianMin;
  int mIncrementalJacobianMinIndexX;
  int mIncrementalJacobianMinIndexY;
  int mIncrementalJacobianMinIndexZ;

  // FFT Stuff

//  ITX_FFT mfft;
  TDLfft mfft;
};

#endif // __FluidTRANSFORMATION_H__
