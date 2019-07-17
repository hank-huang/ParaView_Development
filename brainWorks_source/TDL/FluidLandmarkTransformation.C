///////////////////////////////////////////////////////////////////////////
//
// IntellX, L.L.C., Proprietary
// Do not reproduce without permission in writing.
// Copyright (c) 1998 IntellX, L.L.C.
// All rights reserved
//
// File: FluidLandmarkTransformation.C
//
// FluidLandmarkTransformation class
//
// Author: Keith Doolittle
//	   11/98
//
///////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <math.h>
#include <TDL/FluidLandmarkTransformation.h>
#include <TDL/AffineLandmarkTransformation.h>
#include <TDL/Point.h>
#include <TDL/Timer.h>

const char * const FluidLandmarkTransformation::mClassRevision = "$Id: FluidLandmarkTransformation.C,v 1.1 2004/11/15 04:44:08 joeh Exp $";
const char * const FluidLandmarkTransformation::mClassName = "FluidLandmarkTransformation";
const int FluidLandmarkTransformation::mClassVersion = 0;

//
// LAPACK routines
//
extern "C" {
extern void dsytrf_(char *UPLO, int *N, double *A, int *LDA,
		    int *IPIV, double *work, int *LWORK, int *INFO);

extern void dsytrs_(char *UPLO, int *N, int *NRHS, double *A,
		    int *LDA, int *IPIV, double *B, int *LDB, int *INFO);
};


FluidLandmarkTransformation::FluidLandmarkTransformation() :
	_alpha(0.01),
	_sigma(0.1),
	_beta(0.0002),
	_numt(10),
	_numiter(200)
{
//
}

FluidLandmarkTransformation::~FluidLandmarkTransformation() 
{
// empty
}

int FluidLandmarkTransformation::initialize(
	const Matrix<double> &atlasPnts,
	const Matrix<double> &patientPnts,
        int AtlasSizeX,  int AtlasSizeY,  int AtlasSizeZ,
	double AtlasScaleX, double AtlasScaleY, double AtlasScaleZ,
        int PatientSizeX,  int PatientSizeY,  int PatientSizeZ,
	double PatientScaleX, double PatientScaleY, double PatientScaleZ
	)


{
  if((AtlasScaleX <= 0.0)||(AtlasScaleY <= 0.0)||(AtlasScaleZ <= 0.0)||
     (PatientScaleX <= 0.0)||(PatientScaleY <= 0.0)||(PatientScaleZ <= 0.0)) {
	cerr << 
	  "ERROR: FluidLandmarkTransformation :: invalid pixel sizes" << endl;
	return 1;
	}

  // error checking
  if ((atlasPnts.getNrow() != patientPnts.getNrow()) ||
      (atlasPnts.getNcol() != patientPnts.getNcol())) {
    return 1;
    }

  if ((atlasPnts.getNcol() != 3) ||
      (patientPnts.getNcol() != 3)) {
    return 2;
    }

  if ((atlasPnts.getNrow() < 4) ||
      (patientPnts.getNrow() < 4)) {
    return 3;
    }

  _atlasPoints   = atlasPnts;
  _patientPoints = patientPnts;
  _numPoints     = _atlasPoints.getNrow();
  _nx            = PatientSizeX;
  _ny            = PatientSizeY;
  _nz            = PatientSizeZ;

  _patientSizeX  = PatientSizeX;
  _patientSizeY  = PatientSizeY;
  _patientSizeZ  = PatientSizeZ;
  _atlasSizeX    = AtlasSizeX;
  _atlasSizeY    = AtlasSizeY;
  _atlasSizeZ    = AtlasSizeZ;

  _patientResX   = PatientScaleX;
  _patientResY   = PatientScaleY;
  _patientResZ   = PatientScaleZ;

  _atlasResX     = AtlasScaleX;
  _atlasResY     = AtlasScaleY;
  _atlasResZ     = AtlasScaleZ;

  _TotalHField.setDim(_nz,_ny,_nx*3);

  if(_TotalHField.isEmpty())
	return 5;

//
// hfield is in mm
//
  for(int k=0;k<_nz;k++) 
   for(int j=0;j<_ny;j++)
    for(int i=0;i<_nx;i++) {
	_TotalHField[k][j][i*3+0] = (float)i * _patientResX;
	_TotalHField[k][j][i*3+1] = (float)j * _patientResY;
	_TotalHField[k][j][i*3+2] = (float)k * _patientResZ;
	}

  pathx.setDim(_numt,_numPoints);
  pathy.setDim(_numt,_numPoints);
  pathz.setDim(_numt,_numPoints);

  if(pathx.isEmpty() || pathy.isEmpty() || pathz.isEmpty())
	return 5;

  v1x.setDim(_numPoints);
  v1y.setDim(_numPoints);
  v1z.setDim(_numPoints);
  v2x.setDim(_numPoints);
  v2y.setDim(_numPoints);
  v2z.setDim(_numPoints);
  v3x.setDim(_numt,_numPoints);
  v3y.setDim(_numt,_numPoints);
  v3z.setDim(_numt, _numPoints);

  if(v1x.isEmpty() || v1y.isEmpty() || v1z.isEmpty() ||
     v2x.isEmpty() || v2y.isEmpty() || v2z.isEmpty() ||
     v3x.isEmpty() || v3y.isEmpty() || v3z.isEmpty())
	return 6;

  wtx.setDim(_numPoints);
  wty.setDim(_numPoints);
  wtz.setDim(_numPoints);

  K1.setDim(_numPoints,_numPoints);
  K2.setDim(_numPoints,_numPoints);
  mat.setDim(_numPoints,_numPoints);

  if(wtx.isEmpty() || wty.isEmpty() || wtz.isEmpty() ||
     K1.isEmpty() || K2.isEmpty() || mat.isEmpty()) 
	return 7;

  setOkToCalculate(1);

  return(0);
}




int FluidLandmarkTransformation::calculate()
{

ITXTimer timer1("find path");
ITXTimer timer2("apply affine");
ITXTimer timer3("calc weights");
ITXTimer timer4("update h field");
ITXTimer timer5("scale h field");


cout << "calculate() in FluidLandmarkTransformation" << endl;

cout.precision(7);
cerr.precision(7);

if ( ! isOkToCalculate() )
   return 1;


ShowWorking("Finding path\n");

cout << "Go to findPath()" << endl;

timer1.start();
if(findPath())
	return(1);
timer1.stop();

cout << "Back from findPath()" << endl;

//
// Now apply GLN transform
// (last affine result of findPath())
//

cout << "Applying Affine Transformation to hfield" << endl;

ShowWorking("Applying Affine Transformation to hfield\n");
int i,j,k,l;

Matrix<double> A(4,4),InvA(4,4);

_mAffine.getTransformMatrix(&A);
InvA = MatrixUtils<double>::inv(A);

Vector<double> ipt(4),opt(4);

timer2.start();
cout << "Apply affine : " ;
float *hFieldPtr = _TotalHField.data();
for(k=0;k<_nz;k++)  {
  for(j=0;j<_ny;j++)
    for(i=0;i<_nx;i++) {
      ipt[0] =  hFieldPtr[0];
      ipt[1] =  hFieldPtr[1];
      ipt[2] =  hFieldPtr[2];
      ipt[3] = 1.0;
      
      opt = InvA * ipt;
      
      *hFieldPtr++ = opt[0];
      *hFieldPtr++ = opt[1];
      *hFieldPtr++ = opt[2];
      /*      ipt[0] =  _TotalHField[k][j][3*i+0];
      ipt[1] =  _TotalHField[k][j][3*i+1];
      ipt[2] =  _TotalHField[k][j][3*i+2];
      ipt[3] = 1.0;
      
      opt = InvA * ipt;
      
      //	_mAffine.apply(ipt,&opt);
      
      _TotalHField[k][j][3*i+0] = opt[0];
      _TotalHField[k][j][3*i+1] = opt[1];
      _TotalHField[k][j][3*i+2] = opt[2];
      */
    }
  cout << k;
}
timer2.stop();
cout << endl;

int    step;
float  ux,uy,uz;
double tmp1,dx,dy,dz;

ShowWorking("Updating hfield\n");

float dwork_perc = 1.0/((float)(_numt-1));
float  work_perc = 0.0;

for(step=(int)_numt-2;step>=0;step--,work_perc += dwork_perc) {
  ShowWorking(work_perc);

  cout << "Step " << step << endl;

  timer3.start();
  calculateWeights(
		   pathx.data(step+1), pathy.data(step+1), pathz.data(step+1),
		   pathx.data(step), pathy.data(step), pathz.data(step));
  timer3.stop();

  timer4.start();
  
  float *hFieldPtr = _TotalHField.data();
  double *pathXPtr = pathx.data(step+1, 0);
  double *pathYPtr = pathy.data(step+1, 0);
  double *pathZPtr = pathz.data(step+1, 0);
  
  for(k=0;k<_nz;k++) 
    for(j=0;j<_ny;j++)
      for(i=0;i<_nx;i++) {
	ux = uy = uz = 0.0;
	float hFieldX = hFieldPtr[0];
	float hFieldY = hFieldPtr[1];
	float hFieldZ = hFieldPtr[2];

	for(l=0;l<_numPoints;l++) {
	  dx   = pathXPtr[l] - hFieldX;
	  dy   = pathYPtr[l] - hFieldY;
	  dz   = pathZPtr[l] - hFieldZ;
	  tmp1 = sqrt(dx*dx+dy*dy+dz*dz);
	  tmp1 = exp(-_alpha*tmp1);
	  ux   += (wtx[l]*tmp1);
	  uy   += (wty[l]*tmp1);
	  uz   += (wtz[l]*tmp1);
	}
	
	*hFieldPtr++ = hFieldX + ux;
	*hFieldPtr++ = hFieldY + uy;
	*hFieldPtr++ = hFieldZ + uz;

	/*	
	ux = uy = uz = 0.0;
	for(l=0;l<_numPoints;l++) {
	  dx   = pathx[step+1][l] - _TotalHField[k][j][3*i+0];
	  dy   = pathy[step+1][l] - _TotalHField[k][j][3*i+1];
	  dz   = pathz[step+1][l] - _TotalHField[k][j][3*i+2];
	  tmp1 = sqrt(dx*dx+dy*dy+dz*dz);
	  tmp1 = exp(-_alpha*tmp1);
	  ux   += (wtx[l]*tmp1);
	  uy   += (wty[l]*tmp1);
	  uz   += (wtz[l]*tmp1);
	}
	_TotalHField[k][j][3*i+0] += ux;
	_TotalHField[k][j][3*i+1] += uy;
	_TotalHField[k][j][3*i+2] += uz;
	*/
      }
  timer4.stop();
}

timer5.start();
	// Convert to atlas Voxels.
	for(k=0;k<_nz;k++)
          for(j=0;j<_ny;j++)
             for(i=0;i<_nx;i++) {
	_TotalHField[k][j][3*i+0]  /= _atlasResX;
	_TotalHField[k][j][3*i+1]  /= _atlasResY;
	_TotalHField[k][j][3*i+2]  /= _atlasResZ;
	}
timer5.stop();


setOkToApply(1);


ShowWorking("DONE.\n");
ShowWorking((char*)NULL);

timer1.displayResults();
timer2.displayResults();
timer3.displayResults();
timer4.displayResults();
timer5.displayResults();
return(0);
}



int FluidLandmarkTransformation::findPath()
{
int    i,j;
char   UPLO;
int    INFO,LDA,LWORK,N,NRHS;

Array1D<int>    IPIV(_numPoints);
Array1D<double> work(_numPoints);

Matrix<double> newP(_numPoints,3);
Matrix<double> path(_numPoints,3);

if(_mAffine.initialize(&_patientPoints,&_atlasPoints))
	return(1);


if(_mAffine.calculate())
	return(2);

if(_mAffine.apply(_patientPoints,newP))
	return(3);

for (i=0;i<_numPoints;i++)
  for(j=0;j<_numt;j++){
    pathx[j][i] = ((float)j/((float)_numt-1.0))*
	          (newP[i][0]-_atlasPoints[i][0])+_atlasPoints[i][0];
    pathy[j][i] = ((float)j/((float)_numt-1.0))*
	          (newP[i][1]-_atlasPoints[i][1])+_atlasPoints[i][1];
    pathz[j][i] = ((float)j/((float)_numt-1.0))*
	          (newP[i][2]-_atlasPoints[i][2])+_atlasPoints[i][2];
    }


double tmpx,tmpy,tmpz,tmp;
double dx,dy,dz;

float dwork = 1.0/((float)_numiter);
float work_perc = 0.0;

for(int iter=0;iter<_numiter;iter++,work_perc += dwork) {

cout << "FindPath iter " << iter << " np = " << _numPoints << endl;

	if(iter % 10 == 0)
		ShowWorking(work_perc);
	//
	// Generate affine
	//
	for(i=0;i<_numPoints;i++) {
		path[i][0] = pathx[_numt-1][i];
		path[i][1] = pathy[_numt-1][i];
		path[i][2] = pathz[_numt-1][i];
		}

	if(_mAffine.initialize(&_patientPoints,&path))
		return(4);
	if(_mAffine.calculate())
		return(5);

	if(_mAffine.apply(_patientPoints,newP))
		return(6);

        int t;
	for(t=1;t<_numt;t++) {
	   for(i=0;i<_numPoints;i++) {
	       v1x[i] = 0;
	       v1y[i] = 0;
	       v1z[i] = 0;
	       if(t < (_numt-1)) {
		   v1x[i] = pathx[t+1][i] - pathx[t][i];
		   v1y[i] = pathy[t+1][i] - pathy[t][i];
		   v1z[i] = pathz[t+1][i] - pathz[t][i];
		   }
	       v2x[i] = pathx[t][i] - pathx[t-1][i];
	       v2y[i] = pathy[t][i] - pathy[t-1][i];
	       v2z[i] = pathz[t][i] - pathz[t-1][i];
	       }

	   for(i=0;i<_numPoints;i++)
	     for(j=0;j<_numPoints;j++) {
		double dx = pathx[t][i] - pathx[t][j];
		double dy = pathy[t][i] - pathy[t][j];
		double dz = pathz[t][i] - pathz[t][j];
		double tmp = sqrt(dx*dx+dy*dy+dz*dz);

		K1[i][j] = exp(-_alpha * tmp);

		// *(K1+j+i*np) = exp(-_alpha*tmp);

		dx = pathx[t-1][i] - pathx[t-1][j];
		dy = pathy[t-1][i] - pathy[t-1][j];
		dz = pathz[t-1][i] - pathz[t-1][j];
		tmp = sqrt(dx*dx+dy*dy+dz*dz);

		K2[i][j] = exp(-_alpha * tmp);
		// *(K2+j+i*np) = exp(-_alpha*tmp);
		}

	   UPLO  = 'U';
	   INFO  = 1;
	   LDA   = _numPoints;
	   LWORK = _numPoints;
	   N     = _numPoints;
	   NRHS  = 1;
	
   dsytrf_(&UPLO,&N,K1.data(),&LDA,IPIV.data(),work.data(),&LWORK,&INFO);
   dsytrs_(&UPLO,&N,&NRHS,K1.data(),&LDA,IPIV.data(),v1x.data(),&LDA,&INFO);
   dsytrs_(&UPLO,&N,&NRHS,K1.data(),&LDA,IPIV.data(),v1y.data(),&LDA,&INFO);
   dsytrs_(&UPLO,&N,&NRHS,K1.data(),&LDA,IPIV.data(),v1z.data(),&LDA,&INFO);

	   UPLO  = 'U';
	   INFO  = 1;
	   LDA   = _numPoints;
	   LWORK = _numPoints;
	   N     = _numPoints;
	   NRHS  = 1;

   dsytrf_(&UPLO,&N,K2.data(),&LDA,IPIV.data(),work.data(),&LWORK,&INFO);
   dsytrs_(&UPLO,&N,&NRHS,K2.data(),&LDA,IPIV.data(),v2x.data(),&LDA,&INFO);
   dsytrs_(&UPLO,&N,&NRHS,K2.data(),&LDA,IPIV.data(),v2y.data(),&LDA,&INFO);
   dsytrs_(&UPLO,&N,&NRHS,K2.data(),&LDA,IPIV.data(),v2z.data(),&LDA,&INFO);

	   for(i=0;i<_numPoints;i++) {
		v3x[t][i] = -2*(v1x[i] - v2x[i]);
		v3y[t][i] = -2*(v1y[i] - v2y[i]);
		v3z[t][i] = -2*(v1z[i] - v2z[i]);
		}
	
	   for(i=0;i<_numPoints;i++) {
		tmpx = tmpy = tmpz = 0.0;
		for(j=0;j<_numPoints;j++) {
		    dx = pathx[t][i] - pathx[t][j];
		    dy = pathy[t][i] - pathy[t][j];
		    dz = pathz[t][i] - pathz[t][j];
		    tmp = sqrt(dx*dx+dy*dy+dz*dz);
		    if(i != j) {
			tmpx += (v1x[j]*(_alpha*dx/tmp)*exp(-_alpha*tmp));
			tmpy += (v1y[j]*(_alpha*dy/tmp)*exp(-_alpha*tmp));
			tmpz += (v1z[j]*(_alpha*dz/tmp)*exp(-_alpha*tmp));
			}
		    }

		v3x[t][i] += v1x[i]*tmpx;
		v3y[t][i] += v1y[i]*tmpy;
		v3z[t][i] += v1z[i]*tmpz;

		if(t == _numt-1) {
		    v3x[t][i] += (2.0*(pathx[t][i] - newP[i][0])/_sigma);
		    v3y[t][i] += (2.0*(pathy[t][i] - newP[i][1])/_sigma);
		    v3z[t][i] += (2.0*(pathz[t][i] - newP[i][2])/_sigma);
		    }
		}
	   } // for t = 

	for(t=1;t<_numt;t++)
	  for(i=0;i<_numPoints;i++) {
		pathx[t][i] -= (_beta*v3x[t][i]);
		pathy[t][i] -= (_beta*v3y[t][i]);
		pathz[t][i] -= (_beta*v3z[t][i]);
		}
	} // for iter = 	

return 0;
}





void FluidLandmarkTransformation::calculateWeights(
		double *x1, double *y1, double *z1,
		double *x2, double *y2, double *z2)
{
int i,j;
double tmp,dx,dy,dz;

for(i=0;i<_numPoints;i++)
  for(j=0;j<_numPoints;j++) {
	dx = x1[i] - x1[j];	
	dy = y1[i] - y1[j];	
	dz = z1[i] - z1[j];	
	tmp = sqrt(dx*dx+dy*dy+dz*dz);

	mat[i][j] = exp(-_alpha * tmp);

	// *(mat+j+i*np) = exp(-_alpha*tmp);
	}

for(i=0;i<_numPoints;i++) {
	wtx[i] = x2[i] - x1[i];
	wty[i] = y2[i] - y1[i];
	wtz[i] = z2[i] - z1[i];
	}

char UPLO  = 'U';
int INFO   = 1;
int LDA    = _numPoints;
int LWORK  = 8;
int N      = _numPoints;
int NRHS   = 1;

Array1D<int>    IPIV(_numPoints);
Array1D<double> work(_numPoints);

dsytrf_(&UPLO,&N,mat.data(),&LDA,IPIV.data(),work.data(),&LWORK,&INFO);
dsytrs_(&UPLO,&N,&NRHS,mat.data(),&LDA,IPIV.data(),wtx.data(),&LDA,&INFO);
dsytrs_(&UPLO,&N,&NRHS,mat.data(),&LDA,IPIV.data(),wty.data(),&LDA,&INFO);
dsytrs_(&UPLO,&N,&NRHS,mat.data(),&LDA,IPIV.data(),wtz.data(),&LDA,&INFO);
}

void FluidLandmarkTransformation::getHField(Array3D<float> *hField)
{
*hField = _TotalHField;
}

void FluidLandmarkTransformation::getHField(Array3D<float> *hField,
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


int FluidLandmarkTransformation::saveClassParameters ( ofstream & outStream)
{
  int swapint;
  double swapdouble;
  // Write out size of class name and the class name
  const int nameLength = strlen(mClassName);
  swapint = big_endian_int(nameLength);
  outStream.write((const char *) & swapint, sizeof(nameLength));
  outStream.write(mClassName,                  nameLength);

   // Write out the class version
  const int temp  = mClassVersion;
  swapint = big_endian_int(temp);
  outStream.write((const char *) & swapint, sizeof (temp));

  // Write out size of class revision and the class revision
  const int revisionLength = strlen(mClassRevision);
  swapint = big_endian_int(revisionLength);
  outStream.write((const char *) & swapint, sizeof(revisionLength));
  outStream.write(mClassRevision, revisionLength);

  // write out params
  swapdouble = big_endian_double(_numiter);
  outStream.write((const char *) &swapdouble, sizeof(_numiter));
  swapdouble = big_endian_double(_alpha);
  outStream.write((const char *) &swapdouble, sizeof(_alpha));
  swapdouble = big_endian_double(_beta);
  outStream.write((const char *) &swapdouble,  sizeof(_beta));
  swapdouble = big_endian_double(_sigma);
  outStream.write((const char *) &swapdouble, sizeof(_sigma));
  swapdouble = big_endian_double(_numt);
  outStream.write((const char *) &swapdouble, sizeof(_numt));

  // read in vector field

  _TotalHField.save(outStream);

  // Save parent class's parameters
  if ( FieldTransformation::saveClassParameters( outStream ) )
  {
    return 1;
  };

  return 0;
}

int FluidLandmarkTransformation::loadClassParameters ( ifstream & inStream)
{
  const char * const functionName = "FluidLandmarkTransformation::loadClassParameters";

  // Read in the length of the class name from the file
  int ReadInClassNameLength;
  inStream.read((char *) &ReadInClassNameLength, sizeof (ReadInClassNameLength));
  ReadInClassNameLength =big_endian_int(ReadInClassNameLength);
  if (ReadInClassNameLength < 0 || ReadInClassNameLength > 1000) {
    cerr << "I don't think that this is a FluidLandmarkTransformation Results File." << endl;
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
  // End-of-string so I can use strcmp()
  ReadInClassName[ReadInClassNameLength] = '\0'; 

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
  ReadInClassVersion =big_endian_int(ReadInClassVersion);
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


  // read in params
  inStream.read((char *) &_numiter, sizeof(_numiter));
  inStream.read((char *) &_alpha, sizeof(_alpha));
  inStream.read((char *) &_beta,  sizeof(_beta));
  inStream.read((char *) &_sigma, sizeof(_sigma));
  inStream.read((char *) &_numt, sizeof(_numt));
  
  //byte swap params if neccesary
  _numiter = (big_endian_double(_numiter));
  _alpha= (big_endian_double( _alpha));
  _beta = (big_endian_double( _beta));
  _sigma= (big_endian_double( _sigma));
  _numt= (big_endian_double(_numt));
  
  // read in vector field

  _TotalHField.load(inStream);

  if ( FieldTransformation::loadClassParameters( inStream ) )
  {
    delete ReadInClassName;
    delete ReadInClassRevision;
    return 5;
  };

  // Return with success
  delete ReadInClassName;
  delete ReadInClassRevision;
  return 0;
}


int FluidLandmarkTransformation::load( const char * filename )
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

  // cout << "load() in FluidLandmarkTransformation" << endl;
  inStream.close();

  return 0;
}




int FluidLandmarkTransformation::save( const char * filename )
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

  // cout << "save() in FluidLandmarkTransformation" << endl;
  return 0;
}


//////////////////////////////////////////////////////////////////////////////////
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

ostream & FluidLandmarkTransformation::printClassParameters ( ostream & outStream)
{
  FieldTransformation::printClassParameters ( outStream );
  outStream << mClassName << "::mClassVersion: " << mClassVersion << endl;
  outStream << mClassName << "::mClassRevision: " << mClassRevision << endl;

  return outStream;
}

void FluidLandmarkTransformation::print() 
{
  FluidLandmarkTransformation::printClassParameters(cout);
}




