///////////////////////////////////////////////////////////////////////////
//
// IntellX, L.L.C., Proprietary
// Do not reproduce without permission in writing.
// Copyright (c) 1998 IntellX, L.L.C.
// All rights reserved
//
// File: SphereLandmarkTransformation.C
//
// SphereLandmarkTransformation class
//
// Author: Muge Bakircioglu
//	   08/99
//
///////////////////////////////////////////////////////////////////////////

#include <TDL/SphereLandmarkTransformation.h>
#include <ADL/Matrix_linalg_t.h>
#include <TDL/Point.h>
#include <TDL/SurfaceUtils.h>

const char * const SphereLandmarkTransformation::mClassRevision = "$Id: SphereLandmarkTransformation.C,v 1.1 2004/11/15 04:44:08 joeh Exp $";
const char * const SphereLandmarkTransformation::mClassName ="SphereLandmarkTransformation";
const int SphereLandmarkTransformation::mClassVersion = 0;

//
// LAPACK routines
//
extern "C" {
extern void dsytrf_(char *UPLO, int *N, double *A, int *LDA,
		    int *IPIV, double *work, int *LWORK, int *INFO);

extern void dsytrs_(char *UPLO, int *N, int *NRHS, double *A,
		    int *LDA, int *IPIV, double *B, int *LDB, int *INFO);
};


SphereLandmarkTransformation::SphereLandmarkTransformation() :
  _sigma(0.1),
  _beta(0.0002),
  _numt(10),
  _numiter(200)
{
//
}

SphereLandmarkTransformation::~SphereLandmarkTransformation() {
  // empty
}

int SphereLandmarkTransformation::initialize(
					     const Matrix<double> &atlasPnts,
					     const Matrix<double> &patientPnts,
					     Surface &atlasSurf, 
					     Point &patientNP, 
					     Point & atlasNP) {
  cout<<" initializing SpherelandmarkTransformation "<<endl;   
  // error checking
  if ((atlasPnts.getNrow() != patientPnts.getNrow()) ||
      (atlasPnts.getNcol() != patientPnts.getNcol())) {
    return 1;
  }

  if ((atlasPnts.getNcol() != 3) || (patientPnts.getNcol() != 3)) {
    return 2;
  }

  if ((atlasPnts.getNrow() < 4) || (patientPnts.getNrow() < 4)) {
    return 3;
  }
 
 
  Rot.setDim(3,3);   
 //   _atlasPoints   = atlasPnts;
//    _patientPoints = patientPnts;
  _numVert       = atlasSurf.getNumVert();
  _numPoints     = _atlasPoints.getNrow() +1;
  _atlasSurf     = atlasSurf;

  
	  //add the fixed point to the list of landmarks

  _atlasPoints.setDim(_numPoints,3);
  _patientPoints.setDim(_numPoints,3);
  int i;
  for (i=0; i<_numPoints-1; ++i) {
	  _atlasPoints[i][0]=atlasPnts[i][0];
	  _atlasPoints[i][1]=atlasPnts[i][1];
	  _atlasPoints[i][2]=atlasPnts[i][2];

	  _patientPoints[i][0]=patientPnts[i][0];
	  _patientPoints[i][1]=patientPnts[i][1];
	  _patientPoints[i][2]=patientPnts[i][2];
  }
	 
  _atlasPoints[_numPoints-1][0] = atlasNP.x();
  _atlasPoints[_numPoints-1][1] = atlasNP.y();
  _atlasPoints[_numPoints-1][2] = atlasNP.z();

  _patientPoints[_numPoints-1][0] = patientNP.x();
  _patientPoints[_numPoints-1][1] = patientNP.y();
  _patientPoints[_numPoints-1][2] = patientNP.z();
  
  double  rad;
  
  for (i=0; i<_numPoints; ++i) {
    rad=sqrt(_atlasPoints[i][0]*_atlasPoints[i][0] +
	     _atlasPoints[i][1]*_atlasPoints[i][1] +
	     _atlasPoints[i][2]*_atlasPoints[i][2] );
    _atlasPoints[i][0]/= rad;
    _atlasPoints[i][1]/=rad;
    _atlasPoints[i][2]/=rad;

    rad=sqrt(_patientPoints[i][0]*_patientPoints[i][0]+
	     _patientPoints[i][1]*_patientPoints[i][1]+
	     _patientPoints[i][2]*_patientPoints[i][2]);
    _patientPoints[i][0]/=rad;
    _patientPoints[i][1]/=rad;
    _patientPoints[i][2]/=rad;
  }

	  //store target radius to rescale sphere at the end
  targRad=rad;  

  Surface *sphere;    
  sphere=&_atlasSurf;
  
  SurfaceUtils::projectSurfaceToSphere(*sphere);

 //   SurfaceUtils::alignNP(atlasNP,_atlasPoints);
//    SurfaceUtils::alignNP(patientNP,_patientPoints);
//    SurfaceUtils::alignNP(patientNP,*sphere);

//    sphere->translate(0.0,0.0,1.0);
    
  pathx.setDim(_numt,_numPoints);
  pathy.setDim(_numt,_numPoints);
  pathz.setDim(_numt,_numPoints);
 
  if (pathx.isEmpty() || pathy.isEmpty() || pathz.isEmpty() )
    return 5;

  
  v1x.setDim(_numPoints-1);
  v1y.setDim(_numPoints-1);
  v1z.setDim(_numPoints-1);
  
  v2x.setDim(_numPoints-1);
  v2y.setDim(_numPoints-1);
  v2z.setDim(_numPoints-1);
  
  v3x.setDim(_numt,_numPoints-1);
  v3y.setDim(_numt,_numPoints-1);
  v3z.setDim(_numt,_numPoints-1);

  if(v1x.isEmpty() || v1y.isEmpty()  ||  v1z.isEmpty() ||
     v2x.isEmpty() || v2y.isEmpty()  ||  v2z.isEmpty() ||
     v3x.isEmpty() || v3y.isEmpty()  ||  v3z.isEmpty() )
    return 6;

  wtu.setDim(_numPoints);
  wtv.setDim(_numPoints);

  K1.setDim(_numPoints-1, _numPoints-1);
  K2.setDim(_numPoints-1, _numPoints-1);
  mat.setDim(_numPoints, _numPoints);

  if(wtu.isEmpty() || wtv.isEmpty() ||
     K1.isEmpty() || K2.isEmpty() || mat.isEmpty()) 
    return 7;
  
  setOkToCalculate(1);
  
  cout<<" done initializing "<<endl;
  return(0);	
}



int SphereLandmarkTransformation::calculate(Surface &outSurf) {
  
  Surface *sphere;    
  sphere=&_atlasSurf; 
  
  ITXTimer timer1("find path");
  //ITXTimer timer2("apply affine");
  ITXTimer timer3("calc weights");
  ITXTimer timer4("update h field");
 
  cout << "calculate() in SphereLandmarkTransformation" << endl;

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
  sphere->translate(0.0,0.0,1.0);

  int    step;
  float uu,uv;
  double tmp1;

  ShowWorking("Updating surface\n");
  
  float dwork_perc = 1.0/((float)(_numt-1));
  float  work_perc = 0.0;

  for (step=0; step<=_numt-2; ++step, work_perc += dwork_perc) { 
    ShowWorking(work_perc);

    cout << "Step " << step << endl;
    
    timer3.start();
    
    calculateWeights( pathx.data(step), pathy.data(step), pathz.data(step), 
		      pathx.data(step+1), pathy.data(step+1), pathz.data(step+1) );
    timer3.stop();
    
    timer4.start();
    
    for(int i=0; i<_numVert; i++) {
      uv = uu =  0.0;
      Point path, u_basis, v_basis;
 
      for(int l=0; l<_numPoints; l++) {

	//center back to origin
	
  	float sol_ang = pathx[step][l] * sphere->vertices()[i].x() +
	  pathy[step][l] *  sphere->vertices()[i].y() +
	  (pathz[step][l] - 1) * (sphere->vertices()[i].z() - 1); 
	tmp1=SurfaceUtils::BuildCov(sol_ang,50);
      
	uu   += (wtu[l]*tmp1);
	uv   += (wtv[l]*tmp1);
	
      }
      
      //find the basis for the tangent space for the current point
      SurfaceUtils::findStereoBasis(sphere->vertices()[i], u_basis,v_basis); 
      Point tmp= sphere->vertices()[i] +  Point( uu * u_basis.x() + uv * v_basis.x(),
				     uu * u_basis.y() + uv * v_basis.y(),
				     uu * u_basis.z() + uv * v_basis.z() - 1.0);
	  sphere->changeVert(i,tmp);
    
	  SurfaceUtils::projectSurfaceToSphere(*sphere);
      //project back to sphere
     
	  tmp=sphere->vertices()[i] + Point(0.0,0.0, 1.0);
	  sphere->changeVert(i,tmp);
	} //end i 
    timer4.stop();
  } //end step
  sphere->translate(0.0,0.0,-1.0);
  Point newNP= Point(_patientPoints[_numPoints-1][0],
					 _patientPoints[_numPoints-1][1],
					 _patientPoints[_numPoints-1][2]);
  SurfaceUtils::unalignNP(newNP,*sphere);
  sphere->scale(targRad);
	  

  setOkToApply(1);
  ShowWorking("DONE.\n");
  ShowWorking((char*)NULL);
  
  timer1.displayResults();
  //timer2.displayResults();
  timer3.displayResults();
  timer4.displayResults();

  outSurf=*sphere;
  return (0);
}

int SphereLandmarkTransformation::findPath() {
  int    i,j,t;
  char   UPLO;
  int    INFO,LDA,LWORK,N,NRHS;
  int tempx, tempy,  tempz;

  Array1D<int>    IPIV(_numPoints);
  Array1D<double> work(_numPoints);

  SurfaceUtils::findSphereRotation(_atlasPoints,_patientPoints,&Rot);
  Matrix<double> newP = _atlasPoints * Rot;
  _atlasPoints = newP;  

  Vector<double> B(3);
  Matrix<double> temp1(3,3);
  B[0]=B[1]=B[2]=0;
  Rot.transpose(&temp1);
  _atlasSurf.affine(temp1,B);

  
  Point atlasNP= Point(_atlasPoints[_numPoints-1][0],_atlasPoints[_numPoints-1][1],
					   _atlasPoints[_numPoints-1][2] );
  SurfaceUtils::alignNP(atlasNP,_atlasPoints);
  Point patientNP= Point(_patientPoints[_numPoints-1][0],
						 _patientPoints[_numPoints-1][1],
						 _patientPoints[_numPoints-1][2] );
  SurfaceUtils::alignNP(patientNP, _patientPoints);
  SurfaceUtils::alignNP(patientNP,_atlasSurf);

  
  cout << "FindPath 1" << endl;
  Array1D<Point> path(_numt);
 
  Point temp,targ;
 
  for (i=0; i<_numPoints;i++) {
	  temp.set(_atlasPoints[i][0], _atlasPoints[i][1], _atlasPoints[i][2]);
	  targ.set(_patientPoints[i][0], _patientPoints[i][1], _patientPoints[i][2]);
	  SurfaceUtils::rotatePtToPt(temp,targ,path,_numt);
	  for (t=0; t<_numt; ++t) {
		  pathx[t][i]=path[t].x();
		  pathy[t][i]=path[t].y();
		  pathz[t][i]=path[t].z();
	  }
  }

  double tmpx,tmpy,tmpz,tmp;
  
  float dwork = 1.0/((float)_numiter);
 float work_perc = 0.0;
 
 for(int iter=0;iter<_numiter;iter++,work_perc += dwork) {
   
   cout << "FindPath iter " << iter << " np = " << _numPoints << endl;
   if(iter % 10 == 0)
     ShowWorking(work_perc);

    for(t=1;t<_numt;t++) {
      for(i=0;i<_numPoints-1;i++) {
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
     
      for(i=0;i<_numPoints-1;i++)
		  for(j=0;j<_numPoints-1;j++) {
			  tmp = pathx[t][i] * pathx[t][j] + pathy[t][i] * pathy[t][j] +
				  pathz[t][i] * pathz[t][j];
			  K1[i][j]=(SurfaceUtils::BuildCov(tmp,50));

			  tmp = pathx[t-1][i] * pathx[t-1][j] + pathy[t-1][i] * pathy[t-1][j]
				  + pathz[t-1][i] * pathz[t-1][j];
			  K2[i][j]=(SurfaceUtils::BuildCov(tmp,50));
		  } 
	  

      UPLO  = 'U';
      INFO  = 1;
      LDA   = _numPoints-1;
      LWORK = _numPoints;
      N     = _numPoints-1;
      NRHS  = 1;
	
      dsytrf_(&UPLO,&N,K1.data(),&LDA,IPIV.data(),work.data(),&LWORK,&INFO);
      if (INFO!=0) cout<<" error in K1 "<<iter<<endl;	  
      dsytrs_(&UPLO,&N,&NRHS,K1.data(),&LDA,IPIV.data(),v1x.data(),&LDA,&INFO);
      dsytrs_(&UPLO,&N,&NRHS,K1.data(),&LDA,IPIV.data(),v1y.data(),&LDA,&INFO);
      dsytrs_(&UPLO,&N,&NRHS,K1.data(),&LDA,IPIV.data(),v1z.data(),&LDA,&INFO);

      UPLO  = 'U';
      INFO  = 1;
      LDA   = _numPoints-1;
      LWORK = 8;
      N     = _numPoints-1;
      NRHS  = 1;

      dsytrf_(&UPLO,&N,K2.data(),&LDA,IPIV.data(),work.data(),&LWORK,&INFO);
      if (INFO!=0) cout<<" error in K2 "<<iter<<endl;    
      dsytrs_(&UPLO,&N,&NRHS,K2.data(),&LDA,IPIV.data(),v2x.data(),&LDA,&INFO);
      dsytrs_(&UPLO,&N,&NRHS,K2.data(),&LDA,IPIV.data(),v2y.data(),&LDA,&INFO);
      dsytrs_(&UPLO,&N,&NRHS,K2.data(),&LDA,IPIV.data(),v2z.data(),&LDA,&INFO);
      
     
      for(i=0;i<_numPoints-1;i++) {
		  v3x[t][i] = -2*(v1x[i] - v2x[i]);
		  v3y[t][i] = -2*(v1y[i] - v2y[i]);
		  v3z[t][i] = -2*(v1z[i] - v2z[i]);
      }
	
      for(i=0;i<_numPoints-1;i++) {
		  tmpx = tmpy = tmpz = 0.0;
		  for(j=0;j<_numPoints-1;j++) {
			  tmp = pathx[t][i] * pathx[t][i] + pathy[t][i] * pathy[t][i] +
				  pathz[t][i] * pathz[t][i];
			  double sum=0; 
			  for (int n=1; n<50; ++n) 
				  sum+=(SurfaceUtils::plgndr_der(n,tmp))*(2*n+1)
					  /(n*n*(n+1)*(n+1));
	  
			  tmpx += (v1x[j]*sum*pathx[t][j]);
			  tmpy += (v1y[j]*sum*pathy[t][j]);
			  tmpz += (v1z[j]*sum*pathz[t][j]);
	  
		  } //end for j 
		  
		  v3x[t][i] += v1x[i]*tmpx;
		  v3y[t][i] += v1y[i]*tmpy;
		  v3z[t][i] += v1z[i]*tmpz;

		  if(t == _numt-1) {
			  tmp = (pathx[t][i] * _patientPoints[i][0]) + 
				  (pathy[t][i] * _patientPoints[i][1]) + 
				  (pathz[t][i] * _patientPoints[i][2]);
			  if ( fabs(tmp) < 1) {
				  v3x[t][i] += (2.0*acos(tmp)*(-1/sqrt(1-tmp*tmp))*
								_patientPoints[i][0])/_sigma;
				  v3y[t][i] += (2.0*acos(tmp)*(-1/sqrt(1-tmp*tmp))*
								_patientPoints[i][1])/_sigma;
				  v3z[t][i] += (2.0*acos(tmp)*(-1/sqrt(1-tmp*tmp))*
								_patientPoints[i][2])/_sigma;
				  
	

			  } // end if 
		  }//end t ==
      } // end i  
    } //end t

    for(t=1;t<_numt;t++)
		for(i=0;i<_numPoints-1;i++) {
			pathx[t][i] -= (_beta*v3x[t][i]);
			pathy[t][i] -= (_beta*v3y[t][i]);
			pathz[t][i] -= (_beta*v3z[t][i]);
			double rad = sqrt(pathx[t][i]*pathx[t][i] + 
							  pathy[t][i]*pathy[t][i] + pathz[t][i] * pathz[t][i]);
			pathx[t][i] /= rad;
			pathy[t][i] /= rad;
			pathz[t][i] /= rad;
	
      }
 }// for iter =
 for (t=0; t<_numt; ++t) {
	 pathx[t][_numPoints-1] = pathy[t][_numPoints] =0;
	 pathz[t][_numPoints-1] = 1.0;
 }
 cout<<" paths at the end "<<endl;
 for (i=0; i<_numPoints; ++i) 
	 for (t=0; t<_numt; ++t) {
		 pathz[t][i] +=1;
		 cout<<pathx[t][i]<<"\t"<<pathy[t][i]<<"\t"<<pathz[t][i]<<endl;
	 }
 return 0;
}



void SphereLandmarkTransformation::calculateWeights(double *x1, double *y1, 
						    double *z1, double *x2, 
						    double *y2, double *z2) {
  int i,j;
  double tmp,du,dv;
  Point P,Q; 
 
  for(i=0;i<_numPoints; i++)
    for(j=0;j<_numPoints; j++) {
      tmp=x1[i] * x1[j] + y1[i] * y1[j] + (z1[i]-1) * (z1[j]-1);
      mat[i][j] = SurfaceUtils::BuildCov(tmp,50);
    }
 
  Point u,v;
  
  for(i=0;i<_numPoints-1; i++) {
    double diff_x = x2[i] - x1[i]; 
    double diff_y = y2[i] - y1[i];
    double diff_z = z2[i] - z1[i];
    
    Point P = Point(x1[i], y1[i], z1[i]); 
    Point Q = Point(x2[i], y2[i], z2[i]);
    Point diff = Q - P;
    
    SurfaceUtils::findStereoBasis(P,u,v); 
    wtu[i]= diff.innerProd(u);
    wtv[i]= diff.innerProd(v);
  }
  
  wtu[_numPoints-1]= wtv[_numPoints-1]=0;
 
  char UPLO  = 'U';
  int INFO   = 1;
  int LDA    = _numPoints;
  int LWORK  = 8;
  int N      = _numPoints;
  int NRHS   = 1;


  Array1D<int>    IPIV(_numPoints);
  Array1D<double> work(_numPoints);

  dsytrf_(&UPLO,&N,mat.data(),&LDA,IPIV.data(),work.data(),&LWORK,&INFO);
  cout<<" info "<<INFO<<endl;
  dsytrs_(&UPLO,&N,&NRHS,mat.data(),&LDA,IPIV.data(),wtu.data(),&LDA,&INFO);
  dsytrs_(&UPLO,&N,&NRHS,mat.data(),&LDA,IPIV.data(),wtv.data(),&LDA,&INFO);

}

int SphereLandmarkTransformation::saveClassParameters ( ofstream & outStream) {
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
  // swapdouble = big_endian_double(_alpha);
  //outStream.write((const char *) &swapdouble, sizeof(_alpha));
  swapdouble = big_endian_double(_beta);
  outStream.write((const char *) &swapdouble,  sizeof(_beta));
  swapdouble = big_endian_double(_sigma);
  outStream.write((const char *) &swapdouble, sizeof(_sigma));
  swapdouble = big_endian_double(_numt);
  outStream.write((const char *) &swapdouble, sizeof(_numt));

  // read in vector field
  
  // _TotalHField.save(outStream);
  
  // Save parent class's parameters
  if ( FieldTransformation::saveClassParameters( outStream ) ) {
    return 1;
  };

  return 0;
}

int SphereLandmarkTransformation::loadClassParameters ( ifstream & inStream) {
  const char * const functionName = 
    "SphereLandmarkTransformation::loadClassParameters";

  // Read in the length of the class name from the file
  int ReadInClassNameLength;
  inStream.read((char *) &ReadInClassNameLength, sizeof (ReadInClassNameLength));
  ReadInClassNameLength =big_endian_int(ReadInClassNameLength);
  if (ReadInClassNameLength < 0 || ReadInClassNameLength > 1000) {
    cerr << "I don't think that this is a SphereLandmarkTransformation Results File." << endl;
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
  if (0 != strcmp(ReadInClassName, mClassName)) {
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
    cerr << "Attempted to read in " << mClassName 
	 << " Class parameters from version ";
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
    cerr << "I don't think that this is a " << mClassName 
	 << " Results File." << endl;
    cerr << "Apparent class revision length is " 
	 << ReadInClassRevisionLength << endl;
    // Throw error here!
    return 6;
  }
  // Alloc memory for the class revision and read it in
  char * ReadInClassRevision = new char[ReadInClassRevisionLength + 1];
  if ( ! ReadInClassRevision ){
    cerr << functionName << ": ";
    cerr << "Failed to alloc " << ReadInClassRevisionLength+1 
	 << " bytes." << endl;
    // Throw error here!
    return 4;
  }
  inStream.read((char *) ReadInClassRevision, ReadInClassRevisionLength);
  ReadInClassRevision[ReadInClassRevisionLength] = '\0'; // End-of-string so I can use strcmp()
  // Do nothing with it for now, mainly for debugging outside class


  // read in params
  inStream.read((char *) &_numiter, sizeof(_numiter));
  // inStream.read((char *) &_alpha, sizeof(_alpha));
  inStream.read((char *) &_beta,  sizeof(_beta));
  inStream.read((char *) &_sigma, sizeof(_sigma));
  inStream.read((char *) &_numt, sizeof(_numt));
  
  //byte swap params if neccesary
  _numiter = (big_endian_double(_numiter));
  // _alpha= (big_endian_double( _alpha));
  _beta = (big_endian_double( _beta));
  _sigma= (big_endian_double( _sigma));
  _numt= (big_endian_double(_numt));
  
  // read in vector field

  //  _TotalHField.load(inStream);

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


int SphereLandmarkTransformation::load( const char * filename ) {
  ifstream inStream (filename, ios::in);
  if ( ! inStream ) {
    cerr << "Failed to open file \"" << filename << "\"" << endl;
    return 1;
  }

  if ( loadClassParameters( inStream ) ) {
    return 2;
  }

  // cout << "load() in SphereLandmarkTransformation" << endl;
  inStream.close();

  return 0;
}


int SphereLandmarkTransformation::save( const char * filename ) {
  ofstream outStream (filename, ios::out);

  if ( ! outStream ) {
    cerr << "Failed to open file \"" << filename << "\"" << endl;
    return 1;
  }
  if ( saveClassParameters( outStream ) ) {
    return 2;
  };
  outStream.close();

  // cout << "save() in SphereLandmarkTransformation" << endl;
  return 0;
}


//////////////////////////////////////////////////////////////////////////////////
// Function: printClassParameters
//
// Purpose: Calls parent's print function and then prints parameters 
//          that are owned by this class to the supplied ostream.
//          It is envisioned that the child of this class calls this 
//          function and then prints out it's own parameters.
//
// Inputs: ostream that we want to print to
//
// Outputs: reference to the passed in ostream
//
////////////////////////////////////////////////////////////////////////////////

ostream & SphereLandmarkTransformation::printClassParameters (ostream & 
							      outStream) {
  FieldTransformation::printClassParameters ( outStream );
  outStream << mClassName << "::mClassVersion: " << mClassVersion << endl;
  outStream << mClassName << "::mClassRevision: " << mClassRevision << endl;
  
  return outStream;
}

void SphereLandmarkTransformation::print()  {
  SphereLandmarkTransformation::printClassParameters(cout);
}




