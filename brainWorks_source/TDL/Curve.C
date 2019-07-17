///////////////////////////////////////////////////////////////////////////
// 
// IntellX, L.L.C., Proprietary
// The contents of this file comprise confidential, proprietary and/or
// trade secret information which is the property of IntellX LLC.
// Dissemination, distribution or other disclosure of this information
// is strictly prohibited without the express permission of IntellX LLC.
// Copyright (c) 1999 IntellX, L.L.C.
// All rights reserved
//
// File: CurveMatch.C
//
// Author: Muge Bakircioglu
//
// Purpose:
//
///////////////////////////////////////////////////////////////////////////


#include <TDL/Curve.h>
#include <ADL/Vector.h>
#include <ADL/Matrix.h>

//copy constructor

Curve::Curve(Curve const &C)  {
  *this=C;
}

/******************************************************************/
//assignment operator

Curve & Curve::operator=(Curve const &C) {
 
  mNumPoints=C.mNumPoints;
  mCurvePoints.setDim(C.numSamples());
  mCurvePoints=C.mCurvePoints;
  mSpeed.setDim(C.numSamples());
  mSpeed=C.mSpeed;
  mCurvature.setDim(C.numSamples());
  mCurvature=C.mCurvature;
  mTorsion.setDim(C.numSamples());
  mTorsion=C.mTorsion;
  mArcLength.setDim(C.numSamples());
  mArcLength=C.mArcLength;

  Array1D<int> temp(2);
  temp=(C.mPrefer).getDim();
  mPrefer.setDim(temp);
  mPrefer=C.mPrefer;

  mCost.setDim(temp);
  mCost=C.mCost;
  mMatch.setDim((C.mMatch).getSize());
  mMatch=C.mMatch;
  mInd.setDim((C.mInd).getSize());
  mInd=C.mInd;
  sprintf(mName,"%s",C.mName);
  return *this;
}
/*****************************************************************/
void  Curve::reset() {
  mNumPoints=0;
    mCurvePoints.setDim(0);
    mSpeed.setDim(0);
    mArcLength.setDim(0);
    mCurvature.setDim(0);
    mTorsion.setDim(0);
    mMatch.setDim(0);
    mInd.setDim(0);
    mCost.setDim(0,0);
    mPrefer.setDim(0,0);
}

/*****************************************************************/
Point &  Curve::getPoint(int i) {
	return mCurvePoints[i];
}

/******************************************************************/
//add new point to the curve
void Curve::addPoint(Point  &newPoint, double &newArcLength,
		     double &newSpeed, double &newCurvature,
  		     double &newTorsion) {
  
  //store existing curve
  Array1D<Point> tempPoints=mCurvePoints;
  Array1D<double> tempSpeed=mSpeed;
  Array1D<double> tempCurvature=mCurvature;
  Array1D<double> tempTorsion=mTorsion;
  Array1D<double> tempArclength=mArcLength;
  
//increment curve size and reallocate
  mNumPoints++;
  
  mCurvePoints.setDim(mNumPoints);
  mSpeed.setDim(mNumPoints);
  mCurvature.setDim(mNumPoints);
  mTorsion.setDim(mNumPoints);
  mArcLength.setDim(mNumPoints);
  mTangent.setDim(mNumPoints);
  mNormal.setDim(mNumPoints);
  mBinormal.setDim(mNumPoints);

  //copy old points over
  for (int i=0; i<mNumPoints-1; ++i) {
    mCurvePoints[i]=tempPoints[i];
    mSpeed[i]=tempSpeed[i];
    mCurvature[i]=tempCurvature[i];
    mTorsion[i]=tempTorsion[i];
    mArcLength[i]=tempArclength[i];
  }

  //add new point
  mCurvePoints[mNumPoints-1]=newPoint;
  mArcLength[mNumPoints-1]=newArcLength;
  mSpeed[mNumPoints-1]=newSpeed;
  mTorsion[mNumPoints-1]=newTorsion;
  mCurvature[mNumPoints-1]=newCurvature;

}

/******************************************************************/
void  Curve::resample(int new_numPoints, int orderPoly, 
			 Array2D<double> &result) {
  
  Vector<double> X(mNumPoints);
  X[0]=0;
  Point temp;
  int i;
  
  //points on [0,1] to fit the initial polynomial to 
  for (i=1; i<mNumPoints; ++i) {
    temp=mCurvePoints[i]-mCurvePoints[i-1];
    X[i]=temp.norm()+X[i-1];
  }
  
  //normalized
  for (i=1; i<mNumPoints; ++i)
    X[i]/=X[mNumPoints-1];
  
  //the values of the polynomial at the points above
  Vector<double> Cx_coord(mNumPoints);
  Vector<double> Cy_coord(mNumPoints);
  Vector<double> Cz_coord(mNumPoints);
  
  for (i=0; i<mNumPoints; ++i) {
    Cx_coord[i]=mCurvePoints[i].x();
    Cy_coord[i]=mCurvePoints[i].y();
    Cz_coord[i]=mCurvePoints[i].z();
  }
  
  //solve for polynomial  coefficients
  Vector<double> XpolyCoeff=MatrixUtils<double>::polyfit(X,Cx_coord,orderPoly);
  Vector<double> YpolyCoeff=MatrixUtils<double>::polyfit(X,Cy_coord,orderPoly);
  Vector<double> ZpolyCoeff=MatrixUtils<double>::polyfit(X,Cz_coord,orderPoly);

  //new set of points to evaluate the polynomial at
  X.setDim(new_numPoints);
  for ( i=0; i<new_numPoints; ++i) 
    X[i]=(double) i/(new_numPoints-1);
 
  //allocate space for all the members
  mNumPoints=new_numPoints;
  mCurvePoints.setDim(mNumPoints);
  mSpeed.setDim(mNumPoints);
  mArcLength.setDim(mNumPoints);
  mCurvature.setDim(mNumPoints);
  mTorsion.setDim(mNumPoints);
  mTangent.setDim(mNumPoints);
  mNormal.setDim(mNumPoints);
  mBinormal.setDim(mNumPoints);

  double tempx,tempy,tempz;
  double powerX;
 
  //now evaluate the polynomial at these new points
  for (i=0; i<mNumPoints; ++i) {
    tempx=Vector<double>::computePolynomial(XpolyCoeff,X[i]);
    tempy=Vector<double>::computePolynomial(YpolyCoeff,X[i]);
    tempz=Vector<double>::computePolynomial(ZpolyCoeff,X[i]);
    mCurvePoints[i].set(tempx,tempy,tempz);
  }

  //calculate arclength
  mArcLength[0]=0;
  for (i=1; i<mNumPoints; ++i)
    mArcLength[i]=(mCurvePoints[i]-mCurvePoints[i-1]).norm()+mArcLength[i-1];

  //return the coefficients in a 2D array 

  for (i=0; i<=orderPoly; ++i) {
    result[0][i]=XpolyCoeff[i];
    result[1][i]=YpolyCoeff[i];
    result[2][i]=ZpolyCoeff[i];
  }
 
}

/**************************************************************/
void Curve::calculate_geometry(int new_numPoints,int orderPoly) {
  Vector<double> X(new_numPoints);
  cout<<" calculating geo "<<endl;
  Array2D<double> polyCoeff(3,orderPoly+1);
 
  //resample the curve with the new number of points
  resample(new_numPoints,orderPoly,polyCoeff);
 
  //evaluate evenly spaced points at these coefficients
  //derivatives will be used in curvature and torsion calculations
  int i;
  for(i=0; i<mNumPoints; ++i)
    X[i]=(double)i/(mNumPoints-1);

  Array1D<Point> firstDer(new_numPoints);
  Array1D<Point> secondDer(new_numPoints);
  Array1D<Point> thirdDer(new_numPoints);
  Point temp;
  
  for (i=0; i<mNumPoints; ++i) {
    firstDer[i].set(0.0,0.0,0.0);
    secondDer[i].set(0.0,0.0,0.0);
    thirdDer[i].set(0.0,0.0,0.0);
    for (int j=1; j<=orderPoly; ++j) {
      firstDer[i]+=Point(polyCoeff[0][j]*j*pow(X[i],(j-1)),
			 polyCoeff[1][j]*j*pow(X[i],(j-1)),
			 polyCoeff[2][j]*j*pow(X[i],(j-1)));
      if (j>=2) 
	secondDer[i]+=Point(polyCoeff[0][j]*j*(j-1)*pow(X[i],(j-2)),
			    polyCoeff[1][j]*j*(j-1)*pow(X[i],(j-2)),
			    polyCoeff[2][j]*j*(j-1)*pow(X[i],(j-2)));
      if (j>=3)
	thirdDer[i]+=Point(polyCoeff[0][j]*j*(j-1)*(j-2)*pow(X[i],(j-3)),
			 polyCoeff[1][j]*j*(j-1)*(j-2)*pow(X[i],(j-3)),
		         polyCoeff[2][j]*j*(j-1)*(j-2)*pow(X[i],(j-3)));
      
    }
    
 
    mSpeed[i]=firstDer[i].norm();
    temp=firstDer[i].cross(secondDer[i]);
    mCurvature[i]=temp.norm()/(float)(pow(mSpeed[i],3));
    mTorsion[i]=temp.innerProd(thirdDer[i])/(temp.norm()*temp.norm());

    if (firstDer[i].norm()!=0) {
    mTangent[i]=firstDer[i]/(firstDer[i].norm());
    mBinormal[i]=temp/(temp.norm());
    mNormal[i]=mBinormal[i].cross(mTangent[i]);
    // cout<<i<<"\t"<<mTangent[i]<<"\t"<<mNormal[i]<<"\t"<<mBinormal[i]<<endl;
  
    }
  }
  
  
}

/**************************************************************/
void Curve::regenerate(Curve &Target) {
  Array1D<double> l(numSamples());
  Array1D<Point> newPoints(numSamples());
  double curvMax=0.0;
  double torsMin=HUGE;
  double torsMax=0.0;
  double speedMax=0.0;

  cout<<"regenerating "<<endl;
  int i;
  for (i=0; i<numSamples()-1; ++i) {
    l[i]=mArcLength[i+1]-mArcLength[i];
    //   cout<<mSpeed[i]<<"\t"<<l[i]<<endl;
  }
 
 
 //scale curvature between 0-10, torsion between -10 and 10 
  
  for (i=0; i<numSamples(); ++i) {
    if (mCurvature[i] > curvMax) curvMax=mCurvature[i];
    if (mTorsion[i] > torsMax) torsMax=mTorsion[i];
    if (mTorsion[i] < torsMin) torsMin=mTorsion[i];
  //  if (mSpeed[i] > speedMax) speedMax=mSpeed[i];
  }
  
  
  if (fabs(torsMax) < fabs(torsMin)) torsMax=-torsMin;
   for (i=0; i<numSamples(); ++i) {
    
     mCurvature[i]/=curvMax/10;
     mTorsion[i]/=torsMax/10;
   }
  cout<<" max t in temp. "<<torsMax<<" max curv "<<curvMax<<endl;
  
  curvMax=0.0;
  torsMin=HUGE;
  torsMax=0.0;

  for (i=0; i<Target.numSamples(); ++i) {
    if (Target.mCurvature[i] > curvMax) curvMax=Target.mCurvature[i];
    if (Target.mTorsion[i] > torsMax) torsMax=Target.mTorsion[i];
    else if (Target.mTorsion[i] < torsMin) torsMin=Target.mTorsion[i];
    //  if (mSpeed[i] > speedMax) speedMax=mSpeed[i];
  }
  
  
  if (fabs(torsMax) < fabs(torsMin)) torsMax=-torsMin;
   for (i=0; i<Target.numSamples(); ++i) {
     Target.mCurvature[i]/=curvMax/10;
     Target.mTorsion[i]/=torsMax/10;
   }
  cout<<" max t in targ. "<<torsMax<<" max curv "<<curvMax<<endl;
   
  double kappa,tau;
  Matrix<double> L(3,3); 
  Matrix<double> Z(3,3);
  Matrix<double> id(3,3);
  
  Vector<double> temp(3);
 
   id.eye();
   double beta;
   for ( i=0; i<numSamples()-1; ++i) {
     kappa=(mCurvature[i+1]+mCurvature[i])*l[i]*0.5;
     tau=(mTorsion[i+1]+mTorsion[i])*l[i]*0.5;

     L[0][0]=0;        L[0][1]=kappa;     L[0][2]=0;
     L[1][0]=-kappa;   L[1][1]=0;         L[1][2]=tau;   
     L[2][0]=0;        L[2][1]=-tau;      L[2][2]=0;
     beta=sqrt(kappa*kappa+tau*tau);
     Z=L/beta;
     L=id+Z*sin(beta)+Z*Z*(1-cos(beta));

     temp[0]=mTangent[i].x();
     temp[1]=mTangent[i].y();
     temp[2]=mTangent[i].z();
     temp=L*temp;
     mTangent[i+1].set(temp[0],temp[1],temp[2]);

     temp[0]=mNormal[i].x();
     temp[1]=mNormal[i].y();
     temp[2]=mNormal[i].z();
     temp=L*temp;
     mNormal[i+1].set(temp[0],temp[1],temp[2]);

     temp[0]=mBinormal[i].x();
     temp[1]=mBinormal[i].y();
     temp[2]=mBinormal[i].z();
     temp=L*temp;
     mBinormal[i+1].set(temp[0],temp[1],temp[2]);
    
   }


  newPoints[0]=mCurvePoints[0];
  cout<<newPoints[0]<<endl;
  for (i=0; i<numSamples()-1; ++i) {
    newPoints[i+1]=newPoints[i]+mTangent[i]*l[i]
      +mNormal[i]*mCurvature[i]*l[i]*l[i]/2
      +mBinormal[i]*mTorsion[i]*l[i]*l[i]*l[i]/6;
   
 
    cout<<newPoints[i]<<endl;
  }

}


/************************************************************/
void Curve::setName(const char *nm) {
  if(nm) sprintf(mName,"%s",nm);
  else   sprintf(mName,"<no-name>");
}

/************************************************************/
void Curve::print(ostream &os) const {
  for (int i=0; i<mNumPoints; ++i)
    os<<mCurvePoints[i]<<endl;
}

/************************************************************/
void Curve::print_geo() const {
 
  for (int i=0; i<mNumPoints; ++i)
    cout<<mArcLength[i]<<"\t"<<mSpeed[i]<<"\t"
	<<mCurvature[i]<<"\t"<<mTorsion[i]<<"\t"<<mMatch[i]<<endl;
 
}
/******************************************************************/
void Curve::calculateMatchcost(Curve  &Target, int numSmp,
			       int numLmks,  int  polyOrder,
			       double s_coef, double c_coef,
			       double t_coef, int nghbdSize, 
			       Curve &sampledTarget, int &numNghbs,
			       int &spacing, int &dist) {
 
  //calculate geometries of both Target and Template
  calculate_geometry(numSmp,polyOrder);
  Target.calculate_geometry(numSmp,polyOrder);
 
  //distance between lmks and length of neighborhood
  dist=int (ceil((double)numSmp/numLmks));
  spacing=int ((nghbdSize*numSmp/100)/2);
  int i;
  for(i=0; i<Target.numSamples();  i+=dist) 
    sampledTarget.addPoint(Target.mCurvePoints[i], Target.mArcLength[i], 
			   Target.mSpeed[i], Target.mCurvature[i], 
			   Target.mTorsion[i]);
 
  (sampledTarget.mPrefer).setDim(sampledTarget.numSamples(), spacing*2+1);
  (sampledTarget.mCost).setDim(sampledTarget.numSamples(),spacing*2+1);
  (sampledTarget.mMatch).setDim(sampledTarget.numSamples());
  (sampledTarget.mInd).setDim(sampledTarget.numSamples());
 
  int k;
 
  //calculate cost, if ngh out of range (<0 or >numPoints) assign MAX_COST
  for (i=1; i<sampledTarget.numSamples(); ++i) {
    k=0;
    
    for (int j=i*dist-spacing+1; j<=i*dist+spacing; ++j) {
      if (j<=0 || j>=numSamples()) sampledTarget.mCost[i][k++]=MAX_COST;
      else {
	
// 	sampledTarget.mCost[i][k++]=
// 	  s_coef*pow((sampledTarget.mArcLength[i]-mArcLength[j]),2) +
// 	  c_coef*pow((sampledTarget.mSpeed[i]*sampledTarget.mCurvature[i]-
// 		      mSpeed[j]*mCurvature[j]),2) +
// 	  t_coef*pow((sampledTarget.mSpeed[i]*sampledTarget.mTorsion[i]-
// 		      mSpeed[j]*mTorsion[j]),2);
	sampledTarget.mCost[i][k++]=
	  s_coef*pow((sampledTarget.mArcLength[i]-mArcLength[j]),2) +
	  c_coef*pow((sampledTarget.mCurvature[i]-mCurvature[j]),2) +
	  t_coef*pow((sampledTarget.mTorsion[i]-mTorsion[j]),2);
	
	
      } //end else
      // cout<<sampledTarget.mCost[i][k-1]<<"\t";
      //if (j==i*dist+spacing) cout<<endl;
    } //end j
  }  //end i
   
  numNghbs=k;
  
 
//initialize the prefer array in ascending order
  //this will be sorted in the matching routine below

  for (i=1; i<sampledTarget.numSamples(); ++i) 
    for (int j=0; j<numNghbs; ++j)
      sampledTarget.mPrefer[i][j]=j;
  
}

/*************************************************************************/
double Curve::matching(Curve &Target,int const &numNghbs, 
		     int const &spacing, int const &dist) {

  double cum_cost=0;
  int temp,i,j,k;
  if (Target.numSamples()==0) return 0.0;
  regenerate(Target);
  //sort the prefer array in increasing cost
  
  for (i=1; i<Target.numSamples(); ++i) {
    for (j=0; j<numNghbs-1; ++j) 
      for (k=0; k<numNghbs-1-j;  ++k) {
	
	if (Target.mCost[i][Target.mPrefer[i][k]] >
	    Target.mCost[i][Target.mPrefer[i][k+1]]) {
	
	  temp=Target.mPrefer[i][k];
	  Target.mPrefer[i][k]=Target.mPrefer[i][k+1];
	  Target.mPrefer[i][k+1]=temp;
	  
	}  //endif
      }  //end k 
  }//end i
 
  Target.mMatch[0]=0;
  cout<<"initial matches "<<endl;

  for (i=1; i<Target.numSamples(); ++i) {
   
  //initial match for each point is the one w/ lowest cost
  //0th position in the prefer array
    
    //offset into template array index
    Target.mMatch[i]=Target.mPrefer[i][0]+(i*dist-spacing+1);
    cout<<Target.mMatch[i]<<endl;
    
    //index shows how far down on the prefer list we've progressed
    //initialized to 0
    Target.mInd[i]=0;
  
    //initial cost of match is the sum of all least cost matches
    cum_cost+=Target.mCost[i][Target.mPrefer[i][0]];
  } //end i
  

  cout<<"initial cost "<<cum_cost<<endl;
  double delta_min;
  double diff1, diff2,delta;
  int *curr_ind, *next_ind, *curr_match, *next_match,j_min, k_min,notSorted;

  for (i=1; i<Target.numSamples()-1; ++i) {
    curr_ind=&Target.mInd[i];
    next_ind=&Target.mInd[i+1];
    delta_min=MAX_COST;
    curr_match=&Target.mPrefer[i][(*curr_ind)];
    next_match=&Target.mPrefer[i+1][(*next_ind)];
    
    if (Target.mMatch[i]>=Target.mMatch[i+1]) {
      cout<<"bad match "<<Target.mMatch[i]<<"\t"
	  <<Target.mMatch[i+1];
     
      //go down the preference list of both current and next point
      for (int j=(*curr_ind); j<numNghbs;  ++j) 
	for (int k=(*next_ind); k<numNghbs; ++k) { 
	  
	  //do not repeat the same point
	  if (!(j==*curr_ind &&  k==*next_ind)) {
	   
	    //check to see if diffeomorphic
	    if (Target.mPrefer[i][j]+(i*dist-spacing+1)
		<Target.mPrefer[i+1][k]+((i+1)*dist-spacing+1)) {
 
	      //this is a potential match, calculate change in cost
	      diff1=Target.mCost[i][(Target.mPrefer[i][j])]
		-Target.mCost[i][*curr_match];
	      diff2=Target.mCost[i+1][Target.mPrefer[i+1][k]]
		-Target.mCost[i+1][*next_match]; 
	      delta=diff1+diff2;

	      if (delta<delta_min) {
		j_min=j;
		k_min=k;
		delta_min=delta;
	      } //endif

	    }  //endif diffeomorphic
	  } //endif  
	  
	}  //end k

      *curr_ind=j_min;
      *next_ind=k_min;
      Target.mMatch[i]=Target.mPrefer[i][*curr_ind]+(i*dist-spacing+1);
      Target.mMatch[i+1]=Target.mPrefer[i+1][*next_ind]+((i+1)*dist-spacing+1);
      cum_cost+=delta_min;
      cout<<" new match "<<Target.mMatch[i]<<"\t"<<Target.mMatch[i+1]
	  <<" current cost "<<cum_cost<<endl;
      
      //check the match for previous points
      //if diffeomorphism broken, go back to  where it is broken 
      
      notSorted=Target.isDiffeo();
      if (notSorted!=0) {
      	if (notSorted<i) cout<<"old match broken "<<notSorted<<endl;
	i=notSorted-1;
      }
      
    } //endif not diffeomorphic
 
  } //end for i  
  
  cout<<" final match "<<endl;
  for (i=1; i<Target.numSamples(); ++i)
    cout<<Target.mMatch[i]<<endl;
  
  cout<<" cum_cost "<<cum_cost<<endl;
  return cum_cost;
}

/******************************************************************/
//matching with a window size of 3

void Curve::matching2(Curve &Target,int const &numNghbs, 
		      int const &spacing, int const &dist) {
  double cum_cost=0;
  int temp,i,j,k;
  
  //sort the prefer array in increasing cost
 
  for ( i=1; i<Target.numSamples(); ++i) {
    for (j=0; j<numNghbs-1; ++j) 
      for (k=0; k<numNghbs-1-j;  ++k) {
	
	if (Target.mCost[i][Target.mPrefer[i][k]] >
	    Target.mCost[i][Target.mPrefer[i][k+1]]) {
	
	  temp=Target.mPrefer[i][k];
	  Target.mPrefer[i][k]=Target.mPrefer[i][k+1];
	  Target.mPrefer[i][k+1]=temp;
	
	}  //endif
      }  //end k 
  }//end i
  cout<<" num lmks " <<Target.numSamples()<<endl;  
  cout<<"initial matches "<<endl;

  for (i=1; i<Target.numSamples(); ++i) {
   
  //initial match for each point is the one w/ lowest cost
  //0th position in the prefer array
    
    //offset into template array index
    Target.mMatch[i]=Target.mPrefer[i][0]+(i*dist-spacing+1);
    cout<<Target.mMatch[i]<<endl;
    
    //index shows how far down on the prefer list we've progressed
    //initialized to 0
    Target.mInd[i]=0;
  
    //initial cost of match is the sum of all least cost matches
    cum_cost+=Target.mCost[i][Target.mPrefer[i][0]];
  } //end i
  

  cout<<"initial cost "<<cum_cost<<endl;
  double delta_min;
  double diff1, diff2,diff3,delta;
  int *curr_ind, *next_ind,*far_ind, *curr_match, *next_match,*far_match, 
    j_min, k_min,l_min,notSorted;


  for (i=1; i<Target.numSamples()-2; ++i) {
    curr_ind=&Target.mInd[i];
    next_ind=&Target.mInd[i+1];
    far_ind=&Target.mInd[i+2];

    delta_min=MAX_COST;
    curr_match=&Target.mPrefer[i][(*curr_ind)];
    next_match=&Target.mPrefer[i+1][(*next_ind)];
    far_match=&Target.mPrefer[i+2][(*far_ind)];
 
    if (Target.mMatch[i] >= Target.mMatch[i+1] || 
	Target.mMatch[i-1] >= Target.mMatch[i]) {
    
    cout<<"bad match "<<Target.mMatch[i]<<"\t"
	<<Target.mMatch[i+1]<<"\t"<<Target.mMatch[i+2];
     
      //go down the preference list of both current and next point
    for (int j=(*curr_ind); j<numNghbs;  ++j) 
      for (int k=(*next_ind); k<numNghbs; ++k) 
	for (int l=(*far_ind); l<numNghbs; ++l) { 
	  
	  //do not repeat the same point
	  //if (!(j==*curr_ind &&  k==*next_ind))  {
	   
	  //check to see if diffeomorphic
	  if (Target.mPrefer[i][j]+(i*dist-spacing+1)
		<Target.mPrefer[i+1][k]+((i+1)*dist-spacing+1) &&
		Target.mPrefer[i+1][k]+((i+1)*dist-spacing+1)
		<Target.mPrefer[i+2][l]+((i+2)*dist-spacing+1)) {
 
	      //this is a potential match, calculate change in cost
	      diff1=Target.mCost[i][(Target.mPrefer[i][j])]
		-Target.mCost[i][*curr_match];
	      diff2=Target.mCost[i+1][Target.mPrefer[i+1][k]]
		-Target.mCost[i+1][*next_match];
	      diff3=Target.mCost[i+2][Target.mPrefer[i+2][l]]-
		Target.mCost[i+2][*far_match];
	      delta=diff1+diff2+diff3;

	      if (delta<delta_min) {
		j_min=j;
		k_min=k;
		l_min=l;
		delta_min=delta;
	      } //endif

	    }//endif diffeomorphic
	   // } endif  
	  
	  }//end l

      *curr_ind=j_min;
      *next_ind=k_min;
      *far_ind=l_min;
      Target.mMatch[i]=Target.mPrefer[i][*curr_ind]+(i*dist-spacing+1);
      Target.mMatch[i+1]=Target.mPrefer[i+1][*next_ind]+((i+1)*dist-spacing+1);
      Target.mMatch[i+2]=Target.mPrefer[i+2][*far_ind]+((i+2)*dist-spacing+1);
      cum_cost+=delta_min;
      cout<<" new match "<<Target.mMatch[i]<<"\t"<<Target.mMatch[i+1]
	  <<"\t"<<Target.mMatch[i+2]<<" current cost "<<cum_cost<<endl;
      
      //check the match for previous points
      //if diffeomorphism broken, go back to  where it is broken 
      
      notSorted=Target.isDiffeo();
      if (notSorted!=0) {
      	if (notSorted<i) cout<<"old match broken "<<notSorted<<endl;
	i=notSorted-1;
    }
      
    } //endif not diffeomorphic
 
    } //end for i  
  
  cout<<" final match "<<endl;
  for (i=1; i<Target.numSamples(); ++i)
    cout<<Target.mMatch[i]<<endl;
  
  cout<<" cum_cost "<<cum_cost<<endl;
}

/****************************************************************/
// void Curve::matchingN(Curve &Target,int const &numNghbs, 
// 		      int const &spacing, int const &dist, int const &N) {
//   Array1D<int> currMatch;

//  double cum_cost=0;
//   int temp,i,j,k;
  
//   //sort the prefer array in increasing cost
 
//   for ( i=1; i<Target.numSamples(); ++i) {
//     for (j=0; j<numNghbs-1; ++j) 
//       for (k=0; k<numNghbs-1-j;  ++k) {
	
// 	if (Target.mCost[i][Target.mPrefer[i][k]] >
// 	    Target.mCost[i][Target.mPrefer[i][k+1]]) {
	
// 	  temp=Target.mPrefer[i][k];
// 	  Target.mPrefer[i][k]=Target.mPrefer[i][k+1];
// 	  Target.mPrefer[i][k+1]=temp;
	
// 	}  //endif
//       }  //end k 
//   }//end i
//  cout<<"initial matches "<<endl;

//   for (i=1; i<Target.numSamples(); ++i) {
   
//   //initial match for each point is the one w/ lowest cost
//   //0th position in the prefer array
    
//     //offset into template array index
//     Target.mMatch[i]=Target.mPrefer[i][0]+(i*dist-spacing+1);
//     cout<<Target.mMatch[i]<<endl;
    
//     //index shows how far down on the prefer list we've progressed
//     //initialized to 0
//     Target.mInd[i]=0;
  
//     //initial cost of match is the sum of all least cost matches
//     cum_cost+=Target.mCost[i][Target.mPrefer[i][0]];
//   } //end i
  

//   cout<<"initial cost "<<cum_cost<<endl;
//   int *curr_ind[N], *curr_match[N];
//    for (i=1; i<Target.numSamples()-N; ++i) {
//      for (j=0; j<N; ++j) {
//        curr_ind[j]=&Target.mInd[i];
//        curr_match[j]=&Target.mPrefer[i][*(curr_ind[j])];
//      }
//      delta_min=MAX_COST;

//      //check diffeomorphisms on windows of size N
//      for (j=0; j<N; ++j) {

//        //if not diffeomorphic
//        if (Target.mMatch [j+1] >=Target.mMatch[j]) break;
       
	 
	 

// }

/*************************************************************************/
void Curve::checkMatch(Curve &Target,double s_coef,double c_coef,
		       double t_coef) {
  Array1D<int> currMatch(Target.numSamples());
  Array1D<int> bestMatch(Target.numSamples());

  int i,j,k,l,done=0;
  double cum_cost;
  double min_cost=MAX_COST;

  //initialize the prefer array in ascending order
  //this will be sorted in the matching routine below

  for (i=1; i<Target.numSamples(); ++i) 
    for (int j=0; j<numSamples(); ++j)
      Target.mPrefer[i][j]=j;
  

  for (i=1; i<Target.numSamples(); ++i) {
    Target.mCost[i][0]=MAX_COST;
    for (int j=1; j<numSamples(); ++j) {
      Target.mCost[i][j]=
	s_coef*pow((Target.mArcLength[i]-mArcLength[j]),2) +
	c_coef*pow((Target.mSpeed[i]*Target.mCurvature[i]-
		    mSpeed[j]*mCurvature[j]),2) +
	t_coef*pow((Target.mSpeed[i]*Target.mTorsion[i]-
		    mSpeed[j]*mTorsion[j]),2);
      //  cout<<Target.mCost[i][j]<<"\t";
      //if (j==numSamples()-1) cout<<endl;
    } //end j
  }
  
 int temp;
//sort the prefer array in increasing cost

 for (i=1; i<Target.numSamples(); ++i) {
   for (j=0; j<numSamples()-1; ++j) 
     for (k=0; k<numSamples()-1-j;  ++k) {
     
	if (Target.mCost[i][Target.mPrefer[i][k]] >
	    Target.mCost[i][Target.mPrefer[i][k+1]]) {
	  temp=Target.mPrefer[i][k];
	  Target.mPrefer[i][k]=Target.mPrefer[i][k+1];
	  Target.mPrefer[i][k+1]=temp;
	  
	}  //endif
      }  //end k 
 }//end i



  cout<<"initial matches "<<endl;

  for (i=1; i<Target.numSamples(); ++i) {
   
  //initial match for each point is the one w/ lowest cost
  //0th position in the prefer array
    
    //offset into template array index
    Target.mMatch[i]=Target.mPrefer[i][0];
    cout<<Target.mMatch[i]<<endl;
    
    //index shows how far down on the prefer list we've progressed
    //initialized to 0
    Target.mInd[i]=0;
  
    //initial cost of match is the sum of all least cost matches
    cum_cost+=Target.mCost[i][Target.mMatch[i]];
  } //end i
 cout<<cum_cost<<endl;
 //exit(1);

 Array1D<int> currInd(Target.numSamples());
 currInd=Target.mInd;
// cout<<"currInd "<<currInd<<endl;
 Array1D<int> oldInd(Target.numSamples());
 
 done=0; 
 int count=0;
 while (!done ) {
   count++;
   if (count%100000==0) cout<< count<<endl;
   //initial condition: "fix" diffeomorphism as best possible
   
   for (i=1; i<Target.numSamples()-1; ++i) { 
     if (Target.mMatch[i] >=Target.mMatch[i+1]) { 
       if (currInd[i+1]!=numSamples()-2) { 
 	 currInd[i+1]++;
	 Target.mMatch[i+1]=Target.mPrefer[i+1][currInd[i+1]]; 
	 i--;
       }
     }
   }
  
   
   /*
   for (i=Target.numSamples()-2; i>=1; --i) {
     if (Target.mMatch[i] >=Target.mMatch[i+1]) {
       if  (currInd[i+1]!=numSamples()-2) { 
 	 currInd[i+1]++;
	 Target.mMatch[i+1]=Target.mPrefer[i+1][currInd[i+1]];
	 i++;
       }
     }
   }
   */
   cum_cost=Target.calculateCost();
   //  if (count==1) { cout<<"after sorting "<<Target.mMatch<<endl;  }
   for (i=1; i<Target.numSamples()-1; ++i) 
     //if not diffeomorphic break 
     if (Target.mMatch[i] >=  Target.mMatch[i+1])  break;  

   //if (i==Target.numSamples()-1) cout<<" diffeo "<<cum_cost<<"\t"<<count<<endl;
   //found the diffeomorphism
   if (i==Target.numSamples()-1 && cum_cost<min_cost) { 
     min_cost=cum_cost;
     oldInd=currInd;
     cout<<" found  "<<min_cost<<"   "<<count<<endl; 
   }
     
   // cout<<currInd<<"  "<<count++<<endl;
   
   
   for (k=Target.numSamples()-1; k>=1; --k) {
     
     //if (count==1) { k=i+1;  cout<<"1st "<<k<<endl;  break;} 
     if (currInd[k]!=numSamples()-2)
       break;
   }
   reset:
   if (k!=0) {
     currInd[k]++;  //cout<<"incremented  "<<k<<endl;   
     for (l=k+1; l<Target.numSamples(); ++l)
       currInd[l]=0;
     // if (k==1) {
     double temp_cost=0.0;
     for (l=1; l<Target.numSamples(); ++l)
       temp_cost+=Target.mCost[l][Target.mPrefer[l][currInd[l]]];
     
     if (temp_cost >min_cost ){ 
       if (k==1) { cout<<" exiting " <<currInd[1]<<"\t"<<count<<endl; 
                   done=1; break;} 
       if (currInd[k-1]!=numSamples()-2) {k=k-1; goto reset; }
       else {
	 for (l=k-2; k>=1; --k) 
	   if (currInd[l]!=numSamples()-2)    break;
	 k=l; goto reset;
       } //end else
     }  //end if temp_cost
   } //end if k!=0
   
   if (k==0) { done=1;  break;  }  //cout<<"done k=0 "<<endl; 
  
   for (j=1; j<Target.numSamples(); ++j) {
     Target.mMatch[j]=Target.mPrefer[j][currInd[j]];
     //    cout<<Target.mMatch[j]<<" ";
   }
   cum_cost=Target.calculateCost();
 
    
 }//end while

 //cum_cost=Target.calculateCost();
 //for (i=1; i<Target.numSamples(); ++i)
 //cum_cost+=Target.mCost[i][Target.mMatch[i]];

 cout<<"final cost "<<min_cost<<endl;
// cout<<oldInd<<"\n"<<currInd<<endl;
 for (j=1; j<Target.numSamples(); ++j) {
   Target.mMatch[j]=Target.mPrefer[j][oldInd[j]];
 }

// cout<<Target.mMatch<<endl;

}
 //  currMatch[0]=0;

//   //initial diffeomorphism
//   for (i=1; i<Target.numSamples(); ++i) 
//     currMatch[i]=currMatch[i-1]+1;
    
//   while(!done) {
//     cum_cost=0;
//     for (i=1; i<Target.numSamples(); ++i)
//       cum_cost+=Target.mCost[i][currMatch[i]];
//     if (cum_cost<min_cost) {
//       min_cost=cum_cost;
//       bestMatch=currMatch;
//       cout<<" new cost "<<min_cost<<endl;
//     }
//     for (i=1; i<Target.numSamples(); ++i) {
    
//       //if reached the end of row, increment the previous row 
//       //and reset all the next ones
//       if (currMatch[i]==numSamples()-Target.numSamples()+i) {
// 	//cout<<"reached end of row  "<<i<<endl;
//         if (i==1) {
// 	  done=1;
// 	  cout<<" i=1 quitting "<<endl;
// 	  break;
// 	}
// 	//increment previous row where we haven't reached the end
// 	currMatch[i-1]+=1; //   if i=2  cout<<"incrementing 1 "<<endl;
// 	for (k=i; k<Target.numSamples(); ++k) 
// 	      currMatch[k]=currMatch[k-1]+1;
     
// 	break;  //don't need to check for all the other rows  
//       } //end if 
    
//     } //end for i

//     //if no end is reached increment the last row
//     if (i==Target.numSamples()) { // cout<< "incrementing last" <<endl; 
//     currMatch[Target.numSamples()-1]++; }
//   } //end while

//   cout<<" min cost "<<min_cost<<endl;
//   cout<<bestMatch<<endl;
		  
/*******************************************************************/
double Curve::calculateCost() {
  double sum=0; 
  for (int i=1; i<numSamples(); ++i)
    sum+=mCost[i][mMatch[i]];
  return sum; 
}

/***********************************************************************/  
  int Curve::isDiffeo() {
  int i;
  for(i=1; i<mNumPoints-1; ++i) 
    if (mMatch[i]>=mMatch[i+1]) break;

  if (i!=mNumPoints-1) return i;
  else return 0;
  }

/********************************************************************/
void Curve::getMatch(Curve &Target, Array1D<Point> &temp, Array1D<Point> &targ) {
  cout<<"getting match"<<endl;
  int k;
  for (int i=0; i<Target.numSamples(); ++i) {
    targ[i]=Target.mCurvePoints[i];
    k=Target.mMatch[i];
    //  cout<<k<<endl;
    temp[i]=mCurvePoints[k];
    //cout<<temp[i]<<"\t \t"<<targ[i]<<endl;
  }
  
}

/***********************************************************************/

//print function for cout
ostream & operator<< (ostream &os, Curve const &C) {
 
  C.print(os);
  return os;
}




