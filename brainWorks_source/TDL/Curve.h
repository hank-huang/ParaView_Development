#ifndef __CURVE_H__
#define __CURVE_H__

//File:Curve.h
//Header file for curve class
#include <OS.h>
#include <iostream.h>
#include <stdio.h>
#define ADL_COMPLEX_DATATYPE
#include <ADL/Array1D.h>
#include <ADL/Array2D.h>
#include <TDL/Point.h>
#include <TDL/Surface.h>
#include <ADL/Matrix.h>
#include <ADL/Vector.h> 

#define MAX_COST pow(2,63)
class Curve {

 public:

  //default constructor
  Curve() : mNumPoints(0) {  }

  //copy constructor
  Curve(Curve const &C);

  //destructor
  virtual ~Curve() { }

  //assignment operator
  Curve & operator = (Curve const &C);
  Point &  getPoint(int i);

  //add point to the curve 
  void addPoint(Point  &newPoint, double &newArcLength,
		double &newSpeed, double &newCurvature,
		double &newTorsion);

  //resample curve
  void  resample(int new_numPoints, int orderPoly, Array2D<double> &result);


  //calculate speed,curvature,torsion
  void calculate_geometry(int new_numPoints,int orderPoly);
 
  //return number of points in curve
  int numSamples() const {
    return mNumPoints; 
  }

  double calculateCost();
  void setName(const char *nm);
  void getMatch(Curve &Target, Array1D<Point> &temp, Array1D<Point> &targ);

  //calculate cost between curves
  //user puts in the target curve, number of samples to resample the template 
  //and target curves, number of lmks along the curves, the order of
  //the polynomial to fit to the curves, the coefficients of arclength, 
  //curvature and torsion to be used in the matching and the neighborhood
  //size of each lmk as a percentage of the total number of samples. 
  //The others are parameters that get passed to the matching routine

  void  calculateMatchcost(Curve &Target, int numSamples, 
			   int numLmks, int polyOrder, double s_coef,
			   double c_coef, double t_coef, int ngbdSize,
			   Curve &sampledTarget, int &numNghbs, int &spacing,
			   int &dist);

 void regenerate(Curve &Target);  
  
  //find match between curves
  double matching(Curve &Target, int const &numNghbs, 
	     int const &spacing, int const &dist);
  void matching2(Curve &Target, int const &numNghbs, 
	     int const &spacing, int const &dist);
  
  //brute force solution
  void checkMatch(Curve &Target,double s_coef,double c_coef, double t_coef);

  //check if current match is diffeomorphic
  int isDiffeo();
  void print(ostream &os) const ;
  void reset();

  //print the arclength,speed,curvature,torsion
  void print_geo() const;
 
 private:
  char mName[256];
  int mNumPoints;
  Array1D<Point> mCurvePoints;
  Array1D<double> mArcLength;
  Array1D<double> mSpeed;
  Array1D<double> mCurvature;
  Array1D<double> mTorsion;
  Array1D<int> mMatch;
  Array1D<int> mInd;
  Array2D<double> mCost;
  Array2D<int> mPrefer;

  Array1D<Point> mTangent;
  Array1D<Point> mNormal;
  Array1D<Point> mBinormal;

  
};

//Non-member print function
ostream & operator<< (ostream &os, Curve const &C);

#endif // CURVE_CLASS_H



