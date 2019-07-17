/////////////////////////////////////////////////////////////////////
//
// EMPeaks class definition
//
// Class that performs parametric maximum likelihood estimation
// of point process intensities. The program accepts its input in 
// histogram form. From this data maximum likelihood estimation
// of a set of parameters is performed as specified by the user.
//
// Original program and algorithmic development by Sarang Joshi.
//
// Please keep the header.
//
// Usage:
//
// Array1D<unsigned int> histogram;
// Array1D<double> weights, means, variances;
//
// // fill histogram
//
// EMPeaks peaks(3, histogram);
// peaks.Run();
// weights = peaks.GetWeights();
// means = peaks.GetMeans();
// variances = peaks.GetVariances();
//
// Revision Log
// ------------
// $Log: EMPeaks.h,v $
// Revision 1.3  2009/05/04 15:57:39  mbowers
//
// Calculate and provide public access for EM (expectation) and RSS (difference between built up histogram and actual histogram).
//
// Revision 1.2  2009/04/30 20:00:15  mbowers
// Return the final log sum of Super after running EM.  Keep the value and provide public access to it.
//
// Revision 1.1  2004/11/15 04:44:07  joeh
// first check in of BrainWorks
//
//
// Revision ???? 2000/3/28  afternoon hyena
// Added support for Gamma curves
//
// Revision 1.1  1999/10/01 15:27:49  kem
// Initial revision
//
// Revision 1.2  1999/05/12 17:15:37  RAZORDB
// update
//
// Revision 1.1  1997/08/01 19:47:51  csb
// Initial revision
//
// Revision 1.5  1997/05/22 16:06:48  csb
// Test to show kevin
//
// Revision 1.4  1997/05/19 16:33:43  csb
// 1. Moved ToDo list to EMPeaks.C
//
// Revision 1.3  1997/05/19 16:07:17  csb
// 1. Changed SuperPositionBumps and GaussianBump to return object and not const &.
// 2. Removed mMaxIterations (made it a constant in EMPeaks.C)
// 3. Changed GetBump and GetAllBumps to return objects.
//
// Revision 1.2  1997/05/15 20:14:01  csb
// 1. Added function OneEMIteration().
//
// Revision 1.1  1997/05/14 23:28:40  csb
// Initial revision
//
//
//////////////////////////////////////////////////////////////////////

#ifndef __EMPEAKS__
#define __EMPEAKS__
#include <OS.h>
#include <ADL/Array1D.h>

enum CurveType { GaussianCurve, GammaCurve };

class EMPeaks {
 public:
  // Constructors and Destructor
  EMPeaks(void);
  EMPeaks(unsigned int numpeaks,
          Array1D<unsigned int> histogramarray1D, 
	  CurveType ctype=GaussianCurve);
  EMPeaks(EMPeaks const & P);
  ~EMPeaks(void);

  // User functions to set parameters
  unsigned int GetNumPeaks(void);
  void         SetNumPeaks(unsigned int NumPeaks);
  void         SetHistogram(Array1D<unsigned int> const & histogram);

  void 	       SetInitialMeans(Array1D<double> &initialMeans);
  inline void  SetFixMeans(int fixthem) { mFixMeans = fixthem; }


  // Main user function that must be run for the peaks to be calculated

  // return : 0 (didnt converge), 1 (success), -1 (fail)
  int  Run(Array1D<double> &Xvalues);  
  void RunGamma ();
  void RunGauss ();

  // User functions that return the answers about where the peaks are
  Array1D<double> const & GetMeans(void);      //Gaussian only
  Array1D<double> const & GetVariances(void);  //ditto
  Array1D<double> const & GetWeights(void);    //Gaussian OR gamma

  Array1D<double> const & GetTs ();      //Gamma only
  Array1D<double> const & GetAs ();      //Gamma only
  Array1D<double> const & GetBs ();      //Gamma only

  void SetInitialTs(Array1D<double> &); //Gamma only
  void SetInitialAs(Array1D<double> &);
  void SetInitialBs(Array1D<double> &);

  void setInterpolateData(bool tf) { mbInterpData = tf; }

  // Extra user functions to get pretty curves of the bumps
  Array1D<double> GetBump(int WhichBump, Array1D<double> XValues);
  Array1D<double> GetAllBumps(Array1D<double> XValues);

  //toggle between Gamma curve and Gaussian curve
  void SetToGamma () { myCurveType = GammaCurve; }
  void SetToGaussian () { myCurveType = GaussianCurve; }

  CurveType TheCurveType () { return myCurveType; }

  static Array1D<double> GaussianBump (Array1D<double> const &XValues,
                                double weight,
                                double mean,
                                double variance);
  double eM () { return mEM; }
  double rSS () { return mRSS; }

 private:
  unsigned int    mNumPeaks;           //number of peaks
  Array1D<double> mHistogramArray1D;   //histogram used to generate peaks
  Array1D<double> mWeightArray1D;      //weights (Gaussian or Gamma)
  Array1D<double> mMeanArray1D;        //means (Gaussian only)
  Array1D<double> mVarianceArray1D;    //variances (Gaussian only)

  int		  mFixMeans;           //toggles
  int		  mInitializedMeans;
  int		  mInitializedTs;
  int		  mInitializedAs;
  int		  mInitializedBs;
  // KWD
  bool            mbInterpData;
  // MRB
  double      mEM;
  double      mRSS;

  Array1D<double> mTArray1D;           //T (Gamma only)
  Array1D<double> mAArray1D;           //A (Gamma only)
  Array1D<double> mBArray1D;           //B (Gamma only)
  CurveType myCurveType;               //Gamma or Gaussian?


  // Private member functions
  double          SumArray1D (Array1D<double> const &XValues);
  Array1D<double> GammaBump (Array1D<double> const &,
			     double,
			     double,
			     double,
			     double);
  Array1D<double> SuperPositionBumps (Array1D<double> const &XValues);

  //OneEMIteration for Gaussian curves
  double          OneEMIteration (Array1D<double> & XValues,
                                  Array1D<double> & NewWeights,
                                  Array1D<double> & NewMeans,
                                  Array1D<double> & NewVariances,
                                  double& RSS);
/*
  //OneEMIteration for Gamma curves
  void            OneEMIteration (Array1D<double> &,
                                  Array1D<double> &,
                                  Array1D<double> &,
                                  Array1D<double> &);  
*/

  void interpolateHistogram(Array1D<double>&);
};


#endif /* __EMPEAKS__ */



