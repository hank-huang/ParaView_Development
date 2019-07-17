//////////////////////////////////////////////////////////////////////
//
// Class functions of EMPeaks
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
// To Do
// -----
// 1. Deal with missing data, flag to run engine?
// 2. Check if data all zero
// 3. Set minimum variance?
//    - check if variance goes to zero
// 5. Assume bumps in accending order
// 6. Have constructor that takes Array1D<double>
// 7. Have constructor that takes EMPeaks. But do we really need this?
// 8. Take out printfs
//
// Revision Log
// ------------
// $Log: EMPeaks.C,v $
// Revision 1.4  2010/10/19 14:41:56  mbowers
// Increase the max iterations.
//
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
// Revision 1.1  1999/10/01 15:28:53  kem
// Initial revision
//
// Revision 1.3  1999/07/09 17:46:00  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.2  1999/05/12 17:16:34  RAZORDB
// update
//
// Revision 1.1  1997/08/01 19:48:12  csb
// Initial revision
//
// Revision 1.7  1997/05/22 16:06:09  csb
// Test to show kevin
//
// Revision 1.6  1997/05/19 23:40:24  csb
// 1. Fixed bug in EMPeaks(void) so that zero length histogram in is not assigned
//    to zero.
// 2. Wrote constructors so they call initialzers in a cleaner fashion.
//
// Revision 1.5  1997/05/19 16:33:23  csb
// 1. Move ToDo list from .H to .C
//
// Revision 1.4  1997/05/19 16:30:54  csb
// 1. Added Copy Constructor
//
// Revision 1.3  1997/05/19 16:20:50  csb
// 1. Changed SuperPositionBumps() and GaussianBump() to return by value
// 2. Moved mMaxIterations to be constant variable ConstMinimumVariance.
// 3. Removed MIN_VARIANCE and MAX_ITERATIONS macros and made constant
//    variables ConstMinimumVariance and ConstMaxIterations.
// 4. Changed abs() call to be fabs() call because the arguments are doubles
//
// Revision 1.2  1997/05/15 20:13:12  csb
// 1. Restructured Run() function so that it calls OneEMIteration().  This may
//    have slowed it down because it must new and delete memory for the Array1D's needed
//    in the function.
//
// Revision 1.1  1997/05/14 23:28:51  csb
// Initial revision
//
//
//////////////////////////////////////////////////////////////////////

#include <iostream.h>
#include <math.h>
#include <ADL/Array1D.h>
#include <TDL/EMPeaks.h>

// Constants for algorithm and gaussian generator
//const unsigned int ConstMaxIterations    = 1000; // was 5000, arbitrary
const unsigned int ConstMaxIterations    = 5000; // was 5000, arbitrary
const double       ConstMinimumVariance  = 1e-20; // arbitrary
const double       ConstConvergenceValue = 1e-10; // arbitrary

//
// Constructors
//
EMPeaks::EMPeaks(void)
  :         mNumPeaks( 0 ),
    mHistogramArray1D( 0 ),
       mWeightArray1D( 0 ),
         mMeanArray1D( 0 ),
     mVarianceArray1D( 0 ),
            mFixMeans( 0 ),
    mInitializedMeans( 0 ),
                  mEM( 0.0 ),
                 mRSS( 0.0 )
{
  mbInterpData = false;
};


EMPeaks::EMPeaks(unsigned int numpeaks,
                 Array1D<unsigned int> histogramarray1D,
		 CurveType ct)
  :        mNumPeaks( numpeaks ),
      mWeightArray1D( numpeaks ),
        mMeanArray1D( numpeaks ),
    mVarianceArray1D( numpeaks ),
            mFixMeans( 0 ),
    mInitializedMeans( 0 ),
                  mEM( 0.0 ),
                 mRSS( 0.0 )

{
  mbInterpData = false;

  // Have to initialize the histogram this way because there is no
  // Function to change an Array1D<unsigned int> into a Array1D<double>
  mHistogramArray1D.setDim( histogramarray1D.getNelm() );
  for (int i=0; i<mHistogramArray1D.getNelm(); i++)
    mHistogramArray1D[i] = (double) histogramarray1D[i];

};

// Copy Constructor
// This ensures that the appropriate constructors are called.
EMPeaks::EMPeaks(EMPeaks const & P)
  :         mNumPeaks( P.mNumPeaks ),
    mHistogramArray1D( P.mHistogramArray1D ),
       mWeightArray1D( P.mWeightArray1D ),
         mMeanArray1D( P.mMeanArray1D ),
     mVarianceArray1D( P.mVarianceArray1D ),
	    mFixMeans( P.mFixMeans),
    mInitializedMeans( P.mInitializedMeans),
                  mEM( P.mEM ),
                 mRSS( P.mRSS )
{ 
  mbInterpData = P.mbInterpData;
}


//
// Destructor
//
EMPeaks::~EMPeaks()
{
};


// User (public) functions
// Returns the number of peaks
unsigned int
EMPeaks::GetNumPeaks(void)
{
  return mNumPeaks;
}

void EMPeaks::SetInitialMeans(Array1D<double> &imeans)
{
if(imeans.getNelm() != mMeanArray1D.getNelm()) 
	cerr << "ERROR: EMPeaks::SetInitialMeans invalid number of input means" <<endl;
else	{
	mMeanArray1D      = imeans;
	mInitializedMeans = 1;
	}
}


void
EMPeaks::SetNumPeaks(unsigned int NumPeaks)
{
  mNumPeaks = NumPeaks;

  // Realloc these arrays to new size
  mWeightArray1D.setDim( mNumPeaks );
  mMeanArray1D.setDim( mNumPeaks );
  mVarianceArray1D.setDim( mNumPeaks );
  mInitializedMeans = 0;

  // Set the number of peaks to find
}

void
EMPeaks::SetHistogram(Array1D<unsigned int> const & histogram)
{
  mHistogramArray1D.setDim( histogram.getNelm() );
  for (int i=0; i<mHistogramArray1D.getNelm(); i++)
    mHistogramArray1D[i] = (double) histogram[i];
}

double
EMPeaks::OneEMIteration(Array1D<double> & XValues,
                        Array1D<double> & NewWeights,
                        Array1D<double> & NewMeans,
                        Array1D<double> & NewVariances,
                        double& RSS)
{
  int LengthOfHistogram (mHistogramArray1D.getNelm());

  Array1D<double> TempArray1D(LengthOfHistogram);
  Array1D<double> SingleBump(LengthOfHistogram); // SingleBump is
                                // G_n^k(x) introduced in equation 2
                                // in the design document .
  Array1D<double> Super(LengthOfHistogram); // Super is the superposition of 
                                // the mNumPeaks Gaussian bumps, s^k(x) introduced
                                // in equation 3 of the design document.
  Array1D<double> DataBumpSuper(LengthOfHistogram); // DataBumpSuper is 
                                //  data * bump / superpostigion of bumps,
                                // common term in equations 4, 5, and 6
  double SumDataBumpSuper;
//  double mle = 1.0;
  int i;
  
  Super = SuperPositionBumps ( XValues ); // Uses m{Weight,Mean,Variance}Array1D

    for (int PeakIndex = 0; PeakIndex<mNumPeaks; PeakIndex++) 
    {
      SingleBump = GaussianBump (XValues,1, mMeanArray1D[PeakIndex],
                                 mVarianceArray1D[PeakIndex]);
      //
      // Divide each bump by superposition of bumps and mulitply
      // by data.  If Super[i] == 0 then set SingleBump[i]/Super[i] = 1.
      // 
      for (i=0; i<LengthOfHistogram; i++) {
        double DataToUse = 0;

        DataToUse = mHistogramArray1D[i];

        if (0 != Super[i]) {
          DataBumpSuper[i] = DataToUse *  SingleBump[i] / Super[i];
        } else {
          DataBumpSuper[i] = DataToUse;
        }
      }
      SumDataBumpSuper = SumArray1D( DataBumpSuper );

      if (SumDataBumpSuper > 0) {
        // Calculated new weights, equation 4
        NewWeights[PeakIndex] = mWeightArray1D[ PeakIndex ] * SumDataBumpSuper ;
        // Calculated new means, equation 5

	if(!mFixMeans)
            NewMeans[PeakIndex] = SumArray1D(XValues*DataBumpSuper)/SumDataBumpSuper ;
        // Calculated new variances, equation 6

//KWD
//        TempArray1D = XValues - NewMeans[ PeakIndex ];
	for(int i=0;i<TempArray1D.getNelm();i++)
        	TempArray1D[i] = XValues[i] - NewMeans[ PeakIndex ];


        NewVariances[PeakIndex] = SumArray1D(TempArray1D*TempArray1D*DataBumpSuper ) / SumDataBumpSuper ;
      }
    }

    double logSum = 0.0;
    RSS = 0.0;
    double diff;
    for (int sbIx = 0; sbIx < Super.getNelm(); sbIx++)
    {
      if (sbIx != 128)
      {
        diff = mHistogramArray1D[sbIx] - Super[sbIx];
        RSS += diff * diff;
      }
      logSum += log(Super[sbIx]);
    }
    return logSum;

}


int EMPeaks::Run(Array1D<double> &XValues)
{
  int             i,LengthOfHistogram;
//  Array1D<double> XValues;


  if(0 == SumArray1D(mHistogramArray1D)) 
    return(-1);

  if(0 == mNumPeaks) 
    return(-1);

  // Actual EM algorithm runs here

  // Put the bumps where the data is

/*
  if(mbInterpData) {
      interpolateHistogram(XValues);
  } else {
      XValues.setDim(mHistogramArray1D.getNelm());
      for (i=0; i<XValues.getNelm(); i++) 
            XValues[i] = (double)i;
  }
*/

  LengthOfHistogram = mHistogramArray1D.getNelm();
  

  /////////////////////////////////////////////////////
  // Set up initial conditions for parameter vectors //
  /////////////////////////////////////////////////////

  double sum  = SumArray1D ( mHistogramArray1D );
  double mean = SumArray1D ( XValues * mHistogramArray1D ) / sum;

  Array1D<double> Xmean(LengthOfHistogram);
  for(int xx=0;xx<LengthOfHistogram;xx++)
	Xmean[xx] = XValues[xx] - mean;

  double var  = SumArray1D ( Xmean * Xmean * mHistogramArray1D ) / sum;

  double coverage = 4 * sqrt(var);
 
  mWeightArray1D   = sum / mNumPeaks;
  mVarianceArray1D = 0.5 * (coverage/mNumPeaks) * (coverage/mNumPeaks) ;

  if(!mInitializedMeans) {
      i = 0;
      for (int p = -((int)mNumPeaks-1); p <= ((int)mNumPeaks-1); p=p+2) {
        mMeanArray1D[i] = mean + p * coverage / (2 * mNumPeaks);
        i++;
      }
  }

  Array1D<double> NewWeights(mWeightArray1D);
  Array1D<double> NewMeans(mMeanArray1D);
  Array1D<double> NewVariances(mVarianceArray1D);

  double logSumResult = 0.0;
  double RSS = 0.0;

  int iteration = 0;
  double ConvergedSum = 0;
  int ConvergedFlag = 0;
  while((iteration < ConstMaxIterations)&& !ConvergedFlag) {
    iteration++;

    // Do EM iteration

    logSumResult = OneEMIteration(XValues,NewWeights,NewMeans,NewVariances,
                                                                        RSS);

    ConvergedSum = 0;
    for (int PeakIndex=0; PeakIndex<mNumPeaks; PeakIndex++) {
      ConvergedSum += fabs(mWeightArray1D[PeakIndex] - NewWeights[PeakIndex])
                      + fabs(mMeanArray1D[PeakIndex] - NewMeans[PeakIndex])
                  + fabs(mVarianceArray1D[PeakIndex] - NewVariances[PeakIndex]);
    }

    if ( ConvergedSum < ConstConvergenceValue) {
      ConvergedFlag = 1; 
    }

    // Update new arrays
    mWeightArray1D   = NewWeights;
    mMeanArray1D     = NewMeans;
    mVarianceArray1D = NewVariances;
    
  }

  cout << endl << "Iterations:  " << iteration << endl;
  cout << endl << "Sum of Log Super Result:  " << logSumResult << endl;
  cout <<          "RSS:  " << RSS << endl;
  mEM = logSumResult;
  mRSS = RSS;

  //
  // make sure return values are ordered by ascending means
  //
  for(int ii=0;ii<mNumPeaks;ii++)
    for(int jj=ii;jj<mNumPeaks;jj++)
      if(mMeanArray1D[jj] < mMeanArray1D[ii]) {
         double swp;
         swp = mMeanArray1D[ii]; 
         mMeanArray1D[ii] = mMeanArray1D[jj];
         mMeanArray1D[jj] = swp;

         swp = mWeightArray1D[ii]; 
         mWeightArray1D[ii] = mWeightArray1D[jj];
         mWeightArray1D[jj] = swp;

         swp = mVarianceArray1D[ii];
         mVarianceArray1D[ii] = mVarianceArray1D[jj];
         mVarianceArray1D[jj] = swp;
      }

  if (iteration >= ConstMaxIterations)
    return(0);
  else
    return(1);
}

// Returns means
Array1D<double> const &
EMPeaks::GetMeans(void)
{
  return mMeanArray1D;
}

// Returns variances
Array1D<double> const &
EMPeaks::GetVariances(void)
{
  return mVarianceArray1D;
}

// Returns weights 
Array1D<double> const &
EMPeaks::GetWeights(void)
{
  return mWeightArray1D;
}


// Returns bumps to user
Array1D<double>
EMPeaks::GetBump(int WhichBump, Array1D<double> XValues)
{
  return GaussianBump(XValues,
                      mWeightArray1D[ WhichBump ],
                      mMeanArray1D[ WhichBump ],
                      mVarianceArray1D[ WhichBump ]);
}

Array1D<double>
EMPeaks::GetAllBumps(Array1D<double> XValues)
{
  return SuperPositionBumps(XValues);
}



//
// Private functions to the class
//
double
EMPeaks::SumArray1D(Array1D<double> const &X)
{
  double sum = 0;
  int NumberElements = X.getNelm();

  for (int i=0; i<NumberElements; i++)
  {
    sum += X[i];
  }

  return sum;
}

Array1D<double>
EMPeaks::SuperPositionBumps(Array1D<double> const &X)
{
  int i = 0;
  Array1D<double> result( X.getNelm() );

  result = 0;
  for (i=0; i<mNumPeaks; i++)
  {
    result += GaussianBump(X,
                           mWeightArray1D[i],
                           mMeanArray1D[i],
                           mVarianceArray1D[i]);
  }

  return result;
}

Array1D<double>
EMPeaks::GaussianBump (Array1D<double> const & X,
                       double weight,
                       double mean,
                       double variance)
{

  double PI = atan(1.0) * 4.0;
  int NumberElements = X.getNelm();
  double temp;
  Array1D<double> result(NumberElements);

  if (variance < ConstMinimumVariance) 
    variance = ConstMinimumVariance;

  result = weight / sqrt( 2 * PI * variance);
  for (int i=0; i<NumberElements; i++)
  {
    temp = X[i] - mean;
    result[i] *= exp( -(temp*temp) / 2 / variance );
  }

  return result;
}


//
// Force mHistogramArray1D to be 256 values
//
void EMPeaks::interpolateHistogram(Array1D<double> &Xvalues)
{
   int i,nrange,binl,binr,n,mMinBin,mMaxBin;
   double perc,percl,percr,dbin,fbin,origsum,newsum;
   Array1D<double> nHist;

   origsum = SumArray1D(mHistogramArray1D);

   n = mHistogramArray1D.getNelm();
   for(i=0;i<n;i++)
       if(mHistogramArray1D[i])
	break;

   if(i>=n) // do nothing (zero hist)
	goto no_interp;

   if(i<n) mMinBin = i;

   for(i=n-1;i>mMinBin;i--)
       if(mHistogramArray1D[i])
	break;

   if(i<=mMinBin) mMaxBin = mMinBin;
   else        mMaxBin = i;

   // there are values in bins [mMinBin...mMaxBin]
   nrange = mMaxBin - mMinBin + 1;

   if((nrange < 2)||(nrange >= (n-1))) // too small or full already
	goto no_interp;

   dbin = (double)nrange/(float)n;

   nHist.setDim(n);
   Xvalues.setDim(n);
   
   nHist[0]   = mHistogramArray1D[mMinBin];
   nHist[n-1] = mHistogramArray1D[mMaxBin];

   Xvalues[0]   = (double)mMinBin;
   Xvalues[n-1] = (double)mMaxBin;

   for(i=1;i<n-1;i++) {
         fbin      = (double)i * dbin + (double)mMinBin;
         binl      = (int)fbin;
         binr      = binl + 1;
         percr     = fbin - (double)binl;
         percl     = 1. - percr;

         Xvalues[i] = fbin;
         nHist[i]   = percl * mHistogramArray1D[binl] +
                      percr * mHistogramArray1D[binr];
   }

   // now make sure sums of the histograms are the same
   newsum = SumArray1D(nHist);
   perc = origsum/newsum;

   for(i=0;i<n;i++)
	nHist[i] *= perc;

   newsum = SumArray1D(nHist);
  
   mHistogramArray1D = nHist;

   return;

no_interp:
// 
// No interpretation was performed
//
   n = mHistogramArray1D.getNelm();
   Xvalues.setDim(n);
   for(i=0;i<n;i++)
	Xvalues[i] = (double)i;
}

