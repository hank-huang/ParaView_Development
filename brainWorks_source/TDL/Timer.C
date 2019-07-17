
// $Log: Timer.C,v $
// Revision 1.1  2004/11/15 04:44:08  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:28:58  kem
// Initial revision
//
// Revision 1.3  1999/07/09 17:48:14  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.2  1999/01/20 22:12:11  RAZORDB
// change class "Timer" to "ITXTimer"
//
// Revision 1.1  1998/01/16 01:46:18  csb
// Initial revision
//
// Revision 1.1  1998/01/08 00:15:55  csb
// Initial revision
//

#include <string.h>
#include <TDL/Timer.h>

// Constructor
ITXTimer::ITXTimer()
{
  ITXTimer("");
}

ITXTimer::ITXTimer( const char * timerName)
{
  // cout<< "sizeof(long long) " << sizeof(long long) << endl;
  this->setName( timerName );
  reset();
}

// Destructor
ITXTimer::~ITXTimer() 
{
  delete [] mTimerName;
}

void
ITXTimer::reset()
{
  mStart = 0;
  mStop = 0;
  mTotalTime = 0;
  mSquaredSum = 0;
  mNumSamples = 0;
  mMean = 0;
  mVariance = 0;
};

long long
ITXTimer::sample()
{
  struct timeval t;
  struct timezone tz;
  gettimeofday( &t , &tz);
  return ((long long) t.tv_sec * 1000000 + (long long) t.tv_usec);
}

void
ITXTimer::start()
{
  mStart = sample();
}

void
ITXTimer::stop()
{
  // Stops current timer
  // adds the time onto the total
  mStop = sample();

  mNumSamples++;
  long long diff = mStop - mStart;
  mTotalTime  += diff;
  mSquaredSum += diff * diff;
  mMean = (double) mTotalTime / mNumSamples;
  if (mNumSamples > 1)
  {
    mVariance = ((double) mSquaredSum - (double) mNumSamples * mMean * mMean ) / 
                ((double) mNumSamples - 1 );
  }
  if ( mMean < 0 || mVariance < 0) {
    cout << mNumSamples << " " << diff << " " << mTotalTime << " " << mSquaredSum << " " << mMean << " " << mVariance << endl;
  }
}

void
ITXTimer::displayResults()
{
  if (0 == strcmp("", mTimerName))
  {
    cout << "(no name) : ";
  }
  else
  {
    cout << mTimerName << ": ";
  }
  cout << "Number of samples: " << mNumSamples;
  cout << ", Total Time (sec) " << mTotalTime / 1000000. ;
  cout << ", Mean Time (sec) " << mMean / 1000000. ;
  cout << ", Variance Time (sec) " << mVariance / 1000000. / 1000000. ;
  cout << endl;
}

void
ITXTimer::setName( const char * string)
{
  int length = strlen(string) + 1;
  if (length > 128) length = 128;

  mTimerName = new char[length];  
  if (NULL == mTimerName)
  {
    cerr << "Error mallocing memory for timer name" << endl;
    return;
  }
 else
 {
   mTimerName[0] = '\0';
 }

  strncpy(mTimerName, string, length-1);
  mTimerName[length-1] = '\0';
}

