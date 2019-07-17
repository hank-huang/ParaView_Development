#ifndef TIMER_H
#define TIMER_H

// $Log: Timer.h,v $
// Revision 1.1  2004/11/15 04:44:09  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:27:54  kem
// Initial revision
//
// Revision 1.3  1999/07/09 21:05:59  RAZORDB
// Add method to get the time
//
// Revision 1.2  1999/01/20 22:12:09  RAZORDB
// change class "Timer" to "ITXTimer"
//
// Revision 1.1  1998/01/16 01:46:03  csb
// Initial revision
//
// Revision 1.1  1998/01/08 00:15:29  csb
// Initial revision
//
#include <OS.h>
#include <stream.h>
#include <iostream.h>
#include <sys/time.h>

class ITXTimer {
public:
  ITXTimer();
  ITXTimer( const char * timerName );

  ~ITXTimer();

  void reset( void );

  void start( void );
  void stop( void );

  float getTime() { return mTotalTime/1000000.0; }

  void displayResults( void );

  // Max length is 127 characters
  void setName(const char * name );

private:
  long long sample( void );

private:
  // mStart, mStop, and mTotalTime are in microseconds
  // mSquaredSum is in microseconds squared
  long long mStart;
  long long mStop;

  long long mTotalTime;
  long long mSquaredSum;
  double mMean;
  double mVariance;
  
  int       mNumSamples;
  char *    mTimerName;
};

#endif // TIMER_H


