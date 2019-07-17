///////////////////////////////////////////////////////////////////////////
//
// File: ItXWorkClass.h
//
// Author: Keith Doolittle
//
// Purpose: 
// 
///////////////////////////////////////////////////////////////////////////
#ifndef _ItXWorkClass_h_
#define _ItXWorkClass_h_
#include <OS.h>

class ItXWorkClass {

public:
  ItXWorkClass() : mActive(true)
    { work_proc = NULL; work_data = NULL; 
    sprintf(work_title,"Working"); }

  void SetTitle(char *str)
    { if(str) sprintf(work_title,str); 
    else sprintf(work_title, "Working"); }

  void SetWorkProc( void (*wrk_func)(char*,void*), void *wrk_data)
    { work_proc = wrk_func; work_data = wrk_data; }

  void assignWorkProc(ItXWorkClass *c)
    { work_proc = c->work_proc; work_data = c->work_data; }

  //
  // set dialog to active (by default)
  // or to given parameter value (0 = inactive, other = active)
  // [rst:19990720.1508MST]
  //
  void setWorkingActive(int active = 1)
    { mActive = (active != 0); }

  //
  // set dialog to inactive (no output)
  // [rst:19990720.1508MST]
  //
  void setWorkingInactive()
    { mActive = false; }

  void ShowWorking(float perc_done) const;
  void ShowWorking(char *str) const;

private:
  bool mActive;			// [rst:19990720.1509MST]
  char work_title[1024];
  void (*work_proc)(char *, void *);
  void *work_data;

};

#endif
