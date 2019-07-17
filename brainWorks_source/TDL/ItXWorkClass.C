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
// File: ItXWorkClass.C
//
// Author: Keith Doolittle
//
// Purpose: 
// 
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: ItXWorkClass.C,v 1.1 2004/11/15 04:44:08 joeh Exp $
//
// Revision History
//
// $Log: ItXWorkClass.C,v $
// Revision 1.1  2004/11/15 04:44:08  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:28:55  kem
// Initial revision
//
// Revision 1.2  1999/07/20 22:00:39  RAZORDB
// add empty work dialog capability
//
///////////////////////////////////////////////////////////////////////////

#include <iostream.h>
#include <stdio.h>
#include <stdlib.h>

#include <TDL/ItXWorkClass.h>

void ItXWorkClass::ShowWorking(char *str) const
{
  if (!mActive) return;		// [rst:19990720.1509MST]

  if (work_proc)
    work_proc(str,work_data);
  else if (str)
    cout << str << endl;
}

void ItXWorkClass::ShowWorking(float perc_done) const
{
  if (!mActive) return;		// [rst:19990720.1509MST]

  char tstring[1024];

  if (perc_done < 0.0) {
    if (work_proc)
      work_proc(NULL,work_data);
    else
      cout << "Complete." << endl;
  }
  else
  {
    sprintf(tstring,"%s: %d%% complete",
	    work_title, (int)(100.0 * perc_done));
    if (work_proc) 
      work_proc(tstring, work_data);
    else
      cout << tstring << endl;
  }
}
