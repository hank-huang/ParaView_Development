///////////////////////////////////////////////////////////////////////////
//
// IntellX, L.L.C., Proprietary
// The contents of this file comprise confidential, proprietary and/or
// trade secret information which is the property of IntellX LLC.
// Dissemination, distribution or other disclosure of this information
// is strictly prohibited without the express permission of IntellX LLC.
// Copyright (c) 1997 IntellX, L.L.C.
// All rights reserved
//
// File: ByuSurface.C
//
// Author: Rob Teichman
//
// Purpose: ByuSurface class body
// 
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: ByuSurface.C,v 1.4 2009/02/17 19:56:33 mbowers Exp $
//
// Revision History
//
// $Log: ByuSurface.C,v $
// Revision 1.4  2009/02/17 19:56:33  mbowers
// Go back to old style of byu file...
//
// Revision 1.3  2007/10/29 18:50:02  mbowers
// Read in and write out "true" byu.
//
// Revision 1.2  2006/10/24 13:25:03  mbowers
// Read in with double precision.
//
// Revision 1.1  2004/11/15 04:44:07  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:28:47  kem
// Initial revision
//
// Revision 1.5  1999/09/21 17:10:51  RAZORDB
// debugging changes
//
// Revision 1.4  1998/12/10 19:02:57  RAZORDB
// IDL/AW merge
//
// Revision 1.3  1998/02/06 16:19:32  rst
// Increasing output precision in BYU files
//
// Revision 1.2  1997/11/26 22:34:21  rst
// added check for non-triangulated surfaces
//
// Revision 1.1  1997/08/22 20:36:50  rst
// Initial revision
//
////////////////////////////////////////////////////////////////////////////

#include <iostream.h>
#include <sys/param.h>
#include <stdio.h>
#include <math.h>
#include <TDL/ByuSurface.h>
#include <TDL/SurfaceUtils.h>


//
// handle BYU format files
//

//
// read
//
// read from an input stream
//
// returns true (1) on success and false (0) on failure
//
bool ByuSurface::read(istream& inFile)
{
  int numParts,
    numConn;
  int firstPoly,
    lastPoly;
  int i;

  // first clean out previous surface
  if (!clean()) cerr << "Clean failed." << endl;

  // kwd
  ItXWorkClass::ShowWorking(0.0);

  // read in parameters from first line of file
  inFile >> numParts >> mNumVert >> mNumPoly >> numConn;
  char dummyLine[1024];
  inFile.getline(dummyLine, 1023);

  // read part polygon start and stop
  inFile >> firstPoly >> lastPoly;

  if (mVerbose)
    cerr << "Parts: " << numParts
	 << "\nVertices: " << mNumVert
	 << "\nPolygons: " << mNumPoly
	 << "\nElem in Connectivity Graph: " << numConn
	 << endl;

  // check if parameters can be handled by this program
  // in its current state of development
  if (numParts > 1)
  {
    cerr << "Too many parts, max is " << 1 << endl;

    ItXWorkClass::ShowWorking(-1.0);

    return false;
  }

  if (numParts < 1)
  {
    cerr << "Num parts specified is less than 1.  Assuming 1" << endl;
  }

  // allocate space for arrays
  mVert.setDim(mNumVert);
  mFacet.setDim(mNumPoly,3);

  if (mVerbose) cerr << "space allocated" << endl;

  if (mVerbose)
    cerr << "First Polygon: " << firstPoly
	 << "\nLast Polygon: " << lastPoly
	 << endl;


  //
  // read the list of vertices
  //
  int nevery = mNumVert/5;
  double tx, ty, tz;
  float work_perc = 0.0;

  for (i=0; i<mNumVert; i++)
  {
    if(i % nevery == 0) 
	ItXWorkClass::ShowWorking(work_perc += 0.1);


    inFile >> tx >> ty >> tz;
    mVert[i].set(tx, ty, tz);
  }
  if (mVerbose) cerr << "read vertex list" << endl;

  //
  // read the polygon list
  //
  nevery = mNumPoly/5;
  for (i=0; i<mNumPoly; i++)
  {
    if(i % nevery == 0)
	ItXWorkClass::ShowWorking(work_perc += 0.1);

    inFile >> mFacet[i][0] >> mFacet[i][1] >> mFacet[i][2];
    mFacet[i][2] *= -1;
    if (mFacet[i][2] < 0)
    {
      cerr << "Surface contains non-triangular facets (poly number = " <<
              i << ")!" << endl;

      ItXWorkClass::ShowWorking(-1.0);

      return false;
    }
    if (mFacet[i][0] == 0 || mFacet[i][1] == 0 || mFacet[i][2] == 0)
    {
      cerr << "Illegal vertex index supplied." << endl;
      cerr << "i=" << i << " vertices are " << mFacet[i][0] << " "
	   << mFacet[i][1] << " " << mFacet[i][2] << endl;

      ItXWorkClass::ShowWorking(-1.0);

      return false;
    }
    mFacet[i][0]--; mFacet[i][1]--; mFacet[i][2]--;
  }
  if (mVerbose) cerr << "read polygon list" << endl;

  if (mVerbose)
    cerr << "Read BYU object" << endl;

  ItXWorkClass::ShowWorking(-1.0);

  return true;
}

//
// read
//
// read from a specified file
//
// returns true (1) on success and false (0) on failure
//
bool ByuSurface::read(const char *fileName)
{
  ifstream inFile(fileName);

  if (!inFile)
  {
    cerr << "unable to open file " << fileName << endl;
    return false;
  }

  // call read from stream routine
  return read(inFile);
}


//
// write
//
// write out surface in BYU format
//
// returns true (1) on success and false (0) on failure
//
bool ByuSurface::write(ostream& outFile) const
{
  int i;

  ItXWorkClass::ShowWorking(0.0);

  // print out file header, consisting of
  // numParts numVert numPoly numConn
  outFile << "1 " << mNumVert << " " << mNumPoly << " "
	   << mNumPoly*3 << endl;

  // write out part start/stop polys
  outFile << "1 " << mNumPoly << endl;

  // write out list of vertices
  int op = outFile.precision(6);
  int nevery = mNumVert/5;
  float work_perc = 0.0;

  for (i=0; i<mNumVert; i++)
  {
    if(i % nevery == 0)
	ItXWorkClass::ShowWorking(work_perc += 0.1);

    outFile << mVert[i] << "\n";
  }
  outFile.precision(op);

  // write out polygon list
  nevery = mNumPoly/5;
  for (i=0; i<mNumPoly; i++)
  {
    if(i & nevery == 0)
	ItXWorkClass::ShowWorking(work_perc += 0.1);

    outFile << mFacet[i][0]+1 << " " << mFacet[i][1]+1
	     << " -" << mFacet[i][2]+1 << "\n";
  }

  ItXWorkClass::ShowWorking(-1.0);

  return true;
}

//
// write
//
// write to a specified file
//
// returns true (1) on success and false (0) on failure
//
bool ByuSurface::write(const char *fileName) const
{
  ofstream outFile(fileName);

  if (!outFile)
  {
    cerr << "unable to open file " << fileName << endl;
    return false;
  }

  // call write to stream routine
  return write(outFile);
}


//
// standard assignment operator
//
ByuSurface & ByuSurface::operator = (ByuSurface const &S)
{
  Surface::operator= (S);
  return *this;
}

//
// conversion assignment operator
//
ByuSurface & ByuSurface::operator = (Surface const &S)
{
  Surface::operator= (S);
  return *this;
}

