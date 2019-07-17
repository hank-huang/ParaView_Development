///////////////////////////////////////////////////////////////////////////
//
// IntellX, L.L.C., Proprietary
// Do not reproduce without permission in writing.
// Copyright (c) 1997 IntellX, L.L.C.
// All rights reserved
//
// File: Line.C
//
// Author: Rob Teichman
//
// Purpose: Body of Line Class
// 
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: Line.C,v 1.1 2004/11/15 04:44:08 joeh Exp $
//
// Revision History
//
// $Log: Line.C,v $
// Revision 1.1  2004/11/15 04:44:08  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:28:43  kem
// Initial revision
//
// Revision 1.2  1999/05/12 17:31:12  RAZORDB
// update
//
// Revision 1.1  1997/08/22 20:36:47  rst
// Initial revision
//
// Revision 1.2  1997/08/22 20:24:17  rst
// Cleaning up files
//
// Revision 1.1  1997/08/05 14:25:59  rst
// Initial revision
//
// Revision 1.2  1997/07/15 19:34:51  rst
// testing
//
// Revision 1.1  1997/07/15  19:29:29  rst
// Initial revision
// 
//
///////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <TDL/Line.h>

//
// copy constructor
//
Line::Line (Line const &L) :
  mStartPoint(L.mStartPoint),
  mEndPoint(L.mEndPoint),
  mNumPoints(L.mNumPoints)
{
sprintf(name,"%s",L.name);
}

//
// assigment operator
//
Line & Line::operator= (Line const &L)
{
  mNumPoints = L.mNumPoints;
  mStartPoint = L.mStartPoint;
  mEndPoint = L.mEndPoint;

  return *this;
}


//
// set the name of the line
//
void Line::setName(const char *nm)
{
if(nm) sprintf(name,"%s",nm);
else   sprintf(name,"<no-name>");
}


//
// addSegment
//
// add a segment to the line
//
void Line::addSegment (Point const start, Point const end)
{
  // store existing line
  Array1D<Point> startTemp = mStartPoint;
  Array1D<Point> endTemp = mEndPoint;

  // increment line size and reallocate
  mNumPoints++;

  mStartPoint.setDim(mNumPoints);
  mEndPoint.setDim(mNumPoints);

  // copy old segments over
  for (int i=0; i<mNumPoints-1; i++)
  {
    mStartPoint[i] = startTemp[i];
    mEndPoint[i] = endTemp[i];
  }

  // add new segment
  mStartPoint[mNumPoints-1] = start;
  mEndPoint[mNumPoints-1] = end;

}


//
// addition operators
// combine two lines
//

//
// uniary addition
//
Line & Line::operator+= (Line const &L)
{
  int i;

  // store existing points
  Array1D<Point> startTemp = mStartPoint;
  Array1D<Point> endTemp = mEndPoint;
  int oldSize = mNumPoints;

  // resize array for combined size
  mNumPoints += L.mNumPoints;
  mStartPoint.setDim(mNumPoints);
  mEndPoint.setDim(mNumPoints);

  // copy existing points into new array
  for (i=0; i<oldSize; i++)
  {
    mStartPoint[i] = startTemp[i];
    mEndPoint[i] = endTemp[i];
  }

  // copy new points into array
  for (i=0; i<L.mNumPoints; i++)
  {
    mStartPoint[oldSize+i] = L.mStartPoint[i];
    mEndPoint[oldSize+i] = L.mEndPoint[i];
  }

  return *this;
}

//
// regular addition
//
Line Line::operator+ (Line const &L) const
{
  Line temp(*this);
  temp += L;
  return temp;
}


//
// length
//
// compute cumulative length of line, sum of length of all segments
//
float Line::length() const
{
  float length = 0.0;
  Point vector;

//  cerr << "computing length" << endl;

  for (int i=0; i<mNumPoints; i++)
  {
    vector = mEndPoint[i] - mStartPoint[i];
    length += vector.norm();
  }
  return length;
}


//
// getSegment
//
// return the start and end points of the indicated segment
//
bool Line::getSegment(Point &start, Point &end, int num)
{
  if (num < 0 || num > mNumPoints)
  {
    cerr << "Invalid segment requested";
    return false;
  }

  start = mStartPoint[num];
  end = mEndPoint[num];

  return true;
}



bool Line::setSegment(Point start, Point end, int num)
{
  if (num < 0 || num > mNumPoints)
  {
    cerr << "Invalid segment requested";
    return false;
  }

  mStartPoint[num] = start;
  mEndPoint[num]   = end;

  return true;
}



//
// NOTE: not using Array1D<> istream >> operator
//       since use was not recommended
//
bool Line::read(istream &inFile)
{
inFile >> name;
inFile >> mNumPoints;

mStartPoint.setDim(mNumPoints);
mEndPoint.setDim(mNumPoints);

for(int i=0;i<mNumPoints;i++) {
	double tx,ty,tz;
	inFile >> tx >> ty >> tz;
	mStartPoint[i].set(tx,ty,tz);
	inFile >> tx >> ty >> tz;
	mEndPoint[i].set(tx,ty,tz);
	}
return(true);
}


bool Line::read(const char *fileName)
{
  ifstream inFile(fileName);

  if (!inFile) {
    cerr << "ERROR: Line::read() unable to open file " << fileName << endl;
    return false;
    }

  // call read from stream routine
  return read(inFile);
}


//
// write
//
bool Line::write(ostream &outFile) const
{
outFile << name << endl;
outFile << mNumPoints << endl;
for(int i=0;i<mNumPoints;i++) {
	double tx,ty,tz;

	outFile << mStartPoint[i].x() << " " <<  
		mStartPoint[i].y() << " " << mStartPoint[i].z() << endl;
	outFile << mEndPoint[i].x() << " " << 
		mEndPoint[i].y() << " " << mEndPoint[i].z() << endl;
	}
return(true);
}


//
// write
//
// write to a specified file
//
// returns true (1) on success and false (0) on failure
//
bool Line::write(const char *fileName) const
{
  ofstream outFile(fileName);

  if (!outFile) {
    cerr << "ERROR: Line::write() unable to open file " << fileName << endl;
    return false;
    }

  // call write to stream routine
  return write(outFile);
}




//
// print
//
// print out all the segments in a line (endpoint pairs)
//
void Line::print(ostream &os) const
{
  // first print out number of segments
  os << mNumPoints << "\n";
  // then print out segment start and end points
  for (int i=0; i<mNumPoints; i++)
    os << mStartPoint[i] << " " << mEndPoint[i] << "\n";
}

// Non-member print function
ostream & operator<< (ostream &os, Line const &L)
{
  L.print(os);
  return os;
}
