#ifndef __EMBED_H__
#define __EMBED_H__

///////////////////////////////////////////////////////////////////////////
//
// IntellX, L.L.C., Proprietary
// The contents of this file comprise confidential, proprietary and/or
// trade secret information which is the property of IntellX LLC.
// Dissemination, distribution or other disclosure of this information
// is strictly prohibited without the express permission of IntellX LLC.
// Copyright (c) 1998 IntellX, L.L.C.
// All rights reserved
//
// File: Embed.h
//
// Author: Rob Teichman
//
// Purpose: Class header for embedding a surface in a volume
// 
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: Embed.h,v 1.1 2004/11/15 04:44:07 joeh Exp $
//
// Revision History
//
// $Log: Embed.h,v $
// Revision 1.1  2004/11/15 04:44:07  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:27:56  kem
// Initial revision
//
// Revision 1.4  1999/09/21 21:00:18  RAZORDB
// TDL update
//
// Revision 1.3  1999/07/09 17:46:55  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.2  1998/12/18 18:14:08  RAZORDB
// Fix for 7.2.1 compiler Warning 1356
//
// Revision 1.1  1998/05/18 18:56:17  rst
// Initial revision
//
///////////////////////////////////////////////////////////////////////////
#include <OS.h>
#include <deque.h>
#include <vector.h>
#include <stack.h>
#include <ADL/GreyVolume.h>
#include <TDL/Surface.h>


class Embed : public GreyVolume<unsigned char>
{

public:

  //
  // default constructor and destructor
  //
  Embed() { OVERSAMP = 32; }
  virtual ~Embed() {}

  //
  // new volume constructor
  //
  Embed(int zsize, int ysize, int xsize) :
    GreyVolume<unsigned char>(zsize, ysize, xsize)
    {
      Volume<unsigned char>::operator=(0);
    }

  //
  // existing volume constructor
  //
  Embed(GreyVolume<unsigned char> vol) :
    GreyVolume<unsigned char>(vol) {}

  //
  // new volume and surface constructor
  // automatically do embedding
  //
  Embed(int zsize, int ysize, int xsize,
	Surface &surface,
	Point &seed) :
    GreyVolume<unsigned char>(zsize, ysize, xsize)
    {
      Volume<unsigned char>::operator=(0);
      if (embed(surface, seed) != 0)
	cerr << "Error embedding, called from constructor" << endl;
    }

  //
  // volume and surface constructor
  // automatically do embedding
  //
  Embed(GreyVolume<unsigned char> vol,
	Surface &surface,
	Point &seed) :
    GreyVolume<unsigned char>(vol)
    {
      if (embed(surface, seed) != 0)
	cerr << "Error embedding, called from constructor" << endl;
    }


  //
  // embed
  //
  // this routine takes in a volume and a surface
  // and embeds the surface in the volume, filling it with
  // the fill value.  The default fill value (0) is to create the
  // boundary (buildShell) portion at 127 and fill the interior
  // (seedFill) with 255.  The wx0, wx1 parameters specify the minimum
  // and maximum values to fill within (subvolume).  If set to zero
  // (the default) the it will fill to the volume boundaries.
  //
  // The function returns 0 for success
  //
  int embed(Surface &surface,
	    Point &seed,
	    unsigned char value = 0,
	    int wx0 = 0, int wx1 = 0,
	    int wy0 = 0, int wy1 = 0,
	    int wz0 = 0, int wz1 = 0);

  //
  // buildShell
  //
  // the buildShell routine takes a surface and a volume, and
  // fills in the voxels which intersect the volume with the
  // specified intensity (default of 127).
  //
  // The function returns 0 for success
  //
  int buildShell(Surface &surface,
		 unsigned char value = 127);

  //
  // seedFill
  //
  // the seedFill routine takes a volume and a seed and fills
  // the 6-connected region around the seed (of constant intensity
  // with the specified value.  The default fill value is 255.
  //
  // The function returns 0 for success
  //
  int seedFill(Point seed,
	       unsigned char value = 255,
	       int wx0 = 0, int wx1 = 0,
	       int wy0 = 0, int wy1 = 0,
	       int wz0 = 0, int wz1 = 0);

  //
  // grow
  //
  // takes the volume and turns on all pixels which border
  // a pixel which is already on.  Grows a neigborhood of size n.
  //
  int grow(int n = 1, unsigned char value = 0);


//////////////////////////////////////////////////////////////////////
//
// protected and private
//
//////////////////////////////////////////////////////////////////////

protected:

  //////////////////////////////////////////////////////////////////////
  //
  // seedFill supporting data structures and routines
  //
  //////////////////////////////////////////////////////////////////////

  //
  // data structure for a span
  //
  typedef struct aspan {
    short y, z,			// y,z location of segment
      xl, xr,			// extent in x
      dy;				// traversal direction (in y)
  } Span;

  //
  // data structure for the limiting window
  //
  typedef struct win {
    int x0, x1,
      y0, y1,
      z0, z1;
  } Window;

  //
  // Window for fill
  //
  Window mWin;

  //
  // stack of Span's for traversal in the current z plane
  //
  stack<Span *> mStk;

  // wrapper routines for pushing onto stack
  inline void push(int y, int z, int xl, int xr, int dy);
  inline void pop(int &y, int &z, int &xl, int &xr, int &dy);

  //
  // queue of Span's for traversal in z
  //
  queue<Span *> mQue;

  // wrapper routines for adding to queue
  inline void add(int y, int z, int xl, int xr, int dy);
  inline void get(int &y, int &z, int &xl, int &xr, int &dy);

  //////////////////////////////////////////////////////////////////////
  //
  // buildShell supporting data structures and routines
  //
  //////////////////////////////////////////////////////////////////////

  int OVERSAMP;

  //
  // deta structure for deltas used in rasterizing
  //
  typedef struct _deltas {
    Point dt;
    Point pt;
    int   nstride;
  } Deltas;

  // routines for computing deltas
  void setupDeltas(Deltas *d, Point const &p1, Point const &p2);
  void recalcDeltas(Deltas *d, int newn);

};

#endif // __EMBED_H__
