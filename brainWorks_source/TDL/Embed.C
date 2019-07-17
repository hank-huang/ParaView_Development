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
// File: Embed.C
//
// Author: Rob Teichman
//
// Purpose: Class for embedding a surface in a volume
// 
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: Embed.C,v 1.1 2004/11/15 04:44:07 joeh Exp $
//
// Revision History
//
// $Log: Embed.C,v $
// Revision 1.1  2004/11/15 04:44:07  joeh
// first check in of BrainWorks
//
// Revision 1.1  1999/10/01 15:29:00  kem
// Initial revision
//
// Revision 1.3  1999/09/21 17:10:53  RAZORDB
// debugging changes
//
// Revision 1.2  1999/01/07 16:25:07  RAZORDB
// Minor fixes (rst)
//
// Revision 1.1  1998/05/18 18:57:17  rst
// Initial revision
//
// Revision 1.1  1998/04/02 18:03:44  rst
// Initial revision
//
//
///////////////////////////////////////////////////////////////////////////


#include <TDL/Embed.h>

///////////////////////////////////////////////////////////////////////////
//
//  Private routines
//
///////////////////////////////////////////////////////////////////////////


//
// add new element to list
//
inline void Embed::add(int y, int z, int xl, int xr, int dy)
{
  Span *lp = new Span;

  // only add element to list of there is space and the Span
  // lies within the window
  if (y>=mWin.y0 && y<mWin.y1
      && z>=mWin.z0 && z<mWin.z1
      && xl>=mWin.x0 && xl<=xr && xr<=mWin.x1)
  {
    lp->y = y;
    lp->z = z;
    lp->xl = xl;
    lp->xr = xr;
    lp->dy = dy;
    mQue.push(lp);
  }
  else
  {
    cout << "illegal element attempted to be added to list, ignored" << endl;
  }
}

//
// get next element off of the list
//
inline void Embed::get(int &y, int &z, int &xl, int &xr, int &dy)
{
  Span *lp = mQue.front();

  y = lp->y;
  z = lp->z;
  dy = lp->dy;
  xl = lp->xl;
  xr = lp->xr;

  mQue.pop();
}


//
// push new Span on stack
//
inline void Embed::push(int y, int z, int xl, int xr, int dy)
{
  Span *sp = new Span;

  // only add element to stack if there is space on the stack
  // and the Span lies within the window
  if (y>=mWin.y0 && y<mWin.y1
      && z>=mWin.z0 && z<mWin.z1
      && xl>=mWin.x0 && xl<=xr && xr<=mWin.x1)
  {
    sp->y = y;
    sp->z = z;
    sp->xl = xl;
    sp->xr = xr;
    sp->dy = dy;
    mStk.push(sp);
  }
  else
  {
    cout << "attempted to push out of range segment: ["
	 << xl << "," << xr << "] "
	 << y << " " << z << endl;
  }
}


//
// pop Span off stack
//
inline void Embed::pop(int &y, int &z, int &xl, int &xr, int &dy)
{
  Span *sp = mStk.top();

  y = sp->y;
  z = sp->z;
  dy = sp->dy;
  xl = sp->xl;
  xr = sp->xr;

  mStk.pop();
}


//
// internal routines to support buildShell
//
//

//
// setupDeltas
//
// compute deltas to rasterize between p1 and p2
//
void Embed::setupDeltas(Deltas *d, Point const &p1, Point const &p2)
{
  Point fd;
  float adx, ady, adz, fn;

  d->pt = p1;
//   d->x = p1.x();
//   d->y = p1.y();
//   d->z = p1.z();

  fd = p2 - p1;
//   fdx = p2.x() - p1.x();
//   fdy = p2.y() - p1.y();
//   fdz = p2.z() - p1.z();

  adx = fabs(fd.x());
  ady = fabs(fd.y());
  adz = fabs(fd.z());

  if(adx > ady) {
    if(adx > adz) 
      fn = OVERSAMP * adx;
    else
      fn = OVERSAMP * adz;
  }
  else
  {
    if(ady > adz)
      fn = OVERSAMP * ady;
    else
      fn = OVERSAMP * adz;
  }

/*   d->nstride *= OVERSAMP; */
  d->nstride = (int) (fn + 0.5);
  if (d->nstride == 0)
  {
    cout << "Embed::setupDeltas() : nstride 0, forcing to 1" << endl;
    d->nstride = 1;
  }
  fn = (float)d->nstride;
  d->dt = fd/fn;
}

//
// recalcDeltas
//
// recompute deltas based on new normalizer (newn)
//
void Embed::recalcDeltas(Deltas *d, int newn)
{
  float fn = (float)d->nstride;
  float nfn = (float)newn;

  d->dt = (d->dt * fn)/nfn;
  d->nstride = newn;
}





///////////////////////////////////////////////////////////////////////////
//
//  Public routines
//
///////////////////////////////////////////////////////////////////////////

//
// embed
//
// This function calls the buildShell and seedFill routines
// to embed a surface in a volume
//
int Embed::embed(Surface &surface,
		Point &seed,
		unsigned char value,
		int wx0, int wx1,
		int wy0, int wy1,
		int wz0, int wz1)
{
  // if the default fill value is supplied, then do the
  // buildShell with a value of 127 and the seedFill
  // with a value of 255

  bool defVal = (value == 0);

  if (defVal) value = 127;

  if (surface.getVerbose())
    cout << "\nembed: calling buildShell" << endl;

  if (buildShell(surface, value) != 0)
  {
    cerr << "Error calling buildShell from embed" << endl;
    return 1;
  }

  if (defVal) value = 255;

  if (surface.getVerbose())
    cout << "\nembed: calling seedFill" << endl;

  if (seedFill(seed, value, wx0, wx1, wy0, wy1, wz0, wz1) != 0)
  {
    cerr << "Error calling seedFill from embed" << endl;
    return 1;
  }

//   if (debug)
//     cout << "\nembed: calling grow (x5)" << endl;

//   if (grow(5,25) != 0)
//   {
//     cerr << "Error calling grow" << endl;
//     return 1;
//   }

  return 0;
}
  



//
// seedFill: set the pixel at (x,y,z) and all of its 6-connected neighbors
// with the same pixel value to the new pixel value nv.
// A 6-connected neighbor is a pixel above, below, left, right, in front or
// behind a pixel.
//
//    Based on "A Seed Fill Algorithm" by Paul Heckbert,
//      "Graphics Gems", Academic Press, 1990,
//    Modified for 3D
//
int Embed::seedFill(Point seed,
		    unsigned char nv,	// new pixel value
		    int wx0, int wx1,
		    int wy0, int wy1,
		    int wz0, int wz1)
{
  int l, x1, x2, dy;
  unsigned char ov;		// old pixel value
  int cnt;

  int x = (int)seed.x(),
    y = (int)seed.y(), 
    z = (int)seed.z();

  mWin.x0 = wx0;
  mWin.y0 = wy0;
  mWin.z0 = wz0;
  mWin.x1 = wx1;
  mWin.y1 = wy1;
  mWin.z1 = wz1;
  if (mWin.x1 == 0)
    mWin.x1 = getXsize() - 1;
  if (mWin.y1 == 0)
    mWin.y1 = getYsize() - 1;
  if (mWin.z1 == 0)
    mWin.z1 = getZsize() - 1;

//   if (debug)
//     cout << "Window set to " << mWin.x0 << " " << mWin.y0 << " "
// 	 << mWin.z0 << " by " << mWin.x1 << " " << mWin.y1 << " "
// 	 << mWin.z1 << endl;

  //
  // if the seed is equal to the value to fill with, or the seed
  // is out of bounds, then don't do anything
  //
  if (x<mWin.x0 || x>mWin.x1 ||
      y<mWin.y0 || y>mWin.y1 ||
      z<mWin.z0 || z>mWin.z1)
  {
    cerr << "Seed value not valid (out of image bounds)."
	 << endl;
    cerr << "No filling done." << endl;
    return 1;
  }

  ov = (*this)[z][y][x];
  if (ov == nv)
  {
    cerr << "Seed value not valid (equal to fill value)."
	 << endl;
    cerr << "No filling done." << endl;
    return 1;
  }

//   if (debug)
//     cout << "Replacing " << (int) ov << " with " << (int) nv << endl;

  //
  // the stack is used for filling in the current plane, the list
  // (FIFO) is used for continuing to adjacent planes
  //
  // push starting points onto list
  // from the seed, progressing in each direction in y,
  // for the starting z plane and the shadow of the seed in the
  // adjacent z planes.
  //
  add(y+1,   z, x, x, -1);
  add(y  ,   z, x, x,  1);
  add(y+1, z-1, x, x, -1);
  add(y  , z-1, x, x,  1);
  add(y+1, z+1, x, x, -1);
  add(y  , z+1, x, x,  1);

  cnt = 0;
  while (!mQue.empty())
  {
    get(y, z, x1, x2, dy);
    push(y, z, x1, x2, dy);
    while (!mStk.empty())
    {
      pop(y, z, x1, x2, dy);
      if (!(cnt++ % 1000))
	cout << "." << flush;

//       if (debug) cout << "\nPopped Span: " << y << " "
// 		      << z << " " << x1 << " " << x2 << " "
// 		      << dy << endl;

      //
      // Span of adjacent scan line for x1<=x<=x2 was previously filled,
      // now explore adjacent pixels in this scan line
      //

      // fill from current pixel left
      for (x=x1; x>=mWin.x0 && (*this)[z][y][x] == ov; x--)
      {
//	if (debug) cout << " x" << x << "=" << (int) (*this)[z][y][x] << flush;
	(*this)[z][y][x] = nv;
      }

//      if (debug) cout << " out" << endl;

      if (x <= mWin.x0)
      {
	cerr << "Reached left edge in x, probable leak at y="
	     << y << ", z=" << z << ", x range (" << x1 << "," << x2 << ")" << endl;
	exit(0);
      }

      // if none were filled by above, then skip ahead
      if (x>=x1) goto skip;

      l = x+1;

      //
      // check for possible left handed u-turns
      //
      if (l<x1)
      {
	add(y, z-1, l, x1-1, dy);
	add(y, z+1, l, x1-1, dy);
	push(y-dy, z, l, x1-1, -dy);
// 	if (debug) cout << "\nleft u-turn: " << y-dy << " "
// 			<< z << " " << l << " " << x1-1 << " "
// 			<< -dy << " " << endl;
      }
      x = x1+1;

      //
      // loop through pixels on the right
      //
      do {
	// fill pixels on right
	for (; x<=mWin.x1 && (*this)[z][y][x]==ov; x++)
	{
//	  if (debug) cout << " x" << x << "=" << (int) (*this)[z][y][x];
	  (*this)[z][y][x] = nv;
	}

	if (x >= mWin.x1)
	{
	  cerr << "Reached right edge in x, probable leak at y="
	       << y << ", z=" << z << ", x range (" << x1 << "," << x2 << ")" << endl;
	  exit(0);
	}

	// continue in current direction, next span
	add(y, z-1, l, x-1, dy);
	add(y, z+1, l, x-1, dy);
	push(y+dy, z, l, x-1, dy);
// 	if (debug) cout << "\nPush S turn: " << y+dy << " "
// 			<< z << " " << l << " " << x-1 << " "
// 			<< dy << endl;

	// check for a possible right hand u-turn
	if (x>x2+1) {
	  add(y, z-1, x2+1, x-1, dy);
	  add(y, z+1, x2+1, x-1, dy);
	  push(y-dy, z, x2+1, x-1, -dy);
// 	  if (debug) cout << "\nright u-turn: " << y-dy << " "
// 			  << z << " " << x2+1 << " " << x-1 << " "
// 			  << -dy << endl;
	}

      skip:
	// move toward right edge of previous span, looking for fillable pixels
	for (x++; x<=x2 && (*this)[z][y][x]!=ov; x++);
	l = x;
      } while (x<=x2);
    }
  }
  return 0;
}



//
// buildShell
//
// The buildShell routine takes a surface and a volume and
// turns on all voxels in the volume which are intersected by
// a triangle in the surface
//
// returns 0 on success
//
// this routine, as written, requires Embed to be a
// friend class of surface.  Fix it!!!
//
int Embed::buildShell(Surface &surface,
		      unsigned char value)
{
  int f, fmax, i, j, n,
    n1, n2, n3,
    ix, iy, iz,
    maxx, maxy, maxz;
  float ovlp;
  Point p1, p2, p3, dd;
  Deltas top, left, right;
  float fx,fy,fz,
    adx, ady, adz;

  if (surface.getVerbose())
    cout << "in buildShell..." << endl;

  fmax = surface.mNumPoly;
  maxx = getXsize();
  maxy = getYsize();
  maxz = getZsize();

  if (surface.getVerbose())
    cout << "Starting loop, " << fmax << " polygons..";

  // loop through all facets in the surface
  for(f=0; f<fmax; f++)
  {
    if (f % 1000 == 0)
      cout << "." << flush;

    p1 = surface.mVert[surface.mFacet[f][0]];
    p2 = surface.mVert[surface.mFacet[f][1]];
    p3 = surface.mVert[surface.mFacet[f][2]];

//    cout << " setting up deltas for " << f << "..." << flush;

    setupDeltas(&left,p1,p3);
    setupDeltas(&right,p2,p3);

    if (left.nstride > right.nstride) 
      recalcDeltas(&right, left.nstride);
    else
      recalcDeltas(&left, right.nstride);

//    cout << " rasterizing..." << flush;

    for (i=0; i<left.nstride; i++)
    {

      dd = right.pt - left.pt;
//     ddx = right.x - left.x;
//     ddy = right.y - left.y;
//     ddz = right.z - left.z;
      adx = fabs(dd.x());
      ady = fabs(dd.y());
      adz = fabs(dd.z());
      n1 = (int)(OVERSAMP * adx + 0.5);
      n2 = (int)(OVERSAMP * ady + 0.5);
      n3 = (int)(OVERSAMP * adz + 0.5);
      n = (n1 > n2) ? ((n1 > n3) ? n1 : n3) : ((n2 > n3) ? n2 : n3);

      dd /= (float)n;
      fx = left.pt.x();
      fy = left.pt.y();
      fz = left.pt.z();

      // identify degenerate triangles, degenerate triangles wreak
      // havoc on this algorithm
      if (n < 0)
      {
	cout << "-";
	continue;
      }

      //
      // n == 0 means the triangle is smaller than one OVERSAMP
      // step, fill voxel at the location
      //
      if (n == 0)
      {
	cout << "0";
	ix = (int)(fx + 0.5);
	iy = (int)(fy + 0.5);
	iz = (int)(fz + 0.5);
	(*this)[iz][iy][ix] = value;
	continue;
      }

      for(j=0; j<n; j++, fx+=dd.x(), fy+=dd.y(), fz+=dd.z())
      {
	ix = (int)(fx + 0.5);
	iy = (int)(fy + 0.5);
	iz = (int)(fz + 0.5);

	/* skip triangles which extend outside the volume */
	if (ix < 0 || iy < 0 || iz < 0 ||
	    ix >= maxx || iy >= maxy || iz >= maxz)
	{
	  cout << "!" << flush;
	  continue;
	}

	(*this)[iz][iy][ix] = value;
      }

      // move along edge
      left.pt += left.dt;
      right.pt += right.dt;
    }
  }

  return 0;
}



//
// grow
//
// takes the volume and turns on all pixels which border
// a pixel which is already on
//
int Embed::grow(int n, unsigned char value)
{
  // make sure the n specified is valid
  if (n > 25 || n < 1)
  {
    cerr << "Invalid grow neighborhood specified, using 1" << endl;
    n = 1;
  }

  GreyVolume<unsigned char> orig;
  int rowSize = getXsize(),
    ySize = getYsize(), 
    zSize = getZsize(),
    planeSize = rowSize * getYsize(),
    fullSize = planeSize * getZsize();
  int iz, iy, ix, newval;

  unsigned char *vp;
  unsigned char *op;


  cout << "\nPass    of ";
  cout.width(2);
  cout << n << "\010\010\010\010\010\010";
  for (int iter = 0; iter < n; iter++)
  {
    cout << "\010\010";
    cout.width(2);
    cout << iter + 1 << flush;

    // orig contains a temporary copy of the original volume,
    // the area pointed to by this is actually modified
    orig = *this;

    vp = data();
    op = orig.data();

    // for each voxel in the volume, determine if it borders any non-zero
    // voxels (6-connected neighbors)
    for (iz=0; iz<zSize; iz++)
      for (iy=0; iy<ySize; iy++)
	for (ix=0; ix<rowSize; vp++, op++, ix++)
	{
	  if (*vp == 0)
	  {
	    // sum of 6 neighbors, if non-zero turn this pixel on
	    newval = 0;
	    if (ix > 0)
	      newval += *(op - 1);
	    if (ix < rowSize - 1)
	      newval += *(op + 1);
	    if (iy > 0)
	      newval += *(op - rowSize);
	    if (iy < ySize - 1)
	      newval += *(op + rowSize);
	    if (iz > 0)
	      newval += *(op - planeSize);
	    if (iz < zSize - 1)
	      newval += *(op + planeSize);
	    if (newval != 0)
	      *vp = (value ? value : (newval/6));
	  }
	}
  }

  return 0;
}
