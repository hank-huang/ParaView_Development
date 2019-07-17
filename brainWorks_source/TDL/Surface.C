///////////////////////////////////////////////////////////////////////////
//
// File: Surface.C
//
// Author: Rob Teichman
//
// Purpose: Surface Class Body
//
///////////////////////////////////////////////////////////////////////////
//
// RCS Id: $Id: Surface.C,v 1.4 2011/10/24 15:36:01 mbowers Exp $
//
// Revision History
//
// $Log:
// Revision 1.16  1999/09/21 17:10:56  RAZORDB
// debugging changes
//
// Revision 1.15  1999/07/16 20:04:37  rst
// Add Polygon fliping routines
//
// Revision 1.14  1999/07/09 17:47:52  RAZORDB
// Make Changes to compile under Linux and g++
//
// Revision 1.13  1999/05/12 17:31:23  RAZORDB
// update
//
// Revision 1.12  1999/01/19 17:22:56  RAZORDB
// kwd
//
// Revision 1.9  1998/05/18 19:14:49  rst
// Surface class changes
//
// Revision 1.8  1998/05/12 20:07:46  sarang
// Add functions needed for Schiz.
//
// Revision 1.7  1998/03/05 18:02:28  rst
// inverse apply, += for surfaces
//
// Revision 1.6  1997/12/23 22:16:54  csb
// change prototype in refineTriangle
//
// Revision 1.5  1997/11/26 22:42:19  rst
// updated and added new routines
//
// Revision 1.4  1997/10/23 19:32:34  rst
// fixing and optimizing
//
// Revision 1.3  1997/10/10 19:46:42  rst
// Add neighborhood generation and optimize
//
// Revision 1.2  1997/08/26 15:46:58  kem
// Change invFieldTrans() hField argument to reference
//
// Revision 1.1  1997/08/22 20:36:46  rst
// Initial revision
//
////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include <TDL/Surface.h>
#include <TDL/SurfaceReducer.h> // for MAX_CURV_RES
#include <TDL/TransformationUtils.h>
/*
#include <TDL/ManifoldTransformation.h>
*/

#define MAX_CURV	MAX_CURV_RES

// lapack routines for curvature stuff
extern "C" {
  extern void ssytrf_(char *uplo, int *n, float *a, int *lda, int *ipiv,
		      float *work, int *lwork, int *info);
  extern void ssytrs_(char *uplo, int *n, int *nrhs, float *a, int *lda,
		      int *ipiv, float *b, int *ldb, int *info);
  extern void ssyev_(char *jobz, char *uplo, int *n, float *a, int *lda,
		     float *w, float *work, int *lwork, int *info);
};


//
// copy constructor
//
Surface::Surface (Surface &S) :
	mNumVert(S.mNumVert), mNumPoly(S.mNumPoly),
	mFacet(S.mFacet), mVert(S.mVert),
	mNorm(S.mNorm), mUNorm(S.mUNorm), mNbhd(S.mNbhd),
	mNormDirty(S.mNormDirty), mContDirty(S.mContDirty),
	mCurveDirty(S.mCurveDirty),mCurvature(S.mCurvature),
	mNbhdDepth(S.mNbhdDepth), mCurveType(S.mCurveType), mVerbose(false)
{
  sprintf(filename,"<surface>");
}


//
// geometryChanged()
//
// if topology has changed
// this will invalidate
// surface characteristics like normals,
// curvature, neighborhoods, etc.
//
void Surface::geometryChanged()
{
  verticesChanged();
  freeNeighborhoods();
}

void Surface::verticesChanged()
{
  mNorm.setDim(0);
  mUNorm.setDim(0);
  mNormDirty = true;
  mContDirty = true;
  mCurveDirty = true;
  freeCurvature();
  mCurveType = NotValid;
}


//
// addFacet
//
int Surface::addFacet (int a, int b, int c)
{
  int i,j;

  // first check if facet already exists
  for (i=0; i<mNumPoly; i++)
/*
    if ((mFacet[i][0] == a || mFacet[i][0] == b || mFacet[i][0] == c)
	&& (mFacet[i][1] == a || mFacet[i][1] == b || mFacet[i][1] == c)
	&& (mFacet[i][2] == a || mFacet[i][2] == b || mFacet[i][2] == c))
*/
if (
(mFacet[i][0] == a && mFacet[i][1] == b && mFacet[i][2] == c)||
(mFacet[i][0] == a && mFacet[i][1] == c && mFacet[i][2] == b)||
(mFacet[i][0] == b && mFacet[i][1] == a && mFacet[i][2] == c)||
(mFacet[i][0] == b && mFacet[i][1] == c && mFacet[i][2] == a)||
(mFacet[i][0] == c && mFacet[i][1] == a && mFacet[i][2] == b)||
(mFacet[i][0] == c && mFacet[i][1] == b && mFacet[i][2] == a)
)
    {
      cerr << "Duplicate facet (" << i << ")" << endl;
      return i;
    }


  // see if enough space in mFacet for new facet
  // otherwise, allocate a few more elements

  if(mFacet.getNrow() <= mNumPoly) {
    Array2D<int> tempFacet = mFacet;

    // reallocate space for new poly list
    mFacet.setDim(mNumPoly + 10, 3);

    // copy into new space
    for (i=0; i<mNumPoly; i++)
      for (j=0; j<3; j++)
	mFacet[i][j] = tempFacet[i][j];
  }

  // add new polygon
  mFacet[mNumPoly][0] = a;
  mFacet[mNumPoly][1] = b;
  mFacet[mNumPoly][2] = c;
  mNumPoly++;

  // deal with normals, etc (this is copied from clean())
  geometryChanged();

  return (mNumPoly-1);
}


//
// addVert
//
int Surface::addVert (Point const &p)
{
  int i;

  // first check if new vert is a duplicate of an existing vertex,
  // if so then return its index
  for (i=0; i<mNumVert; i++)
    if (mVert[i] == p)
      return i;

  //
  // if size not big enough, allocate space
  // for additional vertex
  //
  if(mVert.getNelm() <= mNumVert) {
    Array1D<Point> tempVert = mVert;
    mVert.setDim(mNumVert+10);

    // copy into new space
    for (i=0; i<mNumVert; i++)
      mVert[i] = tempVert[i];
  }

  // add new vertex
  mVert[mNumVert] = p;
  mNumVert++;

  geometryChanged();

  return (mNumVert-1);
}


void Surface::flipPolyOrientation(int i)
{
  int tmp;
  tmp = mFacet[i][1];
  mFacet[i][1] = mFacet[i][2];
  mFacet[i][2] = tmp;
  geometryChanged();
}

void Surface::flipOrientation()
{
  int i;
  int tmp;
  for (i=0;i<mNumPoly;i++){
    tmp = mFacet[i][1];
    mFacet[i][1] = mFacet[i][2];
    mFacet[i][2] = tmp;
  }
  geometryChanged();
}



void Surface::SetFilename(const char *nm)
{
  if(nm) sprintf(filename,"%s",nm);
}


//
// change vertex values
//
// return values:
//    0 success
//    1 invalid vertex index supplied
//
int Surface::changeVert (int n, double x, double y, double z)
{
  // make sure n is in range
  if (n<0 || n>mNumVert)
  {
    cerr << "changeVert: illegal vertex index supplied: "
	 << n << endl;
    return 1;
  }

  mVert[n].set(x,y,z);
  verticesChanged();
  return 0;

}

int Surface::changeVert (int n, Point x) {
  // make sure n is in range
  if (n<0 || n>mNumVert)
  {
    cerr << "changeVert: illegal vertex index supplied: "
	 << n << endl;
    return 1;
  }

  mVert[n] = x;
  verticesChanged();
  return 0;
}


//
// assignment operator
//
Surface & Surface::operator= (Surface const &S)
{
  if(this == &S)
	return(*this);

  mNumVert = S.mNumVert;
  mNumPoly = S.mNumPoly;
  mFacet.setDim(S.mFacet.getDim());
  mFacet = S.mFacet;
  mVert.setDim(S.mVert.getNelm());
  mVert = S.mVert;

  if (S.hasNorms())
  {
    mNorm.setDim(S.mNorm.getNelm());
    mNorm = S.mNorm;
    mNormDirty = S.mNormDirty;
  } else {
    mNormDirty = true;
  }
  if (S.hasUNorms()) {
    mUNorm.setDim(S.mUNorm.getNelm());
    mUNorm = S.mUNorm;
  }
  if (S.hasNbhd()) {
    mNbhd.setDim(S.mNbhd.getNelm());
    mNbhd = S.mNbhd;
    mNbhdDepth = S.mNbhdDepth;
  }

  if(S.hasCurvature()) {
    mCurvature.setDim(S.mCurvature.getNelm());
    mCurvature = S.mCurvature;
  }

  mContDirty = true;

  sprintf(filename,"%s", S.filename);

  return *this;
}

//
// addition
// add a surface to the current surface, vertex by vertex
//
Surface & Surface::operator+= (Surface const &S)
{

  // first check to make sure these are compatible (same number of vertices)
  if (S.mNumVert != mNumVert)
  {
    cerr << "operator+=: input Surfaces do not have same structure" << endl;
    return *this;
  }

  for (int i=0; i<mNumVert; i++)
  {
    mVert[i] += S.mVert[i];
  }

  // deal with normals, etc
  verticesChanged();

  return *this;
}


//
// combine
//
// combines two surfaces into a single structure
//
Surface & Surface::combine (Surface const &S)
{
  // save old sizes
  int oNumVert = mNumVert;
  int oNumPoly = mNumPoly;

  // next copy the arrays from this to temp
  Array2D<int> tempFacet = mFacet;
  Array1D<Point> tempVert = mVert;

  // reallocate space in Surface for combined Surface
  mNumVert += S.mNumVert;
  mNumPoly += S.mNumPoly;
  mFacet.setDim(mNumPoly, 3);
  mVert.setDim(mNumVert);

  // copy into new space
  int i, j;
  for (i=0; i<oNumVert; i++)
    mVert[i] = tempVert[i];
  for (i=oNumVert; i<mNumVert; i++)
    mVert[i] = S.mVert[i-oNumVert];
  for (i=0; i<oNumPoly; i++)
    for (j=0; j<3; j++)
      mFacet[i][j] = tempFacet[i][j];
  for (i=oNumPoly; i<mNumPoly; i++)
    for (j=0; j<3; j++)
      mFacet[i][j] = S.mFacet[i-oNumPoly][j] + oNumVert;

  // deal with normals, etc (this is copied from clean())
  geometryChanged();

  return(*this);
}


//
// genNormals
//
// Generates smoothed normals for a surface
//
void Surface::genNormals()
{
  int i;

  if (mVerbose)
    cerr << "in genNormals()" << endl;

  // allocate space and set to zero
  if (mNorm.getNelm() != mNumVert)
    mNorm.setDim(mNumVert);
  for (i=0; i<mNumVert; i++)
    mNorm[i].set(0,0,0);

  // loop through each triangle in the surface and accumulate the
  // contribution from that triangle to each of its vertices
  for (i=0; i<mNumPoly; i++)
  {
    mNorm[mFacet[i][0]] += (mVert[mFacet[i][1]]-mVert[mFacet[i][0]]).cross(
      mVert[mFacet[i][2]]-mVert[mFacet[i][0]]);
    mNorm[mFacet[i][1]] += (mVert[mFacet[i][2]]-mVert[mFacet[i][1]]).cross(
      mVert[mFacet[i][0]]-mVert[mFacet[i][1]]);
    mNorm[mFacet[i][2]] += (mVert[mFacet[i][0]]-mVert[mFacet[i][2]]).cross(
      mVert[mFacet[i][1]]-mVert[mFacet[i][2]]);
  }

  // divide each normal by 6
  // factor of 2 because cross product of vectors is actually
  //   area of a parallelogram, so half that is area of triangle
  // factor of 3 because the area of each triangle is applied in
  //   full to each of its vertices, so each normal has 3 times
  //   the correct area
  for (i=0; i<mNumVert; i++)
  {
    mNorm[i] = mNorm[i]/6.0;
  }

  mNormDirty = false;
}


//
// genUnitNormals
//
// generate unit normals
//
void Surface::genUnitNormals()
{
  double x;

  if (mVerbose)
    cerr << "in genUnitNormals()" << endl;

  // if normals have not been computed, do so
  if (!hasNorms())
    genNormals();

  // create space
  mUNorm.setDim(mNumVert);

  // generate unit normals from normals
  for (int i=0; i<mNumVert; i++)
  {
    if ((x = mNorm[i].norm()) != 0)
      mUNorm[i] = mNorm[i]/x;
    else
    {
//      cerr << "norm = 0 for vertex " << i << endl;
      mUNorm[i] = mNorm[i];
    }

    if (mVerbose)
    {
      x = mUNorm[i].norm();
      if (x > 1.01 || x < 0.99)
	cerr << "Unit normal not mag 1: " << i << " = " << x << endl;
    }
    // cerr << "Unit Normal " << i << ": " << mUNorm[i] << endl;
  }
}



//
// genNeighborhoods
//
// generate neighborhood around each vertex
//
void Surface::genNeighborhoods (int Size)
{
  int i,a,b,c;

  if((mNbhdDepth == Size)||(Size<=0))
    return;

  if (mVerbose) cout << "Building Neighborhood of size " << Size << endl;

  // clear out neighborhood structure
  freeNeighborhoods();

  mNbhd.setDim(mNumVert);

  bool need_depth = (Size != 1);
  for(i=0;i<mNumVert;i++)
    mNbhd[i].setNeedDepth(need_depth);

  mNbhdDepth = Size;

  // loop through all triangles, add each vertex to the other
  // vertices Neighborhoods.

  for (i=0; i<mNumPoly; i++) {
        a = mFacet[i][0];
        b = mFacet[i][1];
        c = mFacet[i][2];

	if(mNbhd[a].addNeighbor(b,1,true))
	   mNbhd[b].addNeighbor(a,1,false); // dont need to check for dup here
	if(mNbhd[a].addNeighbor(c,1,true))
	   mNbhd[c].addNeighbor(a,1,false); // ''
	if(mNbhd[b].addNeighbor(c,1,true))
	   mNbhd[c].addNeighbor(b,1,false); //  ''
  }

  //
  // for size > 1 add neighbors of neighbors
  //
  for (int j=1;j<Size;j++){
	Array1D<Neighborhood> bigNghbd=mNbhd;
	for(i=0;i<mNumVert;i++){
	   if(!(i%50)) ShowWorking((float)i/mNumVert);
	   int numn = mNbhd[i].numNeighbors();

	   for(int k=0;k<numn;k++) {
  		int nei = mNbhd[i].getNeighbor(k);
  		int newnumn = mNbhd[nei].numNeighbors();
		for(int l=0;l<newnumn;l++) {
			int newnei = mNbhd[nei].getNeighbor(l);
			bigNghbd[i].addNeighbor(newnei,j+1,true);
		} //end l
	   }  //end
	} //end i

	mNbhd.setDim(0);
	mNbhd = bigNghbd;
	}
}


//
// euler
//
// compute the Euler Characteristic for the surface
//
int Surface::euler()
{
  int i=0, j=0;
  Array1D<int> used(mNumVert);
  int eulerChar;

  // count number of edges
  int edges=0;
  int vertCount = 0;
  int degenerateTriangles = 0;

  // Tag used vertices
  used = 0;
  for(i=0;i<mNumPoly;i++) {
      used[mFacet[i][0]] = 1;
      used[mFacet[i][1]] = 1;
      used[mFacet[i][2]] = 1;
  }

  // make sure the neighborhood has been computed
  if (!hasNbhd(1))
  {
    if (mVerbose) cerr << "\nGenerating Neighborhoods" << endl;
    genNeighborhoods(1);
  }

  // first loop through all triangles and count degenerate ones
  if (mVerbose) cerr << "Counting degenerate triangles" << endl;
  for (j=0; j<mNumPoly; j++)
    if (mFacet[j][0] == mFacet[j][1] ||
	mFacet[j][1] == mFacet[j][2] ||
	mFacet[j][2] == mFacet[j][0])
      degenerateTriangles++;

  //
  // loop through all vertices
  // for each vertex, go through it's neighbors and count
  // up edges
  // doing this, each edge will get counted twice, once for
  // each vertex endpoint, so divide by two at end
  //
  if (mVerbose) cerr << "counting edges and vertices" << endl;
  int n;

  // make sure not to count neighbors who arent used in a poly
  for(i=0; i<mNumVert; i++) {
    if(used[i] && (n=mNbhd[i].numNeighbors()) > 0) {
      vertCount++;
      edges += n;
    }
  }
  edges /= 2;

  eulerChar = (mNumPoly - degenerateTriangles) + vertCount - edges;

  cerr << "Euler calculation: " << endl;
  cerr << "edges    = " << edges << endl;
  cerr << "vertices = " << vertCount << endl;
  cerr << "polygons = " << mNumPoly << endl;
  cerr << "Euler : " << eulerChar << endl;

  if(degenerateTriangles)
    cerr << "Surface contains " << degenerateTriangles
	 << " degenerate triangles" << endl;

  return eulerChar;
}


//
// clean
//
// Clears out all info about a surface
// should be called by all read routines before reading in
// a new surface
//
bool Surface::clean()
{
  mNumVert = 0;
  mNumPoly = 0;
  mVert.setDim(0);
  mFacet.setDim(0,3);

  geometryChanged();
  return true;
}


//
// affine
//
// Affine transformation, 3x3 matrix (A) and 3x1 vector (B), Ax+B
//
void Surface::affine(Matrix<double> const &A, Vector<double> const &B)
{
  if (A.getNrow() != 3 || A.getNcol() != 3 | B.getNelm() != 3)
  {
    cerr << "Illegal sized matrix or vector in Surface::affine" << endl;
    return;
  }

  for (int i=0; i<mNumVert; i++)
    mVert[i].set(
      mVert[i].x()*A[0][0] + mVert[i].y()*A[0][1] + mVert[i].z()*A[0][2] + B[0],
      mVert[i].x()*A[1][0] + mVert[i].y()*A[1][1] + mVert[i].z()*A[1][2] + B[1],
      mVert[i].x()*A[2][0] + mVert[i].y()*A[2][1] + mVert[i].z()*A[2][2] + B[2]);

  // invalidate normals,etc
  verticesChanged();
}

void Surface::affine(Matrix<float> const &A, Vector<float> const &B)
{
  if (A.getNrow() != 3 || A.getNcol() != 3 | B.getNelm() != 3)
  {
    cerr << "Illegal sized matrix or vector in Surface::affine" << endl;
    return;
  }

  // copy in to doubles and call double affine
  Matrix<double> dA(A.getNrow(), A.getNcol());
  Vector<double> dB(B.getNelm());
  int i,j;

  for (i=0; i<A.getNrow(); i++)
  {
    for (j=0; j<A.getNcol(); j++)
      dA[i][j] = A[i][j];
    dB[i] = B[i];
  }
  affine(dA, dB);
}


//
// scale
//
void Surface::scale(double scaleX, double scaleY, double scaleZ)
{
  // This is not implemented as a call to affine to improve performance

  for (int i=0; i<mNumVert; i++)
    mVert[i].set(mVert[i].x()*scaleX, mVert[i].y()*scaleY, mVert[i].z()*scaleZ);

  verticesChanged();
}

//
// translate
//
// Translate a surface
//
void Surface::translate(double transX, double transY, double transZ)
{
  // not implemented as a call to affine for performance
  Point t(transX, transY, transZ);
  for (int i=0; i<mNumVert; i++)
    mVert[i] += t;

  verticesChanged();
}

//
// rotate
//
// rotate a surface
// input parameters in radians
//
void Surface::rotate(double roll, double pitch, double yaw)
{
  // compute sin and cos combinations needed repeatedly
  double
    cr = cos(roll),  sr = sin(roll),
    cp = cos(pitch), sp = sin(pitch),
    cy = cos(yaw),   sy = sin(yaw);

  double
    srsp = sr * sp,
    crsp = cr * sp;

  // build rotation matrix
  Matrix<double> rot(3,3);

  rot[0][0] = cy * cp;
  rot[0][1] = sy * cp;
  rot[0][2] = -sp;

  rot[1][0] = srsp * cy - cr * sy;
  rot[1][1] = srsp * sy + cr * cy;
  rot[1][2] = sr * cp;

  rot[2][0] = crsp * cy + sr * sy;
  rot[2][1] = crsp * sy - sr * cy;
  rot[2][2] = cr * cp;

  // call affine
  Vector<double> dummy(3);
  dummy = 0;			// set translation to zero

  affine(rot, dummy);		// takes care of mNormDirty
}

//
// volume
//
// Compute volume enclosed by a closed surface using Stoke's Theorm,
// based on code from Sarang Joshi
//
// Returns the volume enclosed by a closed surface.  Returns 0 if
// there is an error (surface not closed)
//
double Surface::volume()
{
 Array1D<int> verts;
        verts.setDim(mNumVert);
        verts=0;

  // Check if surface is closed.  If not, return 0
  if (!isClosed(verts))
  {
    if (mVerbose) cerr << "Surface not closed, unable to compute volume" << endl;
    return 0.0;
  }

  // Check if normals are present.  If not, compute
  if (!hasNorms())
    this->genNormals();

  // Compute volume
  double vol = 0.0;
  for (int i=0; i<mNumVert; i++)
    vol += mVert[i].x() * mNorm[i].x() + mVert[i].y() * mNorm[i].y()
      + mVert[i].z() * mNorm[i].z();

  vol /= 3;

  return vol;
}


//
// isClosed
//
// returns true if the surface is closed.  Algorithm from Sarang Joshi.
//
bool Surface::isClosed( Array1D<int> &vertList)
{
  int i;
  Array1D<int> triangleCount;
  triangleCount.setDim(mNumVert);
  bool failed = false;

  // check if surface is closed by looking at each vertex.  If the number
  // of neighbors of each vertex is equal to the number of triangles
  // it is a part of, then the surface is closed.

  if (!hasNbhd(1))
    genNeighborhoods(1);

  // compute the number of triangles each vertex is part of
  triangleCount = 0;
  for (i=0; i<mNumPoly; i++)
  {
    triangleCount[mFacet[i][0]]++;
    triangleCount[mFacet[i][1]]++;
    triangleCount[mFacet[i][2]]++;
  }

  // check triangle count against num neighbors
  for (i=0; i<mNumVert; i++)
  {
    if (triangleCount[i] != mNbhd[i].numNeighbors())
    {
      vertList[i]=1;
      if (mVerbose)
      {
	cerr << "at vertex " << i << ", #t = " << triangleCount[i]
	     << " but numNeighbors is " << mNbhd[i].numNeighbors() << endl;
      }
      failed = true;
    }
  }

  return !failed;
}

void Surface::getBoundingBox(Point &ul, Point &lr)
{
   double xmin,ymin,zmin;
   double xmax,ymax,zmax;

   xmin = ymin = zmin =  1E20;
   xmax = ymax = zmax = -1E20;
   for(int i=0;i<mNumVert;i++) {
       double x = vertices()[i].x();
       double y = vertices()[i].y();
       double z = vertices()[i].z();
       if(x<xmin) xmin = x; if(x>xmax) xmax = x;
       if(y<ymin) ymin = y; if(y>ymax) ymax = y;
       if(z<zmin) zmin = z; if(z>zmax) zmax = z;
   }
   ul.set(xmin,ymin,zmin);
   lr.set(xmax,ymax,zmax);
}


//
// getSimpleCentroid
//
// Compute the simple centroid of a surface (average sum of vertices)
//
// Takes a pointer to the Point to be filled with the centroid
// Returns 0 on success
//
int Surface::getSimpleCentroid(Point *cent)
{
  double cx = 0, cy = 0, cz = 0;

  for (int i=0; i<mNumVert; i++) {
    cx += mVert[i].x();
    cy += mVert[i].y();
    cz += mVert[i].z();
  }

  double dnum = mNumVert;
  if(dnum!=0)
    cent->set(cx/dnum, cy/dnum, cz/dnum);
  else
    cent->set(0,0,0);

  return 0;
}




//
// getCentroid
//
// Compute the centroid of a surface
//
// Takes a pointer to the Point to be filled with the centroid
// Returns 0 on success
//
int Surface::getCentroid(Point *cent)
{
  // Check if normals are present.  If not, compute
  if (!hasNorms())
    genNormals();

  // Compute centroid using Stokes theorm
  double cx = 0, cy = 0, cz = 0;

  for (int i=0; i<mNumVert; i++)
  {
    cx += mVert[i].x() * mVert[i].x() * mNorm[i].x();
    cy += mVert[i].y() * mVert[i].y() * mNorm[i].y();
    cz += mVert[i].z() * mVert[i].z() * mNorm[i].z();
  }

  double vol = 2 * volume();

  // See if 'volume()' failed
  if(vol == 0.0)
    return(1);

  cent->set(cx/vol, cy/vol, cz/vol);

  return 0;
}


//
// removeDegenerateTriangles
//
// find and remove degenerate triangles from the surface
// (degenerate triangles are those which have two vertices
// the same)
//
// takes as an input parameter the value of "tolerance", two vertices
// with a distance of < tolerance between them are considered the same
//
void Surface::removeDegenerateTriangles(double tolerance)
{
  int i,j;
  int numRemoved = 0;
  Array1D<int> syn(mNumVert);

  // go through all vertices and remove those that are the same
  // (within tolerance)
  if (tolerance > 0)
  {
    syn[0] = 0;
    for (i=1; i<mNumVert; i++)
    {
      syn[i]=i;
      for (j=0; j<i; j++)
      {
	Point diff;
	diff = mVert[i] - mVert[j];
	if (diff.norm() < tolerance)
	{
	  syn[i] = j;
	  numRemoved++;
	}
      }
    }
  }
}


//
// fixNormals
//
// make normals consistent across the entire surface by flipping
// the vertex order for triangles whose normals are in the opposite
// direction of its neighbors
//
void Surface::fixNormals()
{
  cerr << "function not defined." << endl;
}

void Surface::getDimensions(float &xmin, float &ymin, float &zmin,
                            float &xmax, float &ymax, float &zmax)
{
   int i;
   float x,y,z;

   xmin = ymin = zmin =  1E20;
   xmax = ymax = zmax = -1E20;

   for(i=0;i<mNumVert;i++) {
      x = vertices()[i].x();
      y = vertices()[i].y();
      z = vertices()[i].z();
      if(x < xmin) xmin = x; if(x > xmax) xmax = x;
      if(y < ymin) ymin = y; if(y > ymax) ymax = y;
      if(z < zmin) zmin = z; if(z > zmax) zmax = z;
   }
}

//
// findMinMax
//
// Find the minimum and maximum values of a facet in each
// coordinate (bounding box)
//
void Surface::findMinMax(float *xMin, float *xMax,
			 float *yMin, float *yMax,
			 float *zMin, float *zMax, int facet) const
{
  // Find min and max in each dimension
  doFindMinMax(xMin, xMax,
	       mVert[mFacet[facet][0]].x(),
	       mVert[mFacet[facet][1]].x(),
	       mVert[mFacet[facet][2]].x());
  doFindMinMax(yMin, yMax,
	       mVert[mFacet[facet][0]].y(),
	       mVert[mFacet[facet][1]].y(),
	       mVert[mFacet[facet][2]].y());
  doFindMinMax(zMin, zMax,
	       mVert[mFacet[facet][0]].z(),
	       mVert[mFacet[facet][1]].z(),
	       mVert[mFacet[facet][2]].z());
}

void Surface::doFindMinMax (float *min, float *max,
			    double val0, double val1, double val2) const
{
  if (val0 > val1)
  {
    *min = (float)val1;
    *max = (float)val0;
  }
  else
  {
    *min = (float)val0;
    *max = (float)val1;
  }

  if (val2 < *min)
    *min = (float)val2;
  if (val2 > *max)
    *max = (float)val2;

}

//
// genContours
//
// returns the contours of the surface intersecting with
// each of the Sagittal, Coronal, and Axial planes
//
bool Surface::genContours(int xSlice, int ySlice, int zSlice)
{
  int i;
  float xMin, xMax, yMin, yMax, zMin, zMax;
  Point p1, p2;


  int ix = xSlice - m_minYZContour;
  int iy = ySlice - m_minXZContour;
  int iz = zSlice - m_minXYContour;

  if((ix<0)||(iy<0)||(iz<0)||(ix >= m_YZContours.getLen())||
     (iy >= m_XZContours.getLen())||(iz >= m_XYContours.getLen()))
        return(false);

  // clear out the old contours in this element
  m_YZContours[ix].clearLine();
  m_XZContours[iy].clearLine();
  m_XYContours[iz].clearLine();

  //
  // for each triangle, find the min and max in each dimension
  // if the min and max in a dimension suggest there is an intersection,
  // then compute the intersection and add to the contour
  //
  for (i=0; i<mNumPoly; i++) {

    findMinMax(&xMin, &xMax,
	       &yMin, &yMax,
	       &zMin, &zMax, i);

    if(xMin <= xSlice && xMax >= xSlice)
      if(findYZIntersect(&p1, &p2, i, xSlice))
	m_YZContours[ix].addSegment(p1, p2);

    if(yMin <= ySlice && yMax >= ySlice)
      if(findXZIntersect(&p1, &p2, i, ySlice))
	m_XZContours[iy].addSegment(p1, p2);

    if(zMin <= zSlice && zMax >= zSlice)
      if(findXYIntersect(&p1, &p2, i, zSlice))
	m_XYContours[iz].addSegment(p1, p2);

  }

  return true;
}


// get Contour routines
// getSagitalContour, getCoronalContour, getAxialContour
//
// These routines return the contour in the requested plane at the
// specified slice.
// For integer coordinate values if the contour has already been computed,
// it is just returned.  If it has not yet been computed, then it is
// computed (along with slices in other planes) and then returned.
// For non-integer coordinates the contour is computed on the fly and
// returned.
//
// RETURNS: 0 on success
//
void Surface::setupContours()
{
  int   numxy,numxz,numyz;
  int   maxxy,maxxz,maxyz;
  float xMin,yMin,zMin;
  float xMax,yMax,zMax;

  getDimensions(xMin,yMin,zMin,xMax,yMax,zMax);

  m_minXYContour = (int)zMin - 1;
  m_minXZContour = (int)yMin - 1;
  m_minYZContour = (int)xMin - 1;

  maxxy = (int)zMax + 1;
  maxxz = (int)yMax + 1;
  maxyz = (int)xMax + 1;

  numxy = maxxy - m_minXYContour + 1;
  numxz = maxxz - m_minXZContour + 1;
  numyz = maxyz - m_minYZContour + 1;

  m_XYContours.setDim(numxy);
  m_XZContours.setDim(numxz);
  m_YZContours.setDim(numyz);

  mContDirty = false;
}

void Surface::generateContours()
{
  int i,j,k,xSlice,ySlice,zSlice;
  float xMin,yMin,zMin;
  float xMax,yMax,zMax;
  Point p1,p2;

  int numxy = m_XYContours.getLen();
  int numxz = m_XZContours.getLen();
  int numyz = m_YZContours.getLen();

  for(k=0;k<numxy;k++) {
    zSlice = k + m_minXYContour;

    for(j=0;j<numxz;j++) {
      ySlice = j + m_minXZContour;

      for(i=0;i<numyz;i++) {
        xSlice = i + m_minYZContour;

        for (i=0; i<mNumPoly; i++) {

          findMinMax(&xMin, &xMax,
	             &yMin, &yMax,
	             &zMin, &zMax, i);

          if(xMin <= xSlice && xMax >= xSlice)
            if(findYZIntersect(&p1, &p2, i, xSlice))
	      m_YZContours[xSlice].addSegment(p1, p2);

          if(yMin <= ySlice && yMax >= ySlice)
            if(findXZIntersect(&p1, &p2, i, ySlice))
	      m_XZContours[ySlice].addSegment(p1, p2);

          if(zMin <= zSlice && zMax >= zSlice)
            if(findXYIntersect(&p1, &p2, i, zSlice))
	      m_XYContours[zSlice].addSegment(p1, p2);
        }
      }
    }
  }

  mContDirty = false;
cerr << "DONE." << endl;
}



Line *Surface::getYZContour(double x, double y, double z)
{
  int ix;
  int xSlice = (int)floor(x + .5);
  int ySlice = (int)floor(y + .5);
  int zSlice = (int)floor(z + .5);

  if(mContDirty)
     setupContours();

  ix = xSlice - m_minYZContour;
  if(!YZContourExists(ix))
    if(!genContours(xSlice, ySlice, zSlice))
        return(NULL);

  return(&(m_YZContours[ix]));

/*
  // if dirty bit set, clear out all contours
  if (mContDirty)
  {
    m_YZContours.setDim(512);
    m_XZContours.setDim(512);
    m_XYContours.setDim(512);
    mContDirty = false;
  }

  // check if the slice is integer valued
  if (x == xSlice)
  {
    // check if slice needs to be generated
    if (!YZContourExists(xSlice))
    {
      // generate slice, doing others as well
      if (!genContours (xSlice, ySlice, zSlice) )
	cerr << "Error generating contours!!" << endl;
    }
    *cont = m_YZContours[xSlice];
  }
  // if a non-integer coordinate then generate on the fly
  else
  {
    int i;
    double xMin, xMax,
      yMin, yMax,
      zMin, zMax;
    Point p1, p2;
    cont->clearLine();

    for (i=0; i<mNumPoly; i++)
    {
      doFindMinMax(&xMin, &xMax,
		   mVert[mFacet[i][0]].x(),
		   mVert[mFacet[i][1]].x(),
		   mVert[mFacet[i][2]].x());

      if (xMin <= x && xMax >= x)
	if (findYZIntersect(&p1, &p2, i, x))
	  cont->addSegment(p1, p2);
    }
  }

  return 0;
*/
}

Line *Surface::getXZContour(double x, double y, double z)
{
  int iy;
  int xSlice = (int)floor(x + .5);
  int ySlice = (int)floor(y + .5);
  int zSlice = (int)floor(z + .5);

  if(mContDirty)
     setupContours();

  iy = ySlice - m_minXZContour;
  if(!XZContourExists(iy))
    if(!genContours(xSlice, ySlice, zSlice))
        return(NULL);

  return(&(m_XZContours[iy]));

/*
  // if dirty bit set, clear out all contours
  if (mContDirty)
  {
    m_YZContours.setDim(512);
    m_XZContours.setDim(512);
    m_XYContours.setDim(512);
    mContDirty = false;
  }


  // check if the slice is integer valued
  if (y == ySlice)
  {
    // check if slice needs to be generated
    if (!XZContourExists(ySlice))
    {
      // generate slice, doing others as well
      if (!genContours (xSlice, ySlice, zSlice) )
	cerr << "Error generating contours!!" << endl;
    }
    *cont = m_XZContours[ySlice];
  }
  // if a non-integer coordinate then generate on the fly
  else
  {
    int i;
    double xMin, xMax,
      yMin, yMax,
      zMin, zMax;
    Point p1, p2;
    cont->clearLine();

    for (i=0; i<mNumPoly; i++)
    {
      doFindMinMax(&yMin, &yMax,
		   mVert[mFacet[i][0]].y(),
		   mVert[mFacet[i][1]].y(),
		   mVert[mFacet[i][2]].y());

      if (yMin <= y && yMax >= y)
	if (findXZIntersect(&p1, &p2, i, y))
	{
	  cont->addSegment(p1, p2);
//	  cout << "addCoronalSeg: " << p1 << p2 << endl;
	}
    }
  }

  return 0;
*/
}

Line *Surface::getXYContour(double x, double y, double z)
{
  int iz;
  int xSlice = (int)floor(x + .5);
  int ySlice = (int)floor(y + .5);
  int zSlice = (int)floor(z + .5);

  if(mContDirty)
     setupContours();

  iz = zSlice - m_minXYContour;
  if(!XYContourExists(iz))
    if(!genContours(xSlice, ySlice, zSlice))
        return(NULL);

  return(&(m_XYContours[iz]));

/*
  // if dirty bit set, clear out all contours
  if (mContDirty) {
    m_YZContours.setDim(512);
    m_XZContours.setDim(512);
    m_XYContours.setDim(512);
    mContDirty = false;
  }


  // check if the slice is integer valued
  if (z == zSlice)
  {
    // check if slice needs to be generated
    if (!XYContourExists(zSlice)) {

      // generate slice, doing others as well
      if (!genContours (xSlice, ySlice, zSlice) )
	cerr << "Error generating contours!!" << endl;
    }
    *cont = m_XYContours[zSlice];
  }
  // if a non-integer coordinate then generate on the fly
  else
  {
    int i;
    double xMin, xMax,
      yMin, yMax,
      zMin, zMax;
    Point p1, p2;
    cont->clearLine();

    for (i=0; i<mNumPoly; i++)
    {
      doFindMinMax(&zMin, &zMax,
		   mVert[mFacet[i][0]].z(),
		   mVert[mFacet[i][1]].z(),
		   mVert[mFacet[i][2]].z());

      if (zMin <= z && zMax >= z)
	if (findXYIntersect(&p1, &p2, i, z))
	{
	  cont->addSegment(p1, p2);
//	  cout << "addSeg: " << p1 << ", " << p2 << endl;
	}
    }
  }

  return 0;
*/
}


//
// findXYIntersect, findXZIntersect, findYZIntersect
//
// find the intersecting segment endpoints of a triangle with
// one of the primary axis planes.
//
// Returns true if an intersecting segment is found, false if
// non-intersecting
//

//
// findXYIntersect
//
bool Surface::findXYIntersect(Point *p1, Point *p2,
				 int facet, double zVal) const
{
  double factor;
  bool p1Set = false;
  Point tempPoint;

  // check each of the three segments for intersection,
  // no more than two will intersect

  // first pair, 0 and 1
  factor = (zVal - mVert[mFacet[facet][0]].z()) /
    (mVert[mFacet[facet][1]].z() - mVert[mFacet[facet][0]].z());
  if (factor <= 1 && factor >= 0)
  {
    p1->set(factor *
	    (mVert[mFacet[facet][1]].x() - mVert[mFacet[facet][0]].x())
	    + mVert[mFacet[facet][0]].x(),
	    factor *
	    (mVert[mFacet[facet][1]].y() - mVert[mFacet[facet][0]].y())
	    + mVert[mFacet[facet][0]].y(),
	    zVal);
    p1Set = true;
  }

  // second pair, 1 and 2
  factor = (zVal - mVert[mFacet[facet][1]].z()) /
    (mVert[mFacet[facet][2]].z() - mVert[mFacet[facet][1]].z());
  if (factor <= 1 && factor >= 0)
  {
    if (!p1Set)
    {
      p1->set(factor
	      * (mVert[mFacet[facet][2]].x() - mVert[mFacet[facet][1]].x())
	      + mVert[mFacet[facet][1]].x(),
	      factor *
	      (mVert[mFacet[facet][2]].y() - mVert[mFacet[facet][1]].y())
	      + mVert[mFacet[facet][1]].y(),
	      zVal);
      p1Set = true;
    }
    else
    {

      p2->set(factor *
	      (mVert[mFacet[facet][2]].x() - mVert[mFacet[facet][1]].x())
	      + mVert[mFacet[facet][1]].x(),
	      factor *
	      (mVert[mFacet[facet][2]].y() - mVert[mFacet[facet][1]].y())
	      + mVert[mFacet[facet][1]].y(),
	      zVal);
      if (*p2 != *p1)
	return true;
//      else
//	cerr << "degenerate segment, using next vertex" << endl;
    }
  }

  // third pair, 0 and 2
  factor = (zVal - mVert[mFacet[facet][0]].z()) /
    (mVert[mFacet[facet][2]].z() - mVert[mFacet[facet][0]].z());
  if (factor <= 1 && factor >= 0)
    if (!p1Set)
    {
      cerr << "only one intersection!!" << endl;
    }
    else
    {
      p2->set(factor *
	      (mVert[mFacet[facet][2]].x() - mVert[mFacet[facet][0]].x())
	      + mVert[mFacet[facet][0]].x(),
	      factor *
	      (mVert[mFacet[facet][2]].y() - mVert[mFacet[facet][0]].y())
	      + mVert[mFacet[facet][0]].y(),
	      zVal);
      if (*p2 != *p1)
	return true;
//      else
//	cerr << "last combo degenerate, ignoring segment" << endl;
    }

  return false;
}

//
// findXZIntersect
//
bool Surface::findXZIntersect(Point *p1, Point *p2,
				   int facet, double yVal) const
{
  double factor;
  bool p1Set = false;

  // check each of the three segments for intersection,
  // no more than two will intersect

  // first pair, 0 and 1
  factor = (yVal - mVert[mFacet[facet][0]].y()) /
    (mVert[mFacet[facet][1]].y() - mVert[mFacet[facet][0]].y());
  if (factor <= 1 && factor >= 0)
  {
    p1->set(factor *
	    (mVert[mFacet[facet][1]].x() - mVert[mFacet[facet][0]].x())
	    + mVert[mFacet[facet][0]].x(),
	    yVal,
	    factor *
	    (mVert[mFacet[facet][1]].z() - mVert[mFacet[facet][0]].z())
	    + mVert[mFacet[facet][0]].z());
    p1Set = true;
  }

  // second pair, 1 and 2
  factor = (yVal - mVert[mFacet[facet][1]].y()) /
    (mVert[mFacet[facet][2]].y() - mVert[mFacet[facet][1]].y());
  if (factor <= 1 && factor >= 0)
  {
    if (!p1Set)
    {
      p1->set(factor
	      * (mVert[mFacet[facet][2]].x() - mVert[mFacet[facet][1]].x())
	      + mVert[mFacet[facet][1]].x(),
	      yVal,
	      factor *
	      (mVert[mFacet[facet][2]].z() - mVert[mFacet[facet][1]].z())
	      + mVert[mFacet[facet][1]].z());
      p1Set = true;
    }
    else
    {
      p2->set(factor *
	      (mVert[mFacet[facet][2]].x() - mVert[mFacet[facet][1]].x())
	      + mVert[mFacet[facet][1]].x(),
	      yVal,
	      factor *
	      (mVert[mFacet[facet][2]].z() - mVert[mFacet[facet][1]].z())
	      + mVert[mFacet[facet][1]].z());
      if (*p2 != *p1)
	return true;
//      else
//	cerr << "degenerate segment, using next vertex" << endl;
    }
  }

  // third pair, 0 and 2
  factor = (yVal - mVert[mFacet[facet][0]].y()) /
    (mVert[mFacet[facet][2]].y() - mVert[mFacet[facet][0]].y());
  if (factor <= 1 && factor >= 0)
    if (!p1Set)
    {
      cerr << "only one intersection!!" << endl;
    }
    else
    {
      p2->set(factor *
	      (mVert[mFacet[facet][2]].x() - mVert[mFacet[facet][0]].x())
	      + mVert[mFacet[facet][0]].x(),
	      yVal,
	      factor *
	      (mVert[mFacet[facet][2]].z() - mVert[mFacet[facet][0]].z())
	      + mVert[mFacet[facet][0]].z());
      if (*p2 != *p1)
	return true;
//      else
//	cerr << "last combo degenerate, ignoring segment" << endl;
    }

  return false;
}

//
// findYZIntersect
//
bool Surface::findYZIntersect(Point *p1, Point *p2,
				    int facet, double xVal) const
{
  double factor;
  bool p1Set = false;

  // check each of the three segments for intersection,
  // no more than two will intersect

  // first pair, 0 and 1
  factor = (xVal - mVert[mFacet[facet][0]].x()) /
    (mVert[mFacet[facet][1]].x() - mVert[mFacet[facet][0]].x());
  if (factor <= 1 && factor >= 0)
  {
    p1->set(xVal,
	    factor *
	    (mVert[mFacet[facet][1]].y() - mVert[mFacet[facet][0]].y())
	    + mVert[mFacet[facet][0]].y(),
	    factor *
	    (mVert[mFacet[facet][1]].z() - mVert[mFacet[facet][0]].z())
	    + mVert[mFacet[facet][0]].z());
    p1Set = true;
  }

  // second pair, 1 and 2
  factor = (xVal - mVert[mFacet[facet][1]].x()) /
    (mVert[mFacet[facet][2]].x() - mVert[mFacet[facet][1]].x());
  if (factor <= 1 && factor >= 0)
  {
    if (!p1Set)
    {
      p1->set(xVal,
	      factor *
	      (mVert[mFacet[facet][2]].y() - mVert[mFacet[facet][1]].y())
	      + mVert[mFacet[facet][1]].y(),
	      factor *
	      (mVert[mFacet[facet][2]].z() - mVert[mFacet[facet][1]].z())
	      + mVert[mFacet[facet][1]].z());
      p1Set = true;
    }
    else
    {
      p2->set(xVal,
	      factor *
	      (mVert[mFacet[facet][2]].y() - mVert[mFacet[facet][1]].y())
	      + mVert[mFacet[facet][1]].y(),
	      factor *
	      (mVert[mFacet[facet][2]].z() - mVert[mFacet[facet][1]].z())
	      + mVert[mFacet[facet][1]].z());
      if (*p2 != *p1)
	return true;
//      else
//	cerr << "degenerate segment, using next vertex" << endl;
    }
  }

  // third pair, 0 and 2
  factor = (xVal - mVert[mFacet[facet][0]].x()) /
    (mVert[mFacet[facet][2]].x() - mVert[mFacet[facet][0]].x());
  if (factor <= 1 && factor >= 0)
    if (!p1Set)
    {
      cerr << "only one intersection!!" << endl;
    }
    else
    {
      p2->set(xVal,
	      factor *
	      (mVert[mFacet[facet][2]].y() - mVert[mFacet[facet][0]].y())
	      + mVert[mFacet[facet][0]].y(),
	      factor *
	      (mVert[mFacet[facet][2]].z() - mVert[mFacet[facet][0]].z())
	      + mVert[mFacet[facet][0]].z());
      if (*p2 != *p1)
	return true;
//      else
//	cerr << "last combo degenerate, ignoring segment" << endl;
    }

  return false;
}

//
// YZContourExists, coronalControuExists, XYContourExists
//
// functions to determine if the requested slice already has been
// computed
//
bool Surface::YZContourExists(int slice) const
{
  if(!mContDirty && (slice >= 0) && (slice < m_YZContours.getLen()))
    return(m_YZContours[slice].numSegments() == 0 ? false : true);
  else
    return(false);
}

bool Surface::XZContourExists(int slice) const
{
  if(!mContDirty && (slice >= 0) && (slice < m_XZContours.getLen()))
    return(m_XZContours[slice].numSegments() == 0 ? false : true);
  else
    return(false);
}

bool Surface::XYContourExists(int slice) const
{
  if(!mContDirty && (slice >= 0) && (slice < m_XYContours.getLen()))
    return(m_XYContours[slice].numSegments() == 0 ? false : true);
  else
    return(false);
}

//
// temporary functions needed by slicer
// May not be maintained in the future
//
// find minimum x value of the surface
//
double Surface::findXMin() const
{
  int i;
  double min;
  if (mNumVert == 0) {
    return(0);
  }
  min = mVert[0].x();
  for (i=1; i<mNumVert; i++)
    if (mVert[i].x() < min)
      min = mVert[i].x();

  return min;
}

//
// find maximum x value of the surface
//
double Surface::findXMax() const
{
  int i;
  double max;
  if (mNumVert == 0) {
    return(0);
  }
  max = mVert[0].x();
  for (i=1; i<mNumVert; i++)
    if (mVert[i].x() > max)
      max = mVert[i].x();

  return max;
}

//
// find minimum y value of the surface
//
double Surface::findYMin() const
{
  int i;
  double min;
  if (mNumVert == 0) {
    return(0);
  }
  min = mVert[0].y();
  for (i=1; i<mNumVert; i++)
    if (mVert[i].y() < min)
      min = mVert[i].y();

  return min;
}

//
// find maximum y value of the surface
//
double Surface::findYMax() const
{
  int i;
  double max;
  if (mNumVert == 0) {
    return(0);
  }
  max = mVert[0].y();
  for (i=1; i<mNumVert; i++)
    if (mVert[i].y() > max)
      max = mVert[i].y();

  return max;
}

//
// find minimum z value of the surface
//
double Surface::findZMin() const
{
  int i;
  double min;
  if (mNumVert == 0) {
    return(0);
  }
  min = mVert[0].z();
  for (i=1; i<mNumVert; i++)
    if (mVert[i].z() < min)
      min = mVert[i].z();

  return min;
}

//
// find maximum z value of the surface
//
double Surface::findZMax() const
{
  int i;
  double max;
  if (mNumVert == 0) {
    return(0);
  }
  max = mVert[0].z();
  for (i=1; i<mNumVert; i++)
    if (mVert[i].z() > max)
      max = mVert[i].z();

  return max;
}


///////////////////////////////////////////////////////////////////////////
//
// applyFieldTrans
//
// apply an h-field to each vertex in a surface.
// NOTE: this corresponds to mapping backwards, from the patient to the
// atlas!
//
// return values:
//  0  success
//  1  error
///////////////////////////////////////////////////////////////////////////
int Surface::applyFieldTrans(const Array3D<float> &hFields)
{

  // deform Surface
//  cout << "Deforming Surface.." << flush;
  ShowWorking(0.0);

  int i, c, d = 0;
  int nx, ny, nz;

  // first seperate h field into x y and z components for trilinear
  nx = hFields.getXsize()/3;
  ny = hFields.getYsize();
  nz = hFields.getZsize();
  Array3D<float> XhField(nz, ny, nx),
    YhField(nz, ny, nx), ZhField(nz, ny, nx);

  float const * hPtr = hFields.data();
  float * xPtr = XhField.data();
  float * yPtr = YhField.data();
  float * zPtr = ZhField.data();
  for (i=0; i<nx*ny*nz; i++)
  {
    *xPtr++ = *hPtr++;
    *yPtr++ = *hPtr++;
    *zPtr++ = *hPtr++;
  }

  // iterate through all vertices
  Array1D<double> newVert(3);
  for (i=0; i<mNumVert; i++)
  {
    // print status dots
    if ((c = (int)(i*100./mNumVert)) == d)
    {
      //     cout << c << "%.." << flush;
      ShowWorking((float)c/100.0);
      d = d + 10;
    }

    // compute new vertex
//     TransformationUtils::trilinear(&newVert, XhField, YhField, ZhField,
// 	      mVert[i].z(), mVert[i].y(), mVert[i].x(), 0.0f);
    newVert[0] = TransformationUtils::trilinear(XhField,
		 mVert[i].z(), mVert[i].y(), mVert[i].x(), mVert[i].x());
    newVert[1] = TransformationUtils::trilinear(YhField,
		 mVert[i].z(), mVert[i].y(), mVert[i].x(), mVert[i].y());
    newVert[2] = TransformationUtils::trilinear(ZhField,
		 mVert[i].z(), mVert[i].y(), mVert[i].x(), mVert[i].z());

    mVert[i].set(newVert[0], newVert[1], newVert[2]);
  }

  verticesChanged();
  return 0;
}



///////////////////////////////////////////////////////////////////////////
//
// fieldTrans
//
// deforms the surface according to the inverse of the transformation
// specified by the given H field
//
///////////////////////////////////////////////////////////////////////////
bool Surface::invFieldTrans(const Array3D<float> &hFields)
{
  // make sure surface is valid

  // deform Surface
//  cout << "Deforming Surface.." << flush;
  ShowWorking(0.0);

  int c,
    d = 0;

  for (int i=0; i<mNumVert; i++)
  {
    if ((c = (int)(i*100./mNumVert)) == d) {
      ShowWorking((float)c/100.0);
      d = d + 10;
    }
    if (! mVert[i].invFieldTrans(hFields)) {
      cerr << "\n\nError inverting vertex " << i << endl;
      return false;
    }

//     if (mVerbose) {
//       if ((i%100) == 0) {
// 	cerr << "Vertices processed: " << i << endl;
// 	system("date");
//       }
//     }
  }

  verticesChanged();

  //cout << "Complete" << endl;
  return true;
}

///////////////////////////////////////////////////////////////////////////
//
// innerproduct
// Compute and innerproduct of one surface with an orthonormal field.
// The orthonormal field is represented as a surface as well.
//
//
//
double Surface::innerproduct(Surface &basesurf,Surface const &EigenSurface)
{
  // Check if normals are present.  If not, compute it.
  if (!basesurf.hasNorms())
    basesurf.genNormals();

  // Compute innerproduct
  double innerproduct = 0.0;
  double dmu;
  for (int i=0; i<mNumVert; i++) {
    dmu = sqrt((double)(basesurf.mNorm[i].x()*basesurf.mNorm[i].x() +
			basesurf.mNorm[i].y()*basesurf.mNorm[i].y() +
			basesurf.mNorm[i].z()*basesurf.mNorm[i].z()));

    innerproduct +=
      ((mVert[i].x() - basesurf.mVert[i].x())*EigenSurface.mVert[i].x()  +
       (mVert[i].y() - basesurf.mVert[i].y())*EigenSurface.mVert[i].y()  +
       (mVert[i].z() - basesurf.mVert[i].z())*EigenSurface.mVert[i].z() )*dmu;
  }


  return innerproduct;
}


/*
Array1D<double> Surface::ComputeEigenSurf(SurfaceSet &inSurfs,
					  SurfaceSet *EigenSurfaces)
{
  int i,j;
  int numsurfs,numpoints;

  numsurfs = inSurfs.num();
  numpoints = mNumVert;

  // Now Build the Matrix
  // Now Compute the SVD!!
  // Problem ADL only computes the Complete SVD!!
  // That is too big for this.
  // We need the incomplete SVD  I am going to USE the LAPACK
  // routine directlly!!
  double *A;
  double *S;
  double *U;
  double *VT;

  A = (double *)malloc(numsurfs*numpoints*3*sizeof(double));
  S = (double *)malloc(numsurfs*sizeof(double));
  U = (double *)malloc(numsurfs*numpoints*3*sizeof(double));
  VT = (double *)malloc(numsurfs*numpoints*3*sizeof(double));


  double dmu;

  genNormals();

  for (i=0;i<numpoints;i++){
    dmu = sqrt((double) mNorm[i].x()*mNorm[i].x()+
	       mNorm[i].y()*mNorm[i].y()+
	       mNorm[i].z()*mNorm[i].z());
    dmu = sqrt(dmu);
    for (j=0;j<numsurfs;j++){
      *(A+j*3*numpoints+i*3) = ((inSurfs.getSurface(j))->mVert[i].x() - mVert[i].x())*dmu;
      *(A+j*3*numpoints+i*3+1) = ((inSurfs.getSurface(j))->mVert[i].y() - mVert[i].y())*dmu;
      *(A+j*3*numpoints+i*3+2) = ((inSurfs.getSurface(j))->mVert[i].z() - mVert[i].z())*dmu;
    }
  }

  // Now compute the SVD of A;
  char JOBU,JOBVT;
  int INFO;
  int  LDA, LDU, LDVT, LWORK, M, N;
  double *WORK;

  M = numpoints*3;
  N = numsurfs;
  LDA = M;
  LDU = M;
  JOBU = 'S';
  JOBVT = 'N';
  LDVT = N;
  LWORK = 3*M;

  WORK = (double *)malloc(LWORK*sizeof(double));

// dgesvd_((char const &) JOBU, (char const &) JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, &INFO);

  Lapack::xgesvd(JOBU,JOBVT,M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK,&INFO);

  printf("INFO after svd = %d\n",INFO);

  // Now put the Eigen Values in place.

  Array1D<double> EigenValues;

  EigenValues.setDim(numsurfs);

  for(i=0;i<numsurfs;i++){
    EigenValues[i] = (S[i]*S[i])/numsurfs;
  }

  // Now put the vector fields in the Eigen surface set.
  // I have decided to store the Vector Fields instead of
  // the Surfaces. They are equivalent.
  // Get Eigen Surface should return a Surface

  //Delete the old Eigen Surfaces

  EigenSurfaces->reset();
  Surface *tmpsurf;

  for (i=0;i<numsurfs;i++){
    tmpsurf = EigenSurfaces->addSurface(*this);

    for(j=0;j<numpoints;j++){
      dmu = sqrt((double) mNorm[j].x()*mNorm[j].x()+
		 mNorm[j].y()*mNorm[j].y()+
		 mNorm[j].z()*mNorm[j].z());

      dmu = sqrt(dmu);

      double x,y,z;
      x = *(U+i*3*numpoints+j*3)/dmu;
      y = *(U+i*3*numpoints+j*3+1)/dmu;
      z = *(U+i*3*numpoints+j*3+2)/dmu;

      tmpsurf->mVert[j].set(x,y,z);

    }

  }

  free(A);
  free(S);
  free(U);
  free(VT);
  free(WORK);

  // We are done Now!!

  return EigenValues;
}
*/




void Surface::getPlane(double &A, double &B, double &C, double &D,
		       const Point &p, const Point &n)
{
  A = n.x();
  B = n.y();
  C = n.z();
  D = -((n.x())*(p.x())+(n.y())*(p.y())+(n.z())*(p.z()));
}

void Surface::getBasisVectors(Point &b1, Point &b2, Point &b3, const Point &n)
{
  double b1x,b1y,b1z;
  double b2x,b2y,b2z;
  double b3x,b3y,b3z;
  if(fabs(n.x()) > TOLERANCE) {
    b1x = -1.0/n.x()*(n.y()+n.z());
    b1y = b1z = 1.0;
  }
  else if(fabs(n.y()) > TOLERANCE) {
    b1x = b1z = 1.0;
    b1y = -1.0/n.y()*(n.x()+n.z());
  }
  else if(fabs(n.z()) > TOLERANCE) {
    b1x = b1y = 1.0;
    b1z = -1.0/n.z()*(n.x()+n.y());
  }

  double mag = sqrt(b1x*b1x+b1y*b1y+b1z*b1z);
  if(mag==0.0) mag = 1.0;
  b1x /= mag;
  b1y /= mag;
  b1z /= mag;

  b3x = n.x();
  b3y = n.y();
  b3z = n.z();

  b2x = b3y*b1z-b1y*b3z;
  b2y = b3z*b1x-b3x*b1z;
  b2z = b3x*b1y-b3y*b1x;

  mag = sqrt(b2x*b2x+b2y*b2y+b2z*b2z);
  if(mag ==0.0) mag = 1.0;
  b2x /= mag;
  b2y /= mag;
  b2z /= mag;

  b1.set(b1x,b1y,b1z);
  b2.set(b2x,b2y,b2z);
  b3.set(b3x,b3y,b3z);
}


void Surface::setCurvature(Array1D<double> &dat)
{
  if(dat.getNelm() != mNumVert) {
    cerr << "ERROR: Surface::setCurvature invalid number of values" << endl;
    return;
  }
  mCurvature.setDim(mNumVert);
  mCurvature = dat;
}

ItXECode Surface::loadCurvature(char *fname)
{
  int i=0;
  ifstream in;
  in.open(fname,ios::in);
  if(!in) {
    cerr<<"ERROR: could not open file"<<endl;
    return ItXError;
  }

  mCurveType=Surface::UserDefined;

  mCurvature.setDim(mNumVert);

  i = 0;
  double tmp;
  while ( in >> tmp ) {
    i++;
    if ( i == mNumVert+1)
      break;
    mCurvature[i-1] = tmp;
  }
  in.close();

  if (i != (mNumVert) ) {
    cout << "Incorect Read" << endl;
    mCurveDirty=true;
    mCurvature.setDim(0);
    return ItXError;
  }

  mCurveDirty=false;
  return ItXSuccess;
}

ItXECode Surface::saveCurvature(char *fname)
{
  int j=0;
  if(!hasCurvature())
  {
    cerr<<"Error: no curvature to save!"<<endl;
    return ItXError;
  }

  ofstream out;
  out.open(fname,ios::out);

  if(!out) {
    cerr<<"ERROR: could not open file!"<<endl;
    return ItXError;
  }

  for(j=0;j<mNumVert;j++)  {
    out << mCurvature[j] << endl;
  }
  out.close();
  return ItXSuccess;
}


ItXECode Surface::genCurvature(CurvatureType ctype, int ndepth,
		float dx, float dy, float dz)
{
  if((mNumVert <= 0)||(mNumPoly <= 0))  {
    cerr << "ERROR: Surface::genCurvature empty surface" << endl;
    return ItXError;
  }

  ShowWorking("Generating curvature\n");
  if ( ctype != MaxCurvature && ctype != MinCurvature && ctype != MeanCurvature && ctype !=GaussCurvature)
  {
    cerr << "ERROR: invalid curvature type!" << endl;
    return ItXError;
  }

  ShowWorking("Generating curvature\n");

//  if ((hasCurvature()) && mCurveType == ctype) return(ItXSuccess);

  if((ndepth<1)||(ndepth>1000))
    ndepth = 2;

  mCurveType=ctype;
  genNeighborhoods(ndepth);
  if(!hasNbhd(ndepth)) {
    cerr << "ERROR: Surface::genCurvature cannot get neighborhoods" << endl;
    return ItXError;
  }

  if(!hasUNorms())  {
    genUnitNormals();
    if(!hasUNorms())
      return ItXError;
  }

  mCurvature.setDim(mNumVert);
  mCurvature = 0;

  int maxn = 0;

  // sets maxneighborhoods
  for(int m=0;m<mNumVert;m++)
    if(maxn < mNbhd[m].numNeighbors())
      maxn = mNbhd[m].numNeighbors();

//cout << "maximum # neighbors: " << maxn;
  maxn++;

  Array1D<float> h(maxn);
  Array1D<float> u(maxn);
  Array1D<float> v(maxn);
  Array1D<float> uu(maxn);
  Array1D<float> vv(maxn);
  Array1D<float> two_uv(maxn);
  float tmp[3],U[3][3],BU[3],work[3],workc[5];
  float rval;
  int   IPIV[4];
  float C[2][2],W[2];
  char  JOBZ, UPLO;
  int   INFO,  LDA, LDB, LWORK, N, NRHS;
  int   LDC, LWORKC, NC;
  int   i,j;

  for(int pindx=0;pindx<mNumVert;pindx++) {

    double a,b,c,d;
    Point  b1,b2,b3;
    int    num_nei = mNbhd[pindx].numNeighbors();

    double ptx,pty,ptz;
    ptx = mVert[pindx].x() * dx;
    pty = mVert[pindx].y() * dy;
    ptz = mVert[pindx].z() * dz;

    Point pp;
    pp.set(ptx,pty,ptz);

    Point nn = mUNorm[pindx];
    //nn.set(nn.x()*dx,nn.y()*dy,nn.z()*dz);

    Surface::getPlane(a,b,c,d,pp,nn);
    Surface::getBasisVectors(b1,b2,b3,nn);

    double nx = nn.x(); // mUNorm[pindx].x();
    double ny = nn.y(); // mUNorm[pindx].y();
    double nz = nn.z(); // mUNorm[pindx].z();

    for(i=0;i<num_nei;i++) {
      int nei = mNbhd[pindx].getNeighbor(i);
      if(nei < 0) continue;

      // Point pj = mVert[nei];
      double px = mVert[nei].x() * dx;
      double py = mVert[nei].y() * dy;
      double pz = mVert[nei].z() * dz;

      h[i] = a*px + b*py + c*pz + d;
      tmp[0] = px - h[i] * nx - ptx;
      tmp[1] = py - h[i] * ny - pty;
      tmp[2] = pz - h[i] * nz - ptz;

      u[i] = tmp[0]*b1.x()+tmp[1]*b1.y()+tmp[2]*b1.z();
      v[i] = tmp[0]*b2.x()+tmp[1]*b2.y()+tmp[2]*b2.z();

      two_uv[i] = 2.0*u[i]*v[i];
      uu[i] = u[i]*u[i];
      vv[i] = v[i]*v[i];
    }

    for(i=0;i<3;i++) {
      BU[i] = 0;
      for(j=0;j<3;j++)
	U[i][j] = 0;
    }

    for(i=0; i < num_nei; i++) {
      U[0][0] += (uu[i]*uu[i]);
      U[0][1] += (uu[i]*two_uv[i]);
      U[0][2] += (uu[i]*vv[i]);
      U[1][1] += (two_uv[i]*two_uv[i]);
      U[1][2] += (two_uv[i]*vv[i]);
      U[2][2] += (vv[i]*vv[i]);
      BU[0]   += (uu[i]*2.0*h[i]);
      BU[1]   += (two_uv[i]*2.0*h[i]);
      BU[2]   += (vv[i]*2.0*h[i]);
    }

    U[1][0] = U[0][1];
    U[2][0] = U[0][2];
    U[2][1] = U[1][2];

    UPLO  = 'L';
    INFO  = 1;
    LDA   = 3;
    LDB   = 3;
    LWORK = 3;
    N     = 3;
    NRHS  = 1;

    ssytrf_(&UPLO, &N, (float*)U, &LDA, IPIV, work, &LWORK, &INFO);

    if(INFO == 0) {
      ssytrs_(&UPLO, &N, &NRHS, (float*)U, &LDA, IPIV, BU,
	      &LDB , &INFO);
      C[0][0] = BU[0];
      C[0][1] = BU[1];
      C[1][0] = BU[1];
      C[1][1] = BU[2];

      /* get the eigenvalue of the matrix c */
      JOBZ   = 'N';
      LDC    = 2;
      NC     = 2;
      LWORKC = 3*NC-1;
      ssyev_(&JOBZ, &UPLO, &NC, (float*)C, &LDC, W, workc,
	     &LWORKC,&INFO);


      switch(ctype) {
      case Surface::MinCurvature:
	if(fabs(W[0]) < fabs(W[1]))
	  rval = W[0];
	else    rval = W[1];
	break;
      case Surface::MaxCurvature:
	if(fabs(W[0]) > fabs(W[1]))
	  rval = W[0];
	else    rval = W[1];
	break;
      case Surface::MeanCurvature:
	rval = (W[0]+W[1])/2.0;
	break;
      case Surface::GaussCurvature:
      default:
	rval = W[0] * W[1];
	break;
      }
    }
    else    {
      rval = 0.0;
      cerr << "ERROR: systrs failed in genCurvature()" << endl;
    }

    // crop since some degenerate polys give really large curvature
    if(fabs(rval) > MAX_CURV) {
	if(rval < 0) rval = -MAX_CURV;
	else rval = MAX_CURV;
	}

    // neg because we changed sign of polys
    mCurvature[pindx] = -rval;
  }

 skipthis:
  cerr << endl << "DONE" << endl;
  mCurveDirty = false;
  return(ItXSuccess);
}

void Surface::extract() {

    // method creates separate surfaces from each of the components
    // of this surface
    //
    // array to keep track of vertices that have appeared in a surface component
    Array1D<int> vertVisited(this->mNumVert);
    vertVisited = -1;  // "false"

    // what vertices showed up in the latest extracted surface
    Array1D<int>* thisVertVisited = new Array1D<int>(this->mNumVert);

    int vertInLargestSurf = 0;
    int vertCountOfLargestSurf = 0;
    bool done = false;
    while (!done)
    {
      Surface* thisExtractedSurface = new Surface(*this);
      int vertIx = vertInLargestSurf;
      while (vertIx < this->mNumVert && vertVisited[vertIx] != -1)
        vertIx++;

      if (vertVisited[vertIx] == -1)
      {
        // haven't checked this vertex yet...
        thisExtractedSurface->extract(vertIx, thisVertVisited);
        // bigger than current biggest?
        if (thisExtractedSurface->getNumVert() > vertCountOfLargestSurf)
        {
          vertInLargestSurf = vertIx;
          vertCountOfLargestSurf = thisExtractedSurface->getNumVert();
        }
        // check off vertices that have been visited
        for (int visitIx = vertIx; visitIx < this->mNumVert; visitIx++)
          if ((*thisVertVisited)[visitIx] == 1) vertVisited[visitIx] = 1;
      }
      else
        done = true;

      delete thisExtractedSurface;
    }
    this->extract(vertInLargestSurf);
    delete thisVertVisited;
}

void Surface::extract(int ind, Array1D<int>* includedInExtract) {

	 int i;
	genNeighborhoods();

	//array to flag vertices which have been included

	Array1D<int> included;
	included.setDim(mNumVert);
	included=-1;
	included[ind]=1;

	//rest is adapted from regiongrow algorithm
   int count=0;
	int *stack=new int[mNumVert];
    int stacksize = 1;
	//psuh the initial point onto stack
	stack[0]=ind;
    int temp_ind;
	int temp_ngh;
	Neighborhood  ngbhs;
	while(stacksize !=0) {
		temp_ind= stack[stacksize-1];
		ngbhs.clear();
		ngbhs= mNbhd[temp_ind];
        stacksize--;

		for (i=0; i<ngbhs.numNeighbors(); ++i) {
			temp_ngh= ngbhs.getNeighbor(i);
			if (included[temp_ngh] == -1) {
				stack[stacksize++] =temp_ngh;
				included[temp_ngh]=1;

			}
		}
	}
	cout<<"done regiongrow"<<endl;
	delete [] stack;

    if (includedInExtract)
      (*includedInExtract) = included;

    int j;
	count=0;
	Surface newSurf;
    newSurf.mVert.setDim(mNumVert);
	if (included[0]==1) { included[0]=0; count++; newSurf.mVert[0]=mVert[0]; }
	else cout<<"0 not included "<<endl;
	for (i=1; i<mNumVert; ++i) {
		if (included[i]==1) {
			for (j=i-1; j>=0; --j) {
				if (included[j]!=-1) {
					included[i]=included[j]+1;
					break;
				}
			}
				if (j==-1) { included[i]=0;  cout<<"1st vertex   "<<endl; }

			newSurf.mVert[count++]= mVert[i];

		}
	}
	mVert.setDim(count);
	mNumVert=count;
	for (i=0; i<mNumVert; ++i)
		mVert[i] = newSurf.mVert[i];
	cout<<"done deletin unused vertices "<<count<<endl;

//tentetavily allocate space for facets

	newSurf.mFacet.setDim(mNumPoly,3);
	int tempx,tempy,tempz;
	count=0;
	for (i=0; i<mNumPoly; ++i) {
		tempx= mFacet[i][0];
		tempy= mFacet[i][1];
		tempz= mFacet[i][2];

		if (included[tempx]!=-1) {
			newSurf.mFacet[count][0]=included[tempx];
			newSurf.mFacet[count][1]=included[tempy];
			newSurf.mFacet[count][2]=included[tempz];
			count++;
		}
	}

	mNumPoly= count;
	mFacet.setDim(mNumPoly,3);

	for (i=0;  i<mNumPoly; ++i) {
		mFacet[i][0] = newSurf.mFacet[i][0];
		mFacet[i][1] = newSurf.mFacet[i][1];
		mFacet[i][2] = newSurf.mFacet[i][2];
	}

	geometryChanged();
    genNeighborhoods();

}


// Lei 04/17/2003
// keep the triangle indices, just show the extracted ones, however
void Surface::extractKeep(int ind) {

  int i;
  genNeighborhoods();

  //array to flag vertices which have been included

  Array1D<int> included;
  included.setDim(mNumVert);
  included=-1;
  included[ind]=1;

  //rest is adapted from regiongrow algorithm
  int count=0;
  int *stack=new int[mNumVert];
  int stacksize = 1;
  //psuh the initial point onto stack
  stack[0]=ind;
  int temp_ind;
  int temp_ngh;
  Neighborhood  ngbhs;
  while(stacksize !=0) {
    temp_ind= stack[stacksize-1];
    ngbhs.clear();
    ngbhs= mNbhd[temp_ind];
    stacksize--;

    for (i=0; i<ngbhs.numNeighbors(); ++i) {
      temp_ngh= ngbhs.getNeighbor(i);
      if (included[temp_ngh] == -1) {
	stack[stacksize++] =temp_ngh;
	included[temp_ngh]=1;

      }
    }
  }
  cout<<"done regiongrow"<<endl;
  delete [] stack;


  int j;
  count=0;
  Surface newSurf;
  newSurf.mVert.setDim(mNumVert);
  if (included[0]==1) { included[0]=0; count++; newSurf.mVert[0]=mVert[0]; }
  else cout<<"0 not included "<<endl;
  newSurf.mVert[0]=mVert[0];
  for (i=1; i<mNumVert; ++i) {
    if (included[i]==1) {
      for (j=i-1; j>=0; --j) {
	if (included[j]!=-1) {
	  included[i]=included[j]+1;
	  break;
	}
      }
      /*----------------------------<Lei>----------------------------*
       | if (j==-1) { included[i]=0;  cout<<"1st vertex   "<<endl; } |
       | newSurf.mVert[count++]= mVert[i];                           |
       *-------------------------------------------------------------*/
      if (j==-1) { included[i]=0;  cout<<"1st vertex   "<<endl; }
    }
    newSurf.mVert[i]= mVert[i];
  }
  cout<<"done keeping unused vertices "<< i <<endl;
  /*--------<Lei>---------*
   | mVert.setDim(count); |
   | mNumVert=count;      |
   *----------------------*/
  mVert.setDim(mNumVert);
  for (i=0; i<mNumVert; ++i)
    mVert[i] = newSurf.mVert[i];
  cout<<"done keeping unused vertices "<< i <<endl;

  // tentetavily allocate space for facets

  newSurf.mFacet.setDim(mNumPoly,3);
  int tempx,tempy,tempz;
  count=0;
  for (i=0; i<mNumPoly; ++i) {
    tempx= mFacet[i][0];
    tempy= mFacet[i][1];
    tempz= mFacet[i][2];

    if (included[tempx]!=-1) {
      /*-------------------<Lei>-------------------*
       | newSurf.mFacet[count][0]=included[tempx]; |
       | newSurf.mFacet[count][1]=included[tempy]; |
       | newSurf.mFacet[count][2]=included[tempz]; |
       *-------------------------------------------*/
      newSurf.mFacet[count][0]=tempx;
      newSurf.mFacet[count][1]=tempy;
      newSurf.mFacet[count][2]=tempz;
      count++;
    }
  }
  cout<<"done deleting unused polygons "<< count <<endl;

  mNumPoly= count;
  mFacet.setDim(mNumPoly,3);

  for (i=0;  i<mNumPoly; ++i) {
    mFacet[i][0] = newSurf.mFacet[i][0];
    mFacet[i][1] = newSurf.mFacet[i][1];
    mFacet[i][2] = newSurf.mFacet[i][2];
  }

  geometryChanged();
  genNeighborhoods();

}


void Surface::decimate(double tolerance) {

	cout<<"tolerance "<<tolerance<<" numpoly "<<mNumPoly<<endl;
	int i,j,k,m;
	Array1D<int> included;
	included.setDim(mNumVert);
	Array1D<int> replaced;
	replaced.setDim(mNumVert);

	if (!hasNbhd(1))
		genNeighborhoods(1);

	Array1D<unsigned char> edge;
	isEdge(edge);

	for (i=0; i<mNumVert; ++i) {
		included[i]=i;
		replaced[i]=i;

	}

	int temp1,  temp2;
	double length;
	int dummy;
	int count=0;


	for (i=0; i<mNumPoly; ++i) {
		for (j=0;  j<3; ++j) {

		retry:
			//if degenerate triangle (after an edge collapsed to a single point)
		  	if ((mFacet[i][0]==mFacet[i][1]) || (mFacet[i][1]==mFacet[i][2])  || (mFacet[i][0]==mFacet[i][2])) {
				mFacet[i][0]=mFacet[i][1]=mFacet[i][2]=-1;
				j=3;
				break;
  			}

			if (included[mFacet[i][j]]>=0  && included[mFacet[i][(j+1)%3]]>=0) {
				temp1 = mFacet[i][j];
				temp2=  mFacet[i][(j+1)%3];
				Point temp = mVert[temp1]-mVert[temp2];
				length = temp.norm();

				if (length <tolerance  && (!edge[temp1] && !edge[temp2]) ) {

					count=0;
					Array1D<int> comVert;
					comVert.setDim(20);
					for (k=0; k<mNbhd[temp1].numNeighbors(); ++k)
						for (m=0; m<mNbhd[temp2].numNeighbors(); ++m)
							if (mNbhd[temp1].getNeighbor(k) == mNbhd[temp2].getNeighbor(m))
								comVert[count++]=mNbhd[temp1].getNeighbor(k);

					for (m=0; m<count; ++m)
						for (k=0; k<count; ++k)
							if(mNbhd[comVert[k]].isNeighbor(comVert[m])) {
								if (i!=mNumPoly-1) {
									i++;
									j=0;
									goto retry;
								}
								else goto done;
							}

					if (m==count) {
					for (k=0; k<mNbhd[temp2].numNeighbors(); ++k) {
						dummy=mNbhd[temp2].getNeighbor(k);
						if (dummy!=temp1) {
							mNbhd[temp1].addNeighbor(dummy);
							mNbhd[dummy].removeNeighbor(temp2);
							mNbhd[dummy].addNeighbor(temp1);
						}
					}  //end k

					mNbhd[temp2].clear();
					mNbhd[temp1].removeNeighbor(temp2);
					mVert[temp1] = (mVert[temp1]+mVert[temp2])*0.5;
					replaced[temp2]=temp1;
					included[temp2]=-1;
					mFacet[i][0]=mFacet[i][1]=mFacet[i][2]=-1;
					j=3;
					break;
					}
				} //endif length

			} //endif legal vertices
			else {
				if (included[mFacet[i][j]]<0) {
					dummy = mFacet[i][j];
					while (included[replaced[dummy]]<0)
						dummy = replaced[dummy];
					mFacet[i][j] = replaced[dummy];
					goto retry;
				}
				if (included[mFacet[i][(j+1)%3]]<0) {
					dummy=mFacet[i][(j+1)%3];
					while (included[replaced[dummy]]<0)
						dummy=replaced[dummy];
					mFacet[i][(j+1)%3]=replaced[dummy];
					goto retry;
				}

			}  //endelse
		} //end j

	} //end i

	done:
	count = 0;
		//clean up vertex list
	Array1D<Point> newVerts;
	newVerts.setDim(mNumVert);

	if (included[0]>=0) { included[0]=0; count++; newVerts[0]=mVert[0]; }
	for (i=1; i<mNumVert; ++i) {
		if (included[i]>=0) {
			for (j=i-1; j>=0; --j) {
				if (included[j]>=0) {
					included[i]=included[j]+1;
					break;
				}
			}
			if (j==-1)  included[i]=0;
			newVerts[count++]= mVert[i];

		}
	}

	mVert.setDim(count);
	mNumVert=count;
	for (i=0; i<mNumVert; ++i)
		mVert[i] = newVerts[i];


	Array2D<int> newFacets;
	newFacets.setDim(mNumPoly,3);

	count=0;


//clean up poly list
	int tempx,tempy,tempz;
	for (i=0; i<mNumPoly; ++i) {
		tempx= mFacet[i][0];
		tempy= mFacet[i][1];
		tempz= mFacet[i][2];

		if (tempx!=-1) {

			while (included[tempx]<0)
				tempx=replaced[tempx];
			while (included[tempy]<0)
				tempy=replaced[tempy];
			while (included[tempz]<0)
				tempz=replaced[tempz];

			if (included[tempx]!=included[tempy] && included[tempy]!=included[tempz]
				&& included[tempx]!=included[tempz]) {
				newFacets[count][0]=included[tempx];
				newFacets[count][1]=included[tempy];
				newFacets[count][2]=included[tempz];
				count++;
			}
		}
	}

	mNumPoly=count;
	mFacet.setDim(count,3);
	for (i=0; i<count; ++i) {
		mFacet[i][0]=newFacets[i][0];
		mFacet[i][1]=newFacets[i][1];
		mFacet[i][2]=newFacets[i][2];
	}
	geometryChanged();
	extract(mFacet[0][0]);
	geometryChanged();

}

void Surface::isEdge(Array1D<unsigned char> &edge) {

	Array1D<int> triangleCount;
	triangleCount.setDim(mNumVert);
	triangleCount=0;
	int i;

	if (!hasNbhd(1))
		genNeighborhoods(1);

	for (i=0; i<mNumPoly; ++i) {
		triangleCount[mFacet[i][0]]++;
		triangleCount[mFacet[i][1]]++;
		triangleCount[mFacet[i][2]]++;
	}

        edge.setDim(mNumVert);
	for (i=0; i<mNumVert; ++i) {
		if (triangleCount[i] != mNbhd[i].numNeighbors())
		     edge[i]=1;
		else edge[i]=0;
	}

}

int Surface::estimatePolyCount(double tolerance) {

	int i,j;

	Point P, R,temp;
	int x,y;
	Array1D<unsigned char> edge;

	isEdge(edge);

	int count=0;

	for (i=0; i<mNumPoly; ++i) {
		for (j=0; j<3; ++j) {
			x=mFacet[i][j];
			y=mFacet[i][(j+1)%3];
			P=mVert[x];
			R=mVert[y];
			temp=P-R;
			if (temp.norm() <tolerance && !edge[x] && !edge[y])
				count++;

		}
	}
	return (mNumPoly-count);
}



void Surface::generateFlowers(char *file) {
	ofstream fp;
	fp.open(file,ios::out);
	int i,a,b,c;
	struct cell {
		int index;
		struct cell *next;
		int dirty;
	};

	cell **nei=new cell*[mNumVert];
	Array1D<int> neiCount;
	neiCount.setDim(mNumVert);
	neiCount=0;

	for (i=0; i<mNumVert; ++i ) {
		nei[i]=NULL;
		}
	cell *t1=NULL,*t2=NULL,*t3=NULL;
	cell *a_t,*b_t=NULL, *c_t=NULL, *pre=NULL, *old_t=NULL;
	int found=0;


	for (i=0; i<mNumPoly; ++i) {
		a= mFacet[i][0];
		b= mFacet[i][1];
		c= mFacet[i][2];

		t1=NULL;  t2=NULL; t3=NULL;
		a_t=NULL; b_t=NULL;  c_t=NULL;
		pre=NULL;  old_t=NULL;
		found=0;

		t1=nei[a];
		while(t1) {
			if (t1->index == b) {
				b_t=t1;
				found=1;
				if (c_t) break;
			}
			if (t1->index == c) {
				c_t=t1;
				pre=old_t;
				found=1;
				if (b_t) break;
			}
			if (!t1->next)  break;
			else { old_t=t1;  t1=t1->next; }
		}
			//case 0: both not there: add to the end
		if (!found) {
			neiCount[a]+=2;
			t2=new cell;
			t3=new cell;

			t2->index=b;
			t3->index=c;

			t2->next=t3;
			t2->dirty=0;

			t3->next=NULL;
			t3->dirty=1;

			if (t1) {
				t1->next=t2;
				t1->dirty=1;
			}
			else nei[a]=t2;
		}
			//case 1: b is there, c is not

		else if (b_t && !c_t) {
			neiCount[a]++;

			t2=new cell;
			t2->index=c;
			t2->next=b_t->next;
			t2->dirty=1;

			b_t->next=t2;
			b_t->dirty=0;
		}
			//case 12: both there
		else if (b_t && c_t) {
			t3=c_t;
			t2=b_t->next;

				//make sure they are not already next to each other
			if (t2!=t3) {
				while(t3 && !t3->dirty && t3!=b_t) {
					if (t3->next) t3 = t3->next;
					else break;
				}

				if (t3!=b_t ) {
					if (pre) { pre->next=t3->next; pre->dirty=1; }
					else nei[a]=t3->next;

					b_t->next=c_t;
					b_t->dirty=0;
					t3->next=t2;
					t3->dirty=1;
				}
			}
			else { b_t->dirty=0; b_t->next=c_t; }
		}
			//case 2: c there, b is not
		else if (!b_t && c_t) {
			neiCount[a]++;
			t2=new cell;
			t2->index=b;
			t2->dirty=0;
			t2->next=c_t;
				//case 21: b goes to the head
			if (!pre)
				nei[a]=t2;

				//case 22: b goes anwhere else
			else {
				pre->next=t2;
				pre->dirty=1;
			}
		}

		found=0;
		t1=t2=t3=old_t=a_t=b_t=c_t=pre=NULL;

		t1=nei[b];
		while(t1) {

			if (t1->index == c) {
				c_t=t1;
				found=1;
				if (a_t) break;
			}
			 if (t1->index == a) {
				 a_t=t1;
				 pre=old_t;
				 found=1;
				 if (c_t) break;
			 }
			 if (!t1->next)   break;
			 else { old_t=t1;  t1=t1->next; }
		 }
			//case 0: both not there: add to the end
		if (!found) {
			neiCount[b]+=2;
			t2=new cell;
			t3=new cell;

			t2->index=c;
			t3->index=a;

			t2->next=t3;
			t2->dirty=0;

			t3->next=NULL;
			t3->dirty=1;

			if (t1) {
				t1->next=t2;
				t1->dirty=1;
			}
			else   nei[b]=t2;
		}
			//case 1: c is there, a is not

		else if (c_t && !a_t) {
			neiCount[b]++;
			t2=new cell;
			t2->index=a;
			t2->next=c_t->next;
			t2->dirty=1;

			c_t->next=t2;
			c_t->dirty=0;
		}
			//case 12: both there
		else if (c_t && a_t) {
			if (c_t->dirty==0) cout<<"problem c_a clean 2"<<endl;
			t2=c_t->next;
			t3=a_t;


			 if (t2!=t3) {
			 while (t3 && !t3->dirty  &&t3!=c_t) {
				 if (t3->next) t3=t3->next;
				 else break;
			 }

			 if (t3!=c_t) {

				 if (pre) { pre->next=t3->next; pre->dirty=1; }
				 else nei[b]=t3->next;

				 c_t->next=a_t;
				 c_t->dirty=0;
				 t3->next=t2;
				 t3->dirty=1;
			 }
			 }
			 else {c_t->dirty=0; c_t->next=a_t; }

		 }
			 //case 2: a there, c is not
		 else if (!c_t && a_t) {
			 neiCount[b]++;
			 t2=new cell;
			 t2->index=c;
			 t2->dirty=0;
			 t2->next=a_t;
				 //case 21: c goes to the head
			 if (!pre)
				 nei[b]=t2;

				 //case 22: c goes anwhere else
			 else {
				 pre->next=t2;
				 pre->dirty=1;
					 // t2->next=a_t;
			 }
		 }


		 found = 0;
		 t1=t2=t3=old_t=a_t=b_t=c_t=pre=NULL;
		 t1=nei[c];
		 while(t1) {

			 if (t1->index == a) {
				 a_t=t1;
				 found=1;
				 if (b_t) break;

			 }
			 if (t1->index == b) {
				 b_t=t1;
				 pre=old_t;
				 found=1;
				 if (a_t) break;
			 }
			 if (!t1->next ) {  break; }
			 else { old_t=t1;  t1=t1->next;  }
		 }
			 //case 0: both not there: add to the end
		 if (!found) {
			 neiCount[c]+=2;
			 t2=new cell;
			 t3=new cell;

			 t2->index=a;
			 t3->index=b;

			 t2->next=t3;
			 t2->dirty=0;

			 t3->next=NULL;
			 t3->dirty=1;

			if (t1) {
				t1->next=t2;
				t1->dirty=0;
			}
			else nei[c]=t2;
		 }
			 //case 1: a is there, b is not

		 else if (a_t && !b_t) {
			 neiCount[c]++;

			 t2=new cell;
			 t2->index=b;
			 t2->next=a_t->next;
			 t2->dirty=1;

			 a_t->next=t2;
			 a_t->dirty=0;
		 }
		 //case 12: both there
		 else if (a_t && b_t) {

			 t2=a_t->next;
			 t3=b_t;

			 if (t2!=t3) {
				 while(t3 && !t3->dirty && t3!=a_t) {

					 if (t3->next) t3 = t3->next;
					 else break;
				 }

			  if (t3!=a_t) {
				  if (pre) { pre->next=t3->next; pre->dirty=1; }
				  else nei[c]=t3->next;

				  a_t->next=b_t;
				  a_t->dirty=0;
				  t3->next=t2;
				  t3->dirty=1;
			  }
			 }
			 else { a_t->dirty=0; a_t->next=b_t; }

		 }
		 else if (!a_t && b_t) {
			 neiCount[c]++;

			 t2=new cell;
			 t2->index=a;
			 t2->dirty=0;
			 t2->next=b_t;
				 //case 21: a goes to the head
			 if (!pre)
				 nei[c]=t2;

				 //case 22: a goes anwhere else
			 else {
				 pre->next=t2;
				 pre->dirty=1;
					 // t2->next=b_t;
			 }
		 }

	}
	Array1D<int> triangleCount;
	triangleCount.setDim(mNumVert);
	triangleCount=0;


	for (i=0; i<mNumPoly; ++i) {
		triangleCount[mFacet[i][0]]++;
		triangleCount[mFacet[i][1]]++;
		triangleCount[mFacet[i][2]]++;
	}

	Array1D<unsigned char> edge;
	isEdge(edge);

	for (i=0; i<mNumVert; ++i) {
		fp<<i+1<<"\t"<<triangleCount[i]<<"\t";

		cell *t1=nei[i];
		while (t1) {
			fp<<t1->index+1<<" ";
			t1=t1->next;
		}
		if (!edge[i]) fp<<nei[i]->index+1<<endl;
		else fp<<endl;
	}


}

void Surface::cutSurface(Array1D<int> &selected) {
	int i,j,count;

	Array1D<int> included(mNumVert);
	included=-1;

	for (i=0; i<mNumPoly; ++i)
		if (selected[i]) {
			included[mFacet[i][0]]=1;
			included[mFacet[i][1]]=1;
			included[mFacet[i][2]]=1;
		}

	//rest is adapted from extract algorithm


	count=0;
	Surface newSurf;
    newSurf.mVert.setDim(mNumVert);
	if (included[0]==1) { included[0]=0; count++; newSurf.mVert[0]=mVert[0]; }

	for (i=1; i<mNumVert; ++i) {
		if (included[i]==1) {
			for (j=i-1; j>=0; --j) {
				if (included[j]!=-1) {
					included[i]=included[j]+1;
					break;
				}
			}
			if (j==-1)  included[i]=0;

			newSurf.mVert[count++]= mVert[i];

		}
	}
	mVert.setDim(count);
	mNumVert=count;
	for (i=0; i<mNumVert; ++i)
		mVert[i] = newSurf.mVert[i];
	cout<<"done deleting unused vertices "<<count<<endl;


	//tentetavily allocate space for facets

	newSurf.mFacet.setDim(mNumPoly,3);
	int tempx,tempy,tempz;
	count=0;
	for (i=0; i<mNumPoly; ++i) {
		tempx= mFacet[i][0];
		tempy= mFacet[i][1];
		tempz= mFacet[i][2];

		if (included[tempx]!=-1 && included[tempy]!=-1 && included[tempz]!=-1) {
			newSurf.mFacet[count][0]=included[tempx];
			newSurf.mFacet[count][1]=included[tempy];
			newSurf.mFacet[count][2]=included[tempz];
			count++;
		}
	}

	mNumPoly= count;
	mFacet.setDim(mNumPoly,3);

	for (i=0;  i<mNumPoly; ++i) {
		mFacet[i][0] = newSurf.mFacet[i][0];
		mFacet[i][1] = newSurf.mFacet[i][1];
		mFacet[i][2] = newSurf.mFacet[i][2];
	}

	geometryChanged();

}

//find (the index of) closest point on the surface given a point P
int Surface::findClosestPoint(Point P) {

	int i;
	int min_ind;
	double dist, min_dist=HUGE;

	for (i=0; i<mNumVert; ++i) {
		dist=(P-mVert[i]).norm();
		if (dist < min_dist) {
			min_dist = dist;
			min_ind =i;
		}

	}

	return min_ind;
}

void Surface::readSettings(char *fname) {

	cout<<"reading settings "<<endl;
	ifstream fin(fname);
	double x, y,z;
	int dummy;
	int i,count=0;

	mProbes.setDim(100);
	mProbeNames.setDim(100);

	for (i=0; i<100; ++i)
		mProbeNames[i]=new char[255];

	while (fin>>x>>y>>z>>dummy && count<100) {
		mProbes[count].set(x,y,z);
		fin.getline(mProbeNames[count],255);
		count ++;
		cout<<count<<endl;
	}

	fin.close();

	int ind;
	for (i=0; i<count; ++i) {
		ind=findClosestPoint(mProbes[i]);
		mProbes[i]=mVert[ind];
	}
	cout<<"done reading settings "<<endl;
}

void Surface::getSettings(Array1D<Point> &points, Array1D<char *> &names) {
	points=mProbes;
	names=mProbeNames;

}
