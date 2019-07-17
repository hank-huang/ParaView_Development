#ifndef _Surface_h_
#define _Surface_h_
///////////////////////////////////////////////////////////////////////////
//
// File: Surface.h
//
// Author: Rob Teichman
//
// Purpose: Surface Class Definition
//  This class is a base class for representation of surfaces.  It includes
// basic operations performed on surfaces.
//
// All derrived classes must implement the read and write routines.
// 
///////////////////////////////////////////////////////////////////////////
#include <OS.h>
#include <iostream.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/param.h>
#include <ADL/Array1D.h>
#include <ADL/Array2D.h>
#include <ADL/Vector.h>
#include <ADL/Matrix.h>
#include <ADL/Lapack.h>
#include <TDL/Neighborhood.h>
#include <TDL/Point.h>
#include <TDL/Line.h>
#include <TDL/IDLdefines.h>
#include <TDL/ItXWorkClass.h>


#ifndef TOLERANCE
#define TOLERANCE 0.0001                   // for curvature calculations
#endif

//class SurfaceSet;

class Surface : public ItXWorkClass
{

friend class SurfaceUtils;
friend class SurfaceGrow;
friend class SurfaceReducer;
friend class Embed;
friend class SurfaceGen;
friend class SurfaceTracker;

public:
  // for genCurvature()
  enum CurvatureType {
        NotValid,MeanCurvature, MaxCurvature, MinCurvature, GaussCurvature,
        UserDefined};

  // default constructor
  Surface() :
  mNumVert(0), mNumPoly(0), mNormDirty(true), 
    mContDirty(true), mNbhdDepth(0), mVerbose(false) {
      mCurveType=NotValid;   
      sprintf(filename,"<surface>");
    }

  // copy constructor
  Surface(Surface &S);

  // destructor
  virtual ~Surface() {}

 ////  load and save curvature to disk file
  ItXECode loadCurvature(char *);
  ItXECode saveCurvature(char *);

  //
  // virtual functions for file i/o
  //
  virtual bool read (istream &inFile = cin)
  { cerr << "Cannot call read on type Surface" << endl; return false; };

  virtual bool read (const char *fileName)
  { cerr << "Cannot call read on type Surface" << endl; return false; };

  virtual bool write (ostream &outFile = cout) const
  { cerr << "Cannot call write on type Surface" << endl; return false; };

  virtual bool write (const char *fileName) const
  { cerr << "Cannot call write on type Surface" << endl; return false; };

  bool Load(const char *fileName)
	{ SetFilename(fileName); return(read(fileName)); }

  bool Save(const char *fileName)
	{ SetFilename(fileName); return(write(fileName)); }

  void  SetFilename(const char *nm);
  inline const char *Filename() { return(filename); }

  //
  // direct surface modification routines
  //
  int addFacet(int a, int b, int c);
  int addVert(Point const &p);
  int addVert(double x, double y, double z)
  { Point p(x,y,z); return addVert(p); }

  int changeVert(int n, double x, double y, double z);
  int changeVert(int n, Point x);

  //
  // if surface topology or vertex values
  // change, this will invalidate neighborhoods,
  // curvature, normals, etc.
  //

  void verticesChanged();
  void geometryChanged();
  void flipPolyOrientation(int i); // Flip the Orientation of the ith polygon
  void flipOrientation(); // FLip Orientation of the entire surface

  // cleanup routine
  bool clean();

  // assigment operator
  Surface & operator= (Surface const &S);

  // addition operator, adds two surfaces vertex-wise
  Surface & operator+= (Surface const &S);

  // division by scaler
  Surface & operator/= (double const f)
  { this->scale(1.0/f); return *this; }

  // combine two surfaces into one
  Surface & combine (Surface const &S);

  // generate normals
  void genNormals();
  void genUnitNormals();

  // compute neighborhoods of given size at each vertex
  void genNeighborhoods(int Size = 1);
  void freeNeighborhoods()		{ mNbhd.setDim(0); mNbhdDepth = 0; }

  // compute curvature
  ItXECode genCurvature(CurvatureType, int depth = 2, 
			float dx=1.f, float dy=1.f, float dz=1.f);
  void     freeCurvature()		{ mCurvature.setDim(0); }
  void     setCurvature(Array1D<double>&);

  // compute the Euler Characteristic number
  int euler();

  // compute volume enclosed by surface
  double volume();

  // check if surface is closed
  bool isClosed(Array1D<int> &vertlist);

  // compute the centroid of the surface
  int getCentroid(Point *cent);

  int getNeighborhoodDepth()	{ return(mNbhdDepth); }

  // compute the simple centroid of the surface
  // avg of vertices
  int getSimpleCentroid(Point *cent);

  void getBoundingBox(Point &upper_left, Point &lower_right);

  inline int getNumVert()	{ return(mNumVert); }
  inline int getNumPoly()	{ return(mNumPoly); }

  inline Point getVert(int i)  { return mVert[i]; }
  inline int *getFacet(int i)  { return (&mFacet[i][0]); }

  inline const Array1D<Point>  &vertices()  const { return(mVert);      }
  inline const Array1D<Point>  &normals()   const { return(mNorm);      }
  inline const Array1D<Point>  &unormals()  const { return(mUNorm);     }
  inline const Array2D<int>    &facets()    const { return(mFacet);     }
  inline const Array1D<double> &curvature() const { return(mCurvature); }
  inline const Array1D<Neighborhood> &neighborhood() const { return(mNbhd); }

  inline Array1D<Point>  &vertices()  { return(mVert);      }
  inline Array2D<int>    &facets()    { return(mFacet);     }

  inline void setSize(int nvert, int nfacet) {
	mVert.setDim(nvert);
	mFacet.setDim(nfacet,3);
        mNumVert = nvert;
        mNumPoly = nfacet;
	geometryChanged();
	}

  // Get bounding box
  void getDimensions(float &xmin, float &ymin, float &zmin,
                     float &xmax, float &ymax, float &zmax);

  //
  // Affine transformations
  //

  // generic affine transformation (3x3 matrix and 3x1 vector,
  // Ax+B applied to each vertex x in the surface)
  void affine(Matrix<double> const &A, Vector<double> const &B);
  void affine(Matrix<float> const &A, Vector<float> const &B);

  void scale(double scaleFactor)	// scale uniformly in all dimensions
    { this->scale(scaleFactor, scaleFactor, scaleFactor); }

  void scale(double scaleX, double scaleY, double scaleZ);
  void translate(double transX, double transY, double transZ);
  void rotate(double roll, double pitch, double yaw);

  // transformations around centriod
  void rotateCentroid(double roll, double pitch, double yaw);

  //
  // Functions to clean up surfaces with common problems
  //
  void removeDegenerateTriangles(double tolerance = 0);
  void fixNormals();

  // Functions to check if data is valid
  bool isInit() const		// returns true if facet and vert are loaded
  { return ( mNumVert != 0 && mNumPoly != 0); }
  bool hasNorms() const		// returns true if norm are valid
  { return (!mNorm.isEmpty() && !mNormDirty); }
  bool hasUNorms() const	// returns true if uNorm are valid
  { return (!mUNorm.isEmpty() && !mNormDirty); }
  bool hasNbhd(int need_depth = 1) const	// returns true if nbhd has been calculated
  { return ((!mNbhd.isEmpty())&&(mNbhdDepth == need_depth)); }	// at that particular depth
  bool hasCurvature() const	// returns true if curvature is valid
  { return (!mCurvature.isEmpty()); }

  //
  // Functions to get contours
  //
  // force generate all
  void generateContours();

  Line *getXYContour(double x, double y, double z);
  Line *getXZContour(double x, double y, double z);
  Line *getYZContour(double x, double y, double z);

  void setupContours();

  // temporary functions (for contour generation)
  double findXMin() const;
  double findXMax() const;
  double findYMin() const;
  double findYMax() const;
  double findZMin() const;
  double findZMax() const;

  // field transformation
  int applyFieldTrans(const Array3D<float> &hFields);

  // inverse field transformation routine
  bool invFieldTrans(const Array3D<float> &hFields);

  // Function to compute the innerproduct 
  double innerproduct(Surface &basesurf,Surface const &eigensurface);

//  Array1D<double> ComputeEigenSurf(SurfaceSet &inSurfs, SurfaceSet *EigenSurfaces);

  // set/get verbose flag
  inline void setVerbose(int verbose);
  inline bool getVerbose();
  // extract the largest component of the surface
  void extract();
  void extract(int ind, Array1D<int>* includedInExtract = NULL);
  void extractKeep(int ind);	// Lei 02/04/2003
  int estimatePolyCount(double tolerance);
  void decimate(double tolerance);
  void isEdge(Array1D<unsigned char> &edge);
  void generateFlowers(char *file);
  void cutSurface(Array1D<int> &selected);	  
  int findClosestPoint(Point P);
  void readSettings(char *fname);
 void getSettings(Array1D<Point> &points, Array1D<char *> &names);
protected:

  // verbose flag for surface class
  bool mVerbose;

  //const double TOLERANCE=0.0001; // for curvature calculations - now a pound-define
  //const int MAX_PARTS=1;	// For now only handles simple surfaces,
				// not scenes (or groups) of surfaces

  char filename[MAXPATHLEN];	// Filename of surface (if any)

  // Fundamental surface data
  int mNumVert;			// total number of vertices
  int mNumPoly;			// total number of polygons
 CurvatureType mCurveType;
  Array2D<int> mFacet;		// array that describes the connectivity of
				// vertices to form triangles
  Array1D<Point> mVert;		// list of vertices for the surface

  // additional useful surface information, initialized as needed
  Array1D<Point> mNorm;		// list of normals at each vertex
  Array1D<Point> mUNorm;	// list of unit normals at each vertex
  Array1D<Neighborhood> mNbhd;	// contains neighborhood of each vertex
  int     mNbhdDepth;

  // contours
  Array1D<Line> m_XYContours;
  Array1D<Line> m_XZContours;
  Array1D<Line> m_YZContours;
  int           m_minXYContour; 
  int           m_minXZContour; 
  int           m_minYZContour; 

  // curvature
  Array1D<double> mCurvature;

  // status bits
  bool mNormDirty;
  bool mContDirty;
  bool mCurveDirty;

  Array1D<Point> mProbes;
  Array1D<char *> mProbeNames;

private:
  // function to extract contour traces for each plane
  bool genContours (int xSlice, int ySlice, int zSlice);

  // other contour functions
  void findMinMax(float *xMin, float *xMax,
		  float *yMin, float *yMax,
		  float *zMin, float *zMax, int facet) const;
  void doFindMinMax (float *min, float *max,
		     double val0, double val1, double val2) const;
  bool findYZIntersect(Point *p1, Point *p2, int facet, double xVal) const;
  bool findXZIntersect(Point *p1, Point *p2, int facet, double yVal) const;
  bool findXYIntersect(Point *p1, Point *p2, int facet, double zVal) const;

  bool YZContourExists(int xslice) const;
  bool XZContourExists(int yslice) const;
  bool XYContourExists(int zslice) const;

  // curvature fcns
  static void getBasisVectors(Point &b1, Point &b2, Point &b3, const Point &n);
  static void getPlane(double &A, double &B, double &C, double &D,
                        const Point &p, const Point &n);

};


//
// inline bodies
//
inline void Surface::setVerbose(int verbose = 1)
{
  mVerbose = (verbose != 0);
}

inline bool Surface::getVerbose()
{
  return(mVerbose);
}

#endif // __SURFACE_H__
