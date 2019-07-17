///////////////////////////////////////////////////////////////////////////
//
// Surface utilities
//
///////////////////////////////////////////////////////////////////////////

#include <iostream.h>
#include <ADL/Array1D.h>
#include <ADL/Matrix.h>
#include <TDL/Point.h>
#include <TDL/ManifoldTransformation.h>
#include <TDL/SurfaceUtils.h>

//
// unitSphere
//
// create a unit sphere surface of the given level and return it
//
// based on "sphere" by Jon Leech (leech@cs.unc.edu) 3/24/89
//
// generate a triangle mesh approximating a sphere by
// recursive subdivision. First approximation is an octahedron;
// each level of refinement increases the number of triangles by
// a factor of 4.
//
void SurfaceUtils::unitSphere(Surface &surf, int maxlevel)
{
  surf.clean();

  // set up initial vertices of octahedron
  int XPLUS = surf.addVert( 1,  0,  0);
  int XMIN  = surf.addVert(-1,  0,  0);
  int YPLUS = surf.addVert( 0,  1,  0);
  int YMIN  = surf.addVert( 0, -1,  0);
  int ZPLUS = surf.addVert( 0,  0,  1);
  int ZMIN  = surf.addVert( 0,  0, -1);

  // set up initial facets of octahedron
  surf.addFacet(ZPLUS, XPLUS, YPLUS);
  surf.addFacet(ZPLUS, YPLUS, XMIN);
  surf.addFacet(ZPLUS, XMIN , YMIN);
  surf.addFacet(ZPLUS, YMIN , XPLUS);
  surf.addFacet(ZMIN , YPLUS, XPLUS);
  surf.addFacet(ZMIN , XPLUS, YMIN);
  surf.addFacet(ZMIN , YMIN , XMIN);
  surf.addFacet(ZMIN , XMIN , YPLUS);

  //
  // subdivide each starting triangle (maxlevel - 1) times
  //
  // Subdivide each triangle in the old approximation and normalize
  //  the new points thus generated to lie on the surface of the unit
  //  sphere.
  // Each input triangle with vertices labelled [0,1,2] as shown
  //  below will be turned into four new triangles:
  //
  //                      Make new points
  //                          a = (0+1)/2
  //                          b = (1+2)/2
  //                          c = (0+2)/2
  //        2
  //       /\             Normalize a, b, c
  //      /  \
  //    c/____\ b         Construct new triangles
  //    /\    /\              [0,a,c]
  //   /  \  /  \             [c,b,2]
  //  /____\/____\            [a,b,c]
  // 0      a     1           [a,1,b]
  //

  int level, i, oldnum, a, b, c;
  for (level = 1; level<maxlevel; level++) {
    oldnum = surf.mNumPoly;
    for(i=0; i<oldnum; i++) {
      // create new vertices a midpoints of edges
      a = surf.addVert(normalize(midpoint(
	surf.mVert[surf.mFacet[i][0]], surf.mVert[surf.mFacet[i][1]])));
      b = surf.addVert(normalize(midpoint(
	surf.mVert[surf.mFacet[i][1]], surf.mVert[surf.mFacet[i][2]])));
      c = surf.addVert(normalize(midpoint(
	surf.mVert[surf.mFacet[i][0]], surf.mVert[surf.mFacet[i][2]])));

      // add three new triangles
      surf.addFacet(c, b, surf.mFacet[i][2]);
      surf.addFacet(a, b, c);
      surf.addFacet(a, surf.mFacet[i][1], b);
      surf.mFacet[i][1] = a;
      surf.mFacet[i][2] = c;
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// surfaceArea -
//    use Heron's formula for surface area
//
///////////////////////////////////////////////////////////////////////////////

double SurfaceUtils::surfaceArea(Surface &S,
               float xPixelSize,
               float yPixelSize,
               float zPixelSize)
{
    if ((xPixelSize <= 0.0) || (xPixelSize <= 0.0) || (xPixelSize <= 0.0))
    {
        return -1.0;
    }

    int    i, ia, ib, ic;
    double a, b, c, s, pa, area = 0.0;
    double ax, ay, az, bx, by, bz, cx, cy, cz;
    for (i = 0; i < S.getNumPoly(); i++)
    {
        ia = S.facets()[i][0];
        ib = S.facets()[i][1];
        ic = S.facets()[i][2];

        ax = S.vertices()[ia].x() * xPixelSize;
        ay = S.vertices()[ia].y() * yPixelSize;
        az = S.vertices()[ia].z() * zPixelSize;
    
        bx = S.vertices()[ib].x() * xPixelSize;
        by = S.vertices()[ib].y() * yPixelSize;
        bz = S.vertices()[ib].z() * zPixelSize;
    
        cx = S.vertices()[ic].x() * xPixelSize;
        cy = S.vertices()[ic].y() * yPixelSize;
        cz = S.vertices()[ic].z() * zPixelSize;
    
        a = sqrt((bx-cx)*(bx-cx)+(by-cy)*(by-cy)+(bz-cz)*(bz-cz));
        b = sqrt((ax-cx)*(ax-cx)+(ay-cy)*(ay-cy)+(az-cz)*(az-cz));
        c = sqrt((bx-ax)*(bx-ax)+(by-ay)*(by-ay)+(bz-az)*(bz-az));
        s = (a + b + c)/2.0;

        pa = sqrt(s*(s-a)*(s-b)*(s-c));
        area += pa;
    }

    return area;
}
///////////////////////////////////////////////////////////////////////////
//
// maxCoordDist
//
// computes the maximum distance between coordinates (x, y, or z) of
// corresponding vertices of two surfaces which have the same structure
//
// returns 0 on success, 1 on failure
//
///////////////////////////////////////////////////////////////////////////
int SurfaceUtils::maxCoordDist (Surface const &S1,
				Surface const &S2,
				float *xMaxDist,
				float *yMaxDist,
				float *zMaxDist)
{
  float thisDist;
  *xMaxDist = *yMaxDist = *zMaxDist = 0.0;

  if (S1.mNumVert != S2.mNumVert || S1.mNumPoly != S2.mNumPoly)
  {
    cerr << "maxCoordDist: Surfaces do not have same structure" << endl;
    return 1;
  }

  for (int i=0; i<S1.mNumVert; i++)
  {
    thisDist = fabs(S1.mVert[i].x() - S2.mVert[i].x());
    if (thisDist > *xMaxDist)
      *xMaxDist = thisDist;

    thisDist = fabs(S1.mVert[i].y() - S2.mVert[i].y());
    if (thisDist > *yMaxDist)
      *yMaxDist = thisDist;

    thisDist = fabs(S1.mVert[i].z() - S2.mVert[i].z());
    if (thisDist > *zMaxDist)
      *zMaxDist = thisDist;
  }

  return 0;
}


///////////////////////////////////////////////////////////////////////////
//
// maxVertexDist
//
// computes the maximum distance between
// corresponding vertices of two surfaces which have the same structure
//
// returns distance on success, -1 on failure
//
///////////////////////////////////////////////////////////////////////////
float SurfaceUtils::maxVertexDist (Surface const &S1,
				   Surface const &S2)
{
  float thisDist;
  float maxDist = 0.0;

  if (S1.mNumVert != S2.mNumVert || S1.mNumPoly != S2.mNumPoly)
  {
    cerr << "maxVertexDist: Surfaces do not have same structure" << endl;
    return -1;
  }

  for (int i=0; i<S1.mNumVert; i++)
  {
    thisDist = (S1.mVert[i] - S2.mVert[i]).norm();
    if (thisDist > maxDist)
      maxDist = thisDist;
  }

  return maxDist;
}


///////////////////////////////////////////////////////////////////////////
//
// meanDist
//
// computes the mean of the distances between
// corresponding vertices of two surfaces which have the same structure
//
// returns distance on success, -1 on failure
//
///////////////////////////////////////////////////////////////////////////
float SurfaceUtils::meanDist (Surface const &S1,
			      Surface const &S2)
{
  float totDist = 0.0;

  if (S1.mNumVert != S2.mNumVert || S1.mNumPoly != S2.mNumPoly)
  {
    cerr << "meanDist: Surfaces do not have same structure" << endl;
    return -1;
  }

  for (int i=0; i<S1.mNumVert; i++)
  {
    totDist += (S1.mVert[i] - S2.mVert[i]).norm();
  }

  return totDist/S1.mNumVert;
}


///////////////////////////////////////////////////////////////////////////
//
// meanVectorDiff
//
// Computes the mean of the vector differences of
// corresponding vertices of two surfaces which have the same structure
//
// returns 0 on success, -1 on failure
//
///////////////////////////////////////////////////////////////////////////
int SurfaceUtils::meanVectorDiff (Surface const &S1,
				  Surface const &S2,
				  Point *mean)
{
  mean->set(0,0,0);

  if (S1.mNumVert != S2.mNumVert || S1.mNumPoly != S2.mNumPoly)
  {
    cerr << "meanVectorDiff: Surfaces do not have same structure" << endl;
    return -1;
  }

  for (int i=0; i<S1.mNumVert; i++)
  {
    *mean += (S1.mVert[i] - S2.mVert[i]);
  }

  *mean = *mean/S1.mNumVert;
  return 0;

}


///////////////////////////////////////////////////////////////////////////
//
// findRotation
//
// finds the optimal rotation/translation between two surfaces using 
// Froebenus norm
//
// Assumes that the surfaces have the same structure (same triangle list)
//
// returns 0 on success, -1 on failure
//
///////////////////////////////////////////////////////////////////////////
int SurfaceUtils::findRotation (Surface const &surf1,
				Surface const &surf2,
				Matrix<float> *R,
				Vector<float> *trans)
{
  int i;
  Point cent1, cent2;
  

  if (surf1.mNumVert != surf2.mNumVert || surf1.mNumPoly != surf2.mNumPoly)
  {
    cerr << "findRotation: Surfaces do not have same structure" << endl;
    return -1;
  }
  // find centroid of each surface
  cent1.set(0,0,0);
  cent2.set(0,0,0);

  for (i=0; i<surf1.mNumVert; i++)
  {
    cent1 += surf1.mVert[i];
    cent2 += surf2.mVert[i];
  }

  cent1 = cent1/(float)surf1.mNumVert; 
  cent2 = cent2/(float)surf2.mNumVert;

  Matrix<float>       C(3,3);
  C = 0;
  
  for (i=0; i<surf1.mNumVert; i++) 
    {
      C[0][0] +=(surf2.mVert[i].x()-cent2.x())*(surf1.mVert[i].x()-cent1.x());
      C[0][1] +=(surf2.mVert[i].x()-cent2.x())*(surf1.mVert[i].y()-cent1.y());
      C[0][2] +=(surf2.mVert[i].x()-cent2.x())*(surf1.mVert[i].z()-cent1.z());
      C[1][0] +=(surf2.mVert[i].y()-cent2.y())*(surf1.mVert[i].x()-cent1.x());
      C[1][1] +=(surf2.mVert[i].y()-cent2.y())*(surf1.mVert[i].y()-cent1.y());
      C[1][2] +=(surf2.mVert[i].y()-cent2.y())*(surf1.mVert[i].z()-cent1.z());
      C[2][0] +=(surf2.mVert[i].z()-cent2.z())*(surf1.mVert[i].x()-cent1.x());
      C[2][1] +=(surf2.mVert[i].z()-cent2.z())*(surf1.mVert[i].y()-cent1.y());
      C[2][2] +=(surf2.mVert[i].z()-cent2.z())*(surf1.mVert[i].z()-cent1.z());
    }

  Matrix <float> U(3,3),W(3,3),V(3,3),VT(3,3);

  MatrixUtils<float>::svd(C, &U, &W, &V);
  V.transpose(&VT);

  // Determine Rotation Matrix.
  // R = U * Vt
  // --------------------------
  MatrixUtils<float>::multiply(U, VT, R);
  
  // Check For Inversion
  // Note: Hard Coded For 3 Dimensions.
  // ----------------------------------
  if (MatrixUtils<float>::det(*R) < 0.0 ) {
                U[0][2] = -U[0][2];
                U[1][2] = -U[1][2];
                U[2][2] = -U[2][2];
                MatrixUtils<float>::multiply(U, VT, R);
    }

  // Now do the traslation

  
  (*trans)[0] = cent2.x() - ((*R)[0][0]*cent1.x()+(*R)[0][1]*cent1.y()+(*R)[0][2]*cent1.z());
  (*trans)[1] = cent2.y() - ((*R)[1][0]*cent1.x()+(*R)[1][1]*cent1.y()+(*R)[1][2]*cent1.z());
  (*trans)[2] = cent2.z() - ((*R)[2][0]*cent1.x()+(*R)[2][1]*cent1.y()+(*R)[2][2]*cent1.z());

    
  return 0;
}



///////////////////////////////////////////////////////////////////////////
//
// findAffine
//
// finds the optimal Affine/translation between two surfaces using 
//
// Assumes that the surfaces have the same structure (same triangle list)
//
// returns 0 on success, -1 on failure
//
///////////////////////////////////////////////////////////////////////////
int SurfaceUtils::findAffine (Surface const &surf1,
			      Surface const &surf2,
			      Matrix<float> *R,
			      Vector<float> *trans)
{
  int i,j,k;
  int numpoints; 

  if (surf1.mNumVert != surf2.mNumVert || surf1.mNumPoly != surf2.mNumPoly)
  {
    cerr << "findAffine: Surfaces do not have same structure" << endl;
    return -1;
  }

  numpoints = surf1.mNumVert;
	
  Matrix<double> A,B;
  Matrix<double> AA(4,3);
  Matrix<double> BB(4,3);
  Matrix<double> Q(4,4);
  A.setDim(4,numpoints);
  B.setDim(4,numpoints);

  for (i=0;i<numpoints;i++)
  {
    A[0][i] = (surf1.mVert[i]).x();
    A[1][i] = (surf1.mVert[i]).y();
    A[2][i] = (surf1.mVert[i]).z();
    A[3][i] = 1.0;

    B[0][i] = (surf2.mVert[i]).x();
    B[1][i] = (surf2.mVert[i]).y();
    B[2][i] = (surf2.mVert[i]).z();
    B[3][i] = 1.0;
  }
	
  for(j=0;j<4;j++)
    for(k=0;k<4;k++)
    {
      Q[j][k] = 0;
      for(i=0;i<numpoints;i++)
	Q[j][k] = Q[j][k] + A[j][i]*A[k][i];
    }

  // cout << "Q = " << Q;

  for(i=0;i<3;i++)
    for(j=0;j<4;j++)
    {
      AA[j][i] =0;
      for(k=0;k<numpoints;k++)
      {
	AA[j][i] = AA[j][i] + A[j][k]*B[i][k];
      }
    }

  // cout << "AA = " << AA;

  BB = MatrixUtils<double>::inv(Q)*AA;
	
  for (i=0;i<3;i++)
    for(j=0;j<3;j++)
      (*R)[j][i] = BB[i][j];

  (*trans)[0] = BB[3][0];
  (*trans)[1] = BB[3][1];
  (*trans)[2] = BB[3][2];

  return 0;
}



///////////////////////////////////////////////////////////////////////////
//
// applyManiTrans
//
// deform a surface by applying the manifold transformation to each
// vertex in the surface
//
///////////////////////////////////////////////////////////////////////////
int SurfaceUtils::applyManiTrans(ManifoldTransformation &mani,
				 Surface &S)
{
  cout << "Deforming Surface..." << flush;
  int i;
  int c, d = 0;
  double x0, x1, x2, y0, y1, y2;

  for (i=0; i<S.mNumVert; i++)
  {
    if ((c = (int)(i*100./S.mNumVert)) == d)
    {
      cout << c << "%.." << flush;
      d = d + 10;
    }
    x0 = S.mVert[i].x();
    x1 = S.mVert[i].y();
    x2 = S.mVert[i].z();
    mani.transformPoint(x0, x1, x2, y0, y1, y2);
    S.mVert[i].set(y0, y1, y2);
  }

  cout << "Complete" << endl;

  // invalidate normals, curvature, etc.
  S.verticesChanged();
  return 0;    
}



int SurfaceUtils::deformToPoints(Surface &S,int numFixedPoints, Pnt *points)
{
	int i,j;
	
	Array1D<int> ListofVert;
	
	ListofVert.setDim(numFixedPoints);

	// Find the closet vertex to each of the points; 

	double dist; 
	double d;
	Point P1;
	float x,y,z;
	
	for(i=0;i<numFixedPoints;i++){
	dist = 100000;
		for(j=0;j<S.mNumVert;j++){
		x = points[i].x() - S.mVert[j].x();	
		y = points[i].y() - S.mVert[j].y();	
		z = points[i].z() - S.mVert[j].z();	
		P1.set(x,y,z);
		d = P1.norm();
		if (d < dist ) {
			dist=d;
			 ListofVert[i] = j;
			}
		}		
	cout << "min dist = " << d << endl;
	cout << "Vertex " <<  ListofVert[i] << " Found" << endl;
	cout << "Point = " <<  points[i] << endl;
	cout << "Vertex = " <<  S.mVert[ListofVert[i]] << endl;
	}

	// Now set the deformation field at the vertexs.

	S.genUnitNormals();

	if (!S.hasNbhd()) {
    	if (S.getVerbose()) cerr << "\nGenerating Neighborhoods" << endl;
    	S.genNeighborhoods();
  	}


	double *deform,*newdeform;
	double *fixdeform;
	fixdeform = new double[numFixedPoints];
	deform = new double[S.mNumVert];
	newdeform = new double[S.mNumVert];
	for (i=0;i<S.mNumVert;i++) deform[i] = 0;

	for(i=0;i<numFixedPoints;i++){
	j = ListofVert[i];
		x = points[i].x() - S.mVert[j].x();     
                y = points[i].y() - S.mVert[j].y();     
                z = points[i].z() - S.mVert[j].z();     
		P1.set(x,y,z);
	deform[j] = P1.innerProd(S.mUNorm[j]);
	fixdeform[i] = P1.innerProd(S.mUNorm[j]);
		}	

	// Now do the estimation.
	int l;
	float dx;
	int n,k;

	float delta = 0.01;

	for(l=0; l<500; l++) {

	float anorm = 0.0;
	for (i=0;i< S.mNumVert;i++) {
		newdeform[i] = 0;
       		for(j = 0; j < (S.mNbhd[i]).numNeighbors(); j++)
       		{
			n = (S.mNbhd[i]).getNeighbor(j);
			P1 = S.mVert[i]-S.mVert[n];
        		dx =  P1.norm();
        		newdeform[i] += (deform[i]-deform[n])/dx;
			if (i==n) {
			cout <<"Error i==n" << endl;
			}
       		}
	}	
	for (i=0;i< S.mNumVert;i++) {
	float beta_const = 0;
	int notinlist = 1;
		for (k=0;k<numFixedPoints;k++){
		if (ListofVert[k] == i ) notinlist = 0;
		}
		if (notinlist) {
		for(j = 0; j < (S.mNbhd[i]).numNeighbors(); j++)
                {
                        n = (S.mNbhd[i]).getNeighbor(j);
                        P1 = S.mVert[i]-S.mVert[n];
                        dx =  P1.norm();
                        beta_const += (newdeform[i]-newdeform[n])/dx;
                        if (i==n) {
                        cout <<"Error i==n" << endl;
                        }
                }
        deform[i] = deform[i]-delta*beta_const;
	anorm += fabs(beta_const);
		}
        }
	cout << "Anorm = " << anorm << endl;

	}

	for(j=0; j<S.mNumVert; j++){
	S.mVert[j] += (S.mUNorm[j])*((float)deform[j]);
	}

	S.verticesChanged();

	return 0;
}
//
//
// Surface S is assumed to be a unit sphere!!
//
//

void SurfaceUtils::deformSphereToPoints(int numFixedPoints, Pnt
*points,Surface &S)
{
        int i,j;

	// Note that S is assumed to be a unit Sphere!!
	// Not sure what will happen if S is not a Sphere!!
	// It is not checked so make sure you call it with a unit Sphere only!!

	Pnt Center;
	Matrix<double> AA(3,1);
	Matrix<double> P(3,1);
        Vector<double> V(3);
        Vector<double> H(3);
	AA = 0;
	double B=0;
        double t=0;

	for(i=0;i<numFixedPoints;i++){
	V[0] = points[i].x();
        V[1] = points[i].y();
        V[2] = points[i].z();
        P = V;
        AA += P;
        B += V.dot(V);
        }
//	AA = 2.0*AA/(double)numFixedPoints;
	AA *= 2.0;
	AA /= (double)numFixedPoints;
        B = B/(double)numFixedPoints;

        Matrix<double> C(3,3);
        C = 0;
        Vector<double> D(3);
        D = 0;

	V = AA;
        for(i=0;i<numFixedPoints;i++){
        H[0] = points[i].x();
        H[1] = points[i].y();
        H[2] = points[i].z();
        P = H;
        C += AA*AA.transpose()-2.0*(AA*P.transpose())-2.0*(P*AA.transpose()) + 4.0*(P*P.transpose());
        D += 2.0*(H.dot(H))*H -2.0*B*H + B*V - (H.dot(H))*V;
        }

	Center = 0;
        if (MatrixUtils<double>::rank(C) < 3) {
        for(i=0;i<numFixedPoints;i++){
        Center += points[i];
        }
        Center = Center/numFixedPoints;
        } else {
        H = MatrixUtils<double>::inv(C)*D;
        Center.x() = H[0];
        Center.y() = H[1];
        Center.z() = H[2];
        }

	double R=1;

        if (numFixedPoints > 1) {
        R = 0;
        for(i=0;i<numFixedPoints;i++){
	R += Vector<Coord_t>::norm(points[i]-Center)*Vector<Coord_t>::norm(points[i]-Center);
        }
        R = sqrt(R/numFixedPoints);
	}

	S.scale(R);
        S.translate( Center.x(), Center.y(), Center.z());

	if (numFixedPoints < 5) return;

	Matrix<double> K(numFixedPoints,numFixedPoints);

	K = 0;

	int Numeig;
	
	Numeig = 50;

	// Build the Covariance Matrix.
	// Note that it is symetric.
	double A;
	Pnt P1,P2;
	for (i=0;i<numFixedPoints;i++)
		for(j=i;j<numFixedPoints;j++){
	// A is the solid angle between the P_i and P_j
	P1 = points[i]-Center;
	P2 = points[j]-Center;
A = P1.dot(P2)/(Vector<Coord_t>::norm((Vector<Coord_t> const &)P1)*Vector<Coord_t>::norm((Vector<Coord_t> const &)P2));
	K[i][j] = BuildCov((float) A,Numeig);
	K[j][i] = K[i][j];
	}

	cout.flush();
//	cout << "K = " << K << endl;

	Vector<double> Y(numFixedPoints);
	Vector<double> W(numFixedPoints);

	for (i=0;i<numFixedPoints;i++){
	P1 = points[i]-Center;
	Y[i] = Vector<Coord_t>::norm((Vector<Coord_t> const &)P1) - R;
	}

	// Now solve for the Weights

//	cout << "Y = " << Y << endl;

	W = MatrixUtils<double>::inv(K)*Y;

	// Now build the deformation field.
	
//	cout << "W = " << W << endl;

	double *deform;
	
	deform = new double[S.mNumVert];
	
	for (i=0;i<S.mNumVert;i++){
	P1.set(S.mVert[i].x(),S.mVert[i].y(),S.mVert[i].z());
	P1 = P1 - Center;
	deform[i] = 0;
		for(j=0;j<numFixedPoints;j++){
	P2 = points[j]-Center;
A = P1.dot(P2)/(Vector<Coord_t>::norm((Vector<Coord_t> const &)P1)*Vector<Coord_t>::norm((Vector<Coord_t> const &)P2));
		deform[i] += W[j]*BuildCov(A,Numeig);
		}
	}

	// Now deform the Surface!!

	S.genUnitNormals();

        for(j=0; j<S.mNumVert; j++){
        S.mVert[j] += (S.mUNorm[j])*((float)deform[j]);
        }

	S.verticesChanged();
//	S.genUnitNormals();

}	

float SurfaceUtils::BuildCov(float A,int N)
{
// Build A covariance associated with the Laplacian on the 
// Sphere with finite number of Legendre polynomials.

int i;
float c;
	c = 0;

	for (i=1;i<N;i++){
		c += plgndr(i,0,A)*(2*i+1)/((i*(i+1))*(i*(i+1)));
	}
	return c;
}
	

float SurfaceUtils::plgndr(int l,int m,float x)
{
        float fact,pll,pmm,pmmp1,somx2;
        int i,ll;
        void nrerror();

        if (m < 0 || m > l ) 
                {
                cout << "Bad arguments in routine PLGNDR" << endl;
                cout <<  "PLGNDR x = " << x << " M=" << m
		<< " l = " << l << endl;
                exit(0);
                }
	if(x > 1.0) x = 1.0;
	if(x < -1.0) x = -1.0;

        pmm=1.0;
        if (m > 0) {
                somx2=sqrt((1.0-x)*(1.0+x));
                fact=1.0;
                for (i=1;i<=m;i++) {
                        pmm *= -fact*somx2;
                        fact += 2.0;
                }
        }
        if (l == m)
                return pmm;
        else {
                pmmp1=x*(2*m+1)*pmm;
                if (l == (m+1))
                        return pmmp1;
                else {
                        for (ll=(m+2);ll<=l;ll++) {
				pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
                                pmm=pmmp1;
                                pmmp1=pll;
                        }
                        return pll;
                }
        }
}


//compute the derivative of Legendre polynomial usign the recursion
//P'_l(x)=l*P_(l-1)(x) +xP'_(l-1)(x)

float SurfaceUtils::plgndr_der(int l, float x) {
  if (l==0) return 0; 
  return (l*plgndr(l-1,0,x) + x*plgndr_der(l-1,x));
}



void SurfaceUtils::tagEdges(Surface &surf, Array1D<int> &edges)
{
if(!surf.hasNbhd())
        surf.genNeighborhoods();

surf.ShowWorking("locating edge points\n");
Array1D<int> used(surf.mNumVert);
used = 0;
int i,j;
for(i=0;i<surf.mNumPoly;i++) {
        used[surf.mFacet[i][0]]++;
        used[surf.mFacet[i][1]]++;
        used[surf.mFacet[i][2]]++;
        }

edges.setDim(surf.mNumVert);

for(i=0;i<surf.mNumVert;i++) {
        if(surf.mNbhd[i].numNeighbors() != used[i])
                edges[i] = 1;
        else    edges[i] = 0;
        }
}

//
// redistribute
//
// triangle redistribution algorithm from Sarang Joshi
//
ItXECode SurfaceUtils::redistribute(Surface &surf, int redist_iter, double redist_alpha)
{
surf.ShowWorking("Redistributing points\n");

if(!surf.hasUNorms())
        surf.genUnitNormals();

if(!surf.hasNbhd())
        surf.genNeighborhoods();

Array1D<double> newx(surf.mNumVert);
if(newx.isEmpty())
        return(ItXNoMemory);

Array1D<double> newy(surf.mNumVert);
if(newy.isEmpty())
        return(ItXNoMemory);

Array1D<double> newz(surf.mNumVert);
if(newz.isEmpty())
        return(ItXNoMemory);

Array1D<int> isedge(surf.mNumVert);

if(isedge.isEmpty())
        return(ItXNoMemory);

tagEdges(surf,isedge);

float dwork = 1.0/((float)redist_iter);
float work_perc = 0.0;

int j,l;

for(l=0;l<redist_iter;l++,work_perc += dwork) {
   surf.ShowWorking(work_perc);
   for(j=0;j<surf.mNumVert;j++) {
        if(isedge[j]) continue;

        double px  = surf.mVert[j].x();
        double py  = surf.mVert[j].y();
        double pz  = surf.mVert[j].z();

        /* For each point the derivative is
        just average the neighbours minus itself!! */

        int num = surf.mNbhd[j].numNeighbors();

        double tmpx= -num*px;
        double tmpy= -num*py;
        double tmpz= -num*pz;

        for(int i=0;i<num;i++) {
                int nn = surf.mNbhd[j].getNeighbor(i);
                tmpx += surf.mVert[nn].x();
                tmpy += surf.mVert[nn].y();
                tmpz += surf.mVert[nn].z();
                }

        double nrmx = surf.mUNorm[j].x();
        double nrmy = surf.mUNorm[j].y();
        double nrmz = surf.mUNorm[j].z();

        double in   = tmpx*nrmx+tmpy*nrmy+tmpz*nrmz;

        newx[j] = px + redist_alpha*(tmpx - in*nrmx);
        newy[j] = py + redist_alpha*(tmpy - in*nrmy);
        newz[j] = pz + redist_alpha*(tmpz - in*nrmz);
        }

    for(j=0;j<surf.mNumVert;j++)
      if(!isedge[j])
        surf.mVert[j].set(newx[j],newy[j],newz[j]);

    }

surf.verticesChanged();

return(ItXSuccess);
}



ItXECode SurfaceUtils::Ill1(Matrix<float> &r, Matrix<float> &b, float thresh)
{
char JOBU  = 'A';
char JOBVT = 'A';
int  m     = b.getNrow();
int  M     = m;
int  N     = m;
int LDA    = m;

Matrix<float>   S(m,m);
Matrix<float>   U(m,m);
Matrix<float>  VT(m,m);
Matrix<float> INP(m,m);
Array1D<float> WORK(4*m*4*m);

int LDU   = m;
int LDVT  = m;
int LWORK = 5*m;
int INFO;

b.transpose(&INP);

S = 0;
Lapack::xgesvd(JOBU,JOBVT,M,N,INP.data(),LDA,S.data(),U.data(),LDU,
                VT.data(),LDVT,WORK.data(),LWORK,&INFO);

if(INFO != 0)
  return(ItXError);

int i,j;

for(i=0;i<m;i++) {
        float temp = S[i][i];
        S[i][i] = S[0][i];
        S[0][i] = temp;
        }

float S00 = S[0][0];
if(S00 == 0.0) {
        fprintf(stderr,"S00 = ZERO!\n");
        S00 = 1.0;
        }
for(i=1;i<m;i++) {
        if((S[i][i]/S00) < thresh)
                S[i][i] = S[i-1][i-1];
        else    S[i][i] = 1.0/S[i][i];
        }
Matrix<float> tVT(m,m);
Matrix<float> tS(m,m);
Matrix<float> tU(m,m);

for(i=0;i<m;i++)
  for(j=0;j<m;j++) {
        tS[i][j]  =  S[i][j];
        tVT[i][j] = VT[i][j];
        tU[i][j]  =  U[i][j];
        }

Matrix<float> OUT(m,m);
OUT = tVT * tS;
r   = OUT * tU;

return(ItXSuccess);
}


void SurfaceUtils::smooth(Surface &surf, int window_size)
{
if(!surf.hasUNorms())
    surf.genUnitNormals();

if(!surf.hasNbhd())
    surf.genNeighborhoods(window_size);

int i,j,nei,nnei,maxn;

Array1D<Point> newVert(surf.mNumVert);
Array1D<int>   newCount(surf.mNumVert);

maxn = -1;
for(i=0;i<surf.mNumVert;i++) {
   newVert[i].set(0.0,0.0,0.0);
   newCount[i] = 0;
   if((nnei = surf.mNbhd[i].numNeighbors()) > maxn)
        maxn = nnei;
   }

Matrix<float> U(maxn,3);
Matrix<float> H(maxn,1);
Vector<float> h(maxn);
Vector<float> u(maxn);
Vector<float> v(maxn);

H = 0.0;
U = 0.0;

double A,B,C,D;
Point  temp1,temp2,temp3,temp4,temp5,b1,b2;
double ptx,pty,ptz;
double nrmx,nrmy,nrmz;
double b1x,b1y,b1z;
double nptx,npty,nptz;

Point zero;
zero.set(0.0,0.0,0.0);

b2 = zero;


float dwork = 1.0/((float)surf.mNumVert);
float work_perc = 0.0;

for(i=0;i<surf.mNumVert;i++, work_perc += dwork) {
        if(i % 100 == 0) surf.ShowWorking(work_perc);
        ptx  = surf.mVert[i].x();
        pty  = surf.mVert[i].y();
        ptz  = surf.mVert[i].z();
        nrmx = surf.mUNorm[i].x();
        nrmy = surf.mUNorm[i].y();
        nrmz = surf.mUNorm[i].z();

        A = nrmx;
        B = nrmy;
        C = nrmz;
        D = -(ptx*nrmx + pty*nrmy + ptz*nrmz);

        b1x = b1y = b1z = 1.0;
        if(nrmx != 0.0)
                b1x = -(nrmy + nrmz)/nrmx;
        else if(nrmy != 0.0)
                b1y = -(nrmx + nrmz)/nrmy;
        else if(nrmz != 0.0)
                b1z = -(nrmx + nrmy)/nrmz;

        // normalize
        double mag = sqrt(b1x*b1x + b1y*b1y + b1z*b1z);
        if(mag != 0.0) {
                b1x /= mag;
                b1y /= mag;
                b1z /= mag;
                }

        b1.set(b1x,b1y,b1z);

        //veccross(b1,surf.mUNorm[i],b1);
        double erx,ery,erz; 
        erx = (b1.y()*surf.mUNorm[i].z())-(b1.z()*surf.mUNorm[i].y()); 
        ery = (b1.z()*surf.mUNorm[i].x())-(surf.mUNorm[i].z()*b1.x()); 
        erz = (b1.x()*surf.mUNorm[i].y())-(surf.mUNorm[i].x()*b1.y()); 
        b1.set(erx,ery,erz); 

//      b1 = surf.mUNorm[i].cross(b1);

        if((b1 == zero)&&(b2 == zero))
                continue;

//      U = 0.0;
        for(j=0;j<surf.mNbhd[i].numNeighbors();j++) {

                nei = surf.mNbhd[i].getNeighbor(j);
                if((nei < 0)||(nei >= surf.mNumVert)) {
                        cerr << "Invalid Neighbor!" << endl;
                        continue;
                        }

                nptx = surf.mVert[nei].x();
                npty = surf.mVert[nei].y();
                nptz = surf.mVert[nei].z();

                h[j] = A * nptx + B * npty + C * nptz + D;

                temp1 = surf.mUNorm[i] * h[j];
                temp2 = surf.mVert[nei] - temp1;
                temp3 = temp2 - surf.mVert[i];

                u[j] = temp3.innerProd(b1);
                v[j] = temp3.innerProd(b2);

                U[j][0] = u[j]*u[j];
                U[j][1] = u[j]*v[j]*2.0;
                U[j][2] = v[j]*v[j];
                }

        Matrix<float> U_t = U.transpose();
        Matrix<float> dumb3 = U_t * U;

        Matrix<float> dumb4;
        ItXECode rv = Ill1(dumb4,dumb3,0.001);
        if(rv != ItXSuccess)
                continue;

        Matrix<float> dumb5 = dumb4 * U_t * H;

        double C_01 = dumb5[0][0];
        double C_11 = dumb5[1][0];
        double C_02 = dumb5[2][0];

        for(j=0;j<surf.mNbhd[i].numNeighbors();j++) {
                nei = surf.mNbhd[i].getNeighbor(j);
                if(nei < 0) continue;

                temp1 = b1 * u[j];
                temp2 = b2 * v[j];

                double height = 0.5 * (C_01*u[j]*u[j] + 2.0*C_11*u[j]*v[j] +
                                C_02*v[j]*v[j]);

                temp3 = surf.mUNorm[i] * height;

                temp4.set(temp1.x() + temp2.x(),
                          temp1.y() + temp2.y(),
                          temp1.z() + temp2.z());
                temp5.set(temp3.x() + temp4.x(),
                          temp3.y() + temp4.y(),
                          temp3.z() + temp4.z());

                newVert[nei] += (surf.mVert[i] + temp5);
                newCount[nei]++;
                }
        }

for(i=0;i<surf.mNumVert;i++) {
//        double nnum = surf.mNbhd[i].numNeighbors();
        double nnum = newCount[i];
        if(nnum <= 0.0) continue;
        surf.mVert[i] = newVert[i]/nnum;
        }

surf.verticesChanged();
}


///////////////////////////////////////////////////////////////////////////
//
// smoothSurfaceAvg
//
// smooth a Surface by averaging. Used also to reduce any closed surface to
// a sphere.
//
// Assumes that the surfaces have the same structure (same triangle list)
//
// returns 0 on success, -1 on failure
//
///////////////////////////////////////////////////////////////////////////
ItXECode SurfaceUtils::smoothSurfaceAvg (Surface &surf,
                                      Surface &S,int numiter)
{

int i,j,k;

  if (!surf.hasNbhd())
  {
    if (surf.getVerbose()) cerr << "\nGenerating Neighborhoods" << endl;
    surf.genNeighborhoods();
  }

Surface temp;

Point Avg;

// copy surf in to S 
	S = surf;

	// For each vertex move it to the average of it and it's neighbours.
	for (k=0;k<numiter;k++) {
	// Copy the result of the
	temp = S;
	for (i=0;i<S.mNumVert;i++) {
	// Avg = temp.mVert[i];
	   Avg.set(0,0,0);
		for (j=0;j<(S.mNbhd[i]).numNeighbors();j++) {
	// Add the jth neighbour.
		Avg += temp.mVert[(temp.mNbhd[i]).getNeighbor(j)]; 			
	}	
	Avg = Avg / ((temp.mNbhd[i]).numNeighbors());

	S.mVert[i] = Avg;
		}	
	}

	surf.verticesChanged();
	return ItXSuccess;
}

ItXECode SurfaceUtils::normalizeSurface(Surface &surf, Surface &S,double tv) 
{
	Point Centroid;

	// Get the Centroid 

	//surf.genNormals();

	surf.getSimpleCentroid(&Centroid);

	// Compute the "Elpiticity of the surface" 

	// Need to look up the formal def.
	Matrix<double> A(3,3);

	Vector<double> C(3);
	C[0] = -Centroid.x();
	C[1] = -Centroid.y();
	C[2] = -Centroid.z();
	A.eye();
	
	// Move the surface to the origin

	S = surf;

	S.affine(A,C);
	S.genNormals();
	double vol = S.volume();
	if (vol==0) vol=1;
	double xx = 0,xy =0,xz = 0,yy = 0,yz = 0,zz = 0;
	int i;
	for(i=0; i<S.mNumVert; i++) {
	  double xc = (S.mVert[i]).x();
	  double yc = (S.mVert[i]).y();
	  double zc = (S.mVert[i]).z();
	  xx +=  xc*xc*xc*(S.mNorm[i]).x()/3 +
		xc*xc*yc*(S.mNorm[i]).y() + xc*xc*zc*(S.mNorm[i]).z() ;
	  xy +=  xc*xc*yc*(S.mNorm[i]).x()/2 + xc * yc * yc*(S.mNorm[i]).y()/2
		+ xc*yc*zc*(S.mNorm[i]).z();
	  xz +=  xc*xc*zc*(S.mNorm[i]).x()/2 + xc * zc * zc*(S.mNorm[i]).z()/2
		+ xc*zc*yc*(S.mNorm[i]).y();	

	  yy +=  yc*yc*yc*(S.mNorm[i]).y()/3 + xc*yc*yc*(S.mNorm[i]).x()
		+ zc*yc*yc*(S.mNorm[i]).z();
	  yz +=  yc*yc*zc*(S.mNorm[i]).y() + yc * zc * zc*(S.mNorm[i]).z()
		+ yc*zc*xc*(S.mNorm[i]).x();

	  zz +=  zc*zc*zc*(S.mNorm[i]).z()/3 + zc*zc*xc*(S.mNorm[i]).x()
		+ zc*zc*yc*(S.mNorm[i]).y();
	}


	Matrix<double> B(3,3);

	A[0][0] = xx;
	A[0][1] = A[1][0] = xy;
	A[0][2] = A[2][0] = xz;
	A[1][1] = yy;
	A[1][2] = A[2][1] = yz;
	A[2][2] = zz;

        A = A * A.transpose();

//	cout << "Centroid = " << Centroid << endl;


	Matrix<double> U(3,3);
  	Matrix<double> E(3,3);
  	Matrix<double> V(3,3);

//	cout << "A =" << A << endl; 
//	cout << "A(norm) =" << A << endl; 

  	MatrixUtils<double>::svd(A, &U, &E, &V);
  	for (i=0; i<3; i++) {
    		E[i][i] = 1/sqrt(E[i][i]);
    		E[i][i] = sqrt(E[i][i]);
		}
  	B = U * E * V.transpose();

	double dt = MatrixUtils<double>::det(B);
	
	// Take the cube root of the determinent
	dt=exp(log(dt*vol/tv)/3); 

	B = (B/dt);

//	cout << "sqrtA.inv = " << B << endl;

	C = 0;
	S.affine(B,C);

	S.verticesChanged();

	return ItXSuccess;
}	


// Project a surface to a sphere
//

ItXECode SurfaceUtils::projectSurfaceToSphere(Surface &S) 
{

	for (int i=0; i< S.mNumVert; i++)
        {
	S.mVert[i] = S.mVert[i]/(S.mVert[i]).norm();
	}

	S.verticesChanged();

	return ItXSuccess;
}	

ItXECode SurfaceUtils::projectSurfaceToSphere(Surface &S, Surface &out) 
{

	for (int i=0; i< S.mNumVert; i++)
        {
	out.mVert[i] = S.mVert[i]/(S.mVert[i]).norm();
	}

	out.verticesChanged();

	return ItXSuccess;
}	

ItXECode SurfaceUtils::closedSurfaceToSphere(Surface &surf, Surface &S,int numiter)
{
	
	int el;
	
	el = surf.euler();

//  	if ( el !=2 ) {
//  		return ItXError;
//  	}

	Surface surf1;
	Surface surf2;

	surf1 = surf;

	for (int i=0; i < numiter/10; i++) {
        SurfaceUtils::normalizeSurface(surf1, surf2,1.0);
        SurfaceUtils::smoothSurfaceAvg(surf2, surf1, 10);
        }

        SurfaceUtils::normalizeSurface(surf1, surf2,1.0);
        SurfaceUtils::projectSurfaceToSphere(surf2,surf1);

        S = surf1;
	S.verticesChanged();
	
	return ItXSuccess;
	
}

//change rectangular coordinates to spherical
//note theta returns 0 when phi=0,M_PI (at the poles) 

void SurfaceUtils::rectToSph(Point const &P,double &theta, double &phi) {
  
  double rad=sqrt (P.x()*P.x()+P.y()*P.y()+P.z()*P.z());

  if (rad!=0) 
    phi=acos(P.z()/rad);
  //else return;
  
  //if at either one of poles, assign 0 to latitude 
  if (sin(phi)==0)
   {  theta=0;
}  else  {
    double tmpx=P.x()/(rad*sin(phi));
    double tmpy=P.y()/(rad*sin(phi));
    
    if (tmpx>1) tmpx=1;
    if (tmpx<-1) tmpx=-1;

    if (tmpy>1) tmpy=1;
    if (tmpy<-1) tmpy=-1;
    
    if (tmpy>=0)  //1st and 2nd quadrants
      theta=acos(tmpx);
    else {
      if (tmpx>=0)   //4th quadrant
         theta=asin(tmpy)+2*M_PI; 
      else    //3rd quadrant
      theta=M_PI-asin(tmpy);
    }
  }  //end else
}
 
//change spherical coordinates to rectangular

void SurfaceUtils::sphToRect(Point &P, double const &rad, double const &theta,
			     double const &phi) {
  P.set(rad*sin(phi)*cos(theta),rad*sin(phi)*sin(theta),rad*cos(phi));

}

//change rectangular coords to sterographic projection
//note center shifted to (0,0,1) in rectangular

void SurfaceUtils::rectToStereo(Point const &P,double &u, double &v) {
    u=(2*P.x())/(2-P.z());
    v=(2*P.y())/(2-P.z());

}

//change sterographic projection coords to rectangular
//note center is at (0,0,1)
void SurfaceUtils::stereoToRect(Point &P, double const &u, double const &v) {
  P.set(4*u/(u*u+v*v+4),4*v/(u*u+v*v+4),2*(u*u+v*v)/(u*u+v*v+4));   
	
}

//find the basis for the tangent space of sphere at a given point
//in stereroographic coordinates

void SurfaceUtils::findStereoBasis(Point const &P, Point &basis_u, 
				   Point &basis_v) {
  double u,v;
  SurfaceUtils::rectToStereo(P,u,v);
  double f=u*u+v*v;
  basis_u.set( (4*v*v-4*u*u+16)/((f+4)*(f+4)),
	       (-8*u*v)/((f+4)*(f+4)),
	       (16*u)/((f+4)*(f+4))     );
  
  basis_v.set(  (-8*u*v)/((f+4)*(f+4)),
		(4*u*u-4*v*v+16)/((f+4)*(f+4)),
		(16*v)/((f+4)*(f+4))      );

  
  basis_u/=(basis_u.norm());
  basis_v/=(basis_v.norm());
 
}

//rotate point P on sphere to point Q in N steps (N-segment geodesic on
//the sphere between P and Q

void SurfaceUtils::rotatePtToPt(Point const &P, Point const &Q, 
				Array1D<Point> &path, int const &N) {
  
  path[0].set(P.x(),P.y(), P.z());
 
  double sol_ang=P.x()*Q.x()+P.y()*Q.y()+P.z()*Q.z();
  if (sol_ang>1) sol_ang=1.0;
  if (sol_ang<-1) sol_ang=-1.0;
  sol_ang=acos(sol_ang)/(N-1);

  
  //calculate axis of rotation P x Q
  Point rotAx=P.cross(Q);
  rotAx/=rotAx.norm();

  //if 180 degrees apart pick the axis of rotation as +z
  if (rotAx==Point(0.0,0.0,0.0)) rotAx.set(0.0,0.0,1.0);
 
  
  double theta,phi;
  SurfaceUtils::rectToSph(rotAx,theta,phi);

//rotate P so that rotAx is aligned w/ +z-axis

  //step1: rotate by -theta around z-axis
  
  Point temp1;
  
  temp1.set(P.x()*cos(theta)+P.y()*sin(theta),
	   -P.x()*sin(theta)+P.y()*cos(theta),
	   P.z());

  //step 2: rotate by -phi around y-axis
  Point temp2;
  temp2.set(temp1.x()*cos(phi)-temp1.z()*sin(phi),
	    temp1.y(),
	    temp1.x()*sin(phi)+temp1.z()*cos(phi));

  temp1=temp2;
  //now rotate by sol_ang around the z-axis
  
  Array1D<Point> arc1(N), arc2(N);
  int i; 
  for (i=1; i<N; ++i) {
    arc1[i].set( temp1.x()*cos(sol_ang)-temp1.y()*sin(sol_ang),
		temp1.x()*sin(sol_ang)+temp1.y()*cos(sol_ang),
		temp1.z() ); 
    
    temp1=arc1[i];
  }

 
   
  for (i=1; i<N; ++i) {
    //rotate py phi around y-axis to reverse the transformation
  
    arc2[i].set( arc1[i].x()*cos(phi) + arc1[i].z()*sin(phi),
		 arc1[i].y(),
		 arc1[i].x()*-sin(phi) + arc1[i].z()*cos(phi) );

    //rotate by theta around z-axis

    path[i].set( arc2[i].x()*cos(theta)-arc2[i].y()*sin(theta),
		 arc2[i].x()*sin(theta)+arc2[i].y()*cos(theta),
		 arc2[i].z() );

  }
   
}

//align the arbitrary fixed pt on a given sphere w/ the north pole

void SurfaceUtils::alignNP(Point &NP, Surface &surf) {
  
  double rad=sqrt(NP.x()*NP.x()+NP.y()*NP.y()+NP.z()*NP.z());
  NP/=rad;

  double theta,phi;

  SurfaceUtils::rectToSph(NP,theta,phi);

  //step 1: rotate by -theta around z-axis
  Point temp;

  for (int i=0; i<surf.getNumVert(); ++i) {
    temp=(surf.vertices())[i];
    surf.changeVert(i, temp.x()*cos(theta)+temp.y()*sin(theta),
		    -temp.x()*sin(theta)+temp.y()*cos(theta),
		    temp.z());
    
    temp=(surf.vertices())[i];
    //step 2: rotate by -phi around y-axis

    surf.changeVert(i, temp.x()*cos(phi)-temp.z()*sin(phi),
		    temp.y(),
		    temp.x()*sin(phi)+temp.z()*cos(phi) );
  }

  surf.verticesChanged();
}

//take the alignment back after deformation is done

void SurfaceUtils::unalignNP(Point &NP, Surface &surf) {
  
 double rad=sqrt(NP.x()*NP.x()+NP.y()*NP.y()+NP.z()*NP.z());
  NP/=rad;

  double theta,phi;

  SurfaceUtils::rectToSph(NP,theta,phi);

  theta=-theta;
  phi=-phi;
  
  Point temp;

  for (int i=0; i<surf.getNumVert(); ++i) {
    temp=(surf.vertices())[i];
    surf.changeVert(i, temp.x()*cos(phi)-temp.z()*sin(phi),
		    temp.y(),
		    temp.x()*sin(phi)+temp.z()*cos(phi) );
    temp=(surf.vertices())[i];
    surf.changeVert(i, temp.x()*cos(theta)+temp.y()*sin(theta),
		    -temp.x()*sin(theta)+temp.y()*cos(theta),
		    temp.z());
  }
  surf.verticesChanged();
   
}
		    
void SurfaceUtils::alignNP(Point &NP, Matrix<double> &landmks) {
 
 double rad=sqrt(NP.x()*NP.x()+NP.y()*NP.y()+NP.z()*NP.z());
 NP/=rad;
                    
 double theta,phi;
 SurfaceUtils::rectToSph(NP,theta,phi);
  
  //step 1: rotate by -theta around z-axis
  Point temp;
                    
  for (int i=0; i<landmks.getNrow(); ++i) {
    temp.set(landmks[i][0],landmks[i][1],landmks[i][2]);
    landmks[i][0]= temp.x()*cos(theta)+temp.y()*sin(theta);
    landmks[i][1]= -temp.x()*sin(theta)+temp.y()*cos(theta);
    landmks[i][2]=temp.z();
  
    temp.set(landmks[i][0],landmks[i][1],landmks[i][2]);

    //step 2: rotate by -phi around y-axis	  
    landmks[i][0]= temp.x()*cos(phi)-temp.z()*sin(phi);
    landmks[i][1]=temp.y();
    landmks[i][2]=temp.x()*sin(phi)+temp.z()*cos(phi);	
  }

}              

void SurfaceUtils::findSphereRotation(Matrix<double> &A1, Matrix<double>&B1,
Matrix<double> *R) {	
        int m = A1.getNrow();
        int n = A1.getNcol();
        int     i, j;
	//  Matrix<T>       A2(A1.getDim()),B2(B1.getDim());
        Matrix<double>       AT(A1.getNcol(), A1.getNrow());
        
	Matrix<double>       C(n,n);
        Matrix<double>       U(n,n);
        Matrix<double>       V(n,n), VT(n,n);
        Matrix<double>       W(n,n);

	if (( m != B1.getNrow()  )   || ( n != B1.getNcol()  )
	    ||  ( n != R->getNcol()  )   || ( n != R->getNrow()  ) )
	  throwHardError("Matrix::find_transform: bad array dimensions");
	
        // Compute rotation matrix using SVD
        // C = Bt * A
        // ----------
        A1.transpose(&AT);
        MatrixUtils<double>::multiply(AT, B1, &C);
            
        // Decompose C Using SVD.
        // C = U * W * Vt
        // ----------------------
        MatrixUtils<double>::svd(C, &U, &W, &V);
	V.transpose(&VT);
	// Determine Rotation Matrix.
       
	
        // --------------------------
        MatrixUtils<double>::multiply(U, VT, R);  
	// Check For Inversion   
        // Note: Hard Coded For 3 Dimensions.
        // ----------------------------------
	if (MatrixUtils<double>::det(*R) < 0.0 && n == 3) {
	
                U[0][2] = -U[0][2];
                U[1][2] = -U[1][2];
                U[2][2] = -U[2][2];
                MatrixUtils<double>::multiply(U, VT, R);  
    }
}           


// KWD
ItXECode SurfaceUtils::removeSurfacePoints(Surface &S, Array1D<u_char> &isDel)
{

 if(!&S) return(ItXError);

  Surface newS;

  int i,a,b,c,cnt;

  int delv,nv = S.getNumVert();
  int delf,nf = S.getNumPoly();

  Array1D<int> newIndex(nv);

  delv = 0;
  for(i=0;i<nv;i++) {
        newIndex[i] = i - delv;
        if(isDel[i]) delv++;
        }

  delf = 0;
  for(i=0;i<nf;i++) {
        a = S.facets()[i][0];
        b = S.facets()[i][1];
        c = S.facets()[i][2];
        if((isDel[a])||(isDel[b])||(isDel[c]))
                delf++;
  }
 
  if((delv <= 0)&&(delf <= 0))
        return(ItXSuccess);

  if(((nv - delv) <= 0)||((nf - delf) <= 0)) {
	S.clean();
        return(ItXSuccess);
  }

  newS.setSize(nv - delv,nf - delf);

  cnt = 0;
  for(i=0;i<nv;i++)
        if(!isDel[i])
           newS.vertices()[cnt++] = S.vertices()[i];

  cnt = 0;
  for(i=0;i<nf;i++) {
        a = S.facets()[i][0];
        b = S.facets()[i][1];
        c = S.facets()[i][2];
        if((!isDel[a])&&(!isDel[b])&&(!isDel[c])) {
           newS.facets()[cnt][0] = newIndex[a];
           newS.facets()[cnt][1] = newIndex[b];
           newS.facets()[cnt][2] = newIndex[c];
           cnt++;
           }
  }
  S = newS;
  S.geometryChanged();
  return(ItXSuccess);
}


// Lei 01/22/2004
// keep the triangle indices, just show the extracted ones, however
ItXECode SurfaceUtils::keepSurfacePoints(Surface &S, Array1D<u_char> &isKeep)
{

 if(!&S) return(ItXError);

  Surface newS;

  int i,a,b,c,cnt;

// delv and delf are really keepv and keepf
  int delv,nv = S.getNumVert();
  int delf,nf = S.getNumPoly();

  Array1D<int> newIndex(nv);

  delv = 0;
  for(i=0;i<nv;i++) {
        newIndex[i] = i - delv;
        if(isKeep[i]) delv++;
        }

  delf = 0;
  for(i=0;i<nf;i++) {
        a = S.facets()[i][0];
        b = S.facets()[i][1];
        c = S.facets()[i][2];
        if((isKeep[a])||(isKeep[b])||(isKeep[c]))
                delf++;
  }
 
  if((delv <= 0)&&(delf <= 0))
        return(ItXSuccess);

  if(((nv - delv) <= 0)||((nf - delf) <= 0)) {
	S.clean();
        return(ItXSuccess);
  }

  //  newS.setSize(nv - delv,nf - delf);
  // Lei: new surface has the same number of vertices and 
  //      number of facets to KEEP not to DELETE
  newS.setSize(nv,nf - delf);

  cnt = 0;
  for(i=0;i<nv;i++)
    //        if(!isDel[i])
// Lei: comment out the above to use all the original coordinates
           newS.vertices()[cnt++] = S.vertices()[i];

  cnt = 0;
  for(i=0;i<nf;i++) {
        a = S.facets()[i][0];
        b = S.facets()[i][1];
        c = S.facets()[i][2];
        if((!isKeep[a])&&(!isKeep[b])&&(!isKeep[c])) {
	  //           newS.facets()[cnt][0] = newIndex[a];
	  //           newS.facets()[cnt][1] = newIndex[b];
	  //           newS.facets()[cnt][2] = newIndex[c];
	  // Lei: comment out above to use original indices
	  newS.facets()[cnt][0] = a;
	  newS.facets()[cnt][1] = b;
	  newS.facets()[cnt][2] = c;
           cnt++;
           }
  }
  S = newS;
  S.geometryChanged();
  return(ItXSuccess);
}


ItXECode SurfaceUtils::removeSurfacePolys(Surface &S, Array1D<u_char> &isDel)
{
 if(!&S) return(ItXError);

  Surface newS;

  int i,a,b,c,cnt;

  int delv,nv = S.getNumVert();
  int delf,nf = S.getNumPoly();

  Array1D<u_char> ptUsed(nv);

  ptUsed = 0U;
  delf   = 0;
  for(i=0;i<nf;i++) {
        if(isDel[i]) 
           delf++;
        else {
           a = S.facets()[i][0];
           b = S.facets()[i][1];
           c = S.facets()[i][2];
           ptUsed[a] = 1U;
           ptUsed[b] = 1U;
           ptUsed[c] = 1U;
        }
  }

  Array1D<int> newIndex(nv);

  delv = 0;
  for(i=0;i<nv;i++) {
     newIndex[i] = i - delv;
     if(!ptUsed[i])
         delv++;
  }

  if((delv <= 0)&&(delf <= 0))
        return(ItXSuccess);

  if(((nv - delv) <= 0)||((nf - delf) <= 0)) {
	S.clean();
        return(ItXSuccess);
  }

  newS.setSize(nv - delv,nf - delf);

  cnt = 0;
  for(i=0;i<nv;i++)
        if(ptUsed[i])
           newS.vertices()[cnt++] = S.vertices()[i];

  cnt = 0;
  for(i=0;i<nf;i++) 
        if(!isDel[i]) {
           a = S.facets()[i][0];
           b = S.facets()[i][1];
           c = S.facets()[i][2];
           newS.facets()[cnt][0] = newIndex[a];
           newS.facets()[cnt][1] = newIndex[b];
           newS.facets()[cnt][2] = newIndex[c];
           cnt++;
           }

  S = newS;
  S.geometryChanged();
  return(ItXSuccess);
}


ItXECode SurfaceUtils::cleanSurface(Surface &S)
{
  int i,j,k,l,m,n,ai,bi,ci,aj,bj,cj;
  ItXECode r;
  int ndel,nv,nf,maxp,maxn;
  bool done;

  nf = S.getNumPoly();
  nv = S.getNumVert();
  if(nf<=0)
    return(ItXSuccess);

  Array2D<int> hash;

  // max # polys using a given vertex
  maxp = 15; 


retry_hash:
  cerr << "cleanSurface: Building Hash table (" << maxp << ")...";

  hash.setDim(nv,maxp+1);
  if(hash.isEmpty()) {
     cerr << endl << "ERROR: Out of memory in SurfaceUtils::cleanSurface()" << endl;
     return(ItXNoMemory);
  }
  hash = 0;
  maxn = 0;
  //
  // hash polys based on vertices
  //
  for(i=0;i<nf;i++) {
      ai = S.facets()[i][0];
      bi = S.facets()[i][1];
      ci = S.facets()[i][2];

      n  = hash[ai][0];
      if(n >= (maxp-1)) {
         maxp += 10;
         goto retry_hash;
      }
      n++;
      hash[ai][0] = n;
      hash[ai][n] = i;
      if(n>maxn) maxn=n;

      n  = hash[bi][0];
      if(n >= (maxp-1)) {
         maxp += 10;
         goto retry_hash;
      }
      n++;
      hash[bi][0] = n;
      hash[bi][n] = i;
      if(n>maxn) maxn=n;

      n  = hash[ci][0];
      if(n >= (maxp-1)) {
         maxp += 10;
         goto retry_hash;
      }
      n++;
      hash[ci][0] = n;
      hash[ci][n] = i;
      if(n>maxn) maxn=n;

  }
  cerr << "DONE" << endl;
  cerr << "CleanSurface: max of " << maxn << " polys/vertex" << endl;

  Array1D<u_char> delFac(nf);
  delFac = 0U;
  ndel   = 0;

  //
  // Compare each poly at a vertex to all others 
  // that use that vertex. Dups get tagged
  //
  for(i=0;i<nv;i++) {
    n = hash[i][0];
    for(j=1;j<=n;j++) {
        k = hash[i][j];
        if(delFac[k]) 
           continue;
        ai = S.facets()[k][0];
        bi = S.facets()[k][1];
        ci = S.facets()[k][2];
        for(l=j+1;l<=n;l++) {
           m = hash[i][l];
           if(delFac[m]) 
              continue;
           aj = S.facets()[m][0];
           bj = S.facets()[m][1];
           cj = S.facets()[m][2];
           if(((ai==cj)&&(bi==bj)&&(ci==aj))||
              ((ai==aj)&&(bi==cj)&&(ci==bj))||
              ((ai==bj)&&(bi==aj)&&(ci==cj))||
              ((ai==aj)&&(bi==bj)&&(ci==cj))||
              ((ai==bj)&&(bi==cj)&&(ci==aj))||
              ((ai==cj)&&(bi==aj)&&(ci==bj))) {
                delFac[k] = 1U;
                ndel++;
                break;
              }
        }
    }
  }

cerr << "cleanSurface:: Need to delete " << ndel << " facets" << endl;
  
  if(ndel > 0) {
    r = SurfaceUtils::removeSurfacePolys(S,delFac);
    if(r != ItXSuccess) {
       S.geometryChanged();
       return(r);
    }
  }

  nf = S.getNumPoly();
  nv = S.getNumVert();

  if((nf<=0)||(nv<=0))
     return(ItXSuccess);

  Array1D<u_char> delPt(nv);

  delPt  = 1U;
  for(i=0;i<nf;i++) {
     ai = S.facets()[i][0];
     bi = S.facets()[i][1];
     ci = S.facets()[i][2];
     delPt[ai] = 0U;
     delPt[bi] = 0U;
     delPt[ci] = 0U;
  }

  ndel=0;
  for(i=0;i<nv;i++)
     if(delPt[i])
        ndel++;

  cerr << "cleanSurface:: Need to delete " << ndel << " points" << endl;

  if(ndel > 0) {
     r = SurfaceUtils::removeSurfacePoints(S,delPt);
     if(r != ItXSuccess)
        return(r);
  }

  S.geometryChanged();

  return(ItXSuccess);
}

  // find the length of the path specified by the indeces in the array
ItXECode SurfaceUtils::findPathLength(
    Surface &S,
    Array1D<int> &path,
    double &length)
{
    length = 0;
    for (int i = 1; i < path.getLen(); i++)
        length += (S.getVert(path[i - 1]) - S.getVert(path[i])).norm();

    return ItXSuccess;
}




