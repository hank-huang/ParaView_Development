#ifndef _SurfaceGrow_h_
#define _SurfaceGrow_h_

#include <OS.h>

class Surface;
class ItXVolume;

typedef struct _point {
  float x;
  float y;
  float z;
  } point;

typedef struct _triangle {
  int a;
  int b;
  int c;
  } triangle;


typedef struct _nbhd {
  int num;
  int num_alloced;
  int *neigh;
  } neighborhood;


class SurfaceGrow : public ItXWorkClass {
public:
	SurfaceGrow();
       ~SurfaceGrow();
	
	bool initialize(Surface*,const ItXVolume*,int hfact = 20);
	bool apply(  int numiter,
		float bg_val,
		float fg_val,
		float alpha,
		float beta,
		float delta,
		float ath,
		float max_newpoints
		);
	void cleanup();
		
protected:

	float *flag;
	float maxDeform;
	point    *vert;
	point    *newvert;
	point    *normals;
	triangle *face;	
	float    *deform;
	int     **cubetable;
	int      *numtricube;
	int      *tri_compared;
	int       ncubetable;

	neighborhood *vert_nbhd;
	neighborhood *face_nbhd;

	const ItXVolume *img;
	Surface   *surf;

	char work_string[512];

	// # of vertices & normals
	int   numpoints;
	// size of vertices & normals array
	int   numvert_alloced;

	// # of facets
	int   numfaces;

	// size of facets array
	int   numface_alloced;

	// size of search neighborhood in voxels
	float hsizex;
	float hsizey;
	float hsizez;
	float areath;
	int   hfact;
	int   redist_iter;
	int   maxpoly_percube;
	int   max_addpoints;
	int   total_addpoints;

	bool importSurface();
	bool exportSurface();

	void cross(float u1[3], float u2[3], float u[3]);
	void norm(float u[3]);
	void computenormals();
	ItXECode refinetriang();
	int  addPoint(float,float,float);
	int  addFacet(int,int,int);
	void addVertNeighbor(int,int);
	void addFacetNeighbor(int,int);
	void resizeNeighborArray(neighborhood**,int,int);
	void buildVertNeighborhood();

	int  hashindex(int x, int y, int z);
	int  hashnumbers(int,int*,int*,int*);
	int  hashcube(float,float,float);
	void findCubeRegion(
		int inda, int indb, int indc,
                int *outxmin, int *outxmax,
                int *outymin, int *outymax,
                int *outzmin, int *outzmax);
	bool createCubeTable();
	bool calcCubeTable();
	int  checkinter();

	int tri_tri_intersect(
	        float V0[3],float V1[3],float V2[3],
                float U0[3],float U1[3],float U2[3]);

	int coplanar_tri_tri(
                float N[3],float V0[3],float V1[3],float V2[3],
               float U0[3],float U1[3],float U2[3]);

	};

#endif
