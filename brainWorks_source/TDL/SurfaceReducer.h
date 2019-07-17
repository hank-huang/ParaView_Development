#ifndef _SurfaceReducer_h_
#define _SurfaceReducer_h_

#include <OS.h>
#include <TDL/IDLdefines.h>
#include <TDL/ItXWorkClass.h>

class Surface;

#define MIN_ASPECT      	0.05
#define MAX_NEIGHBORS   	50  // max 1 deep neighbors
#define MAX_CURV_NEIGHBORS   	1000 // max 2 deep neighbors
#define TLRNCE          	0.0001
#define MAX_CURV_RES    	20.0
#define MAX_NEW_POLY_POINTS    100


typedef struct _spoint {
        float x,y,z;
        float nx,ny,nz;
        float curvature;
        int num,used,can_move,isedge;
        struct _spoint *neighbors[MAX_NEIGHBORS];
        struct _spoint *next;
        struct _spoint *prev;
        } spoint;

typedef struct _spoly {
        spoint *a,*b,*c;
        struct _spoly *prev;
        struct _spoly *next;
        } spoly;

typedef struct _kpoint {
        float x,y,z;
        } kpoint;


class SurfaceReducer : public ItXWorkClass {
public:
	SurfaceReducer()	{}
       ~SurfaceReducer()	{}

	ItXECode decimate(Surface*,float min_curvature,
                float max_curvature,
		int maxiter, float redist_alpha,
		int redist_iter, int &nremoved);
protected:
	spoint  *spoints,*orig_spoints;
	spoly   *spolys,*orig_spolys;
	spoint  *new_poly[MAX_NEW_POLY_POINTS];
	int	 nnew;
	kpoint  *newvert;
	int      redist_iter;
	float    redist_alpha;
	int	 cur_isccw;
	spoly  **spoly_heap;
	int      nheap;
	int      maxheap;
	float    min_curvature;
	float    max_curvature;

	int  reducePolysByCurvature();
	void ComputeNormals(int);
	void TagEdges();
	void Redistribute();
	int  CreateNewPolyFromPointRemoval(spoint *p);
	void CreateNewPolys();
	int  PolyIsConvex();
	void RemovePoint(spoint*);
	void RemovePoly(spoly*);
	void RemovePointAndPolys(spoint*);
	int  FindOrient(spoint*,spoint*);
	int  generateLinkedArrays(Surface *s);
	int  setDataFromLinkedArrays(Surface *s);
	static void crossvec(float *u1, float *u2, float *u);
	static void crossproduct(spoint *a, spoint *b, spoint *c, float *nor);
	static void normvec(float *u);
	static void AddNeighbor(spoint *p, spoint *n);
	static void AddNeighbor2(spoint *p, spoint *n, spoint **arr);
	static void RemoveNeighbor(spoint *p, spoint *n);
	static int  IsInNeighborhood(spoint *p, spoint *n);
	static float gaussCurvature(spoint *p);
	static int  BadAspectRatio(spoint *a,
			spoint *b, spoint *c);

	};

#endif
