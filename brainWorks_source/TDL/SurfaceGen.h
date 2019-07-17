#ifndef _SurfaceGen_h_
#define _SurfaceGen_h_

#include <OS.h>
#include <TDL/ItXWorkClass.h>

class Surface;
class ItXVolume;

#define TABLE_SIZE 9967

struct hashtype  {
  float *v1;
  float **v;
  struct hashtype *next;
};

typedef struct hashtype hash;

typedef struct {
  int nverts;
  int verts[8];
  int nedges;
  int edges[12];
  int npolys;
  int polys[30];
} CELL_ENTRY;


class SurfaceGen : public ItXWorkClass {
private:

	int VERBOSE;
	int XDIMYDIM;

	int avail;
    	hash *cur;
    	hash *h;


	float *VERTICES;        // a pointer to the x coordinates 
	int NUM_VERTICES;       // number of vertices currently stored 
	int VERT_LIMIT;         // currently allocated space for vertices 
	int VERT_INCR;          // allocate space for this many vertices at a
                                // time 

	hash *htable[TABLE_SIZE];
	hash *manage;
	
 	int iso_surface(u_char *data, int xdim, int ydim, int zdim, 
				float threshold);
 	int iso_surface(u_short *data, int xdim, int ydim, int zdim, 
				float threshold);
 	int iso_surface(short *data, int xdim, int ydim, int zdim, 
				float threshold);
 	int iso_surface(float *data, int xdim, int ydim, int zdim, 
				float threshold);

	void calc_index(int *index, u_char *data, int y1, int z1, 
				int xdim, float thresh);
	void calc_index(int *index, u_short *data, int y1, int z1, 
				int xdim, float thresh);
	void calc_index(int *index, short *data, int y1, int z1, 
				int xdim, float thresh);
	void calc_index(int *index, float *data, int y1, int z1, 
				int xdim, float thresh);

	void get_cell_verts(int *index, u_char *data, int y1, int z1, 
				int xdim, float xtrans, float ytrans, 
				float ztrans, float threshold, 
				float *crossings);
	void get_cell_verts(int *index, u_short *data, int y1, int z1, 
				int xdim, float xtrans, float ytrans, 
				float ztrans, float threshold, 
				float *crossings);
	void get_cell_verts(int *index, short *data, int y1, int z1, 
				int xdim, float xtrans, float ytrans, 
				float ztrans, float threshold, 
				float *crossings);
	void get_cell_verts(int *index, float *data, int y1, int z1, 
				int xdim, float xtrans, float ytrans, 
				float ztrans, float threshold, 
				float *crossings);

	int get_cell_polys(int *index, int xdim, float *crossings);

	int hashit(float *v);
	
	float ***tricompact(float *v,int nv);
		
	hash *newhash();

	void freehash();
	
	int add_polygon(float *p1,float *p2,float *p3);
	int vert_alloc(int size);
	float xtrans,ytrans,ztrans;

	static CELL_ENTRY cell_table[];

public:
	SurfaceGen();
	~SurfaceGen();

	ItXECode apply(ItXVolume *dat,Surface *S, int nx, int ny, int nz,
		          float thresh,float xt = 0,float yt=0, float zt= 0);
	};

#endif
