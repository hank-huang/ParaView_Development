///////////////////////////////////////////////////////////////////////////
// File: SurfaceGen.C
//
// Author: Keith Doolittle 10/97
//
// Purpose: Isosurface generation
// 
// Isosurface extraction algorithm from:
//
//        "Exploiting Triangulated Surface Extraction Using
//         Tetrahedral Decomposition", Andre Gueziec, Robert Hummel,
//        IEEE Transactions on Visualization and Computer Graphics,
//        Vol. 1, No. 4, December 1995.
//
///////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <TDL/Surface.h>
#include <TDL/ItXVolume.h>
#include <TDL/ItXVolumeUtils.h>
#include <TDL/SurfaceGen.h>
#include <string.h>
#include <math.h>
#include <TDL/cell_table.h>

ItXECode SurfaceGen::apply(ItXVolume *vol, Surface *S, 
			   int nx, int ny, int nz, 
			   float thresh,
			   float xt,float yt,float zt)
{
        if(!vol)  return(ItXError);

	int i,j;

        manage = new hash;
        if(!manage)
           return(ItXNoMemory);

	VERTICES = NULL;         // a pointer to the x coordinates
 	NUM_VERTICES = 0;        // number of vertices currently stored
	VERT_LIMIT = 0;          // currently allocated space for vertices
	VERT_INCR = 20000;       // allocate space for this many vertices at a
	VERBOSE = 1;

	avail = 0;
	manage->v1   = NULL;
	manage->v1   = NULL;
        manage->v    = NULL;
        manage->next = NULL;


        for (i=0;i<TABLE_SIZE;i++){
           htable[i] = new hash;
           if(!htable[i]) {
              for(j=0;j<i;j++)
                 delete htable[i];       
              return(ItXNoMemory);
           }
           htable[i]->v1   = NULL;
           htable[i]->v    = NULL;
           htable[i]->next = NULL;
        }

	xtrans = xt;
	ytrans = yt;
	ztrans = zt;

	cout << "Generating Iso Surface" << endl;
	cout << "(nx,ny,nz) = " << nx << "," << ny << "," << nz << endl; 
	cout << "Threshold = " << thresh;

        int rv;
        switch(vol->dataType()) {
        case ItXVolume::UnsignedChar:
	   rv = iso_surface(vol->u_char_data().data(), nx, ny, nz, thresh); 
           break;
        case ItXVolume::UnsignedShort:
	   rv = iso_surface(vol->u_short_data().data(), nx, ny, nz, thresh); 
           break;
        case ItXVolume::Short:
	   rv = iso_surface(vol->short_data().data(), nx, ny, nz, thresh); 
           break;
        case ItXVolume::Float:
	   rv = iso_surface(vol->float_data().data(), nx, ny, nz, thresh); 
           break;
        }

        if(!rv) return(ItXNoMemory);

	cout << "Done Generating iso_surface" << endl;
	cout << "NUM_VERTICES =  " << NUM_VERTICES << endl;

	float ***conn;
	float **vlist;
	int a,b,c;
	int nvert;
        int npol;

	cout << "Coampacting Surface " << endl;

	conn = tricompact(VERTICES, NUM_VERTICES);

        if(!conn) {
           delete [] VERTICES;
           return(ItXNoMemory);
        }

	cout << "Done compacting Surface" << endl;

	nvert = conn[NUM_VERTICES] - conn[0];
	vlist = conn[0];

        npol = NUM_VERTICES/3;

        S->setSize(nvert,npol);

	// Set the vertices

	for(i=0;i<nvert;i++) {
	   float *v = vlist[i];
	   S->vertices()[i].set(v[0],v[1],v[2]);
	}

	// Set the polys

	for(j=0;j<npol;j++) {
	   i = j*3;
	   a = conn[i]-conn[0];
	   b = conn[i+1]-conn[0];
	   c = conn[i+2]-conn[0];
	   S->facets()[j][0] = a;
	   S->facets()[j][1] = b;
	   S->facets()[j][2] = c;
	}

        S->geometryChanged();

	cout << "Numvert = " << S->getNumVert();
	cout << "NumPoly = " << S->getNumPoly();

	avail = 0;

        delete [] conn;
        delete [] vlist;
        delete [] VERTICES;
	
	return(ItXSuccess);
}



SurfaceGen::SurfaceGen()
{
}


SurfaceGen::~SurfaceGen()
{
}


int SurfaceGen::iso_surface(u_char *data, 
			int xdim, int ydim, int zdim, 
			float threshold)
{
    register int x, y, z;
    register int xdim1, ydim1, zdim1;
    int *index;
    int npolys,np;
    float *crossings;
    float *normals;

    zdim1 = zdim - 1;
    ydim1 = ydim - 1;
    xdim1 = xdim - 1;

    XDIMYDIM = xdim * ydim;

    npolys = 0;

    index     = new int[xdim1];
    crossings = new float[xdim1*13*3];
    if(!index || !crossings)
       return(0);

    for(z = 0; z < zdim1; z++)
      for (y = 0; y < ydim1; y++) {
	calc_index(index, data, y, z, xdim, threshold);
	get_cell_verts(index, data, y, z, xdim, xtrans, ytrans, ztrans,
			       threshold, crossings);
	if((np = get_cell_polys(index, xdim, crossings)) < 0) {
            delete [] index;
            delete [] crossings;
            return(0);
        }
        
	npolys += np;
      }

    delete [] index;
    delete [] crossings;

    return(1);
}



int SurfaceGen::iso_surface(u_short *data, 
			int xdim, int ydim, int zdim, 
			float threshold)
{
    register int x, y, z;
    register int xdim1, ydim1, zdim1;
    int *index;
    int npolys,np;
    float *crossings;
    float *normals;

    zdim1 = zdim - 1;
    ydim1 = ydim - 1;
    xdim1 = xdim - 1;

    XDIMYDIM = xdim * ydim;

    npolys = 0;

    index     = new int[xdim1];
    crossings = new float[xdim1*13*3];
    if(!index || !crossings)
       return(0);

    for(z = 0; z < zdim1; z++)
      for (y = 0; y < ydim1; y++) {
	calc_index(index, data, y, z, xdim, threshold);
	get_cell_verts(index, data, y, z, xdim, xtrans, ytrans, ztrans,
			       threshold, crossings);
	if((np = get_cell_polys(index, xdim, crossings)) < 0) {
            delete [] index;
            delete [] crossings;
            return(0);
        }
        
	npolys += np;
      }

    delete [] index;
    delete [] crossings;

    return(1);
}


int SurfaceGen::iso_surface(short *data, 
			int xdim, int ydim, int zdim, 
			float threshold)
{
    register int x, y, z;
    register int xdim1, ydim1, zdim1;
    int *index;
    int npolys,np;
    float *crossings;
    float *normals;

    zdim1 = zdim - 1;
    ydim1 = ydim - 1;
    xdim1 = xdim - 1;

    XDIMYDIM = xdim * ydim;

    npolys = 0;

    index     = new int[xdim1];
    crossings = new float[xdim1*13*3];
    if(!index || !crossings)
       return(0);

    for(z = 0; z < zdim1; z++)
      for (y = 0; y < ydim1; y++) {
	calc_index(index, data, y, z, xdim, threshold);
	get_cell_verts(index, data, y, z, xdim, xtrans, ytrans, ztrans,
			       threshold, crossings);
	if((np = get_cell_polys(index, xdim, crossings)) < 0) {
            delete [] index;
            delete [] crossings;
            return(0);
        }
        
	npolys += np;
      }

    delete [] index;
    delete [] crossings;

    return(1);
}


int SurfaceGen::iso_surface(float *data, 
			int xdim, int ydim, int zdim, 
			float threshold)
{
    register int x, y, z;
    register int xdim1, ydim1, zdim1;
    int *index;
    int npolys,np;
    float *crossings;
    float *normals;

    zdim1 = zdim - 1;
    ydim1 = ydim - 1;
    xdim1 = xdim - 1;

    XDIMYDIM = xdim * ydim;

    npolys = 0;

    index     = new int[xdim1];
    crossings = new float[xdim1*13*3];
    if(!index || !crossings)
       return(0);

    for(z = 0; z < zdim1; z++)
      for (y = 0; y < ydim1; y++) {
	calc_index(index, data, y, z, xdim, threshold);
	get_cell_verts(index, data, y, z, xdim, xtrans, ytrans, ztrans,
			       threshold, crossings);
	if((np = get_cell_polys(index, xdim, crossings)) < 0) {
            delete [] index;
            delete [] crossings;
            return(0);
        }
        
	npolys += np;
      }

    delete [] index;
    delete [] crossings;

    return(1);
}




void SurfaceGen:: calc_index(int *index, u_char *data, 
			     int y1, int z1, int xdim, float thresh)
{
    register u_char *tmp;
    register float threshold = thresh;
    int x1;
    unsigned int i = 0;

    // first compute index of first cube 

    tmp = data + (z1 * XDIMYDIM) + (y1 * xdim) + 0;

    i += (threshold <= tmp[0]);
    i += (threshold <= tmp[1]) * 2;

    tmp += xdim;
    i += (threshold <= tmp[1]) * 4;
    i += (threshold <= tmp[0]) * 8;

    tmp = tmp - xdim + XDIMYDIM;
    i += (threshold <= tmp[0]) * 16;
    i += (threshold <= tmp[1]) * 32;

    tmp += xdim;
    i += (threshold <= tmp[1]) * 64;
    i += (threshold <= tmp[0]) * 128;

    index[0] = i;

    // now compute rest 

    tmp -= xdim + XDIMYDIM;
    for (x1 = 1; x1 < xdim-1; x1++) {
	++tmp;

	// resuse 4 of the bits 
	i = ((i&0x44)<<1) | ((i&0x22)>>1);

	i += (threshold <= tmp[1]) * 2;
	i += (threshold <= tmp[xdim+1]) * 4;
	i += (threshold <= tmp[XDIMYDIM+1]) * 32;
	i += (threshold <= tmp[XDIMYDIM+xdim+1]) * 64;

	index[x1] = i;
    }
}


void SurfaceGen:: calc_index(int *index, u_short *data, 
			     int y1, int z1, int xdim, float thresh)
{
    register u_short *tmp;
    register float threshold = thresh;
    int x1;
    unsigned int i = 0;

    // first compute index of first cube 

    tmp = data + (z1 * XDIMYDIM) + (y1 * xdim) + 0;

    i += (threshold <= tmp[0]);
    i += (threshold <= tmp[1]) * 2;

    tmp += xdim;
    i += (threshold <= tmp[1]) * 4;
    i += (threshold <= tmp[0]) * 8;

    tmp = tmp - xdim + XDIMYDIM;
    i += (threshold <= tmp[0]) * 16;
    i += (threshold <= tmp[1]) * 32;

    tmp += xdim;
    i += (threshold <= tmp[1]) * 64;
    i += (threshold <= tmp[0]) * 128;

    index[0] = i;

    // now compute rest 

    tmp -= xdim + XDIMYDIM;
    for (x1 = 1; x1 < xdim-1; x1++) {
	++tmp;

	// resuse 4 of the bits 
	i = ((i&0x44)<<1) | ((i&0x22)>>1);

	i += (threshold <= tmp[1]) * 2;
	i += (threshold <= tmp[xdim+1]) * 4;
	i += (threshold <= tmp[XDIMYDIM+1]) * 32;
	i += (threshold <= tmp[XDIMYDIM+xdim+1]) * 64;

	index[x1] = i;
    }
}




void SurfaceGen:: calc_index(int *index, short *data, 
			     int y1, int z1, int xdim, float thresh)
{
    register short *tmp;
    register float threshold = thresh;
    int x1;
    unsigned int i = 0;

    // first compute index of first cube 

    tmp = data + (z1 * XDIMYDIM) + (y1 * xdim) + 0;

    i += (threshold <= tmp[0]);
    i += (threshold <= tmp[1]) * 2;

    tmp += xdim;
    i += (threshold <= tmp[1]) * 4;
    i += (threshold <= tmp[0]) * 8;

    tmp = tmp - xdim + XDIMYDIM;
    i += (threshold <= tmp[0]) * 16;
    i += (threshold <= tmp[1]) * 32;

    tmp += xdim;
    i += (threshold <= tmp[1]) * 64;
    i += (threshold <= tmp[0]) * 128;

    index[0] = i;

    // now compute rest 

    tmp -= xdim + XDIMYDIM;
    for (x1 = 1; x1 < xdim-1; x1++) {
	++tmp;

	// resuse 4 of the bits 
	i = ((i&0x44)<<1) | ((i&0x22)>>1);

	i += (threshold <= tmp[1]) * 2;
	i += (threshold <= tmp[xdim+1]) * 4;
	i += (threshold <= tmp[XDIMYDIM+1]) * 32;
	i += (threshold <= tmp[XDIMYDIM+xdim+1]) * 64;

	index[x1] = i;
    }
}


void SurfaceGen:: calc_index(int *index, float *data, 
			     int y1, int z1, int xdim, float thresh)
{
    register float *tmp;
    register float threshold = thresh;
    int x1;
    unsigned int i = 0;

    // first compute index of first cube 

    tmp = data + (z1 * XDIMYDIM) + (y1 * xdim) + 0;

    i += (threshold <= tmp[0]);
    i += (threshold <= tmp[1]) * 2;

    tmp += xdim;
    i += (threshold <= tmp[1]) * 4;
    i += (threshold <= tmp[0]) * 8;

    tmp = tmp - xdim + XDIMYDIM;
    i += (threshold <= tmp[0]) * 16;
    i += (threshold <= tmp[1]) * 32;

    tmp += xdim;
    i += (threshold <= tmp[1]) * 64;
    i += (threshold <= tmp[0]) * 128;

    index[0] = i;

    // now compute rest 

    tmp -= xdim + XDIMYDIM;
    for (x1 = 1; x1 < xdim-1; x1++) {
	++tmp;

	// resuse 4 of the bits 
	i = ((i&0x44)<<1) | ((i&0x22)>>1);

	i += (threshold <= tmp[1]) * 2;
	i += (threshold <= tmp[xdim+1]) * 4;
	i += (threshold <= tmp[XDIMYDIM+1]) * 32;
	i += (threshold <= tmp[XDIMYDIM+xdim+1]) * 64;

	index[x1] = i;
    }
}




#define CROSSINGS(x,a,b) crossings[x*13*3+a*3+b]
#define linterp(a1,a2,a,b1,b2) \
	((float)(((a-a1) * (float)(b2-b1) / (a2-a1)) + (float)b1))








void SurfaceGen:: get_cell_verts(int *index, u_char *data, 
			int y1, int z1, int xdim, 
			float xtrans, float ytrans, float ztrans, 
			float threshold, float *crossings)
{
    int x1, y2, z2;

    y2 = y1 + 1;
    z2 = z1 + 1;
    for (x1 = 0; x1 < xdim-1; x1++) {
	float cx, cy, cz;
	int nedges;
	int crnt_edge;
	int x2 = x1 + 1;
	int i;
	u_char *v1, *v4, *v5, *v8;


	if (!index[x1]) continue;

	v1 = data + z1*XDIMYDIM + y1*xdim + x1;
	v4 = v1 + xdim;
	v5 = v1 + XDIMYDIM;
	v8 = v4 + XDIMYDIM;

	nedges = cell_table[index[x1]].nedges;
	for (i = 0; i < nedges; i++) {
	    crnt_edge = cell_table[index[x1]].edges[i];
	    cx = xtrans; cy = ytrans; cz = ztrans;
	    switch (crnt_edge) {
	    case 1:
	    cx += linterp(v1[0], v1[1], threshold, x1, x2);
	    cy += (float) y1;
	    cz += (float) z1;
	    break;

	    case 2:
	    cy += linterp(v1[1], v4[1], threshold, y1, y2);
	    cx += (float) x2;
	    cz += (float) z1;
	    break;

	    case 3:
	    cx += linterp(v4[0], v4[1], threshold, x1, x2);
	    cy += (float) y2;
	    cz += (float) z1;
	    break;

	    case 4:
	    cy += linterp(v1[0], v4[0], threshold, y1, y2);
	    cx += (float) x1;
	    cz += (float) z1;
	    break;

	    case 5:
	    cx += linterp(v5[0], v5[1], threshold, x1, x2);
	    cy += (float) y1;
	    cz += (float) z2;
	    break;

	    case 6:
	    cy += linterp(v5[1], v8[1], threshold, y1, y2);
	    cx += (float) x2;
	    cz += (float) z2;
	    break;

	    case 7:
	    cx += linterp(v8[0], v8[1], threshold, x1, x2);
	    cy += (float) y2;
	    cz += (float) z2;
	    break;

	    case 8:
	    cy += linterp(v5[0], v8[0], threshold, y1, y2);
	    cx += (float) x1;
	    cz += (float) z2;
	    break;

	    case 9:
	    cz += linterp(v1[0], v5[0], threshold, z1, z2);
	    cy += (float) y1;
	    cx += (float) x1;
	    break;

	    case 10:
	    cz += linterp(v1[1], v5[1], threshold, z1, z2);
	    cy += (float) y1;
	    cx += (float) x2;
	    break;

	    case 11:
	    cz += linterp(v4[0], v8[0], threshold, z1, z2);
	    cy += (float) y2;
	    cx += (float) x1;
	    break;

	    case 12:
	    cz += linterp(v4[1], v8[1], threshold, z1, z2);
	    cy += (float) y2;
	    cx += (float) x2;
	    break;

	    } 
	    CROSSINGS(x1,crnt_edge,0) = cx;
	    CROSSINGS(x1,crnt_edge,1) = cy;
	    CROSSINGS(x1,crnt_edge,2) = cz;
//	    CROSSINGS(x1,crnt_edge,0) = cx+0.5;
//	    CROSSINGS(x1,crnt_edge,1) = cy+0.5;
//	    CROSSINGS(x1,crnt_edge,2) = cz+0.5;
	} 
    }
}


void SurfaceGen:: get_cell_verts(int *index, u_short *data, 
			int y1, int z1, int xdim, 
			float xtrans, float ytrans, float ztrans, 
			float threshold, float *crossings)
{
    int x1, y2, z2;

    y2 = y1 + 1;
    z2 = z1 + 1;
    for (x1 = 0; x1 < xdim-1; x1++) {
	float cx, cy, cz;
	int nedges;
	int crnt_edge;
	int x2 = x1 + 1;
	int i;
	u_short *v1, *v4, *v5, *v8;


	if (!index[x1]) continue;

	v1 = data + z1*XDIMYDIM + y1*xdim + x1;
	v4 = v1 + xdim;
	v5 = v1 + XDIMYDIM;
	v8 = v4 + XDIMYDIM;

	nedges = cell_table[index[x1]].nedges;
	for (i = 0; i < nedges; i++) {
	    crnt_edge = cell_table[index[x1]].edges[i];
	    cx = xtrans; cy = ytrans; cz = ztrans;
	    switch (crnt_edge) {
	    case 1:
	    cx += linterp(v1[0], v1[1], threshold, x1, x2);
	    cy += (float) y1;
	    cz += (float) z1;
	    break;

	    case 2:
	    cy += linterp(v1[1], v4[1], threshold, y1, y2);
	    cx += (float) x2;
	    cz += (float) z1;
	    break;

	    case 3:
	    cx += linterp(v4[0], v4[1], threshold, x1, x2);
	    cy += (float) y2;
	    cz += (float) z1;
	    break;

	    case 4:
	    cy += linterp(v1[0], v4[0], threshold, y1, y2);
	    cx += (float) x1;
	    cz += (float) z1;
	    break;

	    case 5:
	    cx += linterp(v5[0], v5[1], threshold, x1, x2);
	    cy += (float) y1;
	    cz += (float) z2;
	    break;

	    case 6:
	    cy += linterp(v5[1], v8[1], threshold, y1, y2);
	    cx += (float) x2;
	    cz += (float) z2;
	    break;

	    case 7:
	    cx += linterp(v8[0], v8[1], threshold, x1, x2);
	    cy += (float) y2;
	    cz += (float) z2;
	    break;

	    case 8:
	    cy += linterp(v5[0], v8[0], threshold, y1, y2);
	    cx += (float) x1;
	    cz += (float) z2;
	    break;

	    case 9:
	    cz += linterp(v1[0], v5[0], threshold, z1, z2);
	    cy += (float) y1;
	    cx += (float) x1;
	    break;

	    case 10:
	    cz += linterp(v1[1], v5[1], threshold, z1, z2);
	    cy += (float) y1;
	    cx += (float) x2;
	    break;

	    case 11:
	    cz += linterp(v4[0], v8[0], threshold, z1, z2);
	    cy += (float) y2;
	    cx += (float) x1;
	    break;

	    case 12:
	    cz += linterp(v4[1], v8[1], threshold, z1, z2);
	    cy += (float) y2;
	    cx += (float) x2;
	    break;

	    } 
	    CROSSINGS(x1,crnt_edge,0) = cx;
	    CROSSINGS(x1,crnt_edge,1) = cy;
	    CROSSINGS(x1,crnt_edge,2) = cz;
//	    CROSSINGS(x1,crnt_edge,0) = cx + 0.5;
//	    CROSSINGS(x1,crnt_edge,1) = cy + 0.5;
//	    CROSSINGS(x1,crnt_edge,2) = cz + 0.5;
	} 
    }
}


void SurfaceGen:: get_cell_verts(int *index, short *data, 
			int y1, int z1, int xdim, 
			float xtrans, float ytrans, float ztrans, 
			float threshold, float *crossings)
{
    int x1, y2, z2;

    y2 = y1 + 1;
    z2 = z1 + 1;
    for (x1 = 0; x1 < xdim-1; x1++) {
	float cx, cy, cz;
	int nedges;
	int crnt_edge;
	int x2 = x1 + 1;
	int i;
	short *v1, *v4, *v5, *v8;


	if (!index[x1]) continue;

	v1 = data + z1*XDIMYDIM + y1*xdim + x1;
	v4 = v1 + xdim;
	v5 = v1 + XDIMYDIM;
	v8 = v4 + XDIMYDIM;

	nedges = cell_table[index[x1]].nedges;
	for (i = 0; i < nedges; i++) {
	    crnt_edge = cell_table[index[x1]].edges[i];
	    cx = xtrans; cy = ytrans; cz = ztrans;
	    switch (crnt_edge) {
	    case 1:
	    cx += linterp(v1[0], v1[1], threshold, x1, x2);
	    cy += (float) y1;
	    cz += (float) z1;
	    break;

	    case 2:
	    cy += linterp(v1[1], v4[1], threshold, y1, y2);
	    cx += (float) x2;
	    cz += (float) z1;
	    break;

	    case 3:
	    cx += linterp(v4[0], v4[1], threshold, x1, x2);
	    cy += (float) y2;
	    cz += (float) z1;
	    break;

	    case 4:
	    cy += linterp(v1[0], v4[0], threshold, y1, y2);
	    cx += (float) x1;
	    cz += (float) z1;
	    break;

	    case 5:
	    cx += linterp(v5[0], v5[1], threshold, x1, x2);
	    cy += (float) y1;
	    cz += (float) z2;
	    break;

	    case 6:
	    cy += linterp(v5[1], v8[1], threshold, y1, y2);
	    cx += (float) x2;
	    cz += (float) z2;
	    break;

	    case 7:
	    cx += linterp(v8[0], v8[1], threshold, x1, x2);
	    cy += (float) y2;
	    cz += (float) z2;
	    break;

	    case 8:
	    cy += linterp(v5[0], v8[0], threshold, y1, y2);
	    cx += (float) x1;
	    cz += (float) z2;
	    break;

	    case 9:
	    cz += linterp(v1[0], v5[0], threshold, z1, z2);
	    cy += (float) y1;
	    cx += (float) x1;
	    break;

	    case 10:
	    cz += linterp(v1[1], v5[1], threshold, z1, z2);
	    cy += (float) y1;
	    cx += (float) x2;
	    break;

	    case 11:
	    cz += linterp(v4[0], v8[0], threshold, z1, z2);
	    cy += (float) y2;
	    cx += (float) x1;
	    break;

	    case 12:
	    cz += linterp(v4[1], v8[1], threshold, z1, z2);
	    cy += (float) y2;
	    cx += (float) x2;
	    break;

	    } 
	    CROSSINGS(x1,crnt_edge,0) = cx;
	    CROSSINGS(x1,crnt_edge,1) = cy;
	    CROSSINGS(x1,crnt_edge,2) = cz;
//	    CROSSINGS(x1,crnt_edge,0) = cx+0.5;
//	    CROSSINGS(x1,crnt_edge,1) = cy+0.5;
//	    CROSSINGS(x1,crnt_edge,2) = cz+0.5;
	} 
    }
}

void SurfaceGen:: get_cell_verts(int *index, float *data, 
			int y1, int z1, int xdim, 
			float xtrans, float ytrans, float ztrans, 
			float threshold, float *crossings)
{
    int x1, y2, z2;

    y2 = y1 + 1;
    z2 = z1 + 1;
    for (x1 = 0; x1 < xdim-1; x1++) {
	float cx, cy, cz;
	int nedges;
	int crnt_edge;
	int x2 = x1 + 1;
	int i;
	float *v1, *v4, *v5, *v8;


	if (!index[x1]) continue;

	v1 = data + z1*XDIMYDIM + y1*xdim + x1;
	v4 = v1 + xdim;
	v5 = v1 + XDIMYDIM;
	v8 = v4 + XDIMYDIM;

	nedges = cell_table[index[x1]].nedges;
	for (i = 0; i < nedges; i++) {
	    crnt_edge = cell_table[index[x1]].edges[i];
	    cx = xtrans; cy = ytrans; cz = ztrans;
	    switch (crnt_edge) {
	    case 1:
	    cx += linterp(v1[0], v1[1], threshold, x1, x2);
	    cy += (float) y1;
	    cz += (float) z1;
	    break;

	    case 2:
	    cy += linterp(v1[1], v4[1], threshold, y1, y2);
	    cx += (float) x2;
	    cz += (float) z1;
	    break;

	    case 3:
	    cx += linterp(v4[0], v4[1], threshold, x1, x2);
	    cy += (float) y2;
	    cz += (float) z1;
	    break;

	    case 4:
	    cy += linterp(v1[0], v4[0], threshold, y1, y2);
	    cx += (float) x1;
	    cz += (float) z1;
	    break;

	    case 5:
	    cx += linterp(v5[0], v5[1], threshold, x1, x2);
	    cy += (float) y1;
	    cz += (float) z2;
	    break;

	    case 6:
	    cy += linterp(v5[1], v8[1], threshold, y1, y2);
	    cx += (float) x2;
	    cz += (float) z2;
	    break;

	    case 7:
	    cx += linterp(v8[0], v8[1], threshold, x1, x2);
	    cy += (float) y2;
	    cz += (float) z2;
	    break;

	    case 8:
	    cy += linterp(v5[0], v8[0], threshold, y1, y2);
	    cx += (float) x1;
	    cz += (float) z2;
	    break;

	    case 9:
	    cz += linterp(v1[0], v5[0], threshold, z1, z2);
	    cy += (float) y1;
	    cx += (float) x1;
	    break;

	    case 10:
	    cz += linterp(v1[1], v5[1], threshold, z1, z2);
	    cy += (float) y1;
	    cx += (float) x2;
	    break;

	    case 11:
	    cz += linterp(v4[0], v8[0], threshold, z1, z2);
	    cy += (float) y2;
	    cx += (float) x1;
	    break;

	    case 12:
	    cz += linterp(v4[1], v8[1], threshold, z1, z2);
	    cy += (float) y2;
	    cx += (float) x2;
	    break;

	    } 
	    CROSSINGS(x1,crnt_edge,0) = cx;
	    CROSSINGS(x1,crnt_edge,1) = cy;
	    CROSSINGS(x1,crnt_edge,2) = cz;
//	    CROSSINGS(x1,crnt_edge,0) = cx+0.5;
//	    CROSSINGS(x1,crnt_edge,1) = cy+0.5;
//	    CROSSINGS(x1,crnt_edge,2) = cz+0.5;
	} 
    }
}


//
// This subroutine will calculate the polygons 
//
int SurfaceGen::get_cell_polys(int *index, int xdim, float *crossings)
{

    register int num_o_polys, polys = 0;
    register int poly;
    float *p1, *p2, *p3;
    int x1;

    for (x1 = 0; x1 < xdim-1; x1++) {
	if(!index[x1]) continue;
	num_o_polys = cell_table[index[x1]].npolys;
	for (poly = 0; poly < num_o_polys; poly++) {

	    p1 = &CROSSINGS(x1,cell_table[index[x1]].polys[poly*3],0);
	    p2 = &CROSSINGS(x1,cell_table[index[x1]].polys[poly*3 + 1],0);
	    p3 = &CROSSINGS(x1,cell_table[index[x1]].polys[poly*3 + 2],0);
            //
            // Check for degenerate poly
            //
	    if ((p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2]) ||
		(p1[0] == p3[0] && p1[1] == p3[1] && p1[2] == p3[2]) ||
		(p2[0] == p3[0] && p2[1] == p3[1] && p2[2] == p3[2]))  {
		polys--;
		continue;
	    }

	    if(!add_polygon(p1, p2, p3)) 
               return(-1);
	}

	polys += num_o_polys;
    }
    return polys;
}


int SurfaceGen::hashit(float *v)
{
   unsigned long *p = (unsigned long *)v;
   return (((p[0]*283+p[1])*283)+p[2]) % TABLE_SIZE;
}


float ***SurfaceGen::tricompact(float *v,int nv)
{
    int i;

    float ***conn = new float**[nv+1];

    // TODO: vert is NEVER freed

    float  **vert = new float*[nv];

    int nvert = 0;
    int dup = 0;

    if (!conn || !vert) {
        if(conn) delete [] conn;
        if(vert) delete [] vert;
        return(NULL); 
    }

    for(i = 0; i < nv; i++) {
	int h = hashit(v+3*i);
	hash *hh;

	// if vertex already exists - ignore 
	for(hh = htable[h]; hh; hh=hh->next) {
	if (hh->v1) {
	    if (hh->v1[0] == v[3*i+0] &&
		hh->v1[1] == v[3*i+1] &&
		hh->v1[2] == v[3*i+2]) {
		    goto vdup;
	    }
	 }	
	}
	vert[nvert] = v+3*i;
	hh = newhash();
	hh->next = htable[h]; htable[h] = hh;
	hh->v1 = v+3*i; hh->v = vert+nvert;
	conn[i] = vert+nvert;
	nvert++;
	continue;
vdup:
	conn[i] = hh->v;
	dup++;
    }
    conn[i] = vert+nvert;
    if(VERBOSE) printf ("verts = %d dups = %d\n", nv, dup);

    { int vp = 0, np = 0;
    for(i = 0; i < TABLE_SIZE; i++) {
	hash *p;
	if (!htable[i]) continue;
	for(p = htable[i]; p; p = p->next) {
	    vp++;
	}
	np++;
    }
    }
    freehash();
    return(conn);
}



hash *SurfaceGen::newhash() 
{
#define N 1001

    if (!avail) {
        cur = new hash[N];
	cur->next = manage;
	manage = cur;
	avail = N-1;
	cur++;
    }
    avail--; h = cur++;
    return(h);
}



void SurfaceGen::freehash() 
{
    hash *h, *h1;
    for(h = manage; h;) {
	h1 = h->next; delete h; h = h1;
    }
}


//
// This subroutine stores a polygon (triangle) in a list of vertices 
// and connectivity.  This list can then be written out in different 
// file formats. 
//
int SurfaceGen::add_polygon(float *p1,float *p2,float *p3)
{
    unsigned size, offset;
    float *ptr;
    //
    // re-allocate more space for verts if needed
    //
    if (NUM_VERTICES >= (VERT_LIMIT - 3)) {
	VERT_LIMIT += VERT_INCR;

        float *nvert = new float[VERT_LIMIT*3];
        if(!nvert) return(0);

        float *p1 = nvert;
        float *p2 = VERTICES;
        for(int i=0;i<NUM_VERTICES*3;i++)
            *p1++ = *p2++;

        if(VERTICES) delete [] VERTICES;
        VERTICES = nvert;
    }

    // store the vertices 
    ptr = VERTICES + NUM_VERTICES * 3;
    *ptr++ = *p1++;		// x of first vertex 
    *ptr++ = *p1++;		// y of first vertex 
    *ptr++ = *p1++;		// z of first vertex 
    *ptr++ = *p2++;		// x of second vertex
    *ptr++ = *p2++;		// y of second vertex 
    *ptr++ = *p2++;		// z of second vertex 
    *ptr++ = *p3++;		// x of third vertex 
    *ptr++ = *p3++;		// y of third vertex 
    *ptr++ = *p3++;		// z of third vertex

    NUM_VERTICES += 3;
	
    return(1);
}

