///////////////////////////////////////////////////////////////////////////
//
// File: Grow.C
//
// Author: Sarang Joshi		(Algorithm, original C code)
//	   Keith Doolittle	(C++)
//
// Purpose:
//	   Grow surface to a volume
//
///////////////////////////////////////////////////////////////////////////
#include <iostream.h>
#include <ADL/Array1D.h>
#include <TDL/Surface.h>
#include <TDL/ItXVolume.h>
#include <TDL/ItXVolumeUtils.h>
#include <TDL/SurfaceGrow.h>


void SurfaceGrow::cleanup()
{
if(vert) {
	delete [] vert;
	vert = NULL;
	}
if(deform) {
	delete [] deform;
	deform = NULL;
	}
if(newvert) {
	delete [] newvert;
	newvert = NULL;
	}
if(face) {
	delete [] face;
	face = NULL;
	}
if(normals) {
	delete [] normals;
	normals = NULL;
	}
 if (flag) {
	 delete [] flag;
	 flag = NULL;
 }
 int i;
if(vert_nbhd) {
	for(i=0;i<numpoints;i++) 
	    if(vert_nbhd[i].neigh)
		delete [] vert_nbhd[i].neigh;
	delete [] vert_nbhd;
	vert_nbhd = NULL;
	}
if(face_nbhd) {
	for(i=0;i<numpoints;i++) 
	    if(face_nbhd[i].neigh)
		delete [] face_nbhd[i].neigh;
	delete [] face_nbhd;
	face_nbhd = NULL;
	}

if(cubetable) {
	for(i=0;i<ncubetable;i++)
		delete [] cubetable[i];
	delete [] cubetable;
	cubetable = NULL;
	}

if(numtricube) {
	delete [] numtricube;
	numtricube = NULL;
	}

if(tri_compared) {
	delete [] tri_compared;
	tri_compared = NULL;
	}


numfaces        = 0;
numface_alloced = 0;
numpoints       = 0;
numvert_alloced = 0;
}



void SurfaceGrow::findCubeRegion(int inda, int indb, int indc,
                int *outxmin, int *outxmax,
                int *outymin, int *outymax,
                int *outzmin, int *outzmax)
{
int x,y,z,xmin,xmax,ymin,ymax,zmin,zmax;

/* find bounding cube of hash indexes */
xmin = ymin = zmin = 0;
xmax = ymax = zmax = 0;

if(hashnumbers(inda,&x,&y,&z)) {
     xmin = xmax = x;
     ymin = ymax = y;
     zmin = zmax = z;
     }
if(hashnumbers(indb,&x,&y,&z)) {
     if(x<xmin) xmin = x;
     if(x>xmax) xmax = x;
     if(y<ymin) ymin = y;
     if(y>ymax) ymax = y;
     if(z<zmin) zmin = z;
     if(z>zmax) zmax = z;
     }
if(hashnumbers(indc,&x,&y,&z)) {
     if(x<xmin) xmin = x;
     if(x>xmax) xmax = x;
     if(y<ymin) ymin = y;
     if(y>ymax) ymax = y;
     if(z<zmin) zmin = z;
     if(z>zmax) zmax = z;
     }

*outxmin = xmin;
*outxmax = xmax;
*outymin = ymin;
*outymax = ymax;
*outzmin = zmin;
*outzmax = zmax;
}


int SurfaceGrow::hashindex(int x, int y, int z)
{
return(x + hfact*y + hfact*hfact*z);
}

/* single index to 3 coord */
int SurfaceGrow::hashnumbers(int indx, int *xret, int *yret, int *zret)
{
int x,y,z;

if(indx < 0) 
   indx = 0;
else if(indx>=(hfact*hfact*hfact))
   indx = (hfact*hfact*hfact)-1;

z = indx/(hfact*hfact);
indx -= (z*hfact*hfact);
y = indx/hfact;
indx -= (y*hfact);
x = indx;

*xret = x;
*yret = y;
*zret = z;
return(1);
}

/* vertex to single coord */
int SurfaceGrow::hashcube(float x,float y,float z) {
int ind;
int xin,yin,zin;

xin = (int)(x/hsizex);
yin = (int)(y/hsizey);
zin = (int)(z/hsizez);

if (xin < 0 ) xin = 0;
if (xin >= hfact ) xin = hfact-1;
if (yin < 0 ) yin = 0;
if (yin >= hfact ) yin = hfact-1;
if (zin < 0 ) zin = 0;
if (zin >= hfact ) zin = hfact-1;

ind = hashindex(xin,yin,zin);
return(ind);
}

void SurfaceGrow::cross(float u1[3], float u2[3], float u[3])
{
  float amp;

  u[0] = u1[1]*u2[2] - u1[2]*u2[1];
  u[1] = u2[0]*u1[2] - u2[2]*u1[0];
  u[2] = u1[0]*u2[1] - u1[1]*u2[0];
}

/*
 * compute the norm of u
 */
void SurfaceGrow::norm(float u[3])
{
  float amp;

  amp = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);

  if (amp == 0.0 ) {
    u[0] = 0;
    u[1] = 0;
    u[2] = 0;
    amp = 1.0;
  }

  u[0] = u[0]/amp;
  u[1] = u[1]/amp;
  u[2] = u[2]/amp;
}


void SurfaceGrow::computenormals()
{
  int v1,v2,v3;
  int i,j,k;
  float p1[3],p2[3],p3[3],nor[3];
       
  // Initialize

  for(j=0;j<numpoints;j++){
    normals[j].x = 0;
    normals[j].y = 0;
    normals[j].z = 0;
    }

  // Now compute nal at each vertex 
  // and add to the previously computed nal

  for (i=0;i<numfaces;i++) {
    v1 = face[i].a;
    v2 = face[i].b;
    v3 = face[i].c;

    p1[0] = vert[v3].x-vert[v1].x;
    p1[1] = vert[v3].y-vert[v1].y;
    p1[2] = vert[v3].z-vert[v1].z;

    p2[0] = vert[v1].x-vert[v2].x;
    p2[1] = vert[v1].y-vert[v2].y;
    p2[2] = vert[v1].z-vert[v2].z;

    p3[0] = vert[v2].x-vert[v3].x;
    p3[1] = vert[v2].y-vert[v3].y;
    p3[2] = vert[v2].z-vert[v3].z;

    cross(p1,p2,nor);

    normals[v1].x = normals[v1].x-nor[0]/2;
    normals[v1].y = normals[v1].y-nor[1]/2;
    normals[v1].z = normals[v1].z-nor[2]/2;

    cross(p2,p3,nor);

    normals[v2].x = normals[v2].x-nor[0]/2;
    normals[v2].y = normals[v2].y-nor[1]/2;
    normals[v2].z = normals[v2].z-nor[2]/2;

    cross(p3,p1,nor);

    normals[v3].x = normals[v3].x-nor[0]/2;
    normals[v3].y = normals[v3].y-nor[1]/2;
    normals[v3].z = normals[v3].z-nor[2]/2;
    }

  // Normalize to unit normals

  for(j=0;j<numpoints;j++){
    p1[0] = normals[j].x;
    p1[1] = normals[j].y;
    p1[2] = normals[j].z;
    norm(p1);
    normals[j].x = -p1[0];
    normals[j].y = -p1[1];
    normals[j].z = -p1[2];
    }
}


void SurfaceGrow::resizeNeighborArray(neighborhood **lst, 
				      int numcpy, int num2alloc)
{
neighborhood *oldlst = *lst;
neighborhood *newlst = new neighborhood[num2alloc];
if(!newlst) {
	if(oldlst) delete [] oldlst;
	*lst = NULL;
	return;
	}
if(oldlst) {
 for(int i=0;i<numcpy;i++) {
	newlst[i].num = oldlst[i].num;
	newlst[i].num_alloced = oldlst[i].num_alloced;
	newlst[i].neigh = oldlst[i].neigh;
	}
 delete [] oldlst;
 }

*lst = newlst;
}



int SurfaceGrow::addFacet(int a, int b, int c)
{
if(numfaces >= numface_alloced) {
	numface_alloced += numfaces/10;

	triangle *nface = new triangle[numface_alloced];
	if(!nface) return(-1);

        for(int j=0;j<numfaces;j++) {
		nface[j].a = face[j].a;
		nface[j].b = face[j].b;
		nface[j].c = face[j].c;
		}
	delete [] face;
	face = nface;

	if(tri_compared)
	   delete [] tri_compared;
	tri_compared = new int[numface_alloced];
	if(!tri_compared) return(-1);
	}

// set facet

face[numfaces].a = a;
face[numfaces].b = b;
face[numfaces].c = c;

numfaces++;
return(1);
}


int SurfaceGrow::addPoint(float x, float y, float z)
{
		//cout<<" addpoint called "<<endl;
if(numpoints >= numvert_alloced) {
	numvert_alloced += numpoints/10;

	// resize vertex array

	point *nvert = new point[numvert_alloced];
	if(!nvert) return(-1);

	float *newFlag =new float[numvert_alloced];
	if (!newFlag) return(-1);
	
	
	for(int j=0;j<numpoints;j++) {
		nvert[j].x = vert[j].x;
		nvert[j].y = vert[j].y;
		nvert[j].z = vert[j].z;
		newFlag[j] = flag[j];
	}
	delete [] vert;
	vert = nvert;

	delete flag;
	flag = newFlag;
	
	// resize deform array 
	// do not have to initialize 

	if(deform) delete [] deform;
	deform  = new float[numvert_alloced];
	if(!deform) return(-1);

	// resize normals array 
	// do not have to initialize -- see refinetriang()

	if(normals) delete [] normals;
	normals  = new point[numvert_alloced];
	if(!normals) return(-1);

	// resize new vertex array

	if(newvert) delete [] newvert;
	newvert = new point[numvert_alloced];
	if(!newvert) return(-1);

	// resize vertex neighborhood array

	resizeNeighborArray(&vert_nbhd,numpoints,numvert_alloced);
	if(!vert_nbhd) return(-1);

	resizeNeighborArray(&face_nbhd,numpoints,numvert_alloced);
	if(!face_nbhd) return(-1);
	}

// set point

vert[numpoints].x = x;
vert[numpoints].y = y;
vert[numpoints].z = z;

// initialize this points vertex neighbor list

vert_nbhd[numpoints].num          = 0;
vert_nbhd[numpoints].num_alloced  = 0;
vert_nbhd[numpoints].neigh        = NULL;

// initialize this points triangle neighbor list

face_nbhd[numpoints].num          = 0;
face_nbhd[numpoints].num_alloced  = 0;
face_nbhd[numpoints].neigh        = NULL;

numpoints++;
return(numpoints-1);
}


//
// Add vertex 'n' as neighbor of vertex 'i'
//
void SurfaceGrow::addVertNeighbor(int i, int n)
{
if((i<0)||(i>=numpoints)) {
	cerr << "ERROR: SurfaceGrow::addVertNeighbor() invalid index" << endl;
	return;
	}
if(vert_nbhd[i].neigh == NULL) {
	vert_nbhd[i].neigh       = new int[10];
	vert_nbhd[i].num         = 0;
	vert_nbhd[i].num_alloced = 10;
	}
else if(vert_nbhd[i].num >= vert_nbhd[i].num_alloced) {
	vert_nbhd[i].num_alloced += 10;
	int *nneigh = new int[vert_nbhd[i].num_alloced];
	for(int j=0;j<vert_nbhd[i].num;j++)
		nneigh[j] = vert_nbhd[i].neigh[j];
	delete [] vert_nbhd[i].neigh;
	vert_nbhd[i].neigh = nneigh;
	}

// should not have to check for dups

vert_nbhd[i].neigh[vert_nbhd[i].num] = n;
vert_nbhd[i].num++;
}



//
// Add triangle 'n' as a neighbor of vertex 'i'
//
void SurfaceGrow::addFacetNeighbor(int i, int n)
{
if(face_nbhd[i].neigh == NULL) {
	face_nbhd[i].neigh       = new int[10];
	face_nbhd[i].num         = 0;
	face_nbhd[i].num_alloced = 10;
	}
else if(face_nbhd[i].num >= face_nbhd[i].num_alloced) {
	face_nbhd[i].num_alloced += 10;
	int *nneigh = new int[face_nbhd[i].num_alloced];
	for(int j=0;j<face_nbhd[i].num;j++)
		nneigh[j] = face_nbhd[i].neigh[j];
	delete [] face_nbhd[i].neigh;
	face_nbhd[i].neigh = nneigh;
	}

// should not have to check for dups

face_nbhd[i].neigh[face_nbhd[i].num] = n;
face_nbhd[i].num++;
}





ItXECode SurfaceGrow::refinetriang()
{
  int v1,v2,v3;
  int i,j,k,oldnumpoints,l;
  float p1[3],p2[3],p3[3],nor[3],area;
  float tmpx,tmpy,tmpz,in;

  oldnumpoints = numpoints;

  float totarea = 0;

if(total_addpoints >= max_addpoints)
	return(ItXSuccess);

ShowWorking("   refining triangles");

  for (i=0;i<numfaces;i++){

     //
     // Calculate the area of the face 
     //
     v1 = face[i].a;
     v2 = face[i].b;
     v3 = face[i].c;

     p1[0] = vert[v3].x-vert[v1].x;
     p1[1] = vert[v3].y-vert[v1].y;
     p1[2] = vert[v3].z-vert[v1].z;

     p2[0] = vert[v1].x-vert[v2].x;
     p2[1] = vert[v1].y-vert[v2].y;
     p2[2] = vert[v1].z-vert[v2].z;

     p3[0] = vert[v2].x-vert[v3].x;
     p3[1] = vert[v2].y-vert[v3].y;
     p3[2] = vert[v2].z-vert[v3].z;

     cross(p1,p2,nor);

     area = 0.5*sqrt(nor[0]*nor[0]+nor[1]*nor[1]+nor[2]*nor[2]);

     totarea = totarea+ area;

     if(total_addpoints >= max_addpoints) {
	ShowWorking("   stopped refining triangles. Point limit exceeded");
	break;
	}

     // Subdivide the trangle 

     if (area > areath){
	total_addpoints += 3;

        // Add the mid points to the list of vertices 

	// numpoints
	int nv0 = addPoint((vert[v1].x+vert[v2].x)/2,
			   (vert[v1].y+vert[v2].y)/2,
			   (vert[v1].z+vert[v2].z)/2);

	// numpoints + 1
	int nv1 = addPoint((vert[v1].x+vert[v3].x)/2,
			   (vert[v1].y+vert[v3].y)/2,
			   (vert[v1].z+vert[v3].z)/2);

	// numpoints + 2
	int nv2 = addPoint((vert[v2].x+vert[v3].x)/2,
			   (vert[v2].y+vert[v3].y)/2,
			   (vert[v2].z+vert[v3].z)/2);
	
	//
	// check for out of memory
	//
	if((nv0 < 0)||(nv1 < 0)||(nv2 < 0))
		return(ItXNoMemory);

        // Replace the original triangle by new 4 

        face[i].b = nv0;
        face[i].c = nv1;

	int f1 = addFacet(nv0,v2,nv2);
	int f2 = addFacet(nv1,nv2,v3);
	int f3 = addFacet(nv0,nv2,nv1);
	if((f1 < 0)||(f2 < 0)||(f3 < 0))
		return(ItXNoMemory);

        // Add the new point as the neighbour of the three points 

	// first replace old neighbors

        for (k=0;k<vert_nbhd[v1].num;k++) {
           if (vert_nbhd[v1].neigh[k] == v2) 
			vert_nbhd[v1].neigh[k] = nv0;
           if (vert_nbhd[v1].neigh[k] == v3) 
			vert_nbhd[v1].neigh[k] = nv1;
           }

        for (k=0;k<vert_nbhd[v2].num;k++) {
           if (vert_nbhd[v2].neigh[k] == v1) 
			vert_nbhd[v2].neigh[k] = nv0;
           if (vert_nbhd[v2].neigh[k] == v3) 
			vert_nbhd[v2].neigh[k] = nv2;
           }

        for (k=0;k<vert_nbhd[v3].num;k++) {
           if (vert_nbhd[v3].neigh[k] == v1) 
			vert_nbhd[v3].neigh[k] = nv1;
           if (vert_nbhd[v3].neigh[k] == v2) 
			vert_nbhd[v3].neigh[k] = nv2;
           }

	addVertNeighbor(nv0,nv1);
	addVertNeighbor(nv0,nv2);
	addVertNeighbor(nv0,v1);
	addVertNeighbor(nv0,v2);

	addVertNeighbor(nv1,nv0);
	addVertNeighbor(nv1,nv2);
	addVertNeighbor(nv1,v1);
	addVertNeighbor(nv1,v3);

	addVertNeighbor(nv2,nv1);
	addVertNeighbor(nv2,nv0);
	addVertNeighbor(nv2,v2);
	addVertNeighbor(nv2,v3);

        // Now find the adjoining triangles 

        for (k=0;k<numfaces;k++){

         if (face[k].a == v1 && face[k].b == v2 && face[k].c != v3) {
           face[k].b = nv0;
	   if(addFacet(nv0,v2,face[k].c) < 0) return(ItXNoMemory);
	   addVertNeighbor(face[k].c,nv0);
	   addVertNeighbor(nv0,face[k].c);
           }

         if (face[k].b == v1 && face[k].a == v2 && face[k].c != v3) {
           face[k].a = nv0;
	   if(addFacet(v2,nv0,face[k].c) < 0) return(ItXNoMemory);
	   addVertNeighbor(face[k].c,nv0);
	   addVertNeighbor(nv0,face[k].c);
           }

         if (face[k].c == v1 && face[k].b == v2 && face[k].a != v3) {
           face[k].c = nv0;
	   if(addFacet(face[k].a,nv0,v1) < 0) return(ItXNoMemory);
	   addVertNeighbor(face[k].a,nv0);
	   addVertNeighbor(nv0,face[k].a);
           }

         if (face[k].b == v1 && face[k].c == v2 && face[k].a != v3) {
           face[k].c = nv0;
	   if(addFacet(face[k].a,nv0,v2) < 0) return(ItXNoMemory);
	   addVertNeighbor(face[k].a,nv0);
	   addVertNeighbor(nv0,face[k].a);
           }

         if (face[k].a == v1 && face[k].c == v2 && face[k].b != v3) {
           face[k].c = nv0;
	   if(addFacet(face[k].b,v2,nv0) < 0) return(ItXNoMemory);
	   addVertNeighbor(face[k].b,nv0);
	   addVertNeighbor(nv0,face[k].b);
           }

         if (face[k].c == v1 && face[k].a == v2 && face[k].b != v3) {
           face[k].c = nv0;
	   if(addFacet(face[k].b,v1,nv0) < 0) return(ItXNoMemory);
	   addVertNeighbor(face[k].b,nv0);
	   addVertNeighbor(nv0,face[k].b);
           }

         if (face[k].a == v1 && face[k].b == v3 && face[k].c != v2) {
           face[k].b = nv1;
	   if(addFacet(nv1,v3,face[k].c) < 0) return(ItXNoMemory);
	   addVertNeighbor(face[k].c,nv1);
	   addVertNeighbor(nv1,face[k].c);
           }

         if (face[k].b == v1 && face[k].a == v3 && face[k].c != v2) {
           face[k].a = nv1;
	   if(addFacet(v3,nv1,face[k].c) < 0) return(ItXNoMemory);
	   addVertNeighbor(face[k].c,nv1);
	   addVertNeighbor(nv1,face[k].c);
           }

         if (face[k].c == v1 && face[k].b == v3 && face[k].a != v2) {
           face[k].c = nv1;
	   if(addFacet(face[k].a,nv1,v1) < 0) return(ItXNoMemory);
	   addVertNeighbor(face[k].a,nv1);
	   addVertNeighbor(nv1,face[k].a);
           }

         if (face[k].b == v1 && face[k].c == v3 && face[k].a != v2) {
           face[k].c = nv1;
	   if(addFacet(face[k].a,nv1,v3) < 0) return(ItXNoMemory);
	   addVertNeighbor(face[k].a,nv1);
	   addVertNeighbor(nv1,face[k].a);
           }

         if (face[k].a == v1 && face[k].c == v3 && face[k].b != v2) {
           face[k].c = nv1;
	   if(addFacet(face[k].b,v3,nv1) < 0) return(ItXNoMemory);
	   addVertNeighbor(face[k].b,nv1);
	   addVertNeighbor(nv1,face[k].b);
           }

         if (face[k].c == v1 && face[k].a == v3 && face[k].b != v2) {
           face[k].c = nv1;
	   if(addFacet(face[k].b,v1,nv1) < 0) return(ItXNoMemory);
	   addVertNeighbor(face[k].b,nv1);
	   addVertNeighbor(nv1,face[k].b);
           }

         if (face[k].a == v2 && face[k].b == v3 && face[k].c != v1) {
           face[k].b = nv2;
	   if(addFacet(nv2,v3,face[k].c) < 0) return(ItXNoMemory);
	   addVertNeighbor(face[k].c,nv2); // was nv1, typo??
	   addVertNeighbor(nv2,face[k].c);
           }

         if (face[k].b == v2 && face[k].a == v3 && face[k].c != v1) {
           face[k].b = nv2;
	   if(addFacet(nv2,v2,face[k].c) < 0) return(ItXNoMemory);
	   addVertNeighbor(face[k].c,nv2); // was nv1, typo??
	   addVertNeighbor(nv2,face[k].c);
           }

         if (face[k].c == v2 && face[k].b == v3 && face[k].a != v1) {
           face[k].c = nv2;
	   if(addFacet(face[k].a,nv2,v2) < 0) return(ItXNoMemory);
	   addVertNeighbor(face[k].a,nv2);
	   addVertNeighbor(nv2,face[k].a);
           }

         if (face[k].b == v2 && face[k].c == v3 && face[k].a != v1) {
           face[k].c = nv2;
	   if(addFacet(face[k].a,nv2,v3) < 0) return(ItXNoMemory);
	   addVertNeighbor(face[k].a,nv2);
	   addVertNeighbor(nv2,face[k].a);
           }

         if (face[k].a == v2 && face[k].c == v3 && face[k].b != v1) {
           face[k].c = nv2;
	   if(addFacet(face[k].b,v3,nv2) < 0) return(ItXNoMemory);
	   addVertNeighbor(face[k].b,nv2);
	   addVertNeighbor(nv2,face[k].b);
           }

         if (face[k].c == v2 && face[k].a == v3 && face[k].b != v1) {
           face[k].c = nv2;
	   if(addFacet(face[k].b,v2,nv2) < 0) return(ItXNoMemory);
	   addVertNeighbor(face[k].b,nv2);
	   addVertNeighbor(nv2,face[k].b);
           }
         } // for(k..numfaces)

        } // if(area > areath)

     } // for(i..numfaces)


   //
   // build facet neighbor list again
   // TODO: put it above
   //
   for(j=0;j<numpoints;j++)
        face_nbhd[j].num = 0;

   for(i=0;i<numfaces;i++) {
        v1 = face[i].a;
        v2 = face[i].b;
        v3 = face[i].c;
	addFacetNeighbor(v1,i);
	addFacetNeighbor(v2,i);
	addFacetNeighbor(v3,i);
        }
   // 
   // redistribute 
   // 
   for(l=0;l<redist_iter;l++){

        computenormals();

        for(j=0;j<numpoints;j++) {

            // For each point the derivative is
            // just average the neighbours minus itself!!

	    int njnum = vert_nbhd[j].num;

            tmpx = -njnum*vert[j].x;
            tmpy = -njnum*vert[j].y;
            tmpz = -njnum*vert[j].z;

            for(i = 0;i < njnum ;i++){
		int neighi = vert_nbhd[j].neigh[i];
                tmpx += vert[neighi].x;
                tmpy += vert[neighi].y;
                tmpz += vert[neighi].z;
                }

            // Now project the derivative on to the tangentspace 
            // innerproduct the derivative with the normal 
	    float norx = normals[j].x;
	    float nory = normals[j].y;
	    float norz = normals[j].z;

	    in = tmpx*norx + tmpy*nory + tmpz*norz;
            newvert[j].x = vert[j].x + 0.5*(tmpx - in*norx)/njnum;
            newvert[j].y = vert[j].y + 0.5*(tmpy - in*nory)/njnum;
            newvert[j].z = vert[j].z + 0.5*(tmpz - in*norz)/njnum;

            } // for(j..numpoints)

        // copy the newvert in to the old and compute the normals 

        for(j=0;j<numpoints;j++){
            vert[j].x = newvert[j].x;
            vert[j].y = newvert[j].y;
            vert[j].z = newvert[j].z;
            }

        } // for(l..redist_iter)

return(ItXSuccess);
}


bool SurfaceGrow::importSurface()
{
ShowWorking("   importing surface");
cleanup();

if(!surf) {
	cout << "SurfaceGrow::importSurface() invalid surface" << endl;
	return(false);
	}

numpoints = surf->mNumVert;
numfaces  = surf->mNumPoly;
//
// allocate extra facets/verts for adding
// to surface
//
numvert_alloced = numpoints + numpoints/10;
numface_alloced = numfaces  + numfaces/10;

vert      = new point[numvert_alloced];
newvert   = new point[numvert_alloced];
normals   = new point[numvert_alloced];
deform    = new float[numvert_alloced];
face      = new triangle[numface_alloced];
vert_nbhd = new neighborhood[numvert_alloced];
face_nbhd = new neighborhood[numvert_alloced];
tri_compared = new int[numface_alloced];
flag = new float[numface_alloced];
 
if(!vert || !face || !vert_nbhd || !face_nbhd || 
   !normals || !deform || !newvert || !tri_compared || !flag) {
	cerr << "ERROR: SurfaceGrow::importSurface() out of memory" << endl;
	cleanup();
	return(false);
	}
int i;
for(i=0;i<numpoints;i++) {
	vert[i].x = surf->mVert[i].x();
	vert[i].y = surf->mVert[i].y();
	vert[i].z = surf->mVert[i].z();
	vert_nbhd[i].num         = 0;
	vert_nbhd[i].num_alloced = 0;
	vert_nbhd[i].neigh       = NULL;
	face_nbhd[i].num         = 0;
	face_nbhd[i].num_alloced = 0;
	face_nbhd[i].neigh       = NULL;
	}

for(i=0;i<numfaces;i++) {
	face[i].a = surf->mFacet[i][0];
	face[i].b = surf->mFacet[i][1];
	face[i].c = surf->mFacet[i][2];
	}
return(true);
}


bool SurfaceGrow::exportSurface()
{
ShowWorking("   exporting surface");
if(!surf) {
	return(false);
	}
surf->mNumVert = numpoints;
surf->mNumPoly = numfaces;

surf->mVert.setDim(numpoints);
if(surf->mVert.isEmpty()) {
	surf->clean();
	cerr << "ERROR: SurfaceGrow::exportSurface() out of memory" << endl;
	return(false);
	}

surf->mFacet.setDim(numfaces,3);
if(surf->mFacet.isEmpty()) {
	surf->clean();
	cerr << "ERROR: SurfaceGrow::exportSurface() out of memory" << endl;
	return(false);
	}
int i;
for(i=0;i<numpoints;i++) 
	surf->mVert[i].set(vert[i].x,vert[i].y,vert[i].z);

if(surf->hasUNorms()) {
	surf->mUNorm.setDim(numpoints);
	for(i=0;i<numpoints;i++)
		surf->mUNorm[i].set(normals[i].x,normals[i].y,normals[i].z);
	}

if(surf->hasNorms()) {
	surf->mNorm.setDim(numpoints);
	for(i=0;i<numpoints;i++)
		surf->mNorm[i].set(normals[i].x,normals[i].y,normals[i].z);
	}

for(i=0;i<numfaces;i++) {
	surf->mFacet[i][0] = face[i].a;
	surf->mFacet[i][1] = face[i].b;
	surf->mFacet[i][2] = face[i].c;
	}

return(true);
}

bool SurfaceGrow::calcCubeTable()
{
   int i,j,k,ii;
   int a,b,c,n;
   int xmin,ymin,zmin,xmax,ymax,zmax,x,y,z;
   int inda,indb,indc,indx;


   for (j=0;j<hfact*hfact*hfact;j++) { numtricube[j] = 0; }

   for (ii=0;ii<numfaces;ii++) {

   a = face[ii].a;
   b = face[ii].b;
   c = face[ii].c;

   /* get hash indexes for vertices */
   inda=hashcube(vert[a].x,vert[a].y,vert[a].z);
   indb=hashcube(vert[b].x,vert[b].y,vert[b].z);
   indc=hashcube(vert[c].x,vert[c].y,vert[c].z);

   findCubeRegion(inda,indb,indc,&xmin,&xmax,&ymin,&ymax,&zmin,&zmax);

   /* add triangle to each cube in bounding box of neighborhoods */
   for(k=zmin;k<=zmax;k++)
     for(j=ymin;j<=ymax;j++)
       for(i=xmin;i<=xmax;i++) {
           indx = hashindex(i,j,k);
           n = numtricube[indx];
	   if(n >= maxpoly_percube)
		return(false);
	   else {
                cubetable[indx][n] = ii;
                numtricube[indx]++;
		}
           }

   }
return(true);
}


void SurfaceGrow::buildVertNeighborhood()
{
int i,j,addb,addc,addbc;

for(i=0;i<numpoints;i++) {
	vert_nbhd[i].num         = 0;
	vert_nbhd[i].num_alloced = 0;
	vert_nbhd[i].neigh       = NULL;
	}

for(i=0;i<numfaces;i++) {

       // Check if vertex b and c are in nghbd. of a 
       // and check if b is in nghbd of c 
       int a = face[i].a;
       int b = face[i].b;
       int c = face[i].c;

       addb = addc = addbc = 1;
       for(j=0;j<vert_nbhd[a].num;j++) {
           if(vert_nbhd[a].neigh[j] == b) 
		addb = 0;
           if(vert_nbhd[a].neigh[j] == c) 
		addc = 0;
           }

       for (j=0;j<vert_nbhd[b].num;j++) 
           if(vert_nbhd[b].neigh[j] == c) 
		addbc = 0;

       if(addbc) {
	   addVertNeighbor(c,b);
	   addVertNeighbor(b,c);
	   }
       if(addb) {
	   addVertNeighbor(a,b);
	   addVertNeighbor(b,a);
	   }
       if(addc) {
	   addVertNeighbor(a,c);
	   addVertNeighbor(c,a);
	   }	
       }
}



bool SurfaceGrow::initialize(Surface *s,const ItXVolume *i, int h) 
{
if(!s || !i)
	return(false);

ShowWorking("   initializing");

hfact = h;
img   = i;
surf  = s;

hsizex = (float)i->getSizeX()/(float)hfact;
hsizey = (float)i->getSizeY()/(float)hfact;
hsizez = (float)i->getSizeZ()/(float)hfact;

maxpoly_percube = 20;

if(importSurface() == false)
	return(false);

if(createCubeTable() == false)
	return(false);

return(true);
}



bool SurfaceGrow::createCubeTable()
{
int i,j;

sprintf(work_string,"   creating cube table (%d entries per cube)",
			maxpoly_percube);
ShowWorking(work_string);
//
// create cube table
//
if(cubetable) {
	for(i=0;i<ncubetable;i++) 
		delete [] cubetable[i];
	delete [] cubetable;
	}
ncubetable = hfact*hfact*hfact;
cubetable  = new int*[ncubetable];
numtricube = new int[ncubetable];

if(!cubetable || !numtricube)
	return(false);

for(i=0;i<ncubetable;i++) {
	numtricube[i] = 0;
	cubetable[i]  = new int[maxpoly_percube];
	if(!cubetable[i]) {
		// clean up so far
		for(j=0;j<i;j++)
			delete [] cubetable[i];
		delete [] cubetable;
		cubetable = NULL;
		return(false);
		}
	}
return(true);
}


//
// check for new polygons intersecting surface
//
// TODO: use tri_compared[] to prune
//
int SurfaceGrow::checkinter()
{
int i,j,k,l,m;
int ind,inda,indb,indc;
float a1[3],a2[3],a3[3];
float b1[3],b2[3],b3[3];
int intersect;
int tri1,tri2;
int af1,af2,af3;
int f1,f2,f3;
int ii,jj,kk;
int xmin,xmax,ymin,ymax,zmin,zmax;
int nintersect=0;

   for(j=0;j<numfaces;j++)
        tri_compared[j] = -1;

   for (j=0;j<numpoints && flag[j]==HUGE ;j++) {

	//
	// move one point
	//
        newvert[j].x = vert[j].x+deform[j]*normals[j].x;
        newvert[j].y = vert[j].y+deform[j]*normals[j].y;
        newvert[j].z = vert[j].z+deform[j]*normals[j].z;

        intersect = 0;

        for (k=0;k<face_nbhd[j].num && !intersect;k++){

            tri1  = face_nbhd[j].neigh[k];
            af1   = face[tri1].a;
            af2   = face[tri1].b;
            af3   = face[tri1].c;

            a1[0] = newvert[af1].x;
            a1[1] = newvert[af1].y;
            a1[2] = newvert[af1].z;
            a2[0] = newvert[af2].x;
            a2[1] = newvert[af2].y;
            a2[2] = newvert[af2].z;
            a3[0] = newvert[af3].x;
            a3[1] = newvert[af3].y;
            a3[2] = newvert[af3].z;

           inda=hashcube(newvert[af1].x,newvert[af1].y,newvert[af1].z);
           indb=hashcube(newvert[af2].x,newvert[af2].y,newvert[af2].z);
           indc=hashcube(newvert[af3].x,newvert[af3].y,newvert[af3].z);

           findCubeRegion(inda,indb,indc,&xmin,&xmax,&ymin,&ymax,&zmin,&zmax);

           for(kk=zmin;(kk<=zmax) && !intersect;kk++) {
             for(jj=ymin;(jj<=ymax) && !intersect;jj++) {
               for(ii=xmin;(ii<=xmax) && !intersect;ii++) {

                /* check triangle against this neighborhood */

                ind = hashindex(ii,jj,kk);

                for (i=0;i<numtricube[ind] && !intersect;i++) {
                   tri2 = cubetable[ind][i];

		   if(tri_compared[tri2] == tri1)
			continue;

                   f1 = face[tri2].a;
                   f2 = face[tri2].b;
                   f3 = face[tri2].c;
                   if((af1==f1)||(af1==f2)||(af1==f3))
                        continue;
                   if((af2==f1)||(af2==f2)||(af2==f3))
                        continue;
                   if((af3==f1)||(af3==f2)||(af3==f3))
                        continue;
                   b1[0] = newvert[f1].x;
                   b1[1] = newvert[f1].y;
                   b1[2] = newvert[f1].z;
                   b2[0] = newvert[f2].x;
                   b2[1] = newvert[f2].y;
                   b2[2] = newvert[f2].z;
                   b3[0] = newvert[f3].x;
                   b3[1] = newvert[f3].y;
                   b3[2] = newvert[f3].z;
                   if (tri_tri_intersect(a1,a2,a3,b1,b2,b3) == 1 ) {
			nintersect++;
                        intersect = 1;
                        }
                   tri_compared[tri2] = tri1;
                   }
              }
            }
          }
        }

        if(intersect) {
                newvert[j].x = vert[j].x;
                newvert[j].y = vert[j].y;
                newvert[j].z = vert[j].z;
                deform[j] = 0.0;
                }
        }

return(nintersect);
}



bool SurfaceGrow::apply(  
		int numiter,
                float bg_val,
                float fg_val,
                float alpha,
                float beta,
                float delta,
                float ath,
                float max_newpoints
                )
{
if(max_newpoints < 0) max_newpoints = 0;
else if(max_newpoints > 1.0) max_newpoints = 1.0;
total_addpoints = 0;
max_addpoints   = (int)(max_newpoints * numfaces);

//
// see if this has been initialized
//
if(!img || !surf) {
	cerr << "ERROR: SurfaceGrow::apply() not initialized yet" << endl;
	ShowWorking((char*)NULL);
	return(false);
	}

//
// make sure surface hasn't changed since initialization
//
if((surf->mNumVert != numpoints)||(surf->mNumPoly != numfaces)) {
	ShowWorking(" *** Surface changed. Re-initializing");
	if(importSurface() == false) {
			ShowWorking((char*)NULL);
			return(false);
			}
	if(createCubeTable() == false) {
			ShowWorking((char*)NULL);
			return(false);
			}
	}

areath = ath;

computenormals();
buildVertNeighborhood();

 int   i,j,k,l; 
for (j=0; j<numpoints; ++j)
	deform[j]=0.0;

 float alpha_const;
float beta_const;
float val,bv,fv,val2,val3;
 cout<<"numiter  "<<numiter<<endl;
for(k=0;k<numiter;k++) {
	maxDeform = 0;
		
	sprintf(work_string,"Iteration %d : %d polys, %d points",
		k,numfaces,numpoints);
	ShowWorking(work_string);

	if(refinetriang() != ItXSuccess) {
		ShowWorking("ERROR: possibly out of memory");
		break;
		}
	if (k==0) 
		for (int j=0; j<numpoints; ++j)
			deform[j]=0.0;
	computenormals();
	
	//
	// calculate which triangles fall in each of the
	// neighborhood 'cubes'
	//
	// Currently up to 10,000 triangles per cube
	// can be handled (maxpoly_percube = 100*100)
	//
	int ntries = 0;
	while(calcCubeTable() == false) {
		if(ntries++ > 100) {
		  cerr << "ERROR: cannot create neighborhoods" << endl;
		  cerr << "       Bailing out early" << endl;
		  break;
		  }
		//
		// increase estimate of max # cubes per 
		// neighborhood region
		//
		maxpoly_percube += 100;
		if(createCubeTable() == false) {
			cerr << "ERROR: SurfaceGrow() possibly out of memory" 
			     << endl;
			ntries = 101;
			break;
			}
		}
	if(ntries > 100)
		break;
    int ind=0;

	for(j=0; j<numpoints ;  j++) {
		if (k==0) { maxDeform=0.0;  break; }
		if ((float) fabs(deform[j]) > maxDeform)
		{ maxDeform = fabs(deform[j]);  ind=j; }
		//	else if(fabs(flag[j]) > maxDeform && flag[j]!=HUGE)
		//	{  maxDeform = fabs(flag[j]);  ind=j; }
	}
	
	cout<<k<<" maxDeform "<<maxDeform<<endl;
   int count=0;
	for (j=0; j<numpoints; ++j) {
		if (k == 0) { flag[j]=HUGE;  count++; } 
		else {
			//case 1: it has just moved and can move again
			//case 2: it didn't move in the last step, but is ready to move now
			if ((fabs(deform[j]) > 0.0001* maxDeform && flag[j] == HUGE)
				|| ( fabs(flag[j]) > 0.0001* maxDeform && flag[j] != HUGE)  )
			{ flag[j]=HUGE; count++;  }

			//case 3: hasn't moved in the last two steps 
			else if (deform[j]!=0.0) flag[j]=deform[j];  
		}
		deform[j] = 0.0;
	}
	cout<<"count "<<count<<endl;
//	sprintf(work_string," % d %f count  maxdeform ",count, maxDeform);
	// ShowWorking(work_string);
	for(l=0;l<15;l++) {
	   for(i=0;i<numpoints && flag[i]==HUGE; i++) {
	      ItXVolumeUtils::trilinearInterp(*img,
			newvert[i].x,newvert[i].y,newvert[i].z,val,val2,val3);
              // assume non-RGB (ie only 'val' valid)
	      bv = val - bg_val;
	      fv = val - fg_val;
	      alpha_const = alpha * (bv*bv - fv*fv)/1000;
	      beta_const  = 0;
	      for(j=0;j<vert_nbhd[i].num;j++) {
		  int neighij = vert_nbhd[i].neigh[j];
		  float deltx = vert[i].x-vert[neighij].x;
		  float delty = vert[i].y-vert[neighij].y;
		  float deltz = vert[i].z-vert[neighij].z;
		  float dx    = deltx*deltx + delty*delty + deltz*deltz;
		  beta_const += ((deform[i]-deform[neighij])/sqrt(dx));
		  }
	      beta_const *= beta;
              deform[i] += delta*(alpha_const-beta_const);
	      }
	   //
	   // see if moving points will cause new
	   // triangles to intersect surface
	   //
	   int ninter = checkinter();   
	   if(ninter) {
	      sprintf(work_string,"   %d intersections",ninter);
	      ShowWorking(work_string);
	      }
	   else ShowWorking("   no intersections");
	   }
	for(int yy=0;yy<numpoints;yy++) {
		vert[yy].x = newvert[yy].x;
		vert[yy].y = newvert[yy].y;
		vert[yy].z = newvert[yy].z;
		}
	}
//
// 
//
exportSurface();
ShowWorking((char*)NULL);
return(true);
}



SurfaceGrow::SurfaceGrow()
{
vert            = NULL;
newvert         = NULL;
deform          = NULL;
face            = NULL;
normals         = NULL;
vert_nbhd       = NULL;
face_nbhd       = NULL;
cubetable       = NULL;
numtricube      = NULL;
tri_compared    = NULL;
flag            = NULL; 

numpoints       = 0;
numvert_alloced = 0;
numfaces        = 0;
numface_alloced = 0; 
redist_iter     = 30; // was 100
ncubetable      = 0;

}


SurfaceGrow::~SurfaceGrow()
{
cleanup();
}




//******************************************************************
//******************************************************************
//************* Triangle intersection test *************************
//******************************************************************
//******************************************************************


/* Triangle/triangle intersection test routine,
 * by Tomas Moller, 1997.
 * See article "A Fast Triangle-Triangle Intersection Test",
 * Journal of Graphics Tools, 2(2), 1997
 *
 * int tri_tri_intersect(float V0[3],float V1[3],float V2[3],
 *                         float U0[3],float U1[3],float U2[3])
 *
 * parameters: vertices of triangle 1: V0,V1,V2
 *             vertices of triangle 2: U0,U1,U2
 * result    : returns 1 if the triangles intersect, otherwise 0
 *
 */

/* if USE_EPSILON_TEST is true then we do a check:
         if |dv|<EPSILON then dv=0.0;
   else no check is done (which is less robust)
*/
#define USE_EPSILON_TEST TRUE
#define EPSILON 0.000001

/* some macros */
#define CROSS(dest,v1,v2)                      \
              dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
              dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
              dest[2]=v1[0]*v2[1]-v1[1]*v2[0];

#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define SUB(dest,v1,v2)          \
            dest[0]=v1[0]-v2[0]; \
            dest[1]=v1[1]-v2[1]; \
            dest[2]=v1[2]-v2[2];

/* sort so that a<=b */
#define SORT(a,b)       \
             if(a>b)    \
             {          \
               float c; \
               c=a;     \
               a=b;     \
               b=c;     \
             }

#define ISECT(VV0,VV1,VV2,D0,D1,D2,isect0,isect1) \
              isect0=VV0+(VV1-VV0)*D0/(D0-D1);    \
              isect1=VV0+(VV2-VV0)*D0/(D0-D2);

#define COMPUTE_INTERVALS(VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2,isect0,isect1) \
  if(D0D1>0.0f)                                         \
  {                                                     \
    /* here we know that D0D2<=0.0 */                   \
    /* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
    ISECT(VV2,VV0,VV1,D2,D0,D1,isect0,isect1);          \
  }                                                     \
  else if(D0D2>0.0f)                                    \
  {                                                     \
    /* here we know that d0d1<=0.0 */                   \
    ISECT(VV1,VV0,VV2,D1,D0,D2,isect0,isect1);          \
  }                                                     \
  else if(D1*D2>0.0f || D0!=0.0f)                       \
  {                                                     \
    /* here we know that d0d1<=0.0 or that D0!=0.0 */   \
    ISECT(VV0,VV1,VV2,D0,D1,D2,isect0,isect1);          \
  }                                                     \
  else if(D1!=0.0f)                                     \
  {                                                     \
    ISECT(VV1,VV0,VV2,D1,D0,D2,isect0,isect1);          \
  }                                                     \
  else if(D2!=0.0f)                                     \
  {                                                     \
    ISECT(VV2,VV0,VV1,D2,D0,D1,isect0,isect1);          \
  }                                                     \
  else                                                  \
  {                                                     \
    /* triangles are coplanar */                        \
    return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2);      \
  }

/* this edge to edge test is based on Franlin Antonio's gem:
   "Faster Line Segment Intersection", in Graphics Gems III,
   pp. 199-202 */
#define EDGE_EDGE_TEST(V0,U0,U1)                      \
  Bx=U0[i0]-U1[i0];                                   \
  By=U0[i1]-U1[i1];                                   \
  Cx=V0[i0]-U0[i0];                                   \
  Cy=V0[i1]-U0[i1];                                   \
  f=Ay*Bx-Ax*By;                                      \
  d=By*Cx-Bx*Cy;                                      \
  if((f>0 && d>=0 && d<=f) || (f<0 && d<=0 && d>=f))  \
  {                                                   \
    e=Ax*Cy-Ay*Cx;                                    \
    if(f>0)                                           \
    {                                                 \
      if(e>=0 && e<=f) return 1;                      \
    }                                                 \
    else                                              \
    {                                                 \
      if(e<=0 && e>=f) return 1;                      \
    }                                                 \
  }

#define EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2) \
{                                              \
  float Ax,Ay,Bx,By,Cx,Cy,e,d,f;               \
  Ax=V1[i0]-V0[i0];                            \
  Ay=V1[i1]-V0[i1];                            \
  /* test edge U0,U1 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U0,U1);                    \
  /* test edge U1,U2 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U1,U2);                    \
  /* test edge U2,U1 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U2,U0);                    \
}

#define POINT_IN_TRI(V0,U0,U1,U2)           \
{                                           \
  float a,b,c,d0,d1,d2;                     \
  /* is T1 completly inside T2? */          \
  /* check if V0 is inside tri(U0,U1,U2) */ \
  a=U1[i1]-U0[i1];                          \
  b=-(U1[i0]-U0[i0]);                       \
  c=-a*U0[i0]-b*U0[i1];                     \
  d0=a*V0[i0]+b*V0[i1]+c;                   \
                                            \
  a=U2[i1]-U1[i1];                          \
  b=-(U2[i0]-U1[i0]);                       \
  c=-a*U1[i0]-b*U1[i1];                     \
  d1=a*V0[i0]+b*V0[i1]+c;                   \
                                            \
  a=U0[i1]-U2[i1];                          \
  b=-(U0[i0]-U2[i0]);                       \
  c=-a*U2[i0]-b*U2[i1];                     \
  d2=a*V0[i0]+b*V0[i1]+c;                   \
  if(d0*d1>0.0)                             \
  {                                         \
    if(d0*d2>0.0) return 1;                 \
  }                                         \
}





int SurfaceGrow::coplanar_tri_tri(
		float N[3],float V0[3],float V1[3],float V2[3],
               float U0[3],float U1[3],float U2[3])
{
   float A[3];
   short i0,i1;
   /* first project onto an axis-aligned plane, that maximizes the area */
   /* of the triangles, compute indices: i0,i1. */
   A[0]=fabs(N[0]);
   A[1]=fabs(N[1]);
   A[2]=fabs(N[2]);
   if(A[0]>A[1])
   {
      if(A[0]>A[2])
      {
          i0=1;      /* A[0] is greatest */
          i1=2;
      }
      else
      {
          i0=0;      /* A[2] is greatest */
          i1=1;
      }
   }
   else   /* A[0]<=A[1] */
   {
      if(A[2]>A[1])
      {
          i0=0;      /* A[2] is greatest */
          i1=1;
      }
      else
      {
          i0=0;      /* A[1] is greatest */
          i1=2;
      }
    }

    /* test all edges of triangle 1 against the edges of triangle 2 */
    EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2);
    EDGE_AGAINST_TRI_EDGES(V1,V2,U0,U1,U2);
    EDGE_AGAINST_TRI_EDGES(V2,V0,U0,U1,U2);

    /* finally, test if tri1 is totally contained in tri2 or vice versa */
    POINT_IN_TRI(V0,U0,U1,U2);
    POINT_IN_TRI(U0,V0,V1,V2);

    return 0;
}



int SurfaceGrow::tri_tri_intersect(float V0[3],float V1[3],float V2[3],
                      float U0[3],float U1[3],float U2[3])
{
  float E1[3],E2[3];
  float N1[3],N2[3],d1,d2;
  float du0,du1,du2,dv0,dv1,dv2;
  float D[3];
  float isect1[2], isect2[2];
  float du0du1,du0du2,dv0dv1,dv0dv2;
  short index;
  float vp0,vp1,vp2;
  float up0,up1,up2;
  float b,c,max;

  /* compute plane equation of triangle(V0,V1,V2) */
  SUB(E1,V1,V0);
  SUB(E2,V2,V0);
  CROSS(N1,E1,E2);
  d1=-DOT(N1,V0);
  /* plane equation 1: N1.X+d1=0 */

  /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
  du0=DOT(N1,U0)+d1;
  du1=DOT(N1,U1)+d1;
  du2=DOT(N1,U2)+d1;

  /* coplanarity robustness check */
#if USE_EPSILON_TEST==TRUE
  if(fabs(du0)<EPSILON) du0=0.0;
  if(fabs(du1)<EPSILON) du1=0.0;
  if(fabs(du2)<EPSILON) du2=0.0;
#endif
  du0du1=du0*du1;
  du0du2=du0*du2;

  if(du0du1>0.0f && du0du2>0.0f) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */

  /* compute plane of triangle (U0,U1,U2) */
  SUB(E1,U1,U0);
  SUB(E2,U2,U0);
  CROSS(N2,E1,E2);
  d2=-DOT(N2,U0);
  /* plane equation 2: N2.X+d2=0 */

  /* put V0,V1,V2 into plane equation 2 */
  dv0=DOT(N2,V0)+d2;
  dv1=DOT(N2,V1)+d2;
  dv2=DOT(N2,V2)+d2;

#if USE_EPSILON_TEST==TRUE
  if(fabs(dv0)<EPSILON) dv0=0.0;
  if(fabs(dv1)<EPSILON) dv1=0.0;
  if(fabs(dv2)<EPSILON) dv2=0.0;
#endif

  dv0dv1=dv0*dv1;
  dv0dv2=dv0*dv2;

  if(dv0dv1>0.0f && dv0dv2>0.0f) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */

  /* compute direction of intersection line */
  CROSS(D,N1,N2);

  /* compute and index to the largest component of D */
  max=fabs(D[0]);
  index=0;
  b=fabs(D[1]);
  c=fabs(D[2]);
  if(b>max) max=b,index=1;
  if(c>max) max=c,index=2;

        /* this is the simplified projection onto L*/
        vp0=V0[index];
        vp1=V1[index];
        vp2=V2[index];

        up0=U0[index];
        up1=U1[index];
        up2=U2[index];

  /* compute interval for triangle 1 */
  COMPUTE_INTERVALS(vp0,vp1,vp2,dv0,dv1,dv2,dv0dv1,dv0dv2,isect1[0],isect1[1]);

  /* compute interval for triangle 2 */
  COMPUTE_INTERVALS(up0,up1,up2,du0,du1,du2,du0du1,du0du2,isect2[0],isect2[1]);

  SORT(isect1[0],isect1[1]);
  SORT(isect2[0],isect2[1]);

  if(isect1[1]<isect2[0] || isect2[1]<isect1[0]) return 0;
  return 1;
}

