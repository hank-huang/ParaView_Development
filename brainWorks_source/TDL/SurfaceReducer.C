///////////////////////////////////////////////////////////////////////////
//
// File: SurfaceReducer.C
//
// Author: Keith Doolittle
//
// Purpose: Surface polygon reduction class
// 
///////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <TDL/Surface.h>
#include <TDL/SurfaceReducer.h>

#define SPOLY_DALLOC	200
#define VEC(a,b,v) {    \
        v[0] = a->x - b->x;     \
        v[1] = a->y - b->y;     \
        v[2] = a->z - b->z;     \
        }

#define DIST(f,v1) f = (float)sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2])


//
// Fortran matrix factorization sub-routines
//
extern "C" {
extern void ssytrf_(char *uplo, int *n, float *a, int *lda, int *ipiv,
                float *work, int *lwork, int *info);
extern void ssytrs_(char *uplo, int *n, int *nrhs, float *a, int *lda,
                int *ipiv, float *b, int *ldb, int *info);
extern void ssyev_(char *jobz, char *uplo, int *n, float *a, int *lda,
                float *w, float *work, int *lwork, int *info);
        };



void SurfaceReducer::crossvec(float *u1, float *u2, float *u)
{
u[0] = u1[1]*u2[2]-u1[2]*u2[1];
u[1] = u2[0]*u1[2]-u2[2]*u1[0];
u[2] = u1[0]*u2[1]-u1[1]*u2[0];
}

void SurfaceReducer::crossproduct(spoint *a, spoint *b, spoint *c, float *nor)
{
float v1[3],v2[3];
v1[0] = b->x - a->x;
v1[1] = b->y - a->y;
v1[2] = b->z - a->z;
v2[0] = c->x - b->x;
v2[1] = c->y - b->y;
v2[2] = c->z - b->z;
crossvec(v1,v2,nor);
}


void SurfaceReducer::normvec(float *u)
{
float amp;

amp = (float)sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
if(amp==0.0) amp=1.0;

u[0] = u[0]/amp;
u[1] = u[1]/amp;
u[2] = u[2]/amp;
}

void SurfaceReducer::AddNeighbor(spoint *p, spoint *n)
{
int i;

if(p==n) return;
for(i=0;i<p->num;i++)
  if(p->neighbors[i] == n)
        return;

if(p->num >= MAX_NEIGHBORS-1) {
	cerr << "ERROR: surface has too many neighbors" << endl;
	cerr << "       for polygon decimation code!" << endl;
        return;
	}

p->neighbors[p->num++] = n;
}

void SurfaceReducer::AddNeighbor2(spoint *p, spoint *n, spoint **parr)
{
int i;

if(p==n) return;
for(i=0;i<p->num;i++)
  if(parr[i] == n)
        return;

if(p->num >= MAX_CURV_NEIGHBORS-1) {
	cerr << "ERROR: surface has too many neighbors" << endl;
	cerr << "       for polygon decimation code!" << endl;
        return;
	}

parr[p->num++] = n;
}


int SurfaceReducer::IsInNeighborhood(spoint *tst, spoint *pt)
{
int i;

for(i=0;i<pt->num;i++)
  if(pt->neighbors[i] == tst)
        return(1);
return(0);
}


void SurfaceReducer::ComputeNormals(int shw_wrk)
{
if(shw_wrk)
    ShowWorking("Calculating Normals\n");

spoint *pt;
spoly  *pol;
float p1[3],p2[3],p3[3],nor[3];

for(pt=spoints;pt;pt=pt->next)
   pt->nx = pt->ny = pt->nz = 0.0;

for(pol=spolys;pol;pol=pol->next) {
        p1[0] = pol->c->x - pol->a->x;
        p1[1] = pol->c->y - pol->a->y;
        p1[2] = pol->c->z - pol->a->z;

        p2[0] = pol->a->x - pol->b->x;
        p2[1] = pol->a->y - pol->b->y;
        p2[2] = pol->a->z - pol->b->z;

        p3[0] = pol->b->x - pol->c->x;
        p3[1] = pol->b->y - pol->c->y;
        p3[2] = pol->b->z - pol->c->z;

        crossvec(p1,p2,nor);
        pol->a->nx += nor[0];
        pol->a->ny += nor[1];
        pol->a->nz += nor[2];

        crossvec(p2,p3,nor);
        pol->b->nx += nor[0];
        pol->b->ny += nor[1];
        pol->b->nz += nor[2];

        crossvec(p3,p1,nor);
        pol->c->nx += nor[0];
        pol->c->ny += nor[1];
        pol->c->nz += nor[2];
        }

for(pt=spoints;pt;pt=pt->next) {
        p1[0] = pt->nx;
        p1[1] = pt->ny;
        p1[2] = pt->nz;
        normvec(p1);
        pt->nx = p1[0];
        pt->ny = p1[1];
        pt->nz = p1[2];
        }
}


void SurfaceReducer::TagEdges()
{
spoint *pt;
spoly  *pol;

for(pt=spoints;pt;pt=pt->next)
   pt->used = 0;

for(pol=spolys;pol;pol=pol->next) {
        pol->a->used++;
        pol->b->used++;
        pol->c->used++;
        }

/* if # of neighbors == number of polys that use point */
/* then it is not an edge point and can be shifted     */

for(pt=spoints;pt;pt=pt->next)
        pt->isedge = (pt->used != pt->num);
}


void SurfaceReducer::Redistribute()
{
ShowWorking("Redistributing points\n");

float tmpx,tmpy,tmpz;
int   num;
float in;
int   i,j,l;

spoint *pt;
float  px,py,pz;

TagEdges();

for(l=0;l<redist_iter;l++) {

    for(pt=spoints,j=0;pt;pt=pt->next,j++) {
        if(pt->isedge) continue; /* edge point */

        num = pt->num;
        px  = pt->x;
        py  = pt->y;
        pz  = pt->z;

        newvert[j].x = px;
        newvert[j].y = py;
        newvert[j].z = pz;

        /* For each point the derivative is
        just average the neighbours minus itself!! */
        tmpx= -num*px;
        tmpy= -num*py;
        tmpz= -num*pz;
        for(i=0;i<num;i++) {
                tmpx += pt->neighbors[i]->x;
                tmpy += pt->neighbors[i]->y;
                tmpz += pt->neighbors[i]->z;
                }
        in = tmpx*pt->nx+tmpy*pt->ny+tmpz*pt->nz;
        newvert[j].x = px + redist_alpha*(tmpx - in*pt->nx);
        newvert[j].y = py + redist_alpha*(tmpy - in*pt->ny);
        newvert[j].z = pz + redist_alpha*(tmpz - in*pt->nz);

        }

    for(pt=spoints,j=0;pt;pt=pt->next,j++){
        if(pt->isedge) continue;
        pt->x = newvert[j].x;
        pt->y = newvert[j].y;
        pt->z = newvert[j].z;
        }
    ComputeNormals(0);
    }
}


int SurfaceReducer::generateLinkedArrays(Surface *surf)
{
ShowWorking("initializing\n");

spoint *spt,*sprev;
spoly  *spol,*spolprev;

sprev = NULL;

spoints = orig_spoints = new spoint[surf->mNumVert];
if(!spoints) return 1;

spolys  = orig_spolys = new spoly[surf->mNumPoly];
if(!spolys) {
	delete [] spoints;
	spoints = NULL;
	return 1;
	}
int i;
for(i=0,spt=spoints;i<surf->mNumVert;i++,spt++) {
        spt->x    = surf->mVert[i].x();
        spt->y    = surf->mVert[i].y();
        spt->z    = surf->mVert[i].z();
        spt->num  = 0;
        spt->next = NULL;
        if(sprev) {
                spt->prev   = sprev;
                sprev->next = spt;
                }
        else    spt->prev   = NULL;
        sprev = spt;
        }

spolprev=NULL;
for(i=0,spol=spolys;i<surf->mNumPoly;i++,spol++) {
	int a = surf->mFacet[i][0];
	int b = surf->mFacet[i][1];
	int c = surf->mFacet[i][2];

        spol->a = &(spoints[a]);
        spol->b = &(spoints[b]);
        spol->c = &(spoints[c]);
        spol->next = NULL;

        spol->prev = spolprev;
        if(spolprev)
                spolprev->next = spol;

        AddNeighbor(spol->a,spol->b);
        AddNeighbor(spol->a,spol->c);
        AddNeighbor(spol->b,spol->a);
        AddNeighbor(spol->b,spol->c);
        AddNeighbor(spol->c,spol->a);
        AddNeighbor(spol->c,spol->b);
        spolprev = spol;
        }
ComputeNormals(1);
return(0);
}




float SurfaceReducer::gaussCurvature(spoint *p)
{
float harr[MAX_NEIGHBORS];
float uarr[MAX_NEIGHBORS];
float varr[MAX_NEIGHBORS];
float two_uv[MAX_NEIGHBORS];
float uu[MAX_NEIGHBORS];
float vv[MAX_NEIGHBORS];
float a,b,c,d;
float b1x,b1y,b1z;
float b2x,b2y,b2z;
float b3x,b3y,b3z;
float tmp[3],rval;

spoint *neighbors[MAX_CURV_NEIGHBORS];

// Get plane equation for point
a = p->nx;
b = p->ny;
c = p->nz;
d = -((p->nx)*(p->x)+(p->ny)*(p->y)+(p->nz)*(p->z));


// Get basis vector
float mag;

    if(fabs(p->nx) > TLRNCE) {
       b1x = -1.0/p->nx*(p->ny+p->nz);
       b1y = b1z = 1.0;
       goto basis;
    }

    if(fabs(p->ny) > TLRNCE) {
        b1x = b1z = 1.0;
        b1y = -1.0/p->ny*(p->nx+p->nz);
        goto basis;
    }

    if(fabs(p->nz) > TLRNCE) {
        b1x = b1y = 1.0;
        b1z = -1.0/p->nz*(p->nx+p->ny);
    }

basis:
    mag = (float)sqrt(b1x*b1x+b1y*b1y+b1z*b1z);
    if(mag==0.0) mag = 1.0;
    b1x /= mag;
    b1y /= mag;
    b1z /= mag;

    b3x = p->nx;
    b3y = p->ny;
    b3z = p->nz;

    b2x = b3y*b1z-b1y*b3z;
    b2y = b3z*b1x-b3x*b1z;
    b2z = b3x*b1y-b3y*b1x;

    mag = (float)sqrt(b2x*b2x+b2y*b2y+b2z*b2z);
    if(mag ==0.0) mag = 1.0;
    b2x /= mag;
    b2y /= mag;
    b2z /= mag;

// Calculate neighbors of neighbors
int     old_num = p->num;
int     i,j;

// copy level 1 neighbors
for(i=0;i<old_num;i++) 
	neighbors[i] = p->neighbors[i];

spoint *np;

// calculate level 2 neighbors
for(i=0;i<old_num;i++) {
        np = neighbors[i];
        for(j=0;j<np->num;j++)
                AddNeighbor2(p,np->neighbors[j],neighbors);
        }


int num_nghbrs = p->num;

spoint *pj;

for(i=0;i<num_nghbrs; i++) {
        pj = neighbors[i];

        harr[i] = a*pj->x+b*pj->y+c*pj->z+d;
        tmp[0] = pj->x-harr[i]*p->nx-p->x;
        tmp[1] = pj->y-harr[i]*p->ny-p->y;
        tmp[2] = pj->z-harr[i]*p->nz-p->z;

        uarr[i] = tmp[0]*b1x+tmp[1]*b1y+tmp[2]*b1z;
        varr[i] = tmp[0]*b2x+tmp[1]*b2y+tmp[2]*b2z;

        two_uv[i] = 2.0*uarr[i]*varr[i];
        uu[i] = uarr[i]*uarr[i];
        vv[i] = varr[i]*varr[i];
        }

// restore original number of neighbors
p->num = old_num;

float U[3][3], BU[3];
float C[2][2],W[2];
float workc[5],work[3];
int   IPIV[4];
int INFO,  LDA, LDB, LWORK, N, NRHS;
char JOBZ, UPLO;
int LDC, LWORKC, NC;

for(i=0; i < 3; i++) {
    BU[i] = 0;
    for(j=0; j < 3; j++)
        U[i][j] = 0;
    }

for(i=0; i < num_nghbrs; i++) {
    U[0][0] += uu[i]*uu[i];
    U[0][1] += uu[i]*two_uv[i];
    U[0][2] += uu[i]*vv[i];
    U[1][1] += two_uv[i]*two_uv[i];
    U[1][2] += two_uv[i]*vv[i];
    U[2][2] += vv[i]*vv[i];
    BU[0]   += uu[i]*2.0*harr[i];
    BU[1]   += two_uv[i]*2.0*harr[i];
    BU[2]   += vv[i]*2.0*harr[i];
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
        ssytrs_(&UPLO, &N, &NRHS, (float*)U, &LDA, IPIV, BU, &LDB , &INFO);
        C[0][0] = BU[0];
        C[0][1] = BU[1];
        C[1][0] = BU[1];
        C[1][1] = BU[2];

        /* get the eigenvalue of the matrix c */
        JOBZ   = 'N';
        LDC    = 2;
        NC     = 2;
        LWORKC = 3*NC-1;
        ssyev_(&JOBZ, &UPLO, &NC, (float*)C, &LDC, W, workc, &LWORKC,&INFO);

        if(fabs(W[0]) > fabs(W[1]))
                rval = W[0];
        else    rval = W[1];

        if(rval > MAX_CURV_RES)       rval =  MAX_CURV_RES;
        else if(rval < -MAX_CURV_RES) rval = -MAX_CURV_RES;
        }
else    rval = MAX_CURV_RES;

return(rval);
}



int SurfaceReducer::BadAspectRatio(spoint *a, spoint *b, spoint *c)
{
float perim,area,asp;
float ab,ac,bc;
float v1[3],v2[3],v3[3];

VEC(a,b,v1);
VEC(b,c,v2);
VEC(a,c,v3);
DIST(ab,v1);
DIST(bc,v2);
DIST(ac,v3);
perim = ab+ac+bc;
if(perim == 0.0)
        return(1);

VEC(b,a,v1);
VEC(c,a,v2);
crossvec(v1,v2,v3);
DIST(area,v3);
area = area/2.0;

asp = area/perim;
if(asp < MIN_ASPECT)
        return(1);
else    return(0);
}



int SurfaceReducer::CreateNewPolyFromPointRemoval(spoint *pt)
{
int num,i,j;
spoint *last,*next;

if((num=pt->num)<=2)
        return(0);

//---------------------------------------------------
//  This check prevents polys collapsing to a plane
//  (like a 3-sided pyramid collapsing into 2 polys)
//---------------------------------------------------
for(i=0;i<num;i++)
  if(pt->neighbors[i]->num <= 3)
        return(0);

//
// Clear used flag in neighbors
//
for(i=0;i<num;i++)
   pt->neighbors[i]->used = 0;


nnew             = 0;
last             = pt->neighbors[0];
new_poly[nnew++] = last;
last->used       = 1;

for(i=1;i<num;i++) {
        for(j=0;j<num;j++) {
           next = pt->neighbors[j];
           if(!(next->used) && IsInNeighborhood(next,last)) {
                new_poly[nnew++] = next;
                next->used = 1;
                if((nnew >= 3) && BadAspectRatio(new_poly[0],last,next)) {
                        nnew = 0;
                        return(0);
                        }
                last = next;
                break;
                }
           }
        if(j>=num)  {
                //------------------------------------------------
                // Couldnt find closed poly from all neighbors
                // Some wierd poly arrangements cause this
                // Dont remove point, hole is probably non-convex
                //------------------------------------------------
                return(0);
                }
        }
return(1);
}


//
// find if new_poly is convex using turn test (all angles < 90)
//
int SurfaceReducer::PolyIsConvex()
{
spoint *p1,*p2,*p3;
int i;
float d,origv[3],newv[3];

if(nnew < 3)
   return(0);

// get first normal vector
p1 = new_poly[0];
p2 = new_poly[1];
p3 = new_poly[2];
crossproduct(p1,p2,p3,origv);

for(i=3;i<(nnew+2);i++) {
    //
    // now test to see if each normal vector in
    // same direction as first
    //
    p1 = p2;
    p2 = p3;
    p3 = new_poly[(i%nnew)];
    crossproduct(p1,p2,p3,newv);
    d = origv[0]*newv[0] + origv[1]*newv[1] + origv[2]*newv[2];

    if(d < 0.0)
        return(0);
    }

return(1);
}

void SurfaceReducer::RemovePoint(spoint *pt)
{
if(pt->prev == NULL)
     spoints = pt->next;
else pt->prev->next = pt->next;

if(pt->next)
   pt->next->prev = pt->prev;
}



void SurfaceReducer::RemovePoly(spoly *pol)
{
if(pol->prev == NULL)
        spolys = pol->next;
else    pol->prev->next = pol->next;

if(pol->next)
        pol->next->prev = pol->prev;
}


void SurfaceReducer::RemoveNeighbor(spoint *p, spoint *n)
{
int j,k;
for(j=0;j<p->num;j++) {
   if(p->neighbors[j] == n) {
        p->num--;
        for(k=j;k<p->num;k++)
                p->neighbors[k] = p->neighbors[k+1];
        return;
        }
   }
}


void SurfaceReducer::RemovePointAndPolys(spoint *pt)
{
spoly *pol,*npol;
int i;

RemovePoint(pt);
//
// remove point from neighboring points 'neighbor' list
//
for(i=0;i<nnew;i++)
        RemoveNeighbor(new_poly[i],pt);

for(pol=spolys;pol;pol=npol)  {
        npol = pol->next;
        if((pol->a == pt)||(pol->b == pt)||(pol->c == pt))
                RemovePoly(pol);
        }
}


int SurfaceReducer::FindOrient(spoint *a, spoint *b)
{
spoly *pol;

for(pol=spolys;pol;pol=pol->next) {
    if(((pol->a==a)&&(pol->b==b))||
       ((pol->b==a)&&(pol->c==b))||
       ((pol->c==a)&&(pol->a==b)))
                return(1);
    else if(((pol->a==b)&&(pol->b==a))||
            ((pol->b==b)&&(pol->c==a))||
            ((pol->c==b)&&(pol->a==a)))
                return(0);
    }
return(0);
}

void SurfaceReducer::CreateNewPolys()
{
spoint *a,*b,*c;
spoly  *pol;
int     i;

if(nnew <= 2)
    return;

a = new_poly[0];
b = new_poly[1];
for(i=2;i<nnew;i++) {
        c = new_poly[i];

	if(nheap >= maxheap) {
		maxheap += SPOLY_DALLOC;
		spoly **newheap = new spoly*[maxheap];

		if(!newheap) {
			cerr << "ERROR: SurfaceReducer out of memory" << endl;
			}
		else	{
			memcpy((void*)newheap,(const void*)spoly_heap,
				nheap * sizeof(spoly*));
			delete [] spoly_heap;
			spoly_heap = newheap;
			}
		}
	if(nheap >= maxheap) return;

	pol = new spoly;
	spoly_heap[nheap++] = pol;

        if(cur_isccw) {
                pol->a = c;
                pol->b = b;
                pol->c = a;
                }
        else {  pol->a = a;
                pol->b = b;
                pol->c = c;
                }
        AddNeighbor(a,b);
        AddNeighbor(a,c);
        AddNeighbor(b,a);
        AddNeighbor(c,a);
AddNeighbor(b,c);
AddNeighbor(c,b);
        pol->next = spolys;
        pol->prev = NULL;
        if(spolys != NULL)
                spolys->prev = pol;
        spolys = pol;
        b = c;
        }
}



int SurfaceReducer::reducePolysByCurvature()
{
//
// Recalculate curvature
//
spoint *pt,*npt;
for(pt=spoints;pt;pt=pt->next)
        pt->curvature = gaussCurvature(pt);

int r=0;
float curvabs;
ShowWorking("Removing points of low curvature");
for(pt=spoints,npt=NULL;pt;pt=npt)  {
        npt = pt->next;
        curvabs = fabs(pt->curvature);
        if((curvabs >= min_curvature)&&(curvabs <= max_curvature)) {
           if(CreateNewPolyFromPointRemoval(pt) && PolyIsConvex()) {
                        RemovePointAndPolys(pt);
                        cur_isccw = FindOrient(new_poly[0],new_poly[1]);
                        CreateNewPolys();
			r++;
                        if(r++%100==0)
                               ShowWorking(".");
                        }
                }
        }
ShowWorking("\n");
return(r);
}

		
//
//
//
int SurfaceReducer::setDataFromLinkedArrays(Surface *s)
{
ShowWorking("re-creating final surface\n");

spoint *spt;
spoly  *spol;
int i;
int npoints = 0;
for(spt=spoints;spt;spt=spt->next,npoints++)
        spt->used = npoints;

s->mVert.setDim(npoints);
if(s->mVert.isEmpty()) {
	cerr << "ERROR: SurfaceReducer out of memory" << endl;
	s->clean();
	return 1;
	}
s->mNumVert = npoints;

//
// might as well set normals here too
// (these are normalized, so set unit normals also)
//
s->mNorm.setDim(npoints);
s->mUNorm.setDim(npoints);
s->mNormDirty = false;

//
// todo - set neighborhoods correctly
//
s->mNbhd.setDim(0);

if(s->mNorm.isEmpty() || s->mUNorm.isEmpty()) {
	cerr << "ERROR: SurfaceReducer out of memory" << endl;
	s->clean();
	return 1;
	}

for(spt=spoints,i=0;spt;spt=spt->next,i++)  {
	s->mVert[i].set(spt->x,spt->y,spt->z);
	s->mNorm[i].set(spt->nx,spt->ny,spt->nz);
	s->mUNorm[i].set(spt->nx,spt->ny,spt->nz);
	}

int npolys = 0;
for(spol=spolys;spol;spol=spol->next,npolys++);

s->mFacet.setDim(npolys,3);
s->mNumPoly = npolys;

if(s->mFacet.isEmpty()) {
	cerr << "ERROR: SurfaceReducer out of memory" << endl;
	s->clean();
	return 1;
	}

for(spol=spolys,i=0;spol;spol=spol->next,i++) {
        s->mFacet[i][0] = spol->a->used;
        s->mFacet[i][1] = spol->b->used;
        s->mFacet[i][2] = spol->c->used;
        }

delete [] orig_spoints;
delete [] orig_spolys;
return 0;
}


ItXECode SurfaceReducer::decimate(
		Surface *surf,
		float min_curv, 
		float max_curv, 
		int maxiter, 
		float ralpha,
		int riter, 
		int &totnremoved)
{
min_curvature = min_curv;
max_curvature = max_curv;
redist_alpha  = ralpha;
redist_iter   = riter;

if(generateLinkedArrays(surf)) {
	cerr << "ERROR: SurfaceReducer out of memory!" << endl;
	return(ItXNoMemory);
	}

spoly_heap = new spoly*[SPOLY_DALLOC];
nheap = 0;
maxheap = SPOLY_DALLOC;

// Temp vertices for redistribute
newvert = new kpoint[surf->mNumVert];

totnremoved = 0;
int nremoved=1;

Redistribute();
//
// min_curv <= 0.0 means redistribute only
//
if(min_curv > 0.0)
  for(int iter =0; (nremoved > 0) && iter < maxiter; iter++)  {
	char tstr[20];
	sprintf(tstr,"Iteration %d",iter);
	ShowWorking(tstr);
	if(nremoved = reducePolysByCurvature()) {
		ComputeNormals(1);
		Redistribute();
		totnremoved += nremoved;
		}
	}

delete [] newvert;

int rv = setDataFromLinkedArrays(surf);

for(int j=0;j<nheap;j++)
	delete spoly_heap[j];
delete [] spoly_heap;

ShowWorking("DONE.\n");
ShowWorking((char*)NULL);

if(rv)  {
	cerr << "ERROR: SurfaceReducer out of memory!" << endl;
	return(ItXNoMemory);
	}
else 	{
	if(surf) surf->geometryChanged();
	return(ItXSuccess);
	}
}


