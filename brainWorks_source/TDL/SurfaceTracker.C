///////////////////////////////////////////////////////////////////////////
//
// IntellX, L.L.C., Proprietary
// The contents of this file comprise confidential, proprietary and/or
// trade secret information which is the property of IntellX LLC.
// Dissemination, distribution or other disclosure of this information
// is strictly prohibited without the express permission of IntellX LLC.
// Copyright (c) 1999 IntellX, L.L.C.
// All rights reserved
//
// File: SurfaceTracker.C
//
// Author: Keith Doolittle
//
// Purpose: Tracks features on a surface (sulci,gyri,geodesi)
//
///////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include <TDL/Surface.h>
#include <TDL/SurfaceUtils.h>
#include <TDL/SurfaceTracker.h>
#include <ADL/Array1D.h>

typedef struct _snode {
        int            num;
        int            head;
        float          cost;
        } snode;


ItXECode SurfaceTracker::findPath(int a, int b, Array1D<int> &retlines)
{
snode *nodes = (snode*)vnodes;
if(!mSurface||(mNumPoints<=0)||(a<0)||(a>=mNumPoints)||(b<0)||
    (b>=mNumPoints)||!nodes) {
	cerr << "ERROR: SurfaceTracker::findPath() object not initialized" <<
		endl;
	return(ItXError);
	}
int i,j;
for(i=0;i<mNumPoints;i++) {
	nodes[i].cost = 1E10;
	nodes[i].head = -1;
	nodes[i].num  = i;	

	work_space1[i] = 0;
	work_space2[i] = 0;
	mask[i]        = 0;
	}

int *cur_work  = work_space1;
int *next_work = work_space2;
int  counter,present,ngh;
double t;

// set start node
cur_work[0] = a;

// number of nodes to be processed
int stage_count = 1; 
bool found      = false;

nodes[a].cost = 0.0;

for(i=0;i<threshold-1;i++) {
   counter = 0;
   for(j=0;j<stage_count;j++) {
	present = cur_work[j];
	if(present == b) found = true;

	for(int k=0;k<mSurface->mNbhd[present].numNeighbors();k++) {
	    ngh = mSurface->mNbhd[present].getNeighbor(k);
	    if(nodes[ngh].cost > (t=(nodes[present].cost + cost[present][k]))) {
		nodes[ngh].cost = t;
		nodes[ngh].head = present;
		if(mask[ngh] == 0) {
		    next_work[counter++] = ngh;
		    mask[ngh] = 1;
		    }
		}
	    }
	}
   stage_count = counter;
   int *swap_work   = cur_work;
   cur_work    = next_work;
   next_work   = swap_work;

   for(j=0;j<mNumPoints;j++)
	next_work[j] = mask[j] = 0;
   }

if(found) {
   int q,slen;
   for(q=b,slen=0;q>=0;q=nodes[q].head,slen++) {
	if(slen >= mNumPoints) {
		cerr << "WARNING: broken pointer in SurfaceTracker" << endl;
		return(ItXError);
		}
	}
   if(slen <= 1) 
	return(ItXError);

   retlines.setDim(slen);
   for(q=b,slen=0;q>=0;q=nodes[q].head,slen++)
	retlines[slen] = nodes[q].num;

   return(ItXSuccess);
#ifdef notdef
   int indi,indj;
   Line *retline = new Line();

   indi = -1;
   for(q=b;q>=0;q=nodes[q].head) {
	indj = nodes[q].num;
	if((indi>=0)&&(indj>=0))
	   retline->addSegment(mSurface->mVert[indi],
			       mSurface->mVert[indj]);
	indi = indj;	
	}
   return(retline);
#endif
   }
cout << "WARNING: SurfaceTracker() no path found" << endl;
return(ItXError);
}



void SurfaceTracker::setupCost()
{
if(!mSurface || (mNumPoints <= 0))
	return;

if(!mSurface->hasNbhd(1))
	mSurface->genNeighborhoods(1);

cost = new float*[mNumPoints];
if(!cost) {
	cerr << "ERROR: SurfaceTracker::setupCost() out of memory!" << endl;
	clean();
	return;
	}

int i;

// in case of early out for clean
for(i=0;i<mNumPoints;i++) 
	cost[i] = NULL;

//
// in case other entity
// has changed surface or
// curvature
//
if(!mSurface->hasCurvature())
	mSurface->genCurvature(Surface::MaxCurvature);

for(i=0;i<mNumPoints;i++) {
	int nnei = mSurface->mNbhd[i].numNeighbors();
	double curvi = mSurface->curvature()[i];
	cost[i] = new float[nnei];
	if(!cost[i]) {
		cerr << "ERROR: SurfaceTracker::setupCost() out of memory!" 
		     << endl;
		clean();
		}
	for(int j=0;j<nnei;j++) {
		int neigh = mSurface->mNbhd[i].getNeighbor(j);
		double dx = mSurface->mVert[i].x() - mSurface->mVert[neigh].x();
		double dy = mSurface->mVert[i].y() - mSurface->mVert[neigh].y();
		double dz = mSurface->mVert[i].z() - mSurface->mVert[neigh].z();
		double dist = sqrt(dx*dx+dy*dy+dz*dz);
		double curvj,ff;
		switch(mTrackType) {
		    case SurfaceTracker::Sulci:
			curvj = mSurface->curvature()[neigh];
			ff = (mincurv - 0.5 * (curvi + curvj))/0.1;
		// Changed by SCJ
	        //	if(dist == 0.0) cost[i][j] = 0.0;
		//	else		cost[i][j] = 1.0/(dist * ff * ff);
			cost[i][j] = dist * ff * ff;
			break;
		    case SurfaceTracker::Gyri:
			curvj = mSurface->curvature()[neigh];
			ff = (maxcurv - 0.5 * (curvi + curvj))/0.1;
			cost[i][j] = dist * ff * ff;
			break;
		    case SurfaceTracker::Geodesi:
			cost[i][j] = dist;
			break;
		    }
		}
	}
}

void SurfaceTracker::setSurface(Surface *s, TrackType t)
{
if(!s) {
	cerr << "ERROR: SurfaceTracker::setSurface() surface is NULL" << endl;
	return;
	}
clean();

mSurface   = s;
mNumPoints = s->getNumVert();

if(t != SurfaceTracker::SameType)
mTrackType = t;
if(mNumPoints <= 0) 
	return;

work_space1 = new int[mNumPoints];
work_space2 = new int[mNumPoints];
mask        = new int[mNumPoints];
vnodes      = (void*) new snode[mNumPoints];

if(!mSurface->hasCurvature())
	mSurface->genCurvature(Surface::MaxCurvature);

maxcurv = -1E20;
mincurv = 1E20;
for(int i=0;i<mNumPoints;i++)  {
	if(mSurface->curvature()[i] > maxcurv)
		maxcurv = mSurface->curvature()[i];
	if(mSurface->curvature()[i] < mincurv)
		mincurv = mSurface->curvature()[i];
	}

s->genNeighborhoods(1);

setupCost();
}


void SurfaceTracker::setTrackType(TrackType t)
{
if(mTrackType != t) {
	mTrackType = t;
	setupCost();
	}
}


SurfaceTracker::SurfaceTracker()
{
mTrackType  = SurfaceTracker::Sulci;
mSurface    = NULL;
vnodes      = (void*)NULL;
cost        = (float**)NULL;
mask        = (int*)NULL;
work_space1 = (int*)NULL;
work_space2 = (int*)NULL;
mNumPoints  = 0;
threshold   = 1000;
}

SurfaceTracker::~SurfaceTracker()
{
clean();
}

void SurfaceTracker::clean()
{
for(int i=0;i<mNumPoints;i++)
     if(cost[i])
	delete [] cost[i];
if(cost) 	{ delete []          cost;        cost = NULL; }
if(vnodes) 	{ delete [] (snode*)vnodes;     vnodes = NULL; }
if(work_space1) { delete []   work_space1; work_space1 = NULL; }
if(work_space1) { delete []   work_space2; work_space2 = NULL; }
}

