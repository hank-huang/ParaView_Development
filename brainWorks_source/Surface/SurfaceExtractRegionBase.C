///////////////////////////////////////////////////////////////////////////
//
// File: SurfaceExtractRegionBase.C
//
// Author: Michael Bowers
//
// Purpose: Holds data and performs calculations related to extracting
//          a surface based on a region.
//
///////////////////////////////////////////////////////////////////////////
#include <stdlib.h>
#include <iostream.h>
#include <TDL/IDLdefines.h>
#include <TDL/SurfaceUtils.h>
#include <TDL/ItXVolume.h>
#include <SurfaceExtractRegionBase.h>


SurfaceExtractRegionBase::~SurfaceExtractRegionBase()
{
}


SurfaceExtractRegionBase::SurfaceExtractRegionBase(Surface *S)
{
  normalSize     = 1.f;

  bNonZeroSeg      = true;
  bNonZeroImage    = true;
  bSegEquals       = false;
  bSegNotEquals    = false;
  bNeighborNotZero = false;
  bNormalNotZero   = false;
  bNormalFaces     = true;

  if(!S) {
    return;
  }

  surf  = S;

  int l,nv = surf->getNumVert();

  newSurf  = new Surface(*surf);

  // initially, no vertices are removed
  delVert.setDim(nv);
  delVert = 0U;

  // define new vertex mapping (from orig surface to new surface)
  // originally identity

  old2NewIndex.setDim(nv);
  new2OldIndex.setDim(nv);

  for(l=0;l<nv;l++) {
    old2NewIndex[l] = l;
    new2OldIndex[l] = l;
  }

}

void SurfaceExtractRegionBase::reset()
{
  delVert = 0U;
  for(int l=0;l<surf->getNumVert();l++) {
    old2NewIndex[l] = l;
    new2OldIndex[l] = l;
  }

  updateSurface();
}

void SurfaceExtractRegionBase::updateSurface()
{
  *newSurf = *surf;

   // first, define what new mapping will be
   int i,cnt;
   for(i=0,cnt=0;i<surf->getNumVert();i++) {
    if(delVert[i]) { cnt++; old2NewIndex[i] = -1; } 
    else { 
          old2NewIndex[i]     = i - cnt;   
          new2OldIndex[i-cnt] = i;   
        }
   }

   SurfaceUtils::removeSurfacePoints(*newSurf,delVert);

}

void SurfaceExtractRegionBase::calculate(ItXVolume *mriVol, ItXVolume *segVol)
{
   if (!surf || !mriVol || !segVol)
    return;
   
   int i,a,b,c;
   int nv = surf->getNumVert();
   int nf = surf->getNumPoly();

   delVert = 0U;
   for(i=0;i<nv;i++) {
    old2NewIndex[i] = i;
    new2OldIndex[i] = i;
   }

   bool   keep;
   int    ix,iy,iz,ii,jj,kk;
   int    xst,xnd,yst,ynd,zst,znd;
   float  nrmx,nrmy,nrmz,mrival,normalimgval,nfact;
   u_char segval,normalsegval;

   int c1 = 0;
   int c2 = 0;
   int c3 = 0;
   int c4 = 0;
   int c5 = 0;
   int c6 = 0;
   int c7 = 0;

   if((bNormalNotZero||bNormalFaces) && !(surf->hasUNorms()))
    surf->genUnitNormals();

   if(bNormalFaces)
      cerr << "Removing vertices where normal faces value < " <<
    normalValueFaces << endl;

   int nbad = 0;
   for(i=0;i<nv;i++) {
        ix = (int)(surf->vertices()[i].x() + 0.5);
        iy = (int)(surf->vertices()[i].y() + 0.5);
        iz = (int)(surf->vertices()[i].z() + 0.5);

    if((ix<0)||(ix>=nx)||(iy<0)||(iy>=ny)||(iz<0)||(iz>=nz)) {
        nbad++;
        continue;
    }
    
        segval = segVol->u_char_data()[iz][iy][ix];
        mrival = mriVol->pointValue(ix,iy,iz);

    // perform the tests

        keep = true;

    if((ix<minx)||(ix>maxx)||(iy<miny)||(iy>maxy)||(iz<minz)||(iz>maxz))
       keep = false;

    if(keep && bNonZeroSeg && (segval == 0U))
        keep = false, c1++;

    if(keep && bNonZeroImage && (mrival <= 0.f))
        keep = false, c2++;

    if(keep && bSegEquals && (segval != segValueEquals))
        keep = false, c3++;

    if(keep && bSegNotEquals && (segval == segValueNotEquals))
        keep = false, c4++;

    if(keep && bNeighborNotZero) {
        xst = ix-1; if(xst < 0)   xst = 0;
        xnd = ix+1; if(xnd >= nx) xnd = nx-1;
        yst = iy-1; if(yst < 0)   yst = 0;
        ynd = iy+1; if(ynd >= ny) ynd = ny-1;
        zst = iz-1; if(zst < 0)   zst = 0;
        znd = iz+1; if(znd >= nz) znd = nz-1;
        for(kk=zst;(kk<=znd) && keep; kk++)
          for(jj=yst;(jj<=ynd) && keep; jj++)
            for(ii=xst;(ii<=xnd) && keep; ii++)
            if(segVol->u_char_data()[kk][jj][ii] == 0U)
                keep = false,c5++;
        }

    if(keep && (bNormalNotZero || bNormalFaces)) {
        nrmx = surf->unormals()[i].x() * normalSize;
        nrmy = surf->unormals()[i].y() * normalSize;
        nrmz = surf->unormals()[i].z() * normalSize;
        xst  = (int)((float)ix + nrmx + 0.5);
        yst  = (int)((float)iy + nrmy + 0.5);
        zst  = (int)((float)iz + nrmz + 0.5);
        if(xst < 0) xst = 0; if(xst >= nx) xst = nx-1;
        if(yst < 0) yst = 0; if(yst >= ny) yst = ny-1;
        if(zst < 0) zst = 0; if(zst >= nz) zst = nz-1;

        normalimgval = mriVol->pointValue(xst,yst,zst);
        normalsegval = segVol->u_char_data()[zst][yst][xst];

        if((bNormalNotZero && 
                   (normalimgval <= 0.f) && (normalsegval == 0U))||
                   (bNormalFaces && (normalsegval <= normalValueFaces)))
            keep = false,c6++;

    }

    if(keep) delVert[i] = 0U;
    else     delVert[i] = 1U;
   }

   // Now make sure verts are deleted only when all
   // 3 verts of a poly are invalid

   for(i=0;i<nf;i++) {
    a = surf->facets()[i][0];
    b = surf->facets()[i][1];
    c = surf->facets()[i][2];
    if(delVert[a] && delVert[b] && delVert[c])
       delVert[a] = delVert[b] = delVert[c] = 2U;
   }

   for(i=0;i<nv;i++) {
    if(delVert[i] == 2U) {
       delVert[i] = 1U;
       c7++;
    }
    else 
       delVert[i] = 0U;
   }

   this->updateSurface();
}



