
///////////////////////////////////////////////////////////////////////////
//
// File:  SurfaceDistanceBase.C
//
// Author:  Michael Bowers
//
// Purpose:  Non GUI Surface Distance Class
//
///////////////////////////////////////////////////////////////////////////
#include <stdlib.h>
#include <iostream.h>
#include <math.h>
#include <TDL/ItXVolume.h>
#include <OpMessage.h>
#include <TDL/Surface.h>
#include <ADL/Array2D.h>
#include <OcTree.h>
#include <PointGroup.h>
#include <SurfaceDistanceBase.h>
#include <DynArray.h>
#include <time.h>

#define USE_PGROUP

#define MAX_POL_PER_VERT    10

const float SurfaceDistanceBase::InvalidPixel = -0.12345;

#include <AWHistogram.h>

bool SurfaceDistanceBase::calculate(Surface &S,
    ItXVolume &V, ItXVolume &Seg, float white_thresh,
    ItXVolume *rDist, ItXVolume *rInd,
    AWHistogram& hist, AWHistogram& abshist)
{

   char tstr[1024];
   int i,j,k,idx,npt,npol,n,minf;
   int a,b,c,f,g,mini,polyi;
   float x,y,z;
   u_char ch,sval;
   AWPoint p,pC;
   char estr[100];
   double sgn,ptdist,mind,dist,dx,dy,dz,nrmx,nrmy,nrmz,dot;

   DynArray<int> T;
   AWPointList pList;
 
   if(!&S || !&V || !&Seg || (!rDist && !rInd))
      return(false);
  
#ifndef USE_PGROUP
   OcTree ocTree;
   OcTreeNode *onode;
#else
   PointGroup pGroup;
#endif

   hist.SetSize(200,-10.,10.);
   abshist.SetSize(200,-10.,10.);

   int nx = V.getSizeX();
   int ny = V.getSizeY();
   int nz = V.getSizeZ();

   if((Seg.getSizeX() != nx)|| (Seg.getSizeY() != ny)||
      (Seg.getSizeZ() != nz)) {
      cerr << "ERROR: SurfaceDistace::calculate() SEG not ";
      cerr << "same size as IMG" << endl;
      return(false);
  
   }

   float px = V.getPixelDimensionX();
   float py = V.getPixelDimensionY();
   float pz = V.getPixelDimensionZ();

   if(px==0.f) px = 1.f;
   if(py==0.f) py = 1.f;
   if(pz==0.f) pz = 1.f;

   if(rDist) {
     rDist->setSize(nx,ny,nz,ItXVolume::Float);
     rDist->setPixelDimensionX(px);
     rDist->setPixelDimensionY(py);
     rDist->setPixelDimensionZ(pz);
   }

   if(rInd) {
     rInd->setSize(nx,ny,nz,ItXVolume::Int);
     rInd->setPixelDimensionX(px);
     rInd->setPixelDimensionY(py);
     rInd->setPixelDimensionZ(pz);
   }
  
   npt  = S.getNumVert();
   npol = S.getNumPoly();

#ifndef USE_PGROUP
   ocTree.clear();
   ocTree.setSize((float)nx*px+1.f,(float)ny*py+1.f, (float)nz*pz+1.f,3);
#else
   float szx = 2.f * px;
   float szy = 2.f * py;
   float szz = 2.f * pz;
   if(szx < szy) szx = szy;
   if(szx < szz) szx = szz;

   if(!pGroup.setSize(0.f,(float)nx*px+1.f,
                      0.f,(float)ny*py+1.f,
                      0.f,(float)nz*pz+1.f,szx)) {
      return(false);
   }
#endif
   
   //
   // Store centroid of each poly in OcTree
   //
   for(i=0;i<npol;i++) {
      a  = S.facets()[i][0];
      b  = S.facets()[i][1];
      c  = S.facets()[i][2];
      x  = S.vertices()[a].x() * px;
      y  = S.vertices()[a].y() * py;
      z  = S.vertices()[a].z() * pz;
      x += S.vertices()[b].x() * px;
      y += S.vertices()[b].y() * py;
      z += S.vertices()[b].z() * pz;
      x += S.vertices()[c].x() * px;
      y += S.vertices()[c].y() * py;
      z += S.vertices()[c].z() * pz;
      x /= 3.f;
      y /= 3.f;
      z /= 3.f;
      p.set(x,y,z,i);
      p.setValid(true);
#ifndef USE_PGROUP
      ocTree.AddPoint(p);
#else
      pGroup.addPoint(p,i);
#endif
   }

#define KMIN(a,b) (((a)<(b))?(a):(b))
#define KMAX(a,b) (((a)<(b))?(b):(a))

   int ndegen=0;

   Array1D<u_char> degen(npol);
   degen = 0U;
   for(i=0;i<npol;i++) {
      a = S.facets()[i][0];
      b = S.facets()[i][1];
      c = S.facets()[i][2];
      Point e1 = S.vertices()[a] - S.vertices()[b];
      Point e2 = S.vertices()[b] - S.vertices()[c];
      Point e3 = S.vertices()[c] - S.vertices()[a];
      double d1 = e1.innerProd(e1);
      double d2 = e2.innerProd(e2);
      double d3 = e3.innerProd(e3);
      double dmin = KMIN(d1,d2); dmin = KMIN(dmin,d3);
      double dmax = KMAX(d1,d2); dmax = KMAX(dmax,d3);

      if((d1==0.)||(d2==0.)||(d3==0.)||
         (dmin < (dmax/100.)))  {
           degen[i] = 1U;
           ndegen++;
      }
   }

   if(ndegen > 0)
     cerr << "Located " << ndegen << " degenerate polys" << endl;
   else
     cerr << "No degenerate polys" << endl;

int nbadsign = 0;

double tsum = 0.;
double tcnt = 0.;

time_t start_time = time(NULL);

   for(k=0;k<nz;k++)  {
     sprintf(tstr,"Calculating Surface Map: plane %03d/%03d",k,nz);
     OpMessage::GenericWorkProc(tstr);

     for(j=0;j<ny;j++)
       for(i=0;i<nx;i++) {

            if(rDist) rDist->float_data()[k][j][i] = InvalidPixel;
            if(rInd)  rInd->int_data()[k][j][i]    = -1;

            float vval = V.pointValue(i,j,k);
        if(vval == 0.f)
               continue;

            float sval = Seg.pointValue(i,j,k);

        p   = AWPoint::centerPoint(i,j,k);
            p.x = p.x * px;
            p.y = p.y * py;
            p.z = p.z * pz;

#ifndef USE_PGROUP
        onode  = ocTree.FindClosest(p,pC);
        if(!onode || !pC.isValid()) 
        continue;

            // find all polys within 1.5 mm of 'closest' one
            ocTree.FillDataCloseTo(onode,T,1.5f,1.5f,1.5f);
#else
            pGroup.findClosestPoints(p,T);
#endif

        tsum += T.size();
            tcnt++;

        mind = 1E20;
            mini = -1;

        for(g=0;g<T.size();g++) {
                f = T[g];

                if((f<0)||(degen[f]))
                   continue;

        a = S.facets()[f][0];
        b = S.facets()[f][1];
        c = S.facets()[f][2];

                // actually this is squared distance
                distanceToPoly(S,a,b,c,px,py,pz,p,dist,polyi);

        if(fabs(dist) < fabs(mind))  {
           mind = dist;
                   mini = polyi;
                   minf = f;
                }
            }

            if(mini >= 0) {
               sgn = (mind < 0.) ? -1. : 1.;
               mind = sqrt(fabs(mind));

               //
               // fix signs that are > 2mm and
               // have a valid seg tag
               //
               if((white_thresh != -1.)&&(mind > 2.f)&&(sval > 0)) {
                 if(sval >= white_thresh) {
                    if(sgn > 0) { sgn = -1; nbadsign++; }
                 } else {
                    if(sgn < 0) { sgn =  1; nbadsign++; }
                 }
               }

               mind *= sgn;

               hist.Add(mind);
               abshist.Add(fabs(mind));

               if(rDist)
              rDist->float_data()[k][j][i] = mind;
  
               if(rInd)
              rInd->int_data()[k][j][i] = mini;
            }
       }
    }

time_t end_time = time(NULL);

cerr << "Total time: " << (end_time-start_time) << " seconds" << endl;
cerr << "Avg time per voxel: " << (double)(end_time-start_time)/(double)(nx*ny*nz) << " sec" << endl;

cerr << "Avg # nodes searched: " << tsum/tcnt << endl;
cerr << "# mis-matched signs: " << nbadsign << endl;

 if(rDist) {
     OpMessage::GenericWorkProc("Fixing erroneous signs...");
   fixSigns(rDist->float_data(),px,py,pz);
 }

 return(true);
}

//
// The method of calculating inside/outside above has
// a couple of problems for pixels far from surface
// This ensures neighboring voxels have same sign if
// they are not close to surface
//
void SurfaceDistanceBase::fixSigns(Array3D<float> &A, float px, float py, float pz)
{
   Array3D<float> B(A);
   int nx = A.getXsize();
   int ny = A.getYsize();
   int nz = A.getZsize();
   int nfix=1,totcnt=0;
   float v[27];
   float mincut = 1. * sqrt(px*px + py*py + pz*pz);
   float maxcut = 2.f; // > 2mm handled by seg test

   while((nfix>0)&&(totcnt<100)) {
   nfix=0;
   totcnt++;
   for(int k=1;k<nz-1;k++)
     for(int j=1;j<ny-1;j++)
       for(int i=1;i<nx-1;i++) {
         v[0] = A[k][j][i];
         if((v[0] == InvalidPixel)||(fabs(v[0])<mincut)||(fabs(v[0])>=maxcut))
             continue;
  
         v[1] = A[k][j][i+1];
         v[2] = A[k][j][i-1];
         v[3] = A[k][j+1][i];
         v[4] = A[k][j-1][i];
         v[5] = A[k][j+1][i+1];
         v[6] = A[k][j+1][i-1];
         v[7] = A[k][j-1][i+1];
         v[8] = A[k][j-1][i-1];

         v[9]  = A[k+1][j][i];
         v[10] = A[k+1][j][i+1];
         v[11] = A[k+1][j][i-1];
         v[12] = A[k+1][j+1][i];
         v[13] = A[k+1][j-1][i];
         v[14] = A[k+1][j+1][i+1];
         v[15] = A[k+1][j+1][i-1];
         v[16] = A[k+1][j-1][i+1];
         v[17] = A[k+1][j-1][i-1];

         v[18] = A[k-1][j][i];
         v[19] = A[k-1][j][i+1];
         v[20] = A[k-1][j][i-1];
         v[21] = A[k-1][j+1][i];
         v[22] = A[k-1][j-1][i];
         v[23] = A[k-1][j+1][i+1];
         v[24] = A[k-1][j+1][i-1];
         v[25] = A[k-1][j-1][i+1];
         v[26] = A[k-1][j-1][i-1];

         int npos=0; int nneg=0;
         for(int l=1;l<27;l++) {
            if(v[l] == InvalidPixel) continue;
            if(v[l] < 0.f) nneg++; else npos++;
         }
 
         if((nneg==0)&&(npos==0)) continue;

         float dd = (float)(nneg-npos)/(float)(nneg+npos);

         // swap signs if > 60% neighbors are different
         if(((dd < -0.6)&&(v[0]<0.f))||((dd > 0.6)&&(v[0]>0.f))) {
             B[k][j][i] = -B[k][j][i];
             nfix++;
         }
       }
   cout << nfix << " signs changed" << endl;
   A = B;
   }
}

//
// Return closest distance (in mm) from point p to
// poly (a,b,c) in surface S
//
// miniret is vertex # of closest vertex to p
//
// px,py,pz is size (in mm) of a voxel
//
void SurfaceDistanceBase::distanceToPoly(Surface &S,
                            int a, int b, int c,
                            float px, float py, float pz,
                            AWPoint &p, double &distret, int &miniret)
{
  AWPoint Pt[3];

  if(!S.hasUNorms())
     S.genUnitNormals();

  double cdist,dist;

  //
  // Define poly in mm
  //
  Pt[0].x = S.vertices()[a].x() * px;
  Pt[0].y = S.vertices()[a].y() * py;
  Pt[0].z = S.vertices()[a].z() * pz;

  Pt[1].x = S.vertices()[b].x() * px;
  Pt[1].y = S.vertices()[b].y() * py;
  Pt[1].z = S.vertices()[b].z() * pz;

  Pt[2].x = S.vertices()[c].x() * px;
  Pt[2].y = S.vertices()[c].y() * py;
  Pt[2].z = S.vertices()[c].z() * pz;

  cdist = sqPolyDist(p,Pt);

  //
  // check dot sign for inside/outside test
  //
  float nrmxa = S.unormals()[a].x();
  float nrmya = S.unormals()[a].y();
  float nrmza = S.unormals()[a].z();
  float nrmxb = S.unormals()[b].x();
  float nrmyb = S.unormals()[b].y();
  float nrmzb = S.unormals()[b].z();
  float nrmxc = S.unormals()[c].x();
  float nrmyc = S.unormals()[c].y();
  float nrmzc = S.unormals()[c].z();

  double dxa,dya,dza;
  double dxb,dyb,dzb;
  double dxc,dyc,dzc;

  dxa = p.x - Pt[0].x;
  dya = p.y - Pt[0].y;
  dza = p.z - Pt[0].z;
  double dota = dxa*nrmxa + dya*nrmya + dza*nrmza;

  dxb = p.x - Pt[1].x;
  dyb = p.y - Pt[1].y;
  dzb = p.z - Pt[1].z;
  double dotb = dxb*nrmxb + dyb*nrmyb + dzb*nrmzb;

  dxc = p.x - Pt[2].x;
  dyc = p.y - Pt[2].y;
  dzc = p.z - Pt[2].z;
  double dotc = dxc*nrmxc + dyc*nrmyc + dzc*nrmzc;

  if((dota<0.)&&((dotb<0.)||(dotc<0.)))
     cdist = -cdist;
  else if((dotb<0.)&&(dotc<0.))
     cdist = -cdist;

  distret = cdist;

  double da = dxa*dxa + dya*dya + dza*dza;
  double db = dxb*dxb + dyb*dyb + dzb*dzb;
  double dc = dxc*dxc + dyc*dyc + dzc*dzc;
  if(da < db) {
    if(da < dc) miniret = a;
    else        miniret = c;
  } else {
    if(db < dc) miniret = b;
    else        miniret = c;
  }
}

double SurfaceDistanceBase::sqPolyDist(AWPoint &P, AWPoint Pol[3])
{
    AWPoint Pol2P = Pol[0] - P;
    AWPoint e0,e1;
    e0 = Pol[1] - Pol[0];
    e1 = Pol[2] - Pol[0];
    double fA00 = e0.squareNorm();
    double fA01 = e0.dot(e1);
    double fA11 = e1.squareNorm();
    double fB0  = Pol2P.dot(e0);
    double fB1  = Pol2P.dot(e1);
    double fC   = Pol2P.squareNorm();
    double fDet = fabs(fA00*fA11-fA01*fA01);
    double fS   = fA01*fB1-fA11*fB0;
    double fT   = fA01*fB0-fA00*fB1;

    double fSqrDist;

    if(fS + fT <= fDet) {
        if(fS < 0.0f) {
            if(fT < 0.0f) {
                if(fB0 < 0.0f) {
                    fT = 0.0f;
                    if( -fB0 >= fA00) {
                        fS = 1.0f;
                        fSqrDist = fA00+2.0f*fB0+fC;
                    }
                    else {
                        fS = -fB0/fA00;
                        fSqrDist = fB0*fS+fC;
                    }
                }
                else {
                    fS = 0.0f;
                    if(fB1 >= 0.0f) {
                        fT = 0.0f;
                        fSqrDist = fC;
                    }
                    else if(-fB1 >= fA11) {
                        fT = 1.0f;
                        fSqrDist = fA11+2.0f*fB1+fC;
                    }
                    else {
                        fT = -fB1/fA11;
                        fSqrDist = fB1*fT+fC;
                    }
                }
            }
            else {
                fS = 0.0f;
                if(fB1 >= 0.0f) {
                    fT = 0.0f;
                    fSqrDist = fC;
                }
                else if(-fB1 >= fA11) {
                    fT = 1.0f;
                    fSqrDist = fA11+2.0f*fB1+fC;
                }
                else {
                    fT = -fB1/fA11;
                    fSqrDist = fB1*fT+fC;
                }
            }
        }
        else if(fT < 0.0f) {
            fT = 0.0f;
            if(fB0 >= 0.0f) {
                fS = 0.0f;
                fSqrDist = fC;
            }
            else if(-fB0 >= fA00) {
                fS = 1.0f;
                fSqrDist = fA00+2.0f*fB0+fC;
            }
            else {
                fS = -fB0/fA00;
                fSqrDist = fB0*fS+fC;
            }
        }
        else {
            // inside poly
            double fInvDet = 1.0f/fDet;
            fS *= fInvDet;
            fT *= fInvDet;
            fSqrDist = fS*(fA00*fS+fA01*fT+2.0f*fB0) +
                fT*(fA01*fS+fA11*fT+2.0f*fB1)+fC;
        }
    }
    else {
        double fTmp0, fTmp1, fNumer, fDenom;

        if(fS < 0.0f) {
            fTmp0 = fA01 + fB0;
            fTmp1 = fA11 + fB1;
            if(fTmp1 > fTmp0) {
                fNumer = fTmp1 - fTmp0;
                fDenom = fA00-2.0f*fA01+fA11;
                if(fNumer >= fDenom) {
                    fS = 1.0f;
                    fT = 0.0f;
                    fSqrDist = fA00+2.0f*fB0+fC;
                }
                else {
                    fS = fNumer/fDenom;
                    fT = 1.0f - fS;
                    fSqrDist = fS*(fA00*fS+fA01*fT+2.0f*fB0) +
                        fT*(fA01*fS+fA11*fT+2.0f*fB1)+fC;
                }
            }
            else {
                fS = 0.0f;
                if(fTmp1 <= 0.0f) {
                    fT = 1.0f;
                    fSqrDist = fA11+2.0f*fB1+fC;
                }
                else if(fB1 >= 0.0f) {
                    fT = 0.0f;
                    fSqrDist = fC;
                }
                else {
                    fT = -fB1/fA11;
                    fSqrDist = fB1*fT+fC;
                }
            }
        }
        else if(fT < 0.0f) {
            fTmp0 = fA01 + fB1;
            fTmp1 = fA00 + fB0;
            if (fTmp1 > fTmp0) {
                fNumer = fTmp1 - fTmp0;
                fDenom = fA00-2.0f*fA01+fA11;
                if(fNumer >= fDenom) {
                    fT = 1.0f;
                    fS = 0.0f;
                    fSqrDist = fA11+2.0f*fB1+fC;
                }
                else {
                    fT = fNumer/fDenom;
                    fS = 1.0f - fT;
                    fSqrDist = fS*(fA00*fS+fA01*fT+2.0f*fB0) +
                        fT*(fA01*fS+fA11*fT+2.0f*fB1)+fC;
                }
            }
            else {
                fT = 0.0f;
                if(fTmp1 <= 0.0f) {
                    fS = 1.0f;
                    fSqrDist = fA00+2.0f*fB0+fC;
                }
                else if(fB0 >= 0.0f) {
                    fS = 0.0f;
                    fSqrDist = fC;
                }
                else {
                    fS = -fB0/fA00;
                    fSqrDist = fB0*fS+fC;
                }
            }
        }
        else {
            fNumer = fA11 + fB1 - fA01 - fB0;
            if(fNumer <= 0.0f) {
                fS = 0.0f;
                fT = 1.0f;
                fSqrDist = fA11+2.0f*fB1+fC;
            }
            else {
                fDenom = fA00-2.0f*fA01+fA11;
                if(fNumer >= fDenom) {
                    fS = 1.0f;
                    fT = 0.0f;
                    fSqrDist = fA00+2.0f*fB0+fC;
                }
                else {
                    fS = fNumer/fDenom;
                    fT = 1.0f - fS;
                    fSqrDist = fS*(fA00*fS+fA01*fT+2.0f*fB0) +
                        fT*(fA01*fS+fA11*fT+2.0f*fB1)+fC;
                }
            }
        }
    }
    return(fSqrDist);
}



