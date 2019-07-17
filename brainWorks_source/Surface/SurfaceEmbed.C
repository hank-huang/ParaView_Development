//////////////////////////////////////////////////////////////////////////
//
// File: SurfaceBox.C
//
// Author: Keith Doolittle
//
// Purpose:
//
///////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include <OpMessage.h>
#include <AWPoint.h>
#include <TDL/Surface.h>
#include <TDL/ItXVolume.h>
#include <SurfaceEmbed.h>

#define STEPS_PER_VOXEL 1.1


void SurfaceEmbed::embedSurface(ItXVolume &V, Surface &S, float val)
{
  if(!&V || !&S)
     return;

  char tstr[1024];
  int i;

  for(i=0;i<S.getNumPoly();i++) {
        if((i%1000)==999)
           OpMessage::GenericWorkProc("Embedding...");
        
        int a = S.facets()[i][0];
        int b = S.facets()[i][1];
        int c = S.facets()[i][2];

        float ax = S.vertices()[a].x();
        float ay = S.vertices()[a].y();
        float az = S.vertices()[a].z();
        float bx = S.vertices()[b].x();
        float by = S.vertices()[b].y();
        float bz = S.vertices()[b].z();
        float cx = S.vertices()[c].x();
        float cy = S.vertices()[c].y();
        float cz = S.vertices()[c].z();

        AWPoint pa(ax,ay,az);
        AWPoint pb(bx,by,bz);
        AWPoint pc(cx,cy,cz);

        scanTriangle(V,val,pa,pb,pc);
   }
   OpMessage::GenericWorkProc(NULL);
}


void SurfaceEmbed::embedSurface(ItXVolume &V, Surface &S, Array1D<double> &vals)
{
  if(!&V || !&S)
     return;

  char tstr[1024];
  int i,j,k,nx,ny,nz;
  float f,nf;

  nx = V.getSizeX();
  ny = V.getSizeY();
  nz = V.getSizeZ();

  Array3D<float> VV(nz,ny,nx);
  Array3D<float> sumV(nz,ny,nx);
  VV   = 0.f;
  sumV = 0.f;

  for(i=0;i<S.getNumPoly();i++) {
        if((i%1000)==999)
           OpMessage::GenericWorkProc("Embedding...");
        
        int a = S.facets()[i][0];
        int b = S.facets()[i][1];
        int c = S.facets()[i][2];

        float ax = S.vertices()[a].x();
        float ay = S.vertices()[a].y();
        float az = S.vertices()[a].z();
        float bx = S.vertices()[b].x();
        float by = S.vertices()[b].y();
        float bz = S.vertices()[b].z();
        float cx = S.vertices()[c].x();
        float cy = S.vertices()[c].y();
        float cz = S.vertices()[c].z();

        AWPoint pa(ax,ay,az);
        AWPoint pb(bx,by,bz);
        AWPoint pc(cx,cy,cz);

        double v1 = vals[a];
        double v2 = vals[b];
        double v3 = vals[c];

        scanTriangle(VV,sumV,v1,v2,v3,pa,pb,pc);
   }
   OpMessage::GenericWorkProc(NULL);

   for(k=0;k<nz;k++)
     for(j=0;j<ny;j++)
       for(i=0;i<nx;i++) {
           nf = sumV[k][j][i];
           if(nf == 0.f)
                f = 0.f;
           else f = (float)(VV[k][j][i]/nf);
           V.setPointValue(i,j,k,f);
       }
}


void SurfaceEmbed::embedSurface(ItXVolume &V, Surface &S, Array2D<double> &vals)
{
  if(!&V || !&S || (vals.getNcol() != 3))
     return;

  char tstr[1024];
  int i,j,k,nx,ny,nz;
  float fr,fg,fb,nf;

  nx = V.getSizeX();
  ny = V.getSizeY();
  nz = V.getSizeZ();

  Array3D<float> VVr(nz,ny,nx);
  Array3D<float> VVg(nz,ny,nx);
  Array3D<float> VVb(nz,ny,nx);
  Array3D<float> sumV(nz,ny,nx);
  VVr  = 0.f;
  VVg  = 0.f;
  VVb  = 0.f;
  sumV = 0.f;

  for(i=0;i<S.getNumPoly();i++) {
        if((i%1000)==999)
           OpMessage::GenericWorkProc("Embedding...");
        
        int a = S.facets()[i][0];
        int b = S.facets()[i][1];
        int c = S.facets()[i][2];

        float ax = S.vertices()[a].x();
        float ay = S.vertices()[a].y();
        float az = S.vertices()[a].z();
        float bx = S.vertices()[b].x();
        float by = S.vertices()[b].y();
        float bz = S.vertices()[b].z();
        float cx = S.vertices()[c].x();
        float cy = S.vertices()[c].y();
        float cz = S.vertices()[c].z();

        AWPoint pa(ax,ay,az);
        AWPoint pb(bx,by,bz);
        AWPoint pc(cx,cy,cz);

        double v1 = vals[a][0];
        double v2 = vals[b][0];
        double v3 = vals[c][0];

        scanTriangle(VVr,sumV,v1,v2,v3,pa,pb,pc);

        v1 = vals[a][1];
        v2 = vals[b][1];
        v3 = vals[c][1];

        scanTriangle(VVg,sumV,v1,v2,v3,pa,pb,pc);

        v1 = vals[a][2];
        v2 = vals[b][2];
        v3 = vals[c][2];

        scanTriangle(VVb,sumV,v1,v2,v3,pa,pb,pc);
   }
   OpMessage::GenericWorkProc(NULL);

   for(k=0;k<nz;k++)
     for(j=0;j<ny;j++)
       for(i=0;i<nx;i++) {
           nf = sumV[k][j][i]/3.f;
           if(nf == 0.f) {
             fr = 0.f;
             fg = 0.f;
             fb = 0.f;
           } else {
             fr = (float)(VVr[k][j][i]/nf);
             fg = (float)(VVg[k][j][i]/nf);
             fb = (float)(VVb[k][j][i]/nf);
           }
           V.setPointValue(i,j,k,fr,0);
           V.setPointValue(i,j,k,fg,1);
           V.setPointValue(i,j,k,fb,2);
       }
}



void SurfaceEmbed::scanTriangle(ItXVolume &V, float val,
                   AWPoint &a, AWPoint &b, AWPoint &p2)
{
   AWPoint ab,p1,p3,p4;
   int     i,j,nlines,npix,ix,iy,iz;
   int     nx,ny,nz;
   double  d;

   nx = V.getSizeX();
   ny = V.getSizeY();
   nz = V.getSizeZ();

   ab     = a - b;
   nlines = (int)(STEPS_PER_VOXEL * ab.vectNorm() + 1);
   nlines *= 2;

   for(i=0;i<nlines;i++) {
        d  = (double)i/(double)nlines;
        p1 = b + ab * d;

        p3 = p2 - p1;
        npix = (int)(STEPS_PER_VOXEL * p3.vectNorm() + 1);
        for(j=0;j<npix;j++) {
                d = (double)j/(double)npix;
                p4 = p1 + p3 * d;
                ix = (int)(p4.x + 0.5);
                iy = (int)(p4.y + 0.5);
                iz = (int)(p4.z + 0.5);
                if((ix>=0)&&(ix<nx)&&(iy>=0)&&(iy<ny)&&(iz>=0)&&(iz<nz))
                        V.setPointValue(ix,iy,iz,val);
        }
   }
}



void SurfaceEmbed::scanTriangle(Array3D<float> &V, Array3D<float> &sumV,
                   double vala, double valb, double valc,
                   AWPoint &a, AWPoint &b, AWPoint &c)
{
   AWPoint ab,p1,p3;
   int     i,j,nlines,npix,ix,iy,iz;
   int     nx,ny,nz;
   double  d,uu,vv,vw,uv,uw,dd,s,t,tu,tv,tw;
   AWPoint u,v,w,Q;
   float   val;

   // find parametric form of internal point Q then
   // apply (s,t) to vala,valb,valc to
   // get interpolated value

   u = b - a;
   v = c - a;

   uu = u.dot(u);
   vv = v.dot(v);
   uv = u.dot(v);

   tu = valb - vala;
   tv = valc - vala;

   nx = V.getXsize();
   ny = V.getYsize();
   nz = V.getZsize();

   ab     = a - b;
   nlines = (int)(STEPS_PER_VOXEL * ab.vectNorm() + 1);
   nlines *= 2;

   for(i=0;i<nlines;i++) {
        d  = (double)i/(double)nlines;
        p1 = b + ab * d;

        p3 = c - p1;
        npix = (int)(STEPS_PER_VOXEL * p3.vectNorm() + 1);
        for(j=0;j<npix;j++) {
                d = (double)j/(double)npix;
                Q = p1 + p3 * d;

                ix = (int)(Q.x + 0.5);
                iy = (int)(Q.y + 0.5);
                iz = (int)(Q.z + 0.5);
                if((ix>=0)&&(ix<nx)&&(iy>=0)&&(iy<ny)&&(iz>=0)&&(iz<nz)) {
                        w = Q - a; 
                        uw = u.dot(w);
                        vw = v.dot(w);
                        dd = (uv*uv) - (uu*vv);
                        if(dd == 0.)
                           val = vala;
                        else {
                           s = ((uv*vw) - (vv*uw))/dd;
                           t = ((uv*uw) - (uu*vw))/dd;
                           tw = tu * s + tv * t;
                           val = vala + tw;
                        }
                        V[iz][iy][ix]    += val;
                        sumV[iz][iy][ix] += 1.;
                }
        }
   }
}




