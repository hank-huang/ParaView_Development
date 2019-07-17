#include <stdlib.h>
#include <iostream.h>
#include <math.h>
#include <KDApplication.h>
#include <Main.h>
#include <MainUtils.h>
#include <AWHistogram.h>
#include <TDL/Surface.h>
#include <ADL/Array2D.h>
#include <SurfaceDistance.h>
#include <SurfaceThickness.h>
#include <SurfPatches.h>

#define MAX_VOXEL_DIST	50

static bool SmoothFunction(Surface *S, Array1D<double> &D)
{
  int i,j,num,n,npol;
  npol = S->getNumVert();

  if(!S->hasNbhd())
      S->genNeighborhoods();

  if(!S->hasNbhd()) {
      cerr << "ERROR: SmoothFunction: could not calc neighborhoods" << endl;
      return(false);
  }

  bool gotone = true;

  // push function values to non-valued neighbors

  int iter=0;
  while(gotone) {
    if(++iter > 1000) {
       cerr << "ERROR: SmoothFunction : too many ";
       cerr << "iterations (bad nghbds?)" << endl;
       return(false);
    }
    for(i=0;i<npol;i++) {
      if(D[i] <= 0.)
         continue;
      num = S->neighborhood()[i].numNeighbors();
      for(j=0;j<num;j++) {
         n = S->neighborhood()[i].getNeighbor(j);
         if(D[n] <= 0.) {
            D[n] = D[i];
            gotone = true;
         }
      }
    }
  }
  return(true);
}


bool SurfaceThickness::calculate(
                        Surface   *S, 
                        ItXVolume &Seg, 
                        ItXVolume *dist, 
                        ItXVolume *ind, 
                        float patch_sqarea,
                        float max_dist,
                        float grey_value_min,
                        float grey_value_max,
                        float perc_grey,
                        Array1D<double> &thickness)
{
  int nx,ny,nz,i,j,k,p;
  int npol,nvert,npatch;
  double th,t;
  float px,py,pz;
  float fmin,fmax,frange,fthresh;
  float xmin,xmax;
  float ymin,ymax;
  float zmin,zmax;
  float     d,sval;
  int     idx,totcnt;
  int     targcnt,b;
  char    fname[1024];
  char    fpref[1024];
  char    estr[1024];
  AWHistogram hist;
  SurfPatches patches;
  bool hp;

  nx = Seg.getSizeX();
  ny = Seg.getSizeY();
  nz = Seg.getSizeZ();

  px = Seg.getPixelDimensionX();
  py = Seg.getPixelDimensionY();
  pz = Seg.getPixelDimensionZ();

  if(!S) {
    cerr << "ERROR: SurfaceThickness:: invalid surface" << endl;
    return(false);
  }

  if(!dist || !ind) {
    cerr << "ERROR: SurfaceThickness:: invalid distances/indices" << endl;
    return(false);
  }

  npol  = S->getNumPoly();
  nvert = S->getNumVert();

  Array1D<u_char> isPatch(nvert);

  if((npol<1)||(nvert<1)) {
    cerr << "ERROR: SurfaceThickness:: invalid surface" << endl;
    return(false);
  }

  if((dist->getSizeX() != nx)||(dist->getSizeY() != ny)||
     (dist->getSizeZ() != nz)||(ind->getSizeX()  != nx)||
     (ind->getSizeY()  != ny)||(ind->getSizeZ()  != nz)) {
       cerr << "ERROR: SurfaceThickness:: invalid ";
       cerr << "input DISTANCES/INDICES vol(s)"<< endl;
       cerr << "       (size mismatch)" << endl;

cerr << nx << " " << ny << " " << nz << endl;
cerr << dist->getSizeX() << " " << dist->getSizeY() << " " << dist->getSizeZ() << endl;
cerr << ind->getSizeX() << " " << ind->getSizeY() << " " << ind->getSizeZ() << endl;
       return(false);
  }

  // TEST
  if(dist->dataType() != ItXVolume::Float) {
     cerr << "ERROR: SurfaceThickness:: assuming FLOAT Distances. Recompile"<< endl;
     goto bad_thickness;
  }
  if(ind->dataType() != ItXVolume::Int) {
     cerr << "ERROR: SurfaceThickness:: assuming INT indices. Recompile"<< endl;
     goto bad_thickness;
  }

  hp = false;
  if(App->GetYesNo("Do you have an existing PATCH file?"))  {
retry_patch_load:
    fname[0] = 0;
    if(Browser->GetFilename("Select Patches File",fname,1024,"*.txt")) {
     hp = patches.load(fname);
     if(!hp) {
        App->ShowMessage("Patches file did not load. Try again");   
        goto retry_patch_load;
     }
    }
  }

  if(!hp && !patches.calculate(*S,patch_sqarea,px,py,pz)) {
     cerr << "ERROR: could not calculate surface patches" << endl;
     goto bad_thickness;
  }

  npatch = patches.numPatches();
  if(npatch <= 0) {
     cerr << "ERROR: could not calculate surface patches" << endl;
     goto bad_thickness;
  }

  cerr << npatch << " patches created" << endl;

  fpref[0] = 0;
  if(!Browser->GetFilename("Select PREFIX for patch histograms",
                           fpref,1024))
     goto bad_thickness;

  thickness.setDim(nvert);
  thickness = -1.; // <0 means tag invalid dist

  dist->getMinMax(fmin,fmax);
  if(fmax > max_dist)  fmax =  max_dist;
  if(fmin < -max_dist) fmin = -max_dist;
  frange = (fmax - fmin)/100.;
  if(frange == 0.) frange = 1.;

  hist.SetSize(frange,fmin,fmax);

  for(p=0;p<npatch;p++) {
      sprintf(estr,"Calculating Patch Thickness %04d/%04d",p+1,npatch);
      GenericWorkProc(estr);

      hist.reset();
      isPatch = 0U;
      for(i=0;i<patches.getPatch(p).verts.size();i++)
          isPatch[patches.getPatch(p).verts[i]] = 1U;

      totcnt=0;
      for(k=0;k<nz;k++)
        for(j=0;j<ny;j++)
          for(i=0;i<nx;i++) {
             sval = Seg.pointValue(i,j,k);
             if((sval<grey_value_min)||(sval>grey_value_max))
                 continue;

             d=dist->float_data()[k][j][i];
             if((d == SurfaceDistance::InvalidPixel)||(fabs(d) > max_dist))
                continue;

             idx = ind->int_data()[k][j][i];
 
             //if(!patches.getPatch(p).verts.exists(idx))
             //    continue;
             if((idx>=0)&&(idx<nvert)&&(isPatch[idx])) {
               hist.Add((double)d);
               totcnt++;
             }
          }

       if(totcnt == 0) {
          cerr << "*** WARNING: patch " << p;
          cerr << " is not close to ANY GM voxels." << endl;
          cerr << "             Cannot calculate distance for that patch"<<endl;
       } else {

          targcnt = (int)((float)totcnt * perc_grey + 0.5);
          for(t=0.,b=0;b<hist.Bins();b++) {
              t += hist.ValAt(b); // bin count
              if(t >= targcnt) break;
          }
          if(b >= hist.Bins())
             b = hist.Bins()-1;

          th = hist.BinAt(b);

          // set thickness for those verts in patch

          for(b=0;b<patches.getPatch(p).verts.size();b++)
              thickness[patches.getPatch(p).verts[b]] = th;

          sprintf(fname,"%s_%04d.txt",fpref,p);
          if(!hist.Save(fname))
             cerr << "WARNING: Histogram "<< p <<" did not save"<<endl;
       }
  }
  GenericWorkProc(NULL);

/*
  if(!SmoothFunction(S,thickness)) 
     goto bad_thickness;
*/
  return(true);

bad_thickness:

  return(false);
}





