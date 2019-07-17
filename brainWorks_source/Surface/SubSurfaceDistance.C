#include <stdlib.h>
#include <iostream.h>
#include <math.h>
#include <OpMessage.h>
#include <TDL/ItXVolume.h>
#include <TDL/Surface.h>
#include <ADL/Array2D.h>
#include <OcTree.h>
#include <SubSurfaceDistance.h>
#include <SurfaceDistanceBase.h>

#define MAX_POL_PER_VERT	10

//
// (Exactly like SurfaceDistance.C except masks out
//  voxels who arent close to sub-surface `subS')
//
bool SubSurfaceDistance::calculate(Surface &S, Surface &subS,
	ItXVolume &V, ItXVolume &rV)
{
   char tstr[1024];
   int i,j,k,idx,npt,npol,n,cnt,polyi;
   int a,b,c,f,g,mini,nev,minf;
   float x,y,z,subx,suby,subz;
   u_char ch,sval;
   AWPoint p,pC,Pt[3];
   char estr[100];
   double ptdist,mind,dist,dx,dy,dz,sgn,nrmx,nrmy,nrmz,dot;
   Array1D<u_char> isSub;
   DynArray<int> T;
 
   if(!&S || !&subS || !&V || !&rV)
      return(false);
  
   OpMessage::GenericWorkProc("Calculating...");

   OcTree ocTree;
   OcTreeNode *onode;

   int nx = V.getSizeX();
   int ny = V.getSizeY();
   int nz = V.getSizeZ();

   float px = V.getPixelDimensionX();
   float py = V.getPixelDimensionY();
   float pz = V.getPixelDimensionZ();

   rV.setSize(nx,ny,nz,ItXVolume::Float);
   rV.setPixelDimensionX(px);
   rV.setPixelDimensionY(py);
   rV.setPixelDimensionZ(pz);
  
   if(!S.hasUNorms())
	S.genUnitNormals();

   ocTree.clear();
   ocTree.setSize((float)nx*px+1.f,(float)ny*py+1.f, (float)nz*pz+1.f,3);
   
   npt  = S.getNumVert();
   npol = S.getNumPoly();

   isSub.setDim(npt);
   isSub = 0U;
   cnt=0;
   nev = npt/100;
   for(i=0;i<npt;i++) {
      if((i%nev)==0)
      {
         OpMessage::GenericWorkProc("Calculating coincident vertices...");
      }
      x = S.vertices()[i].x();
      y = S.vertices()[i].y();
      z = S.vertices()[i].z();
      for(j=0;j<subS.getNumVert();j++) {
         subx = subS.vertices()[j].x();
         suby = subS.vertices()[j].y();
         subz = subS.vertices()[j].z();
         if((subx==x)&&(suby==y)&&(subz==z)) {
           isSub[i] = 1U;
           cnt++;
           break;
         }
      }
   }
     
   if(cnt == 0) {
      cerr << "ERROR: SubSurfaceDistance : surfaces do not overlap" << endl;
      rV.float_data() = SurfaceDistanceBase::InvalidPixel;
      return(false);
   }

   cerr << "SubSurfaceDistance : " << cnt << " vertices overlap" << endl;

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
      ocTree.AddPoint(p);
   }

   for(k=0;k<nz;k++)  {
     sprintf(tstr,"Calculating Sub-Surface Map: plane %03d/%03d",k,nz);
     OpMessage::GenericWorkProc(tstr);
     for(j=0;j<ny;j++)
       for(i=0;i<nx;i++) {

            rV.float_data()[k][j][i] = SurfaceDistanceBase::InvalidPixel;

	    if(V.pointValue(i,j,k) == 0.f)
               continue;

   	    p   = AWPoint::centerPoint(i,j,k);
            p.x = p.x * px;
            p.y = p.y * py;
            p.z = p.z * pz;

	    onode = ocTree.FindClosest(p,pC);
	
	    if(!onode || !pC.isValid()) 
		continue;

            // find all polys within 1.5 mm of 'closest' one
            ocTree.FillDataCloseTo(onode,T,1.5f,1.5f,1.5f);

            mind = 1E20;
            mini = -1;
            for(g=0;g<T.size();g++) {
                f = T[g];

                a = S.facets()[f][0];
                b = S.facets()[f][1];
                c = S.facets()[f][2];

                // actually this is squared distance
                SurfaceDistanceBase::distanceToPoly(S,a,b,c,
                                   px,py,pz,p,dist,polyi);
                if(fabs(dist) < fabs(mind))  {
                   mind = dist;
                   mini = polyi;
                   minf = f;
                }
            }

            sgn = (mind < 0.) ? -1. : 1.;
            mind = sgn * sqrt(fabs(mind));

            // only record distance if closest vertex is in sub surface
            if((mini>=0) && isSub[mini])
	       rV.float_data()[k][j][i] = mind;
       }
    }

    OpMessage::GenericWorkProc("Fixing erroneous signs...");
    SurfaceDistanceBase::fixSigns(rV.float_data(),px,py,pz);

    OpMessage::GenericWorkProc(NULL);
 return(1);

}

