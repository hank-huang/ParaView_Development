#include <stdlib.h>
#include <iostream.h>
#include <math.h>
#include <OpMessage.h>
#include <TDL/Surface.h>
#include <ADL/Array2D.h>
#include <OcTree.h>
#include <SurfaceDistanceBase.h>
#include <SurfaceSurfaceDistance.h>
#include <DynArray.h>

#define MAX_POL_PER_VERT	10

bool SurfaceSurfaceDistance::calculate(Surface &S1,
	Surface &S2, double maxdist, Array1D<u_char> &isClose)
{
   char tstr[1024];
   int i,j,k,idx,n;
   int a,b,c,f,g,polyi;
   float x,y,z;
   u_char ch,sval;
   AWPoint p,pC;
   char estr[100];
   double sgn,ptdist,mind,dist,dx,dy,dz,nrmx,nrmy,nrmz,dot;
   bool dodist;

   DynArray<int> T;
 
   if(!&S1 || !&S2 || !&isClose)
      return(false);
  
   OcTree ocTree;
   OcTreeNode *onode;

   int npt1  = S1.getNumVert();
   int npol1 = S1.getNumPoly();

   int npt2  = S2.getNumVert();
   int npol2 = S2.getNumPoly();

   isClose.setDim(npt1);

   if(!S2.hasUNorms())
	S2.genUnitNormals();

   // Put surface S2 verts into OcTree

   float nx1 = (float)S1.findXMin();
   float mx1 = (float)S1.findXMax();
   float ny1 = (float)S1.findYMin();
   float my1 = (float)S1.findYMax();
   float nz1 = (float)S1.findZMin();
   float mz1 = (float)S1.findZMax();

   float nx2 = (float)S2.findXMin();
   float mx2 = (float)S2.findXMax();
   float ny2 = (float)S2.findYMin();
   float my2 = (float)S2.findYMax();
   float nz2 = (float)S2.findZMin();
   float mz2 = (float)S2.findZMax();

   float mx = (mx1 < mx2) ? mx2 : mx1;
   float my = (my1 < my2) ? my2 : my1;
   float mz = (mz1 < mz2) ? mz2 : mz1;

//   ocTree.setSize(nx,ny,nz,mx,my,mz,3);

//   ocTree.setSize(nx2,ny2,nz2,mx2,my2,mz2,3);
   ocTree.setSize(mx+1.f,my+1.f,mz+1.f,3);
   
   //
   // Store centroid of each poly in OcTree
   //
   for(i=0;i<npol2;i++) {
      a  = S2.facets()[i][0];
      b  = S2.facets()[i][1];
      c  = S2.facets()[i][2];
      x  = S2.vertices()[a].x();
      y  = S2.vertices()[a].y();
      z  = S2.vertices()[a].z();
      x += S2.vertices()[b].x();
      y += S2.vertices()[b].y();
      z += S2.vertices()[b].z();
      x += S2.vertices()[c].x();
      y += S2.vertices()[c].y();
      z += S2.vertices()[c].z();
      x /= 3.f;
      y /= 3.f;
      z /= 3.f;
      p.set(x,y,z,i);
      p.setValid(true);
      ocTree.AddPoint(p);
   }

#define KMIN(a,b) (((a)<(b))?(a):(b))
#define KMAX(a,b) (((a)<(b))?(b):(a))

   int ndegen=0;

   Array1D<u_char> degen(npol2);
   degen = 0U;
   for(i=0;i<npol2;i++) {
      a = S2.facets()[i][0];
      b = S2.facets()[i][1];
      c = S2.facets()[i][2];
      Point e1 = S2.vertices()[a] - S2.vertices()[b];
      Point e2 = S2.vertices()[b] - S2.vertices()[c];
      Point e3 = S2.vertices()[c] - S2.vertices()[a];
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

   for(k=0;k<npt1;k++)
   {
        if((k%100)==0)
        {
            sprintf(tstr,"Calculating : vertex %03d/%03d",k,npt1);
            OpMessage::GenericWorkProc(tstr);
        }

     p.x = S1.vertices()[k].x();
     p.y = S1.vertices()[k].y();
     p.z = S1.vertices()[k].z();

     onode  = ocTree.FindClosest(p,pC);
	
     if(!onode)  {
cerr << "*";
	continue;
     }
     if(!pC.isValid())  {
cerr << "^";
	continue;
     }

     // find all polys within 1.5 units of 'closest' one
     ocTree.FillDataCloseTo(onode,T,1.5f,1.5f,1.5f);

     mind = 1E20;
     bool got = false;
     for(g=0;g<T.size();g++) {
        f = T[g];
        if(degen[f])
            continue;

	a = S2.facets()[f][0];
	b = S2.facets()[f][1];
	c = S2.facets()[f][2];

        // actually this is squared distance
        SurfaceDistanceBase::distanceToPoly(S2,a,b,c,
                                        1.f,1.f,1.f,p,dist,polyi);

        if(fabs(dist) < fabs(mind))
	   mind = dist;
        }

        mind = sqrt(fabs(mind));
        if(mind <= maxdist)
          isClose[k] = 1U;
        else
          isClose[k] = 0U;
    }
    
    OpMessage::GenericWorkProc(NULL);

 return(true);
}

