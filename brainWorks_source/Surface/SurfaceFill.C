#include <SurfaceFill.h>
#include <Lib/TDL/ItXVolume.h>
#include <Lib/TDL/Surface.h>
#include <SurfaceEmbed.h>
#include <NRFlood.h>

static void flood1D(int y, int z, Array3D<unsigned char>  &D,
		unsigned char fill, unsigned char bound,
		int x1, int x2, bool isin)
{
   int i,l,r;
   unsigned char lv,rv;

   if(x1 == x2) {
      if(isin) D[z][y][x1] = fill;
      return;
   }

   // fill from left
   lv = D[z][y][x1];
   while((x1<=x2)&&(lv != bound)) {
      if(isin) D[z][y][x1] = fill;
      x1++;
      lv = D[z][y][x1];
   }
   if(x1 >= x2) return;
   // pass through boundary
   while((x1<=x2)&&(lv == bound)) {
      x1++;
      lv = D[z][y][x1];
   }

   // fill from right
   rv = D[z][y][x2];
   while((x1<=x2)&&(rv != bound)) {
      if(isin) D[z][y][x2] = fill;
      x2--;
      rv = D[z][y][x2];
   }
   if(x1 >= x2) return;
   // pass through boundary
   while((x1<=x2)&&(rv == bound)) {
      x2--;
      rv = D[z][y][x2];
   }

   flood1D(y,z,D,fill,bound,x1,x2,!isin);
}

static void flood(Array3D<unsigned char> &D, 
		  unsigned char fill, unsigned char bound)
{
   int i,j,k;
   for(k=0;k<D.getZsize();k++)
     for(j=0;j<D.getYsize();j++)
	flood1D(j,k,D,fill,bound,0,D.getXsize()-1,false);
}



bool SurfaceFill::fillSurface(ItXVolume &V, Surface &S, unsigned char fillval)
{
	int i,j,k;
	Array1D<u_char> edges;

	if(V.dataType() != ItXVolume::UnsignedChar) {
	    cerr << "ERROR: fillSurface() only UCHAR volumes allowed" << endl;
	    return(false);
	}
	//
	// Make sure surface is closed (== no edge vertices)
	//
	S.isEdge(edges);
	for(i=0;i<edges.getSize();i++)
		if(edges[i]) 
		   return(false);

        V.u_char_data() = 0U;
	SurfaceEmbed::embedSurface(V,S,250.f);

	flood(V.u_char_data(),150U,250U);
        for(k=0;k<V.getSizeZ();k++)
          for(j=0;j<V.getSizeY();j++)
            for(i=0;i<V.getSizeX();i++)
                 if(V.u_char_data()[k][j][i])
		    V.u_char_data()[k][j][i] = fillval;

	return(true);	
}






