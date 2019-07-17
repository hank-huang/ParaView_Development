///////////////////////////////////////////////////////////////////////////
//
// File: FlatMap.C
//
// Author: Keith Doolittle/Monica Hurdal/Michael Bowers
//
// Purpose:  Provide some support for flat map surfaces
//
///////////////////////////////////////////////////////////////////////////

#include <FlatMap.h>
#include <TDL/ByuSurface.h>

#define EUCLIDIAN       0
#define HGEOM           1
#define SPHERICAL       2

#define BFLATEXIT(a) { sprintf(statusMsg,"%s",a); goto bflat_fail; }

ByuSurface *FlatMap::LoadFlatMap(const char *filename, flatvertex **vertices, char* statusMsg)
{
  FILE *fp;
  char ch,line[512],tmpstr[512];
  int  done,gotcents;
  int  nv,npol,geo,tmp,i,j,nf,nw;
  *(vertices)=NULL;
  ByuSurface *bsurf=NULL;
  double dx,dy,dz;
  Point  pt;
  float tx,ty,x,y,z;
  flatvertex *verts = NULL;

  fp = fopen(filename,"r");
  if(!fp) BFLATEXIT("Cannot open file")

  nv = -1;

  /* get nodecount, geometry, and FLOWERS tag */
  done = 0;
  while(!done) {
        if(feof(fp)) BFLATEXIT("Failed to find FLOWERS tag")
        line[0] = 0;
        fgets(line,1023,fp);

        line[strlen(line)-1] = 0;
        if((line[0] == 0)||(strlen(line) < 1))
                continue;
        if(sscanf(line,"NODECOUNT: %d",&tmp) == 1)
                nv = tmp;
        else if(sscanf(line,"GEOMETRY: %s",tmpstr) == 1) {
                ch = tmpstr[0];
                if((ch=='S')||(ch=='s'))      geo = SPHERICAL;
                else if((ch=='H')||(ch=='h')) geo = HGEOM;
                else                          geo = EUCLIDIAN;
                }
        else if(!strncmp(line,"FLOWERS:",8))
                done = 1;
  }

  if(nv <= 0) BFLATEXIT("failed to find NODECOUNT tag")

  *vertices = new flatvertex[nv];
  verts = *vertices;
  if(!verts) BFLATEXIT("out of memory")

  npol = 0;
  for(i=0;i<nv;i++) {
      if(feof(fp)) 
    BFLATEXIT("invalid FLOWERS section")
      if(fscanf(fp,"%d %d",&tmp,&nf) != 2) { i--; continue; }
      if(nf >= MAX_FLOWER) 
    BFLATEXIT("recompile for MAX_FLOWER > 20")
      verts[i].num = nf;
      verts[i].color = 0;
      for(j=0;j<=nf;j++)
        fscanf(fp,"%d",&verts[i].flower[j]);
      npol += nf;
  }
  npol /= 3;

  bsurf = new ByuSurface();
  if(!bsurf) BFLATEXIT("out of memory");

  // create a new filename
  bsurf->SetFilename(filename); 
  bsurf->setSize(nv,npol);

  gotcents = 0;
  while(!feof(fp)) {
        line[0] = 0;
        fgets(line,1023,fp);
        if((line[0] == 0)||(strlen(line) < 1))
                continue;

        if(!strncmp(line,"CENTERS:",8)) {
           gotcents = 1;
           for(i=0;i<nv;i++) {
              fscanf(fp,"%f %f",&tx,&ty);
              if(geo == SPHERICAL) {
                x = (float)(cos(tx) * sin(ty));
                y = (float)(sin(tx) * sin(ty));
                z = (float)(cos(ty));
              } else { // EUCLIDIAN (or HGEOM?)
                x = tx;
                y = ty;
                z = 0.f;
              }
          bsurf->vertices()[i].set(x,y,z);
           }
        }
        if(!strncmp(line,"CIRCLE-COLORS:",14)) {
           for(i=0;i<nv;i++) fscanf(fp,"%d %d",&tmp,&verts[i].color);
        }
  }
  fclose(fp);
  if(!gotcents) BFLATEXIT("failed to find CENTERS: tag")

  nw = 0;
  for(i=1;i<=nv;i++)
    for(j=0;j<verts[i-1].num;j++)
      if((i < verts[i-1].flower[j]) && (i < verts[i-1].flower[j+1])) {
    if(nw >= npol) {
      cerr << "WARNING: loadFlatMap :: poly out of range" << endl;
      continue;
    }
    bsurf->facets()[nw][0] = i-1;
    bsurf->facets()[nw][1] = verts[i-1].flower[j]-1;
    bsurf->facets()[nw][2] = verts[i-1].flower[j+1]-1;
        nw++;
      }
  bsurf->geometryChanged();

  // Now, make sure flat map is oriented correctly for viewing
  return(bsurf);

bflat_fail:
 
  if(verts) delete [] verts;
  delete bsurf;

  return((ByuSurface*)NULL);
}
