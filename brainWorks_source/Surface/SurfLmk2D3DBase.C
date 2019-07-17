///////////////////////////////////////////////////////////////////////////
//
// File: SurfLmk2D3DBase.C
//
// Author: Keith Doolittle/Michael Bowers
//
// Purpose:  Implement 2D3DLmks
//
///////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <string.h>
#include <sys/param.h>
#include <TDL/Surface.h>
#include <TDL/ByuSurface.h>
#include <TDL/SurfaceUtils.h>
#include <FunctionTracing.h>
#include <SurfLmk2D3DBase.h>

#define LINE_WIDTH  5.0

#define INDICES_TITLE "Indices 1.0"

bool SurfLmk2D3DBase::loadIndicesFromFile(char *fname, char* statusMsg)
{
 int i,n,j;
 unsigned int idx;
 ifstream fp;
 char line[1024];
 char lname[1024];

  if (!surface)
  {
     sprintf(statusMsg, "ERROR: No surface set for SurfLmk2D3DBase");
     return false;
  }

  fp.open(fname);
  if (!fp || fp.fail())
  {
     sprintf(statusMsg, "ERROR: cannot open '%s' for reading", fname);
     return false;
  }

  fp.getline(line,1023);
  if (fp.fail() || strncmp(line, INDICES_TITLE, strlen(INDICES_TITLE)))
  {
      sprintf(statusMsg, "ERROR: '%s' is not a valid indices file", fname);
      return false;
  }

  SurfLmk2D3DBase::clearLandmarks();
  
   while (!fp.eof())
   {
       fp.getline(line, 1023);
       if (sscanf(line, "%ld %n", &idx, &j) < 1)
          continue;
          
       if ((idx >= surface->getNumVert()) || (idx >= fsurface->getNumVert()))
       {
          fp.close();
          clearLandmarks();
          sprintf(statusMsg, "ERROR: the landmark indices don't fit these surfaces");
          return false;
       }
       
       sprintf(lname,"%s",&(line[j]));
       indices.add(idx);

       float x = surface->vertices()[idx].x();
       float y = surface->vertices()[idx].y();
       float z = surface->vertices()[idx].z();

       float x2 = fsurface->vertices()[idx].x();
       float y2 = fsurface->vertices()[idx].y();
       float z2 = fsurface->vertices()[idx].z();

       landmarks.add(x,y,z,lname,1);
       flandmarks.add(x2,y2,z2,lname,1);
   }
   
   fp.close();
   return true;
}

bool SurfLmk2D3DBase::saveIndicesToFile(char *fname, char* statusMsg)
{
  ofstream fp;
  int i,n;
  char lname[1024];
  const Point *P;

  n = indices.size();
  if (n <= 0)
  {
     sprintf(statusMsg, "ERROR: No indices to save");
     return false;
  }

  fp.open(fname);
  if(!fp || fp.fail())
  {
     sprintf(statusMsg, "ERROR: cannot open '%s' for writing", fname);
     return false;
  }

  fp << INDICES_TITLE << endl;

  for (i = 0; i < n; i++)
  {
     sprintf(lname,"<no-name>");
     fp << indices[i] << " " << lname << endl;
  }

  fp.close();
  return true;
}


unsigned int SurfLmk2D3DBase::findLandmark(Surface *S, float x, float y, float z)
{
  unsigned int i,mini;
  double dx,dy,dz,d,mind;

  if(!S) 
      return((unsigned int)-1);
   
  for(i=0;i<S->getNumVert();i++) {
      dx = x - S->vertices()[i].x();
      dy = y - S->vertices()[i].y();
      dz = z - S->vertices()[i].z();
      d  = dx*dx+dy*dy+dz*dz;
      if((i==0)||(d<mind)) {
         mind = d;
         mini = i;
      }
  }
  return(mini);
}

bool SurfLmk2D3DBase::loadLandmarksFromFile(char *fname, int whichSurface, char* statusMsg)
{
  char lname[1024];
  clearLandmarks();
  
  Surface* surface1;
  Surface* surface2;
  Landmark* landmarks1 = &landmarks;
  Landmark* landmarks2 = &flandmarks;

  if (whichSurface == 1)
  {
    surface1 = surface;
    surface2 = fsurface;
  }
  else if (whichSurface == 2)
  {
    surface1 = fsurface;
    surface2 = surface;
    landmarks1 = &flandmarks;
    landmarks2 = &landmarks;
  }
  else
  {
     sprintf(statusMsg, "ERROR: BrainWorks Programming Error.");
     return false;
  }

  if(!surface1 || !surface2)
  {
     sprintf(statusMsg, "ERROR: No surface set for SurfLmk2D3DBase");
     return false;
  }
  
  int i,j;
  const Point *P;

   if(landmarks1->load(fname) != ItXSuccess)
    {
        sprintf(statusMsg, "ERROR: landmarks could not be loaded.");
        return false;
    }

    indices.setSize(landmarks1->num());

    for (i = 0; i < landmarks1->num(); i++)
    {
        if (!(P = landmarks1->point(i)))
        {
            indices[i] = (unsigned int) -1;
            continue;
        }

        if (landmarks1->name(i))
            sprintf(lname, "%s", landmarks1->name(i));
        else
            sprintf(lname, "<no-name>");

        indices[i] = findLandmark(surface1, P->x(), P->y(), P->z());
        if(indices[i] == (unsigned int)-1)
        continue;

        float x1 = P->x();
        float y1 = P->y();
        float z1 = P->z();
        float x2 = surface2->vertices()[indices[i]].x();
        float y2 = surface2->vertices()[indices[i]].y();
        float z2 = surface2->vertices()[indices[i]].z();

        landmarks2->add(x2, y2, z2, lname, 1);
    }
    return true;
}


bool SurfLmk2D3DBase::saveLandmarksToFile(char *fname, int whichSurface, char* statusMsg)
{
VIN("void SurfLmk2D3DBase::saveLandmarks()")

  int i,n;
  const Point *P;
  char lname[1024];
  
  Landmark* theLandmarks = &landmarks;

  if (whichSurface == 1)
    ;
  else if (whichSurface == 2)
    theLandmarks = &flandmarks;
  else
  {
     sprintf(statusMsg, "ERROR: BrainWorks Programming Error.");
     return false;
  }
  
  n = theLandmarks->num();
  if (n <= 0)
  {
     sprintf(statusMsg, "ERROR: No landmarks to save.");
     return false;
  }
    
  if (theLandmarks->save(fname) != ItXSuccess)
  {
    sprintf(statusMsg, "ERROR: landmarks did not save.");
    return false;
  }
  
  return true;

VOUT("void SurfLmk2D3DBase::saveLandmarks()")
}


void SurfLmk2D3DBase::clearLandmarks()
{
VIN("void SurfLmk2D3DBase::clearLandmarks()")

  landmarks.clear();
  flandmarks.clear();
  indices.clear();
  
VOUT("void SurfLmk2D3DBase::clearLandmarks()")
}

//
// Same as SurfaceViewer(), only 1 view (surf1)
//
SurfLmk2D3DBase::SurfLmk2D3DBase(Surface *surf, Surface *fsurf) :
                surface(surf), fsurface(fsurf)
{
VIN("SurfLmk2D3DBase::SurfLmk2D3DBase(LoadedSurface *surf)")
VOUT("SurfLmk2D3DBase::SurfLmk2D3DBase(LoadedSurface *surf)")
}


SurfLmk2D3DBase::~SurfLmk2D3DBase()
{
VIN("SurfLmk2D3DBase::~SurfLmk2D3DBase()")
VOUT("SurfLmk2D3DBase::~SurfLmk2D3DBase()")
}


