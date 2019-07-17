///////////////////////////////////////////////////////////////////////////
//
// File: LoadedSurface.C
//
// Author: Keith Doolittle
//
// Purpose:
//
///////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <iostream.h>
#include <stdlib.h>
#include <math.h>
#include <KDApplication.h>
#include <KDShell.h>
#include <TDL/Surface.h>
#include <AppColors.h>
#include <LoadedSurface.h>
#include <GLSurface.h>
#include <SurfaceColor.h>

LoadedSurface::LoadedSurface(Surface *bsurf)
{
VIN("LoadedSurface::LoadedSurface(Surface *bsurf)")
surf = bsurf;
if(!bsurf) 
     App->ShowMessage("WARNING: invalid surface");

glsurf = new GLSurface();

ref_count = 0;
surfaceChanged();
VOUT("LoadedSurface::LoadedSurface(Surface *bsurf)")
}


LoadedSurface::~LoadedSurface()
{
VIN("LoadedSurface::~LoadedSurface()")
if(surf)   delete surf;
if(glsurf) delete glsurf;
VOUT("LoadedSurface::~LoadedSurface()")
}

void LoadedSurface::surfaceChanged()
{
VIN("void LoadedSurface::surfaceChanged()")
glsurf->setSurface(surf,glsurf->hasStrips());
VOUT("void LoadedSurface::surfaceChanged()")
}

const char *LoadedSurface::Filename()
{
VIN("(and Leaving) const char *LoadedSurface::Filename()")
if(surf) return(surf->Filename());
else return("<surface>");
}


