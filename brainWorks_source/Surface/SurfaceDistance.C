#include <stdlib.h>
#include <iostream.h>
#include <math.h>
#include <KDApplication.h>
#include <Main.h>
#include <OpMessage.h>
#include <TDL/ItXVolume.h>
#include <TDL/Surface.h>
#include <ADL/Array2D.h>
#include <OcTree.h>
#include <PointGroup.h>
#include <SurfaceDistance.h>
#include <DynArray.h>
#include <time.h>

#define USE_PGROUP

#define MAX_POL_PER_VERT	10

#include <AWHistogram.h>
#include <PlotImage.h>

bool SurfaceDistance::calculate(Surface &S,
	ItXVolume &V, ItXVolume &Seg, float white_thresh,
    ItXVolume *rDist, ItXVolume *rInd)
{
   AWHistogram hist;
   AWHistogram abshist;

   if (!SurfaceDistanceBase::calculate(S, V, Seg, white_thresh,
        rDist, rInd, hist, abshist))
   {
        OpMessage::GenericWorkProc(NULL);
        return false;
   }
    
    OpMessage::GenericWorkProc(NULL);
    
    PlotImage *pp = new PlotImage("Histogram",400,400);
    pp->addCurve(hist,"Histogram");
    pp->addCurve(abshist,"Abs-Histogram");

 return(true);
}



