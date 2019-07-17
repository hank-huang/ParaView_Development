///////////////////////////////////////////////////////////////////////////
//
// File: CurveMatchBase.C
//
// Author: Muge Bakircioglu / Michael Bowers
//
// Purpose:
//
///////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <iostream.h>
#include <sys/param.h>
#include <ADL/Array1D.h>
#include <CurveMatchBase.h>
#include <TDL/ByuSurface.h>
#include <TDL/SurfaceUtils.h>

Curve CurveMatchBase::origTarget;
Curve CurveMatchBase::origTemplate;
Curve CurveMatchBase::newTarg;
Array1D<Point> CurveMatchBase::tempPoints, CurveMatchBase::targPoints;
int CurveMatchBase::numNghbs, CurveMatchBase::spacing, CurveMatchBase::dist;
     
void CurveMatchBase::addToArray(Array1D<Point> &P, Point A)
{
  for(int i=0;i<P.getNelm();i++)
    P[i] += A;
}

CurveMatchBase::CurveMatchBase()
{
}

CurveMatchBase::~CurveMatchBase()
{
}

bool CurveMatchBase::readFile(  int         curveno,
                                const char* fileName,
                                Array1D<Point> &P,
                                      char* statusMsg)
{
    double x,y,z;
   
    ifstream fp(fileName, ios::in);
    if (!fp)
    {
        sprintf(statusMsg, "ERROR: can\'t open file %s", fileName);
        return false;
    }
    
    int count = 0;
    int numPoints;
    fp >> numPoints;
    P.setDim(numPoints);
    while (fp >> x >> y >> z)
    {
       P[count++].set(x, y, z);
    }
    double dummy = 0.0;
    for (int i = 0; i < numPoints; ++i)
    {
        if (curveno == 1)
            _target.addPoint(P[i], dummy, dummy, dummy, dummy);
        else 
            _template.addPoint(P[i], dummy, dummy, dummy, dummy);
    }
    //store original curves in case we run matching multiple times
    //we want to run it on the original curve, not the resampled one
    //this especially matters if we change the order of polynomial we
    //are fitting the curves, cannot go back to a higher order poly if 
    //original curves are not stored
        
    if (curveno==1) origTarget=_target;
    else origTemplate=_template;
    
    return true;
}

void CurveMatchBase::setCurvesApart(Array1D<Point> &P, int curveno) 
{
  
  double diff;
     //if both are loaded, set one apart for display purposes
     //first find the direction of least variation, then set P
     //apart in that direction
     if (_template.numSamples()!=0 &&  _target.numSamples()!=0) {
       double x_min=HUGE;  double x_max=-HUGE;
       double y_min=HUGE;  double y_max=-HUGE;
       double z_min=HUGE;  double z_max=-HUGE;
       for (int i=0; i<P.getNelm(); ++i) {
     if (P[i].x() < x_min) x_min=P[i].x(); 
     else  if (P[i].x() > x_max) x_max=P[i].x();

     if (P[i].y() < y_min) y_min=P[i].y();
     else if (P[i].y() > y_max) y_max=P[i].y();

     if (P[i].z() < z_min) z_min=P[i].z();
     else  if (P[i].z() > z_max) z_max=P[i].z();
    
       }
       double x_diff=x_max-x_min;
       double y_diff=y_max-y_min;
       double z_diff=z_max-z_min;

       int min;
       if (x_min < y_min) {
     if (x_min < z_min){
       min=0;
       diff=(_target.getPoint(0)).x()-(_template.getPoint(0)).x(); 
     }
     else {
       min=2;
       diff=(_target.getPoint(0)).z()-(_template.getPoint(0)).z();
     }
       }
       else {
     if (y_min < z_min) {
       min=1;
       diff=(_target.getPoint(0)).y()-(_template.getPoint(0)).y();
     }
     else {
       min=2;
       diff=(_target.getPoint(0)).z()-(_template.getPoint(0)).z();
     }
       }
       _minFlag=min;
       if (fabs(diff) <50 ) { 
     
     if ((diff < 0 && curveno==1) || (diff >0 && curveno==2))  
       diff=-50-diff;
     else diff=50-diff;
    
     switch(min) {
     case 0: //KWD P+=Point(diff, 0.0,0.0); break;
            addToArray(P,Point(diff,0.0,0.0));
        break;
     case 1: //KWD P+=Point(0.0, diff, 0.0); break;
            addToArray(P,Point(0.0,diff,0.0));
        break;
     case 2: //KWD P+=Point(0.0, 0.0, diff); break;
            addToArray(P,Point(0.0,0.0,diff));
        break;
     default: break;
     }
        
       }
     }
     else { return; }
     
     _offset=diff;
}

bool CurveMatchBase::checkMatchingParameters(int nsamples,
                                         int nlmks,
                                         int polyOrder,
                                         int nghbdsize,
                                         char* statusMsg)
{
    if (nsamples <= 0)
    {
        sprintf(statusMsg, "ERROR: number of samples must be positive");
        return false;
    }
    if (nlmks <= 0)
    {
        sprintf(statusMsg, "ERROR: number of lmks must be positive");
        return false;
    }
    if (polyOrder <= 0)
    {
        sprintf(statusMsg, "ERROR: polynomial order must be positive");
        return false;
    }
    if (nghbdsize <= 0)
    {
        sprintf(statusMsg, "ERROR: neighborhood size  must be positive");
        return false;
    }
    
    return true;
}


int CurveMatchBase::match(int nsamples, int nlmks, int polyOrder, int nghbdsize, double a_coeff, double c_coeff, double t_coeff)
{
   
    //use the original curves stored
    _target = origTarget;
    _template = origTemplate;
    newTarg.reset();
    Array1D<Point> storeMatch(2);
    int i;
    
    char statusMsg[1024];
    if (!checkMatchingParameters(nsamples, nlmks, polyOrder, nghbdsize, statusMsg))
        return 0;

    _template.calculateMatchcost(_target, nsamples, nlmks, polyOrder,
                                 a_coeff, c_coeff, t_coeff, nghbdsize,
                                 newTarg, numNghbs, spacing, dist);
    double matchCost = _template.matching(newTarg, numNghbs, spacing, dist);
    tempPoints.setDim(newTarg.numSamples());
    targPoints.setDim(newTarg.numSamples());

    _template.getMatch(newTarg, tempPoints, targPoints);
    
    Array1D<Point> P(nsamples);
    Array1D<Point> Q(nsamples);

    for (i = 0; i < nsamples; ++i)
    {
        P[i] = _target.getPoint(i);
        Q[i] = _template.getPoint(i);
    }

    for (i = 0; i < newTarg.numSamples(); ++i)
    {
        storeMatch[0] = tempPoints[i];
        storeMatch[1] = targPoints[i];
    }
    return 1;

}

int CurveMatchBase::saveMatch(  const char* curve1LmkFile,
                                const char* curve2LmkFile,
                                      bool  saveAsBWLmks,
                                      char* statusMsg)
{

    // check that there are points to output
    if (tempPoints.getSize() == 0)
    {
        sprintf(statusMsg, "ERROR: There is no match to save\n");
        return 0;
    }
    
    ofstream fpCurve1(curve1LmkFile, ios::out);
    if (!fpCurve1)
    {
        sprintf(statusMsg, "ERROR:  can't save to Curve 1 Landmark Output File:  %s.",
                curve1LmkFile);
        return 0;
    }
    
    ofstream fpCurve2(curve2LmkFile, ios::out);
    if (!fpCurve2)
    {
        sprintf(statusMsg, "ERROR:  can't save to Curve 2 Landmark Output File:  %s.",
                curve2LmkFile);
        return 0;
    }

    if (saveAsBWLmks)
    {   
        int i;
        fpCurve1 << "Landmarks-1.0" << endl;
        fpCurve2 << "Landmarks-1.0" << endl;
        fpCurve1 << newTarg.numSamples() << endl;
        fpCurve2 << newTarg.numSamples() << endl;
        for (i = 0; i < newTarg.numSamples(); i++)
        {
            fpCurve1 << "\"L" << i + 1 << "\"" << endl;
            fpCurve2 << "\"L" << i + 1 << "\"" <<endl;
            fpCurve1 << targPoints[i] << " 1 1" << endl;
            fpCurve2 << tempPoints[i] << " 1 1" << endl;
        }
    }
    else
    {
        int i;
        for (i = 0; i < newTarg.numSamples(); ++i)
        {
            fpCurve1 << targPoints[i] << endl;
            fpCurve2 << tempPoints[i] << endl;
        }
        fpCurve1.close();
        fpCurve2.close();
    }
}
