///////////////////////////////////////////////////////////////////////////
//
// File: Reducer.C
//
// Author: Keith Doolittle
//
// Purpose:
//
///////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <iostream.h>
#include <KDApplication.h>
#include <KDShell.h>
#include <KDWorkingChild.h>
#include <TDL/Surface.h>
#include <TDL/SurfaceReducer.h>
#include <AppIcons.h>
#include <IconList.h>
#include <Main.h>
#include <Reducer.h>

#define MIN_CURV	0.0
#define MAX_CURV	MAX_CURV_RES
#define MIN_ALPHA	0.0
#define MAX_ALPHA	1.0

static void OKCallbackCB(Widget wid, XtPointer cld, XtPointer cad)
{
VIN("static void OKCallbackCB(Widget wid, XtPointer cld, XtPointer cad)")
	if(cld) ((Reducer*)cld)->Callback(1);
	
VOUT("static void OKCallbackCB(Widget wid, XtPointer cld, XtPointer cad)")
}

static void CancelCallbackCB(Widget wid, XtPointer cld, XtPointer cad)
{
VIN("static void CancelCallbackCB(Widget wid, XtPointer cld, XtPointer cad)")
if(cld) ((Reducer*)cld)->Callback(-1);
VOUT("static void CancelCallbackCB(Widget wid, XtPointer cld, XtPointer cad)")
}

static void WorkingCB(char *str, void *data)
{
VIN("static void WorkingCB(char *str, void *data)")
if(data) ((Reducer*)data)->Working(str);
VOUT("static void WorkingCB(char *str, void *data)")
}

static void ChildKilledCB(int val, void *dat)
{
VIN("static void ChildKilledCB(int val, void *dat)")
if(dat) ((Reducer*)dat)->Callback(val);
VOUT("static void ChildKilledCB(int val, void *dat)")
}

static void SurfaceChangedCB(Icon *ic, void *dat)
{
  if(dat) ((Reducer*)dat)->surfaceChanged();
}

void Reducer::surfaceChanged()
{
  LoadedSurface *lsurf;
  Surface *surf;
  Icon *ic;
  double f,fmin,fmax;
  int   i;
  char estr[100];

  ic = surface_list->getSelectedIcon();
  if(!ic || (!ic->Data()))
	return;

  lsurf = (LoadedSurface*)ic->Data();
  surf = lsurf->surface();
  if(!surf) return;

  if(!surf->hasCurvature())
      surf->genCurvature(Surface::MaxCurvature);

  fmin =  1E20;
  fmax = -1E20;
  for(i=0;i<surf->curvature().getNelm();i++) {
      f = fabs(surf->curvature()[i]);
      if(f<fmin) fmin = f;
      if(f>fmax) fmax = f;
  }

  sprintf(estr,"%.4f",fmin);
  App->SetTextString(min_curv_wid,estr);
  sprintf(estr,"%.4f",fmax);
  App->SetTextString(max_curv_wid,estr);
}


Reducer::Reducer()
{
VIN("Reducer::Reducer()")
int nitems = 7;
int w = App->BigFontWidth() * 20;
int h = nitems * (App->BigFontHeight() + 15);
child = NULL;
int x = App->Width()/2  - w/2;
int y = App->Height()/2 - h/2;
char estr[100];
float mincurv,maxcurv;

u_long flgs = KDSHELL_TITLEBAR;

shell = new KDShell("Polygon Reduction",
                x,y,w,h,flgs, DECIMATE_ICON,
                MainWindow->Shell(),LABELCOLOR);

int ddy = 100/nitems;
int yst = 0;
int ynd = ddy;
ATTACHMENTS_TBLR(yst,ynd,0,100)
surface_list = new IconList(shell->Main(),"Surface",IXSurfaceIcon,
		NULL,0,arglist,c);
surface_list->AddCallback(SurfaceChangedCB,(void*)this);

yst  = ynd;
ynd += ddy;
ATTACHMENTS_TBLR(yst,ynd,0,100)
min_curv_wid = App->LabelText(shell->Main(),"Abs Min Curvature","0.1",10,
                        App->SmallFontList(),arglist,c);

yst  = ynd;
ynd += ddy;
ATTACHMENTS_TBLR(yst,ynd,0,100)
max_curv_wid = App->LabelText(shell->Main(),"Abs Max Curvature","0.1",10,
                        App->SmallFontList(),arglist,c);

yst  = ynd;
ynd += ddy;
ATTACHMENTS_TBLR(yst,ynd,0,100)
Widget iter_wid = App->LabelText(shell->Main(),"Maximum Iterations","1",10,
                        App->SmallFontList(),arglist,c);

yst  = ynd;
ynd += ddy;
ATTACHMENTS_TBLR(yst,ynd,0,100)
Widget alpha_wid = App->LabelText(shell->Main(),"Redistribute Alpha","0.005",10,
                        App->SmallFontList(),arglist,c);

yst  = ynd;
ynd += ddy;
ATTACHMENTS_TBLR(yst,ynd,0,100)
Widget riter_wid = App->LabelText(shell->Main(),"Redistribute Iterations",
			"50",10, App->SmallFontList(),arglist,c);

yst  = ynd;
ynd += ddy;
yst += 2;
ynd -= 2;
ATTACHMENTS_TBLR(yst,ynd,15,45)
Widget btn = App->Button(shell->Main(),"Apply",
                        App->MedFontList(),OKCallbackCB,
                        (XtPointer)this,arglist,c);

ATTACHMENTS_TBLR(yst,ynd,55,85)
btn = App->Button(shell->Main(),"CANCEL",
                        App->MedFontList(),CancelCallbackCB,
                        (XtPointer)this,arglist,c);
XmChangeColor(btn,CANCELCOLOR);

surfaceChanged();

retry_this:

shell->Show();
done = 0;
while(!done) App->ProcessEvents();

shell->Hide();

if(done < 0) {
	delete this;
	return;
	}

Icon *ic = surface_list->getSelectedIcon();
if(!ic) {
	App->ShowMessage("ERROR: selected surface is invalid");
	goto retry_this;
	}
LoadedSurface *lsurf = (LoadedSurface*)ic->Data();
if(!lsurf) {
	App->ShowMessage("ERROR: selected surface is invalid");
	goto retry_this;
	}
Surface *surf = lsurf->surface();
if(!surf) {
	App->ShowMessage("ERROR: selected surface is invalid");
	goto retry_this;
	}

char estring[1024];
float min_curv = App->GetFloat(min_curv_wid);
float max_curv = App->GetFloat(max_curv_wid);

if((min_curv<MIN_CURV)||(min_curv > MAX_CURV)) {
	sprintf(estring,"ERROR: min curvature must be between %.2f and %.2f",
		MIN_CURV,MAX_CURV);
	App->ShowMessage(estring);
	goto retry_this;
	}

if((max_curv<MIN_CURV)||(max_curv > MAX_CURV)) {
	sprintf(estring,"ERROR: max curvature must be between %.2f and %.2f",
		MIN_CURV,MAX_CURV);
	App->ShowMessage(estring);
	goto retry_this;
	}

int niter = App->GetInt(iter_wid);
if((niter<0)||(niter > 1000)) {
	App->ShowMessage("ERROR: max iterations must be between 1 and 1000");
	goto retry_this;
	}

float alpha = App->GetFloat(alpha_wid);
if((alpha<MIN_ALPHA)||(alpha > MAX_ALPHA)) {
	sprintf(estring,"ERROR: redist alpha must be between %.2f and %.2f",
		MIN_ALPHA,MAX_ALPHA);
	App->ShowMessage(estring);
	goto retry_this;
	}

int rniter = App->GetInt(riter_wid);
if((rniter<0)||(rniter > 1000)) {
	App->ShowMessage("ERROR: redist iterations must be between 1 and 1000");
	goto retry_this;
	}


sprintf(estring,"Reducing %s",App->FilenameTail(surf->Filename()));

child = new KDWorkingChild(estring,GEN_XFORM_ICON);
child->SetNotify(ChildKilledCB,-1,(void*)this);

int nremoved;

SurfaceReducer r;

lsurf->reference();


r.SetWorkProc(WorkingCB,(void*)this);
ItXECode rv = r.decimate(surf,min_curv,max_curv,niter,alpha,rniter,nremoved);

r.SetWorkProc(NULL,NULL);

lsurf->unreference();
lsurf->surfaceChanged();

if(rv != ItXSuccess) App->ShowMessage("ERROR: decimate failed!");
else 	{
	sprintf(estring,"%s : %d points removed",
		App->FilenameTail(surf->Filename()),nremoved);
	App->ShowMessage(estring);
	}

delete this;
VOUT("Reducer::Reducer()")
}


Reducer::~Reducer()
{
VIN("Reducer::~Reducer()")
if(surface_list) delete surface_list;
if(child) delete child;
if(shell) delete shell;
VOUT("Reducer::~Reducer()")
}
