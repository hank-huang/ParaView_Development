///////////////////////////////////////////////////////////////////////////
//
// File: Redistributer.C
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
#include <TDL/SurfaceUtils.h>
#include <TDL/ByuSurface.h>
#include <AppIcons.h>
#include <LoadedSurface.h>
#include <IconList.h>
#include <Main.h>
#include <Redistributer.h>

#define MIN_CURV	0.0
#define MAX_CURV	0.9
#define MIN_ALPHA	0.0
#define MAX_ALPHA	1.0

static void OKCallbackCB(Widget wid, XtPointer cld, XtPointer cad)
{
VIN("static void OKCallbackCB(Widget wid, XtPointer cld, XtPointer cad)")
if(cld) ((Redistributer*)cld)->Callback(1);
VOUT("static void OKCallbackCB(Widget wid, XtPointer cld, XtPointer cad)")
}

static void CancelCallbackCB(Widget wid, XtPointer cld, XtPointer cad)
{
VIN("static void CancelCallbackCB(Widget wid, XtPointer cld, XtPointer cad)")
if(cld) ((Redistributer*)cld)->Callback(-1);
VOUT("static void CancelCallbackCB(Widget wid, XtPointer cld, XtPointer cad)")
}

static void WorkingCB(char *str, void *data)
{
VIN("static void WorkingCB(char *str, void *data)")
if(data) ((Redistributer*)data)->Working(str);
VOUT("static void WorkingCB(char *str, void *data)")
}

static void ChildKilledCB(int val, void *dat)
{
VIN("static void ChildKilledCB(int val, void *dat)")
if(dat) ((Redistributer*)dat)->Callback(val);
VOUT("static void ChildKilledCB(int val, void *dat)")
}

Redistributer::Redistributer()
{
VIN("Redistributer::Redistributer()")
child = NULL;
 
 int nitems = 4;
int w = App->BigFontWidth() * 20;
int h = nitems * (App->BigFontHeight() + 15);

int x = App->Width()/2  - w/2;
int y = App->Height()/2 - h/2;

u_long flgs = KDSHELL_TITLEBAR;

shell = new KDShell("Polygon Redistribution",
                x,y,w,h,flgs, DECIMATE_ICON,
                MainWindow->Shell(),LABELCOLOR);

int ddy = 100/nitems;
int yst = 0;
int ynd = ddy;
ATTACHMENTS_TBLR(yst,ynd,0,100)
surface_list = new IconList(shell->Main(),"Surface",IXSurfaceIcon,
			NULL,0,arglist,c);

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

sprintf(estring,"Redistributing %s",App->FilenameTail(surf->Filename()));

child = new KDWorkingChild(estring,GEN_XFORM_ICON);
child->SetNotify(ChildKilledCB,-1,(void*)this);

lsurf->reference();
surf->SetWorkProc(WorkingCB,(void*)this);

ItXECode rv = SurfaceUtils::redistribute(*surf,rniter,alpha);
surf->SetWorkProc(NULL,NULL);
lsurf->unreference();

if(rv != ItXSuccess) App->ShowMessage("ERROR: decimate failed!");
else lsurf->surfaceChanged();

delete this;
VOUT("Redistributer::Redistributer()")
}


Redistributer::~Redistributer()
{
VIN("Redistributer::~Redistributer()")
if(surface_list)
    delete surface_list;
if(child)
    delete child;
if(shell)
    delete shell;
VOUT("Redistributer::~Redistributer()")
}
 
