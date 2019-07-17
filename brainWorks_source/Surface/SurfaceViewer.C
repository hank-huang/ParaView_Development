///////////////////////////////////////////////////////////////////////////
//
// File: SurfaceViewer.C
//
// Author: Keith Doolittle
//
// Purpose:
//
///////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <sys/param.h>
#include <KDApplication.h>
#include <Main.h>
#include <TDL/Surface.h>
#include <TDL/ByuSurface.h>
#include <TDL/SurfaceUtils.h>
#include <SurfaceBox.h>
#include <LoadedSurface.h>
#include <SurfaceViewer.h>
#include <KDShell.h>
#include <Icon.h>
#include <IconList.h>
#include <GLPrimitive.h>

int select_surface_done;
DEFINECB(SurfaceViewer,AddSurfaceCB,addSurf)

static void SelectSurfaceCB(Widget wid, XtPointer cld, XtPointer cad) 
{
VIN("static void SelectSurfaceCB(Widget wid, XtPointer cld, XtPointer cad) ")
  // who did this???? select_surface_done = 1;
  select_surface_done = (int)cld;
VOUT("static void SelectSurfaceCB(Widget wid, XtPointer cld, XtPointer cad) ")
}

static void PickCB(u_long dat, void *s, int indx)
{
VIN("static void PickCB(u_long dat, void *s, int indx)")
  if(dat) ((SurfaceViewer*)dat)->pickPoint((Surface*)s,indx);
VOUT("static void PickCB(u_long dat, void *s, int indx)")
}
	
SurfaceViewer::SurfaceViewer(LoadedSurface *lsurf) :
  KDChildWindow(App->FilenameTail(lsurf->Filename ()),50,50) 
{
VIN("SurfaceViewer::SurfaceViewer(LoadedSurface *lsurf) :")
  destroy_func = NULL;
		
  destroy_data = NULL;
  surfBox = NULL;
  
  if(!lsurf || !lsurf->glsurface()) {
      App->ShowMessage("ERROR: invalid or NULL surface. Cannot view");
      delete this;
      return;
    }

  c = 0;
  Set(XmNleftAttachment,XmATTACH_FORM);
  Set(XmNrightAttachment,XmATTACH_FORM);
  Set(XmNtopAttachment,XmATTACH_FORM);
  Set(XmNbottomAttachment,XmATTACH_POSITION);
  Set(XmNbottomPosition,95);
  surfBox = new SurfaceBox(this->Form(),NULL,arglist,c);
	
  ATTACHMENTS_TBLR(95,100,0,10)
  Widget  btn=   App->Button(this->Form(),"Add",App->MedFontList(),
			     AddSurfaceCB,(XtPointer)this,arglist,c);

//keith
  c = 0;
  Set(XmNleftAttachment,XmATTACH_WIDGET);
  Set(XmNleftWidget,btn);
  Set(XmNrightAttachment,XmATTACH_POSITION);
  Set(XmNrightPosition,100);
  Set(XmNtopAttachment,XmATTACH_POSITION);
  Set(XmNtopPosition,95);
  Set(XmNbottomAttachment,XmATTACH_POSITION);
  Set(XmNbottomPosition,100);
  Set(XmNresizeWidth,False);
  Set(XmNresizeHeight,False);
  Set(XmNrecomputeSize,False);
  Set(XmNalignment,XmALIGNMENT_BEGINNING);
  m_text = App->Label(this->Form(),"",App->SmallFontList(),arglist,c);

  int mainw = App->Width()/3;
  int mainh = App->Height()/3;
  SetSize(mainw,mainh);
  surfBox->SetPickStyle(PickPoint,PickCB,(u_long)this);
  Show();

  doAddSurf(lsurf);
VOUT("SurfaceViewer::SurfaceViewer(LoadedSurface *lsurf) :")
}


void SurfaceViewer::addSurf()
{
VIN("void SurfaceViewer::addSurf()")
  LoadedSurface *S = SelectLoadedSurface("Select Surface to Add");
  if(!S || !surfBox) 
     return;

  doAddSurf(S);
VOUT("void SurfaceViewer::addSurf()")
}

void SurfaceViewer::doAddSurf(LoadedSurface *S)
{
VIN("void SurfaceViewer::doAddSurf(LoadedSurface *S)")
  surfBox->addSurf(S);
VOUT("void SurfaceViewer::doAddSurf(LoadedSurface *S)")
}


void SurfaceViewer::pickPoint(Surface *S,int indx) 
{
VIN("void SurfaceViewer::pickPoint(Surface *S,int indx) ")
  char estr[1024];
  char cstr[100];

  if(S && (S->hasCurvature()) && (indx>=0) && (indx < S->getNumVert()))
       sprintf(cstr,"curvature = %.3f",S->curvature()[indx]);
  else cstr[0] = 0;


  sprintf(estr,"Vertex %d (%.2lf, %.2lf, %.2lf) %s",indx,
           S->vertices()[indx].x(),
           S->vertices()[indx].y(),
           S->vertices()[indx].z(), cstr);

  App->SetLabelText(m_text,estr,App->SmallFontTag());
VOUT("void SurfaceViewer::pickPoint(Surface *S,int indx) ")
}

SurfaceViewer::~SurfaceViewer()
{
VIN("SurfaceViewer::~SurfaceViewer()")
  if(destroy_func)
    destroy_func(destroy_data);

  if(surfBox) 
	delete surfBox;
VOUT("SurfaceViewer::~SurfaceViewer()")
}






