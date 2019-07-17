///////////////////////////////////////////////////////////////////////////
//
// File: FlatTracker.C
//
// Author: Keith Doolittle
//
// Purpose:
//
///////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <KDApplication.h>
#include <AppColors.h>
#include <Xm/Form.h>
#include <Xm/ScrolledW.h>
#include <Xm/RowColumn.h>
#include <TDL/Surface.h>
#include <TDL/SurfaceUtils.h>
#include <Main.h>
#include <MainUtils.h>
#include <FlatTracker.h>
#include <LoadedSurface.h>
#include <GLSurface.h>
#include <GLPrimitive.h>
#include <SurfaceBox.h>
#include <IconList.h>
#include <TDL/SurfaceTracker.h>
#include <SurfaceColor.h>

#define CURVE_R 0.2
#define CURVE_G 0.0
#define CURVE_B 0.9

#define T1_BACK   GRAY40
#define T1_LABEL  GRAY50
#define T1_BUTTON GRAY50
#define T1_TOGGLE GRAY50

DEFINECBP(FlatTracker,SetViewS,setView,FlatTracker::FLAT_S)
DEFINECBP(FlatTracker,SetViewF,setView,FlatTracker::FLAT_F)
DEFINECBP(FlatTracker,SetViewSF,setView,FlatTracker::FLAT_SF)
DEFINECB(FlatTracker,LoadCurvatureCB,loadCurvature)
DEFINECB(FlatTracker,ColorSelectedCB,setSelectedColor)
DEFINECB(FlatTracker,WidthSelectedCB,setSelectedWidth)
DEFINECB(FlatTracker,LoadPathsCB,loadPaths)
DEFINECB(FlatTracker,SavePathsCB,savePaths)
DEFINECB(FlatTracker,Use3DVertsToggleCB,use3DVertsToggle)
DEFINECB(FlatTracker,Base1ToggleCB,base1Toggle)
DEFINECB(FlatTracker,ManualTrackCB,manualTrack)
DEFINECB(FlatTracker,DeletePathCB,deletePath)
DEFINECB(FlatTracker,ClearPathsCB,doClearPaths)
DEFINECB(FlatTracker,LineWidthCB,setLineWidth)
DEFINECBP(FlatTracker,SelectButtonCB,selectButton,wid)

void PickPointCB(u_long d, void *, int idx)
{
VIN("void PickPointCB(u_long d, void *, int idx)")
  if(d) ((FlatTracker*)d)->pickPoint(idx);
VOUT("void PickPointCB(u_long d, void *, int idx)")
}

static char *type_names[3] = {
 "Geodesic",
 "Sulcus",
 "Gyrus"
};

void SetTypeCB(char *str, void *dat)
{
VIN("void SetTypeCB(char *str, void *dat)")
  SurfaceTracker::TrackType t;

  if(!strcmp(str,"Geodesic"))    t = SurfaceTracker::Geodesi;
  else if(!strcmp(str,"Sulcus")) t = SurfaceTracker::Sulci;
  else if(!strcmp(str,"Gyrus"))  t = SurfaceTracker::Gyri;
  else return;

  if(dat) ((FlatTracker*)dat)->setTrackType(t);
VOUT("void SetTypeCB(char *str, void *dat)")
}


FlatTracker::FlatTracker(LoadedSurface *S3, LoadedSurface *Sf) :
          KDChildWindow(App->FilenameTail(Sf->Filename()),
                        App->Width(),App->Height())
{
VIN("FlatTracker::FlatTracker(LoadedSurface *S3, LoadedSurface *Sf) :")
  m_surfbox       = NULL;
  m_flatbox       = NULL;
  m_trackTypeList = NULL;
  m_surfCurves    = NULL;
  m_flatCurves    = NULL;
  m_bUse3D        = true;
  m_trackType     = SurfaceTracker::Geodesi;
  m_startIndex    = -1;
  m_endIndex      = -1;
  m_numPaths      = 0;
  m_paths         = NULL;
  m_lineWidth     = 3;
  m_curSurf       = NULL;
  m_curFlat       = NULL;

  if(!S3 || !Sf || !(S3->surface()) || !(Sf->surface())) {
     App->ShowMessage("ERROR: one or more of the surfaces is invalid/empty");
     delete this;
     return;
  }

  if(S3->surface()->getNumVert() != Sf->surface()->getNumVert()) {
     App->ShowMessage("ERROR: the surfaces dont have the same # vertices");
     delete this;
     return;
  }

  if(S3->glsurface()->isFlat()) {
     App->ShowMessage("ERROR: the 3D surface is a Flat Map");
     delete this;
     return;
  }

  if(!Sf->glsurface()->isFlat()) {
     App->ShowMessage("ERROR: the Flat Map is a 3D surface");
     delete this;
     return;
  }

  m_lsurface  = S3;
  m_surface   = S3->surface();
  m_lflat     = Sf;
  m_flat      = Sf->surface();
  m_numPoints = m_surface->getNumVert();

  m_lsurface->reference();
  m_lflat->reference();

  // initially show both surface and flat map
  m_viewType = FLAT_SF;

  Arg arglist[30];
  Widget   wid;
  Cardinal c;
  int      bh;

  ATTACHMENTS_TBLR(1,94,1,39)
  m_surfbox = new SurfaceBox(this->Form(),m_lsurface,arglist,c);

  ATTACHMENTS_TBLR(1,94,40,79)
  m_flatbox = new SurfaceBox(this->Form(),m_lflat,arglist,c);

  ATTACHMENTS_TBLR(95,99,1,9)
  wid = App->Button(this->Form(),"S",App->MedFontList(),SetViewS,
              (XtPointer)this,arglist,c);
  XmChangeColor(wid,T1_BUTTON);

  ATTACHMENTS_TBLR(95,99,10,18)
  wid = App->Button(this->Form(),"F",App->MedFontList(),SetViewF,
              (XtPointer)this,arglist,c);
  XmChangeColor(wid,T1_BUTTON);

  ATTACHMENTS_TBLR(95,99,19,30)
  wid = App->Button(this->Form(),"SF",App->MedFontList(),SetViewSF,
              (XtPointer)this,arglist,c);
  XmChangeColor(wid,T1_BUTTON);

  ATTACHMENTS_TBLR(95,99,30,90)
  m_infoLabel = App->Label(this->Form(),"",App->SmallFontList(),arglist,c);

  // right hand side menu 
  bh  = App->BigFontHeight() + 3;
  ATTACHMENTS_TLRH(0,80,100,bh)
  wid = App->Label(this->Form(),"Type",App->MedFontList(),arglist,c);
  App->RecursiveColor(wid,T1_LABEL);

  ATTACHMENTS_WLRH(wid,80,100,bh*2)
  m_trackTypeList = new IconList(this->Form(),"",IXNoIcon,
                        type_names,3,arglist,c);
  m_trackTypeList->AddExtraCallback(SetTypeCB,(void*)this);
  wid = m_trackTypeList->Form();
  XmChangeColor(wid,T1_TOGGLE);

  ATTACHMENTS_WLRH(wid,80,100,bh)
  wid = App->Button(this->Form(),"Load Curvature",App->MedFontList(),
                    LoadCurvatureCB,(XtPointer)this,arglist,c);
  XmChangeColor(wid,T1_BUTTON);

  ATTACHMENTS_WLRH(wid,80,100,bh)
  wid = App->Button(this->Form(),"Load Paths",App->MedFontList(),
                    LoadPathsCB,(XtPointer)this,arglist,c);
  XmChangeColor(wid,T1_BUTTON);

  ATTACHMENTS_WLRH(wid,80,100,bh)
  wid = App->Button(this->Form(),"Save Paths",App->MedFontList(),
                    SavePathsCB,(XtPointer)this,arglist,c);
  XmChangeColor(wid,T1_BUTTON);

  ATTACHMENTS_WLRH(wid,80,100,bh)
  wid = App->Button(this->Form(),"Delete Current",App->MedFontList(),
                    DeletePathCB,(XtPointer)this,arglist,c);
  XmChangeColor(wid,T1_BUTTON);

  ATTACHMENTS_WLRH(wid,80,100,bh)
  wid = App->Button(this->Form(),"Clear All Paths",App->MedFontList(),
                    ClearPathsCB,(XtPointer)this,arglist,c);
  XmChangeColor(wid,T1_BUTTON);

  ATTACHMENTS_WLRH(wid,80,100,bh)
  wid = App->Button(this->Form(),"Set Line Width",App->MedFontList(),
                    LineWidthCB,(XtPointer)this,arglist,c);
  XmChangeColor(wid,T1_BUTTON);

  ATTACHMENTS_WLRH(wid,80,100,bh)
  wid = App->Button(this->Form(),"Set Color",App->MedFontList(),
                   ColorSelectedCB,(XtPointer)this,arglist,c);
  XmChangeColor(wid,T1_TOGGLE);

  ATTACHMENTS_WLRH(wid,80,100,bh)
  wid = App->Button(this->Form(),"Set Width",App->MedFontList(),
                   WidthSelectedCB,(XtPointer)this,arglist,c);
  XmChangeColor(wid,T1_TOGGLE);

  ATTACHMENTS_WLRH(wid,80,100,bh)
  wid = App->Toggle(this->Form(),"1-Based Vertices",App->MedFontList(),m_bBase1,
                   Base1ToggleCB,(XtPointer)this,arglist,c);
  m_base1Toggle = wid;
  XmChangeColor(wid,T1_TOGGLE);

  ATTACHMENTS_WLRH(wid,80,100,bh)
  wid = App->Toggle(this->Form(),"Use 3D Verts",App->MedFontList(),m_bUse3D,
                   Use3DVertsToggleCB,(XtPointer)this,arglist,c);
  m_vertToggle = wid;
  XmChangeColor(wid,T1_TOGGLE);

  ATTACHMENTS_WLRH(wid,80,100,bh)
  wid = App->Label(this->Form(),"Manual Path",App->MedFontList(),arglist,c);
  XmChangeColor(wid,T1_LABEL);

  ATTACHMENTS_WLRH(wid,80,100,bh)
  wid = App->LabelText(this->Form(),"Start Vert","0",10,App->SmallFontList(),
                   arglist,c);
  m_startText = wid;
  XmChangeColor(XtParent(wid),T1_LABEL);

  ATTACHMENTS_WLRH(wid,80,100,bh)
  wid = App->LabelText(this->Form(),"End Vert","0",10,App->SmallFontList(),
                   arglist,c);
  m_endText = wid;
  XmChangeColor(XtParent(wid),T1_LABEL);

  ATTACHMENTS_WLRH(wid,80,100,bh)
  wid = App->Button(this->Form(),"Track",App->MedFontList(),
                    ManualTrackCB,(XtPointer)this,arglist,c);
  XmChangeColor(wid,T1_BUTTON);

  c = 0;
  Set(XmNleftAttachment,XmATTACH_POSITION);
  Set(XmNleftPosition,80);
  Set(XmNrightAttachment,XmATTACH_FORM);
  Set(XmNtopAttachment,XmATTACH_WIDGET);
  Set(XmNtopWidget,wid);
  Set(XmNbottomAttachment,XmATTACH_FORM);
  Set(XmNscrollBarDisplayPolicy,XmAUTOMATIC);
  Set(XmNscrollBarPlacement,XmBOTTOM_RIGHT);
  Set(XmNscrollingPolicy,XmAUTOMATIC);
  Set(XmNresizeWidth,False);
  Set(XmNspacing,0);
  wid = XmCreateScrolledWindow(this->Form(),"ScrolledWindow",arglist,c);
  XtManageChild(wid);
  App->RecursiveColor(wid,T1_BACK);

  Dimension ww,ww2;
  Widget    cw,sb;

  cw = XtNameToWidget(wid,"ClipWindow");
  sb = XtNameToWidget(wid,"VertScrollBar");

  c = 0;
  Set(XmNwidth,&ww);
  XtGetValues(cw,arglist,c);

  c = 0;
  Set(XmNwidth,&ww2);
  XtGetValues(sb,arglist,c);

  ww -= ww2;
  ww -= 5;

  if(ww < 100) ww = 100;
  if(ww > 500) ww = 100;

  ATTACHMENTS_TBLR(0,100,0,100)
  Set(XmNwidth,ww);
  Set(XmNresizeWidth,False);
  Set(XmNorientation,XmVERTICAL);
  Set(XmNisHomogeneous,False);
  Set(XmNrubberPositioning,False);
  Set(XmNpacking,XmPACK_COLUMN);
  Set(XmNnumColumns,1);
  Set(XmNspacing,1);
  m_curveRC = XmCreateRowColumn(wid,"CurveRC",arglist,c);
  XtManageChild(m_curveRC);
  XmChangeColor(m_curveRC,T1_BACK);

  m_surfbox->SetPickStyle(PickPoint,PickPointCB,(u_long)this);
  m_flatbox->SetPickStyle(PickPoint,PickPointCB,(u_long)this);
  setTrackType(m_trackType); // this will initialize m_tracker

  if(!m_surface->hasCurvature())
     m_surface->genCurvature(Surface::MeanCurvature);

  Array1D<double> D = m_surface->curvature();
  m_flat->setCurvature(D);
VOUT("FlatTracker::FlatTracker(LoadedSurface *S3, LoadedSurface *Sf) :")
}

FlatTracker::~FlatTracker()
{
VIN("FlatTracker::~FlatTracker()")
  if(m_lsurface) m_lsurface->unreference();
  if(m_lflat)    m_lflat->unreference();
  cleanup();
VOUT("FlatTracker::~FlatTracker()")
}

void FlatTracker::cleanup()
{
VIN("void FlatTracker::cleanup()")
  if(m_surfCurves) {
     GLPrimitive::clearList(m_surfCurves);
     m_surfCurves = NULL;
  }

  if(m_flatCurves) {
     GLPrimitive::clearList(m_flatCurves);
     m_flatCurves = NULL;
  }

  if(m_trackTypeList) {
     delete m_trackTypeList;
     m_trackTypeList = NULL;
  }

  if(m_flatbox) {
     delete m_flatbox;
     m_flatbox = NULL;
  }

  if(m_surfbox) {
     delete m_surfbox;
     m_surfbox = NULL;
  }
VOUT("void FlatTracker::cleanup()")
}


void FlatTracker::pickPoint(int idx)
{
  if((idx<0)||(idx>=m_numPoints))
     return;

  if(m_startIndex < 0)
     setStart(idx);
  else
     setEnd(idx);
}

void FlatTracker::setStart(int idx)
{
  if(m_surfStart)
     m_surfCurves = GLPrimitive::removeFromList(m_surfCurves,m_surfStart);
  if(m_flatStart)
     m_flatCurves = GLPrimitive::removeFromList(m_flatCurves,m_flatStart);

  m_startIndex = idx;

  float x,y,z,sx,sy,sz;

  x = (float)m_surface->vertices()[idx].x();
  y = (float)m_surface->vertices()[idx].y();
  z = (float)m_surface->vertices()[idx].z();

  sx = (float)m_flat->vertices()[idx].x();
  sy = (float)m_flat->vertices()[idx].y();
  sz = (float)m_flat->vertices()[idx].z();

  m_surfStart = new GLPrimitive(x,y,z,1.f);
  m_surfStart->setColor(1.f,0.f,0.f);
  m_surfCurves = GLPrimitive::addToList(m_surfCurves,m_surfStart);

  m_flatStart = new GLPrimitive(sx,sy,sz,1.f);
  m_flatStart->setColor(1.f,0.f,0.f);
  m_flatCurves = GLPrimitive::addToList(m_flatCurves,m_flatStart);

  m_surfbox->setPrimitives(m_surfCurves);
  m_flatbox->setPrimitives(m_flatCurves);
  redraw();

  //
  // show info on label
  //
  if(m_bBase1)
     idx++;

  char tstr[1024];
  sprintf(tstr,"Vertex %04d (3D: %.2f, %.2f, %.2f) (2D: %.2f, %.2f)",
       idx,x,y,z,sx,sy);
  App->SetLabelText(m_infoLabel,tstr,App->SmallFontTag());

  sprintf(tstr,"%d",idx);
  App->SetTextString(m_startText,tstr);
}


void FlatTracker::setEnd(int idx)
{
  Array1D<int> inds;
  float x,y,z,sx,sy;

  m_endIndex = idx;
  //
  // show info on label
  //
  x = (float)m_surface->vertices()[idx].x();
  y = (float)m_surface->vertices()[idx].y();
  z = (float)m_surface->vertices()[idx].z();

  sx = (float)m_flat->vertices()[idx].x();
  sy = (float)m_flat->vertices()[idx].y();

  if(m_bBase1)
     idx++;

  char tstr[1024];
  sprintf(tstr,"Vertex %04d (3D: %.2f, %.2f, %.2f) (2D: %.2f, %.2f)",
       idx,x,y,z,sx,sy);
  App->SetLabelText(m_infoLabel,tstr,App->SmallFontTag());

  sprintf(tstr,"%d",idx);
  App->SetTextString(m_endText,tstr);

  //
  // calculate path between m_startIndex and m_endIndex
  //
  if(m_tracker.findPath(m_startIndex,m_endIndex,inds) != ItXSuccess) {
     App->ShowMessage("ERROR: tracking failed. See stdout");
     m_startIndex = -1;
     m_endIndex   = -1;
     return;
  }

  addPath(inds,"<no-name>");
  m_startIndex = -1;
  m_endIndex   = -1;
}


void FlatTracker::manualTrack()
{
  int sidx = App->GetInt(m_startText);
  int eidx = App->GetInt(m_endText);

  if((sidx<0)||(sidx>=m_numPoints)) {
     App->ShowMessage("ERROR: starting vertex is out of range");
     return;
  }
  if((eidx<0)||(eidx>=m_numPoints)) {
     App->ShowMessage("ERROR: ending vertex is out of range");
     return;
  }
  if(eidx == sidx) {
     App->ShowMessage("ERROR: start and end vertices are the same");
     return;
  }

  if(m_bBase1) {
     sidx--;
     eidx--;
  }

  setStart(sidx);
  setEnd(eidx);
}

GLPrimitive *FlatTracker::findPath(GLPrimitive *lst, Widget w)
{
  GLPrimitive *p;

  for(p=lst;p;p=p->Next())
      if((Widget)p->getHook()==w)
         return(p);

  return(NULL);
}

void FlatTracker::setSelectedColor()
{
  if(!m_curSurf || !m_curFlat) {
     App->ShowMessage("Please select a track first");
     return;
  }

  float r,g,b;
  m_curSurf->getColor(r,g,b);
  if(GetColor(&r,&g,&b)) {
     m_curSurf->setColor(r,g,b);
     m_curFlat->setColor(r,g,b);
     XmChangeColor((Widget)m_curSurf->getHook(),App->RGBToPixel(r,g,b));
     redraw();
  }
}

void FlatTracker::setSelectedWidth()
{
  if(!m_curSurf || !m_curFlat) {
     App->ShowMessage("Please select a track first");
     return;
  }

  float l,nl;
  char tstr[100];
  m_curSurf->getLineWidth(l);
retry_lw:
  sprintf(tstr,"%.2f",l);
  if(App->GetText("New Line Width:",tstr,99)) {
       if((sscanf(tstr,"%f",&nl)!=1)||(nl<=0.)||(nl>50)) {
           App->ShowMessage("Please enter a line width > 0 and < 50");
           goto retry_lw;
       }
       m_curSurf->setLineWidth(nl);
       m_curFlat->setLineWidth(nl);
       redraw();
  }
}


void FlatTracker::selectButton(Widget swid)
{
VIN("void FlatTracker::selectButton(Widget swid)")
  Widget       wid,frm;
  GLPrimitive *p;

  // swid is the button, want the parent form
  //frm = XtParent(swid);

  m_curSurf = findPath(m_surfCurves,swid);
  m_curFlat = findPath(m_flatCurves,swid);

  // make sure both valid
  if(!m_curSurf && m_curFlat)
      m_curFlat = NULL;
  else if(!m_curFlat && m_curSurf)
      m_curSurf = NULL;

  redraw();
VOUT("void FlatTracker::selectButton(Widget swid)")
}

void FlatTracker::addPath(Array1D<int> &P, char *name)
{
VIN("void FlatTracker::addPath(Array1D<int> &P, char *name)")
  int i;
  Arg arglist[20];
  Cardinal c;

  Dimension bh = App->BigFontHeight() + 8;
  //
  // add a button to the row/col
  //
  ATTACHMENTS_TLRH(0,0,100,bh)
  Widget form = XmCreateForm(m_curveRC,"ButtonForm",arglist,c);
  Set(XmNresizeHeight,False);
  Set(XmNresizeWidth,False);
  XtManageChild(form);
  XmChangeColor(form,T1_BACK);

  ATTACHMENTS_TBLR(0,100,0,15)
  Widget btn = App->Button(form,"",App->SmallFontList(),
                           SelectButtonCB,(void*)this,arglist,c);
  XmChangeColor(btn,T1_BUTTON);

  ATTACHMENTS_TBLR(0,100,15,100)
  Widget txt = App->Text(form,name,128,
                         NULL,NULL,App->SmallFontList(),arglist,c);
  //
  // now add new path to internal list
  //
  Array1D<int> *newp = new Array1D<int>[m_numPaths+1];

  if(!newp) {
     App->ShowMessage("ERROR: Out of memory. Save everything and re-start");
     return;
  }

  for(i=0;i<m_numPaths;i++)
      newp[i] = m_paths[i];

  if(m_paths) delete [] m_paths;
  m_paths = newp;
  m_paths[m_numPaths] = P;
  // 
  // Remove starting spheres from views
  //
  if(m_surfStart) {
     m_surfCurves = GLPrimitive::removeFromList(m_surfCurves,m_surfStart);
     m_surfStart  = NULL;
  }
  if(m_flatStart) {
     m_flatCurves = GLPrimitive::removeFromList(m_flatCurves,m_flatStart);
     m_flatStart  = NULL;
  }
  //
  // add the GLPrimtive curve to the views
  //
  int idx,np;
  double sx,sy,sz,fx,fy,fz;
  Array1D<Point> sP,fP;
  GLPrimitive *prim;

  np = P.getNelm();
  sP.setDim(np);
  fP.setDim(np);

  if(sP.isEmpty() || fP.isEmpty()) {
     App->ShowMessage("ERROR: out of memory or zero-length path");
     return;
  }

  for(i=0;i<P.getNelm();i++) {
      idx = P[i];
      if((idx<0)||(idx>=m_numPoints))
         continue;
      sx = m_surface->vertices()[idx].x(); 
      sy = m_surface->vertices()[idx].y(); 
      sz = m_surface->vertices()[idx].z(); 
      sP[i].set(sx,sy,sz);

      fx = m_flat->vertices()[idx].x(); 
      fy = m_flat->vertices()[idx].y(); 
      fz = m_flat->vertices()[idx].z(); 
      fP[i].set(fx,fy,fz);
  }

  prim = new GLPrimitive(sP);
  prim->setColor(CURVE_R,CURVE_G,CURVE_B);
  prim->setLineWidth(m_lineWidth);
  prim->setHook((void*)btn);
  prim->setID(m_numPaths);
  m_surfCurves = GLPrimitive::addToList(m_surfCurves,prim);
  m_surfbox->setPrimitives(m_surfCurves);

  XmChangeColor(btn,App->RGBToPixel(CURVE_R,CURVE_G,CURVE_B));

  prim = new GLPrimitive(fP);
  prim->setColor(CURVE_R,CURVE_G,CURVE_B);
  prim->setLineWidth(m_lineWidth);
  prim->setHook((void*)btn);
  prim->setID(m_numPaths);
  m_flatCurves = GLPrimitive::addToList(m_flatCurves,prim);
  m_flatbox->setPrimitives(m_flatCurves);

  m_numPaths++;
  redraw();
VOUT("void FlatTracker::addPath(Array1D<int> &P, char *name)")
}

void FlatTracker::setTrackType(SurfaceTracker::TrackType t)
{
VIN("void FlatTracker::setTrackType(SurfaceTracker::TrackType t)")
 m_trackType = t;

  if(m_bUse3D)
     m_tracker.setSurface(m_surface,t);
  else
     m_tracker.setSurface(m_flat,t);
VOUT("void FlatTracker::setTrackType(SurfaceTracker::TrackType t)")
}

void FlatTracker::setView(FlatTrackViewType t)
{
VIN("void FlatTracker::setView(FlatTrackViewType t)")
  if(m_viewType == t)
     return;

  Arg arglist[10];
  Cardinal c;

  switch(t) {
  case FLAT_S:
     XtUnmapWidget(m_flatbox->Main());
     ATTACHMENTS_TBLR(1,94,1,79)
     XtSetValues(m_surfbox->Main(),arglist,c);
     XtMapWidget(m_surfbox->Main());
     break;
  case FLAT_F:
     XtUnmapWidget(m_surfbox->Main());
     ATTACHMENTS_TBLR(1,94,1,79)
     XtSetValues(m_flatbox->Main(),arglist,c);
     XtMapWidget(m_flatbox->Main());
     break;
  case FLAT_SF:
     ATTACHMENTS_TBLR(1,94,1,39)
     XtSetValues(m_surfbox->Main(),arglist,c);
     ATTACHMENTS_TBLR(1,94,40,79)
     XtSetValues(m_flatbox->Main(),arglist,c);
     XtMapWidget(m_surfbox->Main());
     XtMapWidget(m_flatbox->Main());
     break;
  default:
     return;
  }

  m_viewType = t;
VOUT("void FlatTracker::setView(FlatTrackViewType t)")
}

void FlatTracker::redraw()
{
VIN("void FlatTracker::redraw()")
  switch(m_viewType) {
  case FLAT_S:
      m_surfbox->redraw();
      break;
  case FLAT_F:
      m_flatbox->redraw();
      break;
  case FLAT_SF:
      m_surfbox->redraw();
      m_flatbox->redraw();
      break;
  }
VOUT("void FlatTracker::redraw()")
}

void FlatTracker::loadCurvature()
{
VIN("void FlatTracker::loadCurvature()")
  char fname[1024];
  fname[0] = 0;
  if(Browser->GetFilename("Select Curvature",fname,1024)) {
     if(m_surface->loadCurvature(fname) != ItXSuccess) {
        App->ShowMessage("ERROR: curvature failed to load. See stdout");
        return;
     }
     Array1D<double> D = m_surface->curvature();
     m_flat->setCurvature(D);
     m_lsurface->surfaceChanged();
     m_lflat->surfaceChanged();

     SurfaceColor col;
     col.colorSurface(m_lsurface->glsurface());
     m_lflat->glsurface()->setColorSurface(*(m_lsurface->glsurface()));
     redraw();
  }
VOUT("void FlatTracker::loadCurvature()")
}

void FlatTracker::base1Toggle()
{
  m_bBase1 = (App->GetToggleValue(m_base1Toggle)==True);
}

void FlatTracker::use3DVertsToggle()
{
VIN("void FlatTracker::use3DVertsToggle()")
  m_bUse3D = (App->GetToggleValue(m_vertToggle)==True);
  setTrackType(m_trackType);
VOUT("void FlatTracker::use3DVertsToggle()")
}

void FlatTracker::setLineWidth()
{
VIN("void FlatTracker::setLineWidth()")
  char tstr[100];
  int  lw;

retry_lw:
  sprintf(tstr,"%d",m_lineWidth);
  if(!App->GetText("Enter Line Width",tstr,99))
      return;

  if((sscanf(tstr,"%d",&lw) != 1)||(lw<0)||(lw>20)) {
    App->ShowMessage("ERROR: line width must be between 1..20");
    goto retry_lw;
  }

  if(lw == m_lineWidth)
    return;

  m_lineWidth = lw;

  GLPrimitive *p;
  for(p=m_surfCurves;p;p=p->Next())
      p->setLineWidth(m_lineWidth);

  for(p=m_flatCurves;p;p=p->Next())
      p->setLineWidth(m_lineWidth);

  redraw();
VOUT("void FlatTracker::setLineWidth()")
}

void FlatTracker::clearPaths()
{
VIN("void FlatTracker::clearPaths()")
  GLPrimitive *p;
  for(p=m_surfCurves;p;p=p->Next())
      XtDestroyWidget(XtParent((Widget)p->getHook()));

  GLPrimitive::clearList(m_surfCurves); 
  GLPrimitive::clearList(m_flatCurves); 
  m_surfCurves = NULL;
  m_flatCurves = NULL;
  m_surfStart  = NULL;
  m_flatStart  = NULL;

  m_curSurf = NULL;
  m_curFlat = NULL;

  m_surfbox->setPrimitives(NULL);
  m_flatbox->setPrimitives(NULL);
  redraw();

  for(int i=0;i<m_numPaths;i++)
      m_paths[i].setDim(0);

  if(m_paths) delete [] m_paths;

  m_paths    = NULL;
  m_numPaths = 0;
VOUT("void FlatTracker::clearPaths()")
}


void FlatTracker::doClearPaths()
{
VIN("void FlatTracker::doClearPaths()")
  if((m_numPaths > 0)&&(App->GetYesNo("Delete ALL paths. Are you sure?")))
     clearPaths();
VOUT("void FlatTracker::doClearPaths()")
}

void FlatTracker::deletePath()
{
VIN("void FlatTracker::deletePath()")
  if(m_numPaths < 1) {
     App->ShowMessage("No paths to delete");
     return;
  }

  if(!m_curSurf || !m_curFlat) {
     App->ShowMessage("Please select a path to delete");
     return;
  }

  // destroy buttons/text
  XtDestroyWidget(XtParent((Widget)m_curSurf->getHook()));

  int i,idx = m_curSurf->getID();
  if((idx<0)||(idx>=m_numPaths))
      return;
  
  m_paths[idx].setDim(0);
  for(i=idx;i<(m_numPaths-1);i++)
      m_paths[i] = m_paths[i+1];

  m_numPaths--;

  m_surfCurves = GLPrimitive::removeFromList(m_surfCurves,m_curSurf);
  m_flatCurves = GLPrimitive::removeFromList(m_flatCurves,m_curFlat);

  m_surfbox->setPrimitives(m_surfCurves);
  m_flatbox->setPrimitives(m_flatCurves);
  m_curSurf = NULL;
  m_curFlat = NULL;
  redraw();
VOUT("void FlatTracker::deletePath()")
}


//
// Parse a line for a number followed by a string (which can contain ws)
//
void FlatTracker::parseName(char *line, int &np, char *cname)
{
VIN("void FlatTracker::parseName(char *line, int &np, char *cname)")
  int i,l;
  sscanf(line,"%d",&np);
  l = strlen(line);

  for(i=0;i<l;i++)
    if((line[i]<'0')||(line[i]>'9'))
       break;

  while((i<l)&&(line[i] == ' ')) i++;

  if((i>=l)||(strlen(&(line[i])) < 2))
     sprintf(cname,"<no-name>");
  else
     sprintf(cname,"%s",&(line[i]));
VOUT("void FlatTracker::parseName(char *line, int &np, char *cname)")
}


void FlatTracker::loadPaths()
{
  char fname[1024];
  fname[0] = 0;
  if(!Browser->GetFilename("Load Curves From",fname,1024))
      return;

  ifstream fp(fname);
  if(!fp) {
     App->ShowMessage("ERROR: cannot open curves file");
     return;
  }

  clearPaths();

  int i,j,nl,np,elm,adder;
  char line[1024];
  char cname[200];
  Array1D<int> P;

  adder = (m_bBase1) ? 1 : 0;
 
  fp >> nl;
  if((nl<=0)||(nl>1024)) {
     App->ShowMessage("ERROR: this does not appear to be a curv file");
     return;
  }

  for(i=0;i<nl;i++) {
     line[0] = 0;
     while(strlen(line)<2)
        fp.getline(line,1024);
     parseName(line,np,cname);
     if((np<=0)||(np>=m_numPoints)) {
        App->ShowMessage("ERROR: this does not appear to be a curv file");
        clearPaths();
        return;
     }
     P.setDim(np);
     for(j=0;j<np;j++) {
        fp >> elm;
        elm -= adder;
        if((elm<0)||(elm>=m_numPoints)) {
          App->ShowMessage("ERROR: this path file is not for these surfaces");
          clearPaths();
          return;
        }
        P[j] = elm;
     }
     addPath(P,cname);
  }
  
  GLPrimitive *p,*p1,*p2;
  float r,g,b,l;
  for(i=0;i<nl;i++) {

    p1 = p2 = NULL;
    for(p=m_surfCurves;p;p=p->Next()) 
      if(p->getID() == i) 
       {  p1 = p; break; }

    for(p=m_flatCurves;p;p=p->Next()) 
      if(p->getID() == i) 
       {  p2 = p; break; }

    if(p1 && p2) {
       fp >> r >> g >> b >> l;
       if(!fp.fail() && !fp.eof()) {
          p1->setColor(r,g,b);
          p2->setColor(r,g,b);
          p1->setLineWidth(l);
          p2->setLineWidth(l);
          XmChangeColor((Widget)p1->getHook(),App->RGBToPixel(r,g,b));
       }
    }
  }
  fp.close();
  redraw();
}

void FlatTracker::savePaths()
{
  if(m_numPaths <= 0) {
     App->ShowMessage("ERROR: no paths to save");
     return;
  }

  char fname[1024];
  fname[0] = 0;
  if(!Browser->GetFilename("Save Curves To",fname,1024))
      return;

  ofstream fp(fname);
  if(!fp) {
     App->ShowMessage("ERROR: cannot open curve file for writing");
     return;
  }
  int i,j,np,adder;
  char cname[1024];

  adder = (m_bBase1) ? 1 : 0;

  fp << m_numPaths << endl;;
  for(i=0;i<m_numPaths;i++) {
      np = m_paths[i].getNelm(); 
      getName(i,cname);
      fp << np << " " << cname << endl;
      for(j=0;j<np;j++) {
          fp << (m_paths[i][j]+adder) << " ";
          if(j % 10 == 9) fp << endl;
      }
      fp << endl;
  }

  GLPrimitive *p;
  float r,g,b,l;
  for(i=0;i<m_numPaths;i++) {
    for(p=m_surfCurves;p;p=p->Next()) 
      if(p->getID() == i) 
        break;
    if(p) {
      p->getColor(r,g,b);
      p->getLineWidth(l);
      fp << r << " " << g << " " << b << " " << l << endl;
    } else {
      fp << "1.0 0.0 0.0 5.0" << endl;
    }
  }

  fp.close();
}

void FlatTracker::getName(int idx, char *cname)
{
VIN("void FlatTracker::getName(int idx, char *cname)")
  GLPrimitive *p;
  Widget       wid,txt;

  for(p=m_surfCurves;p;p=p->Next()) 
    if(p->getID() == idx) {
       wid = XtParent((Widget)p->getHook());  
       txt = XtNameToWidget(wid,"TextWidget"); 
       if(!txt) 
         sprintf(cname,"<no-name>");
       else
         App->GetString(txt,cname,200);
       return;
    }

  sprintf(cname,"<no-name>");
VOUT("void FlatTracker::getName(int idx, char *cname)")
}




