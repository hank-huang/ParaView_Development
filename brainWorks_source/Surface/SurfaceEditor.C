///////////////////////////////////////////////////////////////////////////
//
// File: SurfaceEditor.C
//
// Author: Keith Doolittle
//
// 0Purpose:
//
///////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <sys/param.h>
#include <KDApplication.h>
#include <TDL/Surface.h>
#include <TDL/ByuSurface.h>
#include <SurfaceBox.h>
#include <LoadedSurface.h>
#include <KDShell.h>
#include <Icon.h>
#include <IconList.h>
#include <GLPrimitive.h>
#include <SurfaceEditor.h>
// Lei 11/17/2003: added two below
#include <Main.h>
#include <MainUtils.h>

#define NCUT_TYPES 2
static char *cut_names[NCUT_TYPES] = {
	"Geodesic",
	"Window Line",
};

static enum SurfaceEditor::CutType cut_types[NCUT_TYPES] = {
	SurfaceEditor::Geodesic,
	SurfaceEditor::WindowLine,
};


#define NGEO_TYPES 3
static char *geo_names[NGEO_TYPES] = {
	"Sulci",
	"Gyri",
	"Geodesic",
};

static enum SurfaceTracker::TrackType geo_types[NGEO_TYPES] = {
	SurfaceTracker::Sulci,
	SurfaceTracker::Gyri,
	SurfaceTracker::Geodesi,
};


DEFINECB(SurfaceEditor,CutCB,cut)
DEFINECBP(SurfaceEditor,DelRedCB,delPoints,1)
DEFINECBP(SurfaceEditor,DelBlueCB,delPoints,0)
DEFINECB(SurfaceEditor,SaveCB,save)
DEFINECBP(SurfaceEditor,ClearPathCB,clearPath,0)
DEFINECB(SurfaceEditor,ClosePathCB,closePath)
DEFINECB(SurfaceEditor,DelSegCB,deleteSegment)
// Lei 04/17/2003: add two below 
DEFINECBP(SurfaceEditor,HideRedCB,hidePoints,1)
DEFINECBP(SurfaceEditor,HideBlueCB,hidePoints,0)
// Lei 11/11/2003: add below 
DEFINECB(SurfaceEditor,LoadPathCB,loadPath)
DEFINECB(SurfaceEditor,SavePathCB,savePath)
// Lei 12/03/2003: add below
DEFINECBP(SurfaceEditor,ColorPathRedCB,colorPath,1) // red
DEFINECBP(SurfaceEditor,ColorPathBlueCB,colorPath,0) // blue

static void PickPointCB(u_long dat, void *s, int indx)
{
VIN("static void PickPointCB(u_long dat, void *s, int indx)")
  if(dat) ((SurfaceEditor*)dat)->pickPoint(indx);
VOUT("static void PickPointCB(u_long dat, void *s, int indx)")
}

static void RedrawCB(u_long dat)
{
VIN("static void RedrawCB(u_long dat)")
  if(dat) ((SurfaceEditor*)dat)->surfaceBoxRedraw();
VOUT("static void RedrawCB(u_long dat)")
}

static void StartClickCB(Widget w, XEvent *ev, String *parms, Cardinal *np)
{
VIN("static void StartClickCB(Widget w, XEvent *ev, String *parms, Cardinal *np)")
  if(*np != 1) return;

  u_long ll;
  sscanf(parms[0],"%ld",&ll);

  int x = ((XButtonEvent*)ev)->x;
  int y = ((XButtonEvent*)ev)->y;

  if(ll) ((SurfaceEditor*)ll)->mouseDown(w,x,y);
VOUT("static void StartClickCB(Widget w, XEvent *ev, String *parms, Cardinal *np)")
}

static void MotionCB(Widget w, XEvent *ev, String *parms, Cardinal *np)
{
VIN("static void MotionCB(Widget w, XEvent *ev, String *parms, Cardinal *np)")
  if(*np != 1) return;

  u_long ll;
  sscanf(parms[0],"%ld",&ll);

  int x = ((XMotionEvent*)ev)->x;
  int y = ((XMotionEvent*)ev)->y;

  if(ll) ((SurfaceEditor*)ll)->mouseMotion(w,x,y);
VOUT("static void MotionCB(Widget w, XEvent *ev, String *parms, Cardinal *np)")
}


static void EndClickCB(Widget w, XEvent *ev, String *parms, Cardinal *np)
{
VIN("static void EndClickCB(Widget w, XEvent *ev, String *parms, Cardinal *np)")
  if(*np != 1) return;

  u_long ll;
  sscanf(parms[0],"%ld",&ll);

  int x = ((XButtonEvent*)ev)->x;
  int y = ((XButtonEvent*)ev)->y;

  if(ll) ((SurfaceEditor*)ll)->mouseUp(w,x,y);
VOUT("static void EndClickCB(Widget w, XEvent *ev, String *parms, Cardinal *np)")
}


static void CutCB(char *nm, void *dat)
{
VIN("static void CutCB(char *nm, void *dat)")
  SurfaceEditor *ed = (SurfaceEditor*)dat;
  if(!ed) return;

  for(int i=0;i<NCUT_TYPES;i++)
      if(!strcmp(nm,cut_names[i]))
	ed->setCutType(cut_types[i]);
VOUT("static void CutCB(char *nm, void *dat)")
}


static void GeoCB(char *nm, void *dat)
{
VIN("static void GeoCB(char *nm, void *dat)")
  SurfaceEditor *ed = (SurfaceEditor*)dat;
  if(!ed) return;

  for(int i=0;i<NGEO_TYPES;i++)
      if(!strcmp(nm,geo_names[i]))
	ed->setTrackType(geo_types[i]);
VOUT("static void GeoCB(char *nm, void *dat)")
}



static int seactions_added = 0;
static XtActionsRec seactions[] = {
  { "seMouseDown",   StartClickCB },
  { "seMouseMotion", MotionCB },
  { "seMouseUp",     EndClickCB },
};


SurfaceEditor::SurfaceEditor(LoadedSurface *lsurf) :
	KDChildWindow(App->FilenameTail(lsurf->Filename()),
		App->Width(),App->Height()) // Lei 11/01/2003: use full size
{
VIN("SurfaceEditor::SurfaceEditor(LoadedSurface *lsurf) :")
   if(!seactions_added) {
	XtAppAddActions(App->AppContext(),seactions,XtNumber(seactions));
	seactions_added = 1;
	}

   m_loadedSurface = NULL;
   m_surface       = NULL;
   m_icon          = NULL;
   m_surfBox       = NULL;
   m_cutList       = NULL;
   m_bModified     = false;
   m_segments      = NULL;

   m_button1x      = 20;
   m_button1y      = 20;
   m_button2x      = 70;
   m_button2y      = 90;

   if(!lsurf || !lsurf->surface() || (lsurf->surface()->getNumVert() < 1)) {
	App->ShowMessage("ERROR: Cannot edit this surface (empty)");
	return;
	}

   m_surface       = new Surface(*(lsurf->surface()));
   m_loadedSurface = new LoadedSurface(m_surface);
   m_icon          = new Icon(m_loadedSurface);

   m_isRed.setDim(m_surface->getNumVert());
   m_isRed = 0U;

   Cardinal c;
   Arg arglist[30];
   char tstr[1024];

   c = 0;
   ATTACHMENTS_TBLR(2,98,2,75)
   m_surfBox = new SurfaceBox(this->Form(),
			      m_loadedSurface,arglist,c);
   m_surfBox->SetRedrawFunction(RedrawCB,(u_long)this);

   m_bdim = 20;
   c = 0;
   Set(XmNwidth,m_bdim)
   Set(XmNheight,m_bdim)
   Set(XmNresizeWidth,False)
   Set(XmNresizeHeight,False)
   Set(XmNmappedWhenManaged,False)
   Set(XmNtopAttachment,XmATTACH_NONE)
   Set(XmNbottomAttachment,XmATTACH_NONE)
   Set(XmNrightAttachment,XmATTACH_NONE)
   Set(XmNleftAttachment,XmATTACH_NONE)
   m_button1 = App->Label(this->Form(),"",App->MedFontList(),arglist,c);
   m_button2 = App->Label(this->Form(),"",App->MedFontList(),arglist,c);
   XmChangeColor(m_button1,RED);
   XmChangeColor(m_button2,RED);

   // Now, override mouse controls on the buttons
   sprintf(tstr,"<Btn1Down>:seMouseDown(%ld)\n\
		 <Btn1Motion>:seMouseMotion(%ld)\n\
		 <Btn1Up>:seMouseUp(%ld)",
		(u_long)this,(u_long)this,(u_long)this);

   XtTranslations ft = XtParseTranslationTable(tstr);
   XtOverrideTranslations(m_button1,ft);
   XtOverrideTranslations(m_button2,ft);

   int lh  = App->BigFontHeight() + 8;
   Widget btn;

   ATTACHMENTS_TLRH(0,77,98,lh)
   m_cutList = new IconList(this->Form(),"Cut Type",IXNoIcon,
		            cut_names,NCUT_TYPES,arglist,c);
   m_cutList->AddExtraCallback(CutCB,(void*)this);
   btn = m_cutList->Form();


   ATTACHMENTS_WLRH(btn,77,98,lh)
   m_geoList = new IconList(this->Form(),"Geo. Type",IXNoIcon,
		            geo_names,NGEO_TYPES,arglist,c);
   m_geoList->AddExtraCallback(GeoCB,(void*)this);
   btn = m_geoList->Form();


// Lei 11/17/2003: below two
   ATTACHMENTS_WLRH(btn,77,98,lh)
   btn = App->Button(this->Form(),"Load Path",App->MedFontList(),
			LoadPathCB,(XtPointer)this,arglist,c);
   ATTACHMENTS_WLRH(btn,77,98,lh)
   btn = App->Button(this->Form(),"Save Path",App->MedFontList(),
			SavePathCB,(XtPointer)this,arglist,c);

   ATTACHMENTS_WLRH(btn,77,98,lh)
   btn = App->Button(this->Form(),"Clear Path",App->MedFontList(),
			ClearPathCB,(XtPointer)this,arglist,c);
   m_W1 = btn;

   ATTACHMENTS_WLRH(btn,77,98,lh)
   btn = App->Button(this->Form(),"Close Path",App->MedFontList(),
			ClosePathCB,(XtPointer)this,arglist,c);
   m_W2 = btn;

   ATTACHMENTS_WLRH(btn,77,98,lh)
   btn = App->Button(this->Form(),"Del Last Segment",App->MedFontList(),
			DelSegCB,(XtPointer)this,arglist,c);
   m_W3 = btn;
   m_deleteButton = btn;

   ATTACHMENTS_WLRH(btn,77,98,lh)
   btn = App->Button(this->Form(),"Show Cut",App->MedFontList(),
			CutCB,(XtPointer)this,arglist,c);
// Lei 12/03/2003: add two below
   ATTACHMENTS_WLRH(btn,77,98,lh)
   btn = App->Button(this->Form(),"Color Path/BLUE",App->MedFontList(),
			ColorPathBlueCB,(XtPointer)this,arglist,c);
   ATTACHMENTS_WLRH(btn,77,98,lh)
   btn = App->Button(this->Form(),"Color Path/RED",App->MedFontList(),
			ColorPathRedCB,(XtPointer)this,arglist,c);
// Lei 12/03/2003: add two above

   ATTACHMENTS_WLRH(btn,77,98,lh)
   btn = App->Button(this->Form(),"Delete BLUE",App->MedFontList(),
			DelBlueCB,(XtPointer)this,arglist,c);
   ATTACHMENTS_WLRH(btn,77,98,lh)
   btn = App->Button(this->Form(),"Delete RED",App->MedFontList(),
			DelRedCB,(XtPointer)this,arglist,c);

// Lei 04/17/2003: the below two
   ATTACHMENTS_WLRH(btn,77,98,lh)
   btn = App->Button(this->Form(),"Hide BLUE",App->MedFontList(),
			HideBlueCB,(XtPointer)this,arglist,c);
   ATTACHMENTS_WLRH(btn,77,98,lh)
   btn = App->Button(this->Form(),"Hide RED",App->MedFontList(),
			HideRedCB,(XtPointer)this,arglist,c);
   ATTACHMENTS_WLRH(btn,77,98,lh)
// Lei 04/17/2003: the above two

// Lei 11/17/2003: changed "Save" to "Save Surface"
   ATTACHMENTS_WLRH(btn,77,98,lh)
   btn = App->Button(this->Form(),"Save Surface",App->MedFontList(),
			SaveCB,(XtPointer)this,arglist,c);

   m_cutType = cut_types[NCUT_TYPES-1];
   setCutType(cut_types[0]);

   m_tracker.setTrackType(geo_types[0]);

   colorSurface(0);
VOUT("SurfaceEditor::SurfaceEditor(LoadedSurface *lsurf) :")
}


void SurfaceEditor::setButtonPosition(Widget b, int bs, int bx, int by)
{
VIN("void SurfaceEditor::setButtonPosition(Widget b, int bs, int bx, int by)")
  int x,y;

  switch(bs) {
  case 0: // top
     x = m_boxx + bx - m_bdim/2;
     y = m_boxy - m_bdim;
     break;
  case 1: // bottom
     x = m_boxx + bx - m_bdim/2;
     y = m_boxy + m_boxh; 
     break;
  case 2: // left
     x = m_boxx - m_bdim;
     y = m_boxy + by - m_bdim/2;
     break;
  case 3: // right
     x = m_boxx + m_boxw;
     y = m_boxy + by - m_bdim/2;
     break;
  default:
     return;
  }

  App->SetPosition(b,x,y);
VOUT("void SurfaceEditor::setButtonPosition(Widget b, int bs, int bx, int by)")
}


SurfaceEditor::~SurfaceEditor()
{
VIN("SurfaceEditor::~SurfaceEditor()")
   if(m_surfBox) 
	delete m_surfBox;
   if(m_cutList)
	delete m_cutList;
   if(m_geoList)
	delete m_geoList;
VOUT("SurfaceEditor::~SurfaceEditor()")
}

bool SurfaceEditor::AllowClose()
{
VIN("bool SurfaceEditor::AllowClose()")
  if(m_bModified && m_icon && App->GetYesNo("Save Edited Surface First?"))
     m_icon->Save();

  colorSurface(0);

VOUT("bool SurfaceEditor::AllowClose()")
  return(true);
}


void SurfaceEditor::setCutType(CutType t)
{
VIN("void SurfaceEditor::setCutType(CutType t)")
  if((t == m_cutType)||(!m_surfBox))
	return;

  m_cutType = t;
  switch(m_cutType) {
  case Geodesic:
       XtSetSensitive(m_W1,true);
       XtSetSensitive(m_W2,true);
       XtSetSensitive(m_W3,true);
       App->SetVisible(m_button1,false);
       App->SetVisible(m_button2,false);
       clearPath(1);
       m_surfBox->SetPickStyle(PickPoint,PickPointCB,(u_long)this);
       m_tracker.setSurface(m_surface,SurfaceTracker::SameType);
       break;
  case WindowLine:
       XtSetSensitive(m_W1,false);
       XtSetSensitive(m_W2,false);
       XtSetSensitive(m_W3,false);
       App->GetPosition(m_surfBox->Main(),m_boxx,m_boxy);
       m_boxw = m_surfBox->getWidth();
       m_boxh = m_surfBox->getHeight();
       m_button1x = m_boxw/2; m_button1y = 0;        m_bside1 = 0;
       m_button2x = m_boxw/2; m_button2y = m_boxh-1; m_bside2 = 1;
       setButtonPosition(m_button1,m_bside1,m_button1x,m_button1y);
       setButtonPosition(m_button2,m_bside2,m_button2x,m_button2y);
       App->SetVisible(m_button1,true);
       App->SetVisible(m_button2,true);
       clearPath(1);
       m_surfBox->SetPickStyle(PickNone);
       break;
  }
  m_surfBox->redraw(); 
VOUT("void SurfaceEditor::setCutType(CutType t)")
}

//
// This is only called when m_cutType == WindowLine.
// Prepare for moving buttons around
//
void SurfaceEditor::mouseDown(Widget w, int x, int y)
{
VIN("void SurfaceEditor::mouseDown(Widget w, int x, int y)")
  App->GetPosition(m_surfBox->Main(),m_boxx,m_boxy);
  m_boxw = m_surfBox->getWidth();
  m_boxh = m_surfBox->getHeight();
//  App->GetDimensions(m_surfBox->Main(),m_boxw,m_boxh);
VOUT("void SurfaceEditor::mouseDown(Widget w, int x, int y)")
}


//
// Move the button (either m_button1 or m_button2) around
// the sides of the m_surfBox main window
//
void SurfaceEditor::mouseMotion(Widget w, int x, int y)
{
VIN("void SurfaceEditor::mouseMotion(Widget w, int x, int y)")
  bool restrictx;
  int  *xpos,*ypos,*sid;

  if(w == m_button1) {
	xpos = &m_button1x;	
	ypos = &m_button1y;	
        sid  = &m_bside1;
  } else if(w == m_button2) {
	xpos = &m_button2x;	
	ypos = &m_button2y;	
        sid  = &m_bside2;
  } else return;

  switch(*sid) {
  case 0: // top
     *xpos = *xpos + x;
     if(*xpos < 0) { // wrap around to left
	*ypos = -(*xpos);
        *xpos = 0;
	*sid  = 2;
     } else if(*xpos >= m_boxw) { // wrap around to right
	*ypos = *xpos - m_boxw;
        *xpos = m_boxw - 1;
	*sid  = 3;
     }
     break;
  case 1: // bottom
     *xpos = *xpos + x;
     if(*xpos < 0) { // wrap around to left
	*ypos = m_boxh + (*xpos) - 1;
        *xpos = 0;
	*sid  = 2;
     } else if(*xpos >= m_boxw) { // wrap around to right
	*ypos = m_boxh - (*xpos - m_boxw);
        *xpos = m_boxw - 1;
	*sid  = 3;
     }
     break;
  case 2: // left
     *ypos = *ypos + y;
     if(*ypos < 0) { // wrap around to top
        *xpos = - (*ypos); 
	*ypos = 0;
        *sid  = 0;
     } else if(*ypos >= m_boxh) { // wrap around to bottom
        *xpos = *ypos - m_boxh;
        *ypos = m_boxh - 1;
        *sid  = 1;
     }
     break;
  case 3: // right
     *ypos = *ypos + y;
     if(*ypos < 0) { // wrap around to top
        *xpos = m_boxw + (*ypos) - 1; 
	*ypos = 0;
        *sid  = 0;
     } else if(*ypos >= m_boxh) { // wrap around to bottom
        *xpos = m_boxw - (*ypos - m_boxh);
        *ypos = m_boxh - 1;
        *sid  = 1;
     }
     break;
  default:
     return;
  }

  if(w == m_button1)
	setButtonPosition(m_button1,m_bside1,m_button1x,m_button1y);
  else if(w == m_button2)
	setButtonPosition(m_button2,m_bside2,m_button2x,m_button2y);
VOUT("void SurfaceEditor::mouseMotion(Widget w, int x, int y)")
}


void SurfaceEditor::mouseUp(Widget w, int x, int y)
{
VIN("void SurfaceEditor::mouseUp(Widget w, int x, int y)")
  m_surfBox->redraw();
VOUT("void SurfaceEditor::mouseUp(Widget w, int x, int y)")
}


// Cut actually only colors the surface
void SurfaceEditor::cut()
{
VIN("void SurfaceEditor::cut()")
   if(!m_surfBox || !m_surface)
	return;

   switch(m_cutType) {
   case Geodesic:
	labelByGeodesic();
	break;
   case WindowLine:
	labelByWindowLine();
	break;
   default:
	return;
   }
   colorSurface(1);
VOUT("void SurfaceEditor::cut()")
}

void SurfaceEditor::delPoints(int isred)
{
VIN("void SurfaceEditor::delPoints(int isred)")
  if(!m_surfBox || !m_surface)
	return;

  Surface newS;

  int i,a,b,c,cnt;

  int delv,nv = m_surface->getNumVert();
  int delf,nf = m_surface->getNumPoly();

  Array1D<int> newIndex(nv);

  delv = 0;
  for(i=0;i<nv;i++) {
        newIndex[i] = i - delv; // KWD
	if(m_isRed[i] == isred)
		delv++;
        }

  delf = 0;
  for(i=0;i<nf;i++) {
	a = m_surface->facets()[i][0];
	b = m_surface->facets()[i][1];
	c = m_surface->facets()[i][2];
	if((m_isRed[a] == isred)||(m_isRed[b] == isred)||(m_isRed[c] == isred))
		delf++;
  }
  
  if((delv <= 0)||(delf <= 0))
	return;

  if(((nv - delv) <= 0)||((nf - delf) <= 0)) {
	App->ShowMessage("Doing this will delete the entire surface",1);
	return;
  }
 
  newS.setSize(nv - delv,nf - delf);

  cnt = 0;
  for(i=0;i<nv;i++)
	if(m_isRed[i] != isred) 
	   newS.vertices()[cnt++] = m_surface->vertices()[i];

  cnt = 0;
  for(i=0;i<nf;i++) {
	a = m_surface->facets()[i][0];
	b = m_surface->facets()[i][1];
	c = m_surface->facets()[i][2];
	if((m_isRed[a] != isred)&&(m_isRed[b] != isred)&&(m_isRed[c] != isred)) {
		a = newIndex[a];
		b = newIndex[b];
		c = newIndex[c];
		if((a<0)||(a>=newS.getNumVert()))
			cerr << "A BAD" << endl;
		if((b<0)||(b>=newS.getNumVert()))
			cerr << "B BAD" << endl;
		if((c<0)||(c>=newS.getNumVert()))
			cerr << "C BAD" << endl;
		newS.facets()[cnt][0] = a;
		newS.facets()[cnt][1] = b;
		newS.facets()[cnt][2] = c;
		cnt++;
		}
  }

  *m_surface = newS;
  m_loadedSurface->surfaceChanged();
  m_surface->freeCurvature();

  m_isRed.setDim(m_surface->getNumVert());
  m_isRed = 0U;
  colorSurface(0);
  m_bModified = true;

  if(m_cutType == SurfaceEditor::Geodesic) {
     clearPath(1);
     m_tracker.setSurface(m_surface,SurfaceTracker::SameType);
  }
VOUT("void SurfaceEditor::delPoints(int isred)")
}


void SurfaceEditor::save()
{
VIN("void SurfaceEditor::save()")
  if(m_icon) {
	m_icon->Save();
	m_bModified = false;
	}
VOUT("void SurfaceEditor::save()")
}


void SurfaceEditor::surfaceBoxRedraw()
{
VIN("void SurfaceEditor::surfaceBoxRedraw()")
  int w,h;

  if(!m_surfBox)
	return;

  switch(m_cutType) {
  case WindowLine:
	m_surfBox->renderLine(m_button1x,m_button1y,
			      m_button2x,m_button2y,3,255U,0U,0U);
	break;
  case Geodesic:
	// do nothing
	break;
  default:
	break;
  }
VOUT("void SurfaceEditor::surfaceBoxRedraw()")
}


void SurfaceEditor::colorSurface(int isrb)
{
VIN("void SurfaceEditor::colorSurface(int isrb)")
  if(!m_surfBox) return;

  GLSurface *s = m_loadedSurface->glsurface();
  if(!s) return;

  int i,nv = s->getNumVert();
  if(nv <= 0) return;

  if(isrb) {
	for(int i=0;i<nv;i++)
		if(m_isRed[i])
		     s->setVertexColor(i,1.f,0.f,0.f);
		else s->setVertexColor(i,0.f,0.f,1.f);
  } else {
	s->setColorSurface(0.8f,0.8f,0.8f);
  }

  m_surfBox->redraw();
VOUT("void SurfaceEditor::colorSurface(int isrb)")
}


void SurfaceEditor::labelByGeodesic()
{
VIN("void SurfaceEditor::labelByGeodesic()")
  int i,j,a,b,c;
  int np = m_Path.getNelm();
  int i1 = m_Path[0];
  int i2 = m_Path[np-1];
  int i3;

  if(i1 != i2) {
	App->ShowMessage("Path is not closed",1);
	return;
	}

  int nv = m_surface->getNumVert();

  Array1D<u_char> inPath(nv);
  inPath = 0U;

  // Tag vertices that are along the path
  //
  // Make sure no path vertices are duplicated
  // (this is from a path back-tracking along itself)
  // which will cause problems with coloring below

  for(i=0;i<np;i++) {
        int v = m_Path[i];
        for(j=0;j<np;j++)
          if((j != i)&&(m_Path[j]==v))
              break;
     
        if(j >= np)
	   inPath[v] = 1U;
  }

  m_isRed.setDim(nv);
  m_isRed = 0U;

  // label vertices on 'right' side of path as red
  // (vertices in a triangle whose 2 other vertices
  //  are on the path)
  for(i=0;i<m_surface->getNumPoly();i++) {
	a = m_surface->facets()[i][0];
	b = m_surface->facets()[i][1];
	c = m_surface->facets()[i][2];
	if(!inPath[a] && !inPath[b] && !inPath[c])
		continue;

	for(j=0;j<np-1;j++) {
	    i1 = m_Path[j];
	    i2 = m_Path[j+1];
	    if((i1==a)&&(i2==b) && !inPath[c]) {
	        m_isRed[c] = 1U;
                break;
            }
	    else if((i1==b)&&(i2==c) && !inPath[a]) {
	        m_isRed[a] = 1U;
                break;
            }
	    else if((i1==c)&&(i2==a && !inPath[b])) {
	        m_isRed[b] = 1U;
                break;
            }
	}
  }

  //
  // now 'flood-fill' neighbor vertices
  // to fill red region without crossing path
  //
  int iter=0;
  int hit=1;
  while(hit) {
    if(++iter > 1000) break;
    hit = 0;
    for(i=0;i<m_surface->getNumPoly();i++) {
	a = m_surface->facets()[i][0];
	b = m_surface->facets()[i][1];
	c = m_surface->facets()[i][2];
	if(inPath[a] || inPath[b] || inPath[c])
		continue;

	if(m_isRed[a]&&(!m_isRed[b] || !m_isRed[c]))      
		{ m_isRed[b] = m_isRed[c] = 1U; hit = 1; }
	else if(m_isRed[b]&&(!m_isRed[a] || !m_isRed[c]))      
		{ m_isRed[a] = m_isRed[c] = 1U; hit = 1; }
	else if(m_isRed[c]&&(!m_isRed[a] || !m_isRed[b]))      
		{ m_isRed[a] = m_isRed[b] = 1U; hit = 1; }
    }
  }
VOUT("void SurfaceEditor::labelByGeodesic()")
}


//
// Color each surface vertex by which side of the
// window line it projects to
//
void SurfaceEditor::labelByWindowLine()
{
VIN("void SurfaceEditor::labelByWindowLine()")
   float  x0 = m_button1x; 
   float  y0 = m_button1y; 
   float  x1 = m_button2x; 
   float  y1 = m_button2y; 
   float  wx,wy,wz;
   float  p,sx,sy,sz;

   // this is a hack until SurfaceBox fixed

   m_isRed.setDim(m_surface->getNumVert());
   for(int i=0;i<m_surface->getNumVert();i++) {
	wx = m_surface->vertices()[i].x();
	wy = m_surface->vertices()[i].y();
	wz = m_surface->vertices()[i].z();
	m_surfBox->project(wx,wy,wz,sx,sy,sz);
	//
	// use dot product to find side of 2D line 
	// point projects to
	//
	p = (sy - y0)*(x1 - x0) - (sx - x0)*(y1 - y0);
	if(p < 0) m_isRed[i] = 1;
	else      m_isRed[i] = 0;
   }
VOUT("void SurfaceEditor::labelByWindowLine()")
}


void SurfaceEditor::clearPath(int force)
{
VIN("void SurfaceEditor::clearPath(int force)")
  if(force || App->GetYesNo("Delete entire path?")) {
	m_Path.setDim(0);
        GLPrimitive::clearList(m_segments);
        m_segStartIndex = -1;
        m_segments = NULL;
        m_surfBox->setPrimitives(m_segments);
        colorSurface(0);
	m_surfBox->redraw();
  }
VOUT("void SurfaceEditor::clearPath(int force)")
}


void SurfaceEditor::closePath()
{
VIN("void SurfaceEditor::closePath()")
  if(m_Path.getNelm() <= 1)
	return;

  int i1 = m_Path[0];
  int i2 = m_Path[m_Path.getNelm()-1];

  if(i1 == i2) {
	App->ShowMessage("Path is already closed");
	return;
	}
  addSegment(i2,i1);
VOUT("void SurfaceEditor::closePath()")
}


void SurfaceEditor::addSegment(int i1, int i2)
{
VIN("void SurfaceEditor::addSegment(int i1, int i2)")
  Array1D<int> p,newp;
  //
  // Tell tracker to track from i2 --> i1 since the path
  // it returns goes backwards 
  //
  if((m_tracker.findPath(i2,i1,p) != ItXSuccess)||(p.getNelm() <= 0)) {
      App->ShowMessage("could not calculate path (1)",1);
      return;
  }

  int nold = m_Path.getNelm()-1;
  int nnew = p.getNelm();
  int i,n  = nold + nnew;

  if(n <= 0) {
      App->ShowMessage("could not calculate path (2)",1);
      return;
  }
  m_segStartIndex = nold;

  //	
  // Append new path 'p' to list of old verts 'm_Path'
  //
  newp.setDim(n);
  for(i=0;i<nold;i++)
	newp[i] = m_Path[i];
  for(;i<n;i++)
	newp[i] = p[i - nold];
  m_Path.setDim(newp.getNelm());
  m_Path = newp;

  //
  // now add line segment to surface primitives
  //
  Array1D<Point> pts(nnew);
  for(i=0;i<nnew;i++)
    pts[i] = m_surface->vertices()[p[i]];

  GLPrimitive *seg = new GLPrimitive(pts);
  if(!seg) return;
  seg->setColor(0.f,0.f,1.f);
  seg->setLineWidth(3.f);
  m_segments = GLPrimitive::addToList(m_segments,seg);
  m_surfBox->setPrimitives(m_segments); 
  m_surfBox->redraw();
  XtSetSensitive(m_deleteButton,True);
VOUT("void SurfaceEditor::addSegment(int i1, int i2)")
}


void SurfaceEditor::deleteSegment()
{
VIN("void SurfaceEditor::deleteSegment()")
  if((m_segStartIndex <0)||(m_segStartIndex == (m_Path.getNelm()-1)))
	return;

  int i1 = m_segStartIndex;

  Array1D<int> newPath(m_segStartIndex);
  for(int i=0;i<m_segStartIndex;i++)
	newPath[i] = m_Path[i];


  // As it happens, the last primitive added to list is the front
 
  m_segments = GLPrimitive::removeFromList(m_segments,m_segments);
  m_surfBox->setPrimitives(m_segments);

  m_Path.setDim(m_segStartIndex);
  m_Path = newPath;
  m_surfBox->redraw();
  m_segStartIndex = m_Path.getNelm() - 1;

  XtSetSensitive(m_deleteButton,False);
VOUT("void SurfaceEditor::deleteSegment()")
}


void SurfaceEditor::pickPoint(int indx)
{
VIN("void SurfaceEditor::pickPoint(int indx)")
  if(!m_surfBox || !m_surface)
	return;

  // first click?
  if(m_segStartIndex < 0) {
       float x = m_surface->vertices()[indx].x();
       float y = m_surface->vertices()[indx].y();
       float z = m_surface->vertices()[indx].z();
       float r = m_loadedSurface->glsurface()->getMaxExtent() * 0.005;
       GLPrimitive *sph = new GLPrimitive(x,y,z,r);
       sph->setColor(1.f,0.f,0.f);
       m_segments = GLPrimitive::addToList(m_segments,sph);
       m_surfBox->setPrimitives(m_segments);
       m_segStartIndex = 0;
       m_Path.setDim(1);
       m_Path[0] = indx;
  } else {
       int i1 = m_Path[m_Path.getNelm()-1];
       int i2 = indx;
       addSegment(i1,i2);
  }
VOUT("void SurfaceEditor::pickPoint(int indx)")
}


// Lei 04/17/2003: added this one below
// keep the triangle indices, just show the extracted ones, however
void SurfaceEditor::hidePoints(int isred)
{
  if(!m_surfBox || !m_surface)
	return;

  Surface newS;

  int i,a,b,c,cnt;

// delv and delf are really hidev and hidef
  int delv,nv = m_surface->getNumVert();
  int delf,nf = m_surface->getNumPoly();

  Array1D<int> newIndex(nv);

  delv = 0;
  for(i=0;i<nv;i++) {
        newIndex[i] = i - delv; // KWD
	if(m_isRed[i] == isred)
		delv++;
        }

  delf = 0;
  for(i=0;i<nf;i++) {
	a = m_surface->facets()[i][0];
	b = m_surface->facets()[i][1];
	c = m_surface->facets()[i][2];
	if((m_isRed[a] == isred)||(m_isRed[b] == isred)||(m_isRed[c] == isred))
		delf++;
  }
  
  if((delv <= 0)||(delf <= 0))
	return;

  if(((nv - delv) <= 0)||((nf - delf) <= 0)) {
	App->ShowMessage("Doing this will hide the entire surface",1);
	return;
  }
 
  //  newS.setSize(nv - delv,nf - delf);
  // Lei: new surface has the same number of vertices and 
  //      number of facets to KEEP not to DELETE
  newS.setSize(nv,nf-delf);

  cnt = 0;
  for(i=0;i<nv;i++)
// 	if(m_isRed[i] != isred) 
// Lei: comment out the above to use all the original coordinates
	   newS.vertices()[cnt++] = m_surface->vertices()[i];

  cnt = 0;
  for(i=0;i<nf;i++) {
	a = m_surface->facets()[i][0];
	b = m_surface->facets()[i][1];
	c = m_surface->facets()[i][2];
	if((m_isRed[a] != isred)&&(m_isRed[b] != isred)&&(m_isRed[c] != isred)) {
	  //		a = newIndex[a];
	  //		b = newIndex[b];
	  //		c = newIndex[c];
	  // Lei: comment out above to use original indices
		if((a<0)||(a>=newS.getNumVert()))
			cerr << "A BAD" << endl;
		if((b<0)||(b>=newS.getNumVert()))
			cerr << "B BAD" << endl;
		if((c<0)||(c>=newS.getNumVert()))
			cerr << "C BAD" << endl;
		newS.facets()[cnt][0] = a;
		newS.facets()[cnt][1] = b;
		newS.facets()[cnt][2] = c;
		cnt++;
		}
  }

  *m_surface = newS;
  m_loadedSurface->surfaceChanged();
  m_surface->freeCurvature();

  m_isRed.setDim(m_surface->getNumVert());
  m_isRed = 0U;
  colorSurface(0);
  m_bModified = true;

  if(m_cutType == SurfaceEditor::Geodesic) {
     clearPath(1);
     m_tracker.setSurface(m_surface,SurfaceTracker::SameType);
  }
}


// Lei 11/18/2003: added this one below
void SurfaceEditor::loadPath()
{
VIN("void SurfaceEditor::loadPath()")
 char fname[1024];
 clearPath(1);
retry_load:
 fname[0] = 0;
 if(!Browser->GetFilename("Load Path",fname,1024,"*.idx"))
    return;

 FILE *fp=fopen(fname,"r");
 if(!fp) {
    App->ShowMessage("ERROR: could not open file!",1);
    goto retry_load;
 }

 int nl,np,i; 
 fscanf(fp,"%d",&nl);
 if(nl != 1) {
    fclose(fp);
    App->ShowMessage("ERROR: invalid path file",1);
    return;
 }
 char line[255];
 // this line ensures eating newlines
 while(!feof(fp) && (!fgets(line,255,fp) || (strlen(line)<3)));
 line[strlen(line)-1] = 0;
 sscanf(line,"%d",&np);
 m_Path.setDim(np);
 Array1D<Point> pts(np);
 for(i=0;i<np;i++) {
    fscanf(fp,"%d",&m_Path[i]);
    pts[i] = m_surface->vertices()[m_Path[i]];
 }
 fclose(fp);

 // the first point, borrowed from pickPoint()
 int indx=m_Path[0];
 float x = m_surface->vertices()[indx].x();
 float y = m_surface->vertices()[indx].y();
 float z = m_surface->vertices()[indx].z();
 float r = m_loadedSurface->glsurface()->getMaxExtent() * 0.005;
 GLPrimitive *sph = new GLPrimitive(x,y,z,r);
 if(!sph) return;
 sph->setColor(1.f,0.f,0.f);
 m_segments = GLPrimitive::addToList(m_segments,sph);
 m_surfBox->setPrimitives(m_segments);
 m_segStartIndex = 0;

 // the other points, borrowed from addSegment()
 GLPrimitive *seg = new GLPrimitive(pts);
 if(!seg) return;
 seg->setColor(0.f,0.f,1.f);
 seg->setLineWidth(3.f);
 m_segments = GLPrimitive::addToList(m_segments,seg);
 m_surfBox->setPrimitives(m_segments); 
 m_surfBox->redraw();
 XtSetSensitive(m_deleteButton,True);

VOUT("void SurfaceEditor::loadPath()")
}


// Lei 11/17/2003: added this one below
void SurfaceEditor::savePath()
{
VIN("void SurfaceEditor::savePath()")
   if( m_Path.getNelm() <= 0 ) {
      App->ShowMessage("No curves to save",1);
      return;
   }
 
 char fname[1024];
retry_save:
 fname[0] = 0;
 if(Browser->GetFilename("Save Path to",fname,1024)) {
    ofstream fp;
    fp.open(fname,ios::out);
	
    if(!fp || fp.fail()) {
       App->ShowMessage("ERROR: could not open file!",1);
       goto retry_save;
    }
    fp << "1" << endl
       << m_Path.getNelm() << " boundary path" << endl;
    m_Path.print(fp);
    fp << "0.00 0.00 1.00 3.00" << endl;
    fp.close();
 }
VOUT("void SurfaceEditor::savePath()")
}


// Lei 12/03/2003: added this one below
// Re-color entire path to be part of BLUE or RED surface
void SurfaceEditor::colorPath(int isred)
{
VIN("void SurfaceEditor::colorPath(int isred)")
 int i;
 int np=m_Path.getNelm();
 Array1D<Point> pts(np);
 for(i=0;i<np;i++) {
    pts[i] = m_surface->vertices()[m_Path[i]];
    if (isred)
       m_isRed[m_Path[i]] = 1U;
    else
       m_isRed[m_Path[i]] = 0U;
 }
 
 GLPrimitive *seg = new GLPrimitive(pts);
 if(!seg) return;
 if (isred)
    seg->setColor(1.f,0.f,0.f);
 else 
    seg->setColor(0.f,0.f,1.f);
 seg->setLineWidth(3.f);
 m_segments = GLPrimitive::addToList(m_segments,seg);
 m_surfBox->setPrimitives(m_segments); 
 m_surfBox->redraw();
 XtSetSensitive(m_deleteButton,True);

 m_surfBox->setPrimitives(m_segments);
 m_surfBox->redraw();
VOUT("void SurfaceEditor::colorPath(int isred)")
}




