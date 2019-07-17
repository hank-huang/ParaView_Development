//////////////////////////////////////////////////////////////////////////
//
// File: SurfaceBox.C
//
// Author: Keith Doolittle
//
// Purpose:
//
///////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include <KDApplication.h>
#include <Main.h>
#include <Xm/Form.h>
#include <Xm/DrawingA.h>
#include <Xm/RowColumn.h>
#include <Xm/Scale.h>
#include <KDShell.h>
#include <SurfaceBox.h>
#include <SurfaceColor.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GLSurface.h>
#include <GLPrimitive.h>

#define NEAR_ZPLANE     0.1
#define FAR_ZPLANE      10024.0

#ifndef DEFCOLOR_R
#define DEFCOLOR_R      0.8
#endif
#ifndef DEFCOLOR_G
#define DEFCOLOR_G      0.8
#endif
#ifndef DEFCOLOR_B
#define DEFCOLOR_B      0.9
#endif

#ifdef GNU_COMPILER
#define POLYGON_OFFSET	GL_POLYGON_OFFSET_FILL
#else
#define POLYGON_OFFSET	GL_POLYGON_OFFSET_EXT
#endif

#define BUFSIZE 4096*10
bool   SurfaceBox::m_bInitialized = false;
u_int  SurfaceBox::m_selectBuffer[BUFSIZE];

//
// Callback declarations
//
static void SetDrawStyleSolidCB(Widget, XtPointer, XtPointer);
static void SetDrawStyleBlendCB(Widget, XtPointer, XtPointer);
static void SetDrawStyleHiddenLineCB(Widget, XtPointer, XtPointer);
static void SetDrawStyleWireframeCB(Widget, XtPointer, XtPointer);
static void HomeCB(Widget,XtPointer,XtPointer);
static void SetHomeCB(Widget,XtPointer,XtPointer);
static void SetBlendAlphasCB(Widget,XtPointer,XtPointer);
static void GoToVertexCB(Widget,XtPointer,XtPointer);
static void GoToPointCB(Widget,XtPointer,XtPointer);
static void ViewAllCB(Widget,XtPointer,XtPointer);
static void ConfigureCB(Widget,XEvent*,String*,Cardinal*);
static void StartClickCB(Widget,XEvent*,String*,Cardinal*);
static void MotionCB(Widget,XEvent*,String*,Cardinal*);
static void EndClickCB(Widget,XEvent*,String*,Cardinal*);
static void ExposeCB(Widget,XEvent*,String*,Cardinal*);
static void ResizeCB(Widget,XEvent*,String*,Cardinal*);

    void GoToVertex();
//
// Xt Actions
//
static XtActionsRec sbactions[] = {
        { "configure",    ConfigureCB  },
        { "sbStartClick", StartClickCB },
        { "sbMotion",     MotionCB     },
        { "sbEndClick",   EndClickCB   },
        { "sbExpose",     ExposeCB     },
        { "sbResize",     ResizeCB     }
        };

//
// Menus
//
static menu_item style_items[] = {
        { KDButton, "Solid",(void *) SetDrawStyleSolidCB,NULL,NULL },
        { KDButton, "Wireframe",(void *) SetDrawStyleWireframeCB,NULL,NULL },
        { KDButton, "Hidden Line",(void *)SetDrawStyleHiddenLineCB,NULL,NULL},
        { KDButton, "Blend",(void *)SetDrawStyleBlendCB, NULL,NULL}
        };

static menu_item surf_items[] = {
        { KDSubMenu, "Draw Style",          (void*)style_items,
                                            (void*)XtNumber(style_items),NULL },
        { KDSeparator,"",                    NULL,NULL,NULL},
        { KDButton,   "View All",   (void *) ViewAllCB, NULL,NULL},
        { KDButton,   "Go to Home position",(void *) HomeCB, NULL,NULL},
        { KDButton,   "Set Home position", (void *)  SetHomeCB, NULL,NULL},
        { KDSeparator,"",                    NULL,NULL,NULL},
        { KDButton,   "Set Blend Alphas", (void *)  SetBlendAlphasCB, NULL,NULL},
        { KDSeparator,"",                    NULL,NULL,NULL},
        { KDButton,   "Go To Vertex", (void *)  GoToVertexCB, NULL,NULL},
// TODO        { KDButton,   "Go To Point", (void *)  GoToPointCB, NULL,NULL},
        };


//------------------------------------------------------------------
//  Callbacks for mouse motion, resizing, expose events, etc.
//------------------------------------------------------------------

DEFINECB(SurfaceBox, AlphaOKCB,AlphaOK)
DEFINECB(SurfaceBox, AlphaApplyCB,AlphaApply)
DEFINECB(SurfaceBox, HomeCB,   Home)
DEFINECB(SurfaceBox, SetHomeCB,SetHome)
DEFINECB(SurfaceBox, SetBlendAlphasCB,SetBlendAlphas)
DEFINECB(SurfaceBox, GoToVertexCB,GoToVertex)
DEFINECB(SurfaceBox, GoToPointCB,GoToPoint)
DEFINECB(SurfaceBox, ViewAllCB,ViewAll)
DEFINECBP(SurfaceBox,SetDrawStyleBlendCB,SetDrawStyle,SbDrawStyleBlend)
DEFINECBP(SurfaceBox,SetDrawStyleSolidCB,SetDrawStyle,SbDrawStyleSolid)
DEFINECBP(SurfaceBox,SetDrawStyleWireframeCB,SetDrawStyle,SbDrawStyleWireframe)
DEFINECBP(SurfaceBox,SetDrawStyleHiddenLineCB,SetDrawStyle,SbDrawStyleHiddenLine)

static void ConfigureCB(Widget wid, XEvent *ev, String *p, Cardinal *np)
{
VIN("static void ConfigureCB(Widget wid, XEvent *ev, String *p, Cardinal *np)")
  if(*np < 1) return;

  u_long ll;
  sscanf(p[0],"%lu",&ll);

  int w = ((XConfigureEvent*)ev)->width;
  int h = ((XConfigureEvent*)ev)->height;
  if(ll) ((SurfaceBox*)ll)->Configure(w,h);
VOUT("static void ConfigureCB(Widget wid, XEvent *ev, String *p, Cardinal *np)")
}

static void StartClickCB(Widget, XEvent *ev, String *p, Cardinal *np)
{
VIN("static void StartClickCB(Widget, XEvent *ev, String *p, Cardinal *np)")
  if(*np != 2) return;

  u_long ll,meta;
  sscanf(p[0],"%ld",&ll);
  sscanf(p[1],"%ld",&meta);

  int x = ((XButtonEvent*)ev)->x;
  int y = ((XButtonEvent*)ev)->y;

  if(ll) {
    // meta values: left - 0; shift left - 1; shift middle - 2; middle - 3; 
//    if((meta == 3)||(meta == 2)||(meta == 1)) // Lei 10/31/2003
    if((meta == 3)||(meta == 1)) // Lei 10/31/2003
    //    if((meta == 3)||(meta == 2))
       ((SurfaceBox*)ll)->select(x,y,meta);
    else
       ((SurfaceBox*)ll)->startClick(x,y,meta);
  }
VOUT("static void StartClickCB(Widget, XEvent *ev, String *p, Cardinal *np)")
}

static void MotionCB(Widget, XEvent *ev, String *p, Cardinal *np)
{
VIN("static void MotionCB(Widget, XEvent *ev, String *p, Cardinal *np)")
  if(*np != 2) return;

  u_long ll,meta;
  sscanf(p[0],"%ld",&ll);
  sscanf(p[1],"%ld",&meta);

  int x = ((XMotionEvent*)ev)->x;
  int y = ((XMotionEvent*)ev)->y;

  if(ll) ((SurfaceBox*)ll)->motion(x,y,meta);
VOUT("static void MotionCB(Widget, XEvent *ev, String *p, Cardinal *np)")
}


static void EndClickCB(Widget, XEvent *ev, String *p, Cardinal *np)
{
VIN("static void EndClickCB(Widget, XEvent *ev, String *p, Cardinal *np)")
  if(*np != 2) return;

  u_long ll,meta;
  sscanf(p[0],"%ld",&ll);
  sscanf(p[1],"%ld",&meta);

  int x = ((XButtonEvent*)ev)->x;
  int y = ((XButtonEvent*)ev)->y;

  if(ll) ((SurfaceBox*)ll)->endClick(x,y,meta);
VOUT("static void EndClickCB(Widget, XEvent *ev, String *p, Cardinal *np)")
}


static void ExposeCB(Widget, XEvent *ev, String *p, Cardinal *np)
{
VIN("static void ExposeCB(Widget, XEvent *ev, String *p, Cardinal *np)")
  if(*np != 1) return;

  // skip this if more expose events to follow
  if(((XExposeEvent*)ev)->count > 0) return;

  u_long ll;
  sscanf(p[0],"%ld",&ll);

  if(ll) ((SurfaceBox*)ll)->redraw();
VOUT("static void ExposeCB(Widget, XEvent *ev, String *p, Cardinal *np)")
}


static void ResizeCB(Widget, XEvent *ev, String *p, Cardinal *np)
{
VIN("static void ResizeCB(Widget, XEvent *ev, String *p, Cardinal *np)")
  if(*np != 1) return;

  int w = ((XConfigureEvent*)ev)->width;
  int h = ((XConfigureEvent*)ev)->height;

  u_long ll;
  sscanf(p[0],"%ld",&ll);

  if(ll) ((SurfaceBox*)ll)->resize(w,h);
VOUT("static void ResizeCB(Widget, XEvent *ev, String *p, Cardinal *np)")
}


//------------------------------------------------------------------
// Private procedures
//------------------------------------------------------------------

bool SurfaceBox::setCurrent()
{
VIN("bool SurfaceBox::setCurrent()")
  if(glx_context) {
     GLwDrawingAreaMakeCurrent(render_widget,glx_context);
     VOUT("bool SurfaceBox::setCurrent()")
     return(true);
  } 
VOUT("bool SurfaceBox::setCurrent()")
  return(false); 
}

void SurfaceBox::calculateCenter()
{
VIN("void SurfaceBox::calculateCenter()")
  double sumx = 0.;
  double sumy = 0.;
  double sumz = 0.;
  double nsum = 0.;

  for(int i=0;i<numSurfaces();i++) {
     if(!m_Surfaces[i] || !m_Surfaces[i]->glsurface()) 
        continue;

     sumx += m_Surfaces[i]->glsurface()->getCentX();
     sumy += m_Surfaces[i]->glsurface()->getCentY();
     sumz += m_Surfaces[i]->glsurface()->getCentZ();
     nsum += 1.;
  }

  if(nsum > 0)
     center.set(sumx/nsum,sumy/nsum,sumz/nsum);
  else {
     // No surfaces, use primitives
     GLPrimitive *p;
     float xx,yy,zz,ww,hh,dd;
     for(p=m_primitives;p;p=p->Next()) {
        p->getExtents(ww,hh,dd,xx,yy,zz);
        sumx += xx;
        sumy += yy;
        sumz += zz;
        nsum += 1;
     }
     if(nsum > 0)
          center.set(sumx/nsum,sumy/nsum,sumz/nsum);
     else center.set(0.,0.,0.);
  }
VOUT("void SurfaceBox::calculateCenter()")
}

void SurfaceBox::resize(int w, int h)
{
VIN("void SurfaceBox::resize(int w, int h)")
  XWindowAttributes watt;
  XGetWindowAttributes(App->AppDisplay(),XtWindow(render_widget),&watt);
  if(watt.map_state != IsViewable)
	return;
 
  m_width  = w;
  m_height = h;
  setCamera(m_eye.x(),m_eye.y(),m_eye.z());
VOUT("void SurfaceBox::resize(int w, int h)")
}

void SurfaceBox::Configure(int w, int h)
{
VIN("void SurfaceBox::Configure(int w, int h)")
  XWindowAttributes watt;
  XGetWindowAttributes(App->AppDisplay(),XtWindow(render_widget),&watt);
  if(watt.map_state != IsViewable)
	return;

  m_width  = w; 
  m_height = h; 

  Arg arglist[2];
  Cardinal c = 0;
  Set(XmNwidth,m_width);
  Set(XmNheight,m_height);
  XtSetValues(render_widget,arglist,c);

  rotate_per_pixel = 360.f/(float)m_width;
VOUT("void SurfaceBox::Configure(int w, int h)")
}

void SurfaceBox::ViewAll()
{
VIN("void SurfaceBox::ViewAll()")
  int   i;
  float x,y,z,d,maxd;

  maxd = 0.f;
  for(i=0;i<m_numSurfaces;i++) {
	if(!m_Surfaces[i] || !m_Surfaces[i]->glsurface()) 
            continue;

	x = center.x() - m_Surfaces[i]->glsurface()->getCentX();
	y = center.y() - m_Surfaces[i]->glsurface()->getCentY();
	z = center.z() - m_Surfaces[i]->glsurface()->getCentZ();
	d = sqrt(x*x+y*y+z*z) + m_Surfaces[i]->glsurface()->getMaxExtent();
	if(d > maxd) maxd = d;
  }
  if(maxd == 0.f) maxd = center.z();

  setCamera(0.f,0.f,-1.5f * maxd);
  redraw();
VOUT("void SurfaceBox::ViewAll()")
}


void SurfaceBox::SetBlendAlphas()
{
VIN("void SurfaceBox::SetBlendAlphas()")
  if(m_numSurfaces <= 0)
     return;

  Arg arglist[30];
  Cardinal c;
  Widget wid;
  XmString str;
  int i;
  int nitems = m_numSurfaces + 1;
  int bh = App->BigFontHeight() + 15;
  int w  = 20 * App->MedFontWidth(); 
  int h  = nitems * bh;
  int alpha,ty,by,dy = 100/nitems;
  int x = App->Width()/2 - w/2;
  int y = App->Height()/2 - h/2;

  if(m_alphaShell) {
     m_alphaShell->Hide();
     delete m_alphaShell;
  }
  m_alphaShell = new KDShell("Set Blend Alphas", x,y,w,h,KDSHELL_TITLEBAR,0,
			MainWindow->Shell(),SHELLCOLOR);

  ty = 0;
  for(i=0;i<m_numSurfaces;i++) {
     if(m_Surfaces[i])
        alpha = (int)(m_Surfaces[i]->glsurface()->getAlpha() * 100.f);

     str = XmStringCreate(&(m_names[i][0]),App->SmallFontTag());

     by = ty + dy;
     ATTACHMENTS_TBLR(ty,by,0,100)
     Set(XmNminimum,0);
     Set(XmNmaximum,100);
     Set(XmNvalue,alpha);
     Set(XmNorientation,XmHORIZONTAL);
     Set(XmNfontList,App->SmallFontList());
     Set(XmNtitleString,str);
     m_sliders[i] = XmCreateScale(m_alphaShell->Main(),"ShellAlpha",arglist,c);
     XtManageChild(m_sliders[i]);
     App->RecursiveColor(m_sliders[i],SHELLCOLOR);
     ty = by;
  } 

  by = ty + dy;
  ATTACHMENTS_TBLR(ty,by,20,40)
  wid = App->Button(m_alphaShell->Main(),"OK",App->MedFontList(),
		AlphaOKCB,(XtPointer)this, arglist,c);
  ATTACHMENTS_TBLR(ty,by,60,80)
  wid = App->Button(m_alphaShell->Main(),"APPLY",App->MedFontList(),
		AlphaApplyCB, (XtPointer)this, arglist,c);
  
  m_alphaShell->Show();
VOUT("void SurfaceBox::SetBlendAlphas()")
}

void SurfaceBox::AlphaOK()
{
VIN("void SurfaceBox::AlphaOK()")
  AlphaApply();
  if(m_alphaShell) {
     m_alphaShell->Hide();
     delete m_alphaShell;
     m_alphaShell = NULL;
  }
VOUT("void SurfaceBox::AlphaOK()")
}

void SurfaceBox::AlphaApply()
{
VIN("void SurfaceBox::AlphaApply()")
  int i,alpha;

  if(!m_alphaShell)
     return;

  for(i=0;i<m_numSurfaces;i++) {
      XmScaleGetValue(m_sliders[i],&alpha);
      m_Surfaces[i]->glsurface()->setAlpha((float)alpha/100.f);
      }
  redraw();
VOUT("void SurfaceBox::AlphaApply()")
}


void SurfaceBox::GoToPoint()
{
VIN("void SurfaceBox::GoToPoint()")

VOUT("void SurfaceBox::GoToPoint()")
}

void SurfaceBox::GoToVertex()
{
VIN("void SurfaceBox::GoToVertex()")
  if(m_numSurfaces <= 0)
     return;

  if(!pick_func)
  {
    App->ShowMessage("No pick function in this view.");
    return;
  }

  char tstr[1024];
  sprintf(tstr,"1");

retry_bin:
  int vertex = 0;
  if(!App->GetText("Enter Surface 1 Vertex to Mark (-1 to Cancel):",tstr,20))
      return;

  if((sscanf(tstr,"%d",&vertex) != 1)) {
      App->ShowMessage("Please an integer vertex value");
      goto retry_bin;
  } else if (vertex == -1) return;
  else if (vertex < 0) {
      App->ShowMessage("Invalid Vertex value");
      goto retry_bin;
  } else if (vertex > m_Surfaces[0]->surface()->getNumVert() - 1) {
      char estr[1024];
      sprintf(estr, "Max vertex value for surface 1 is %d.",
              m_Surfaces[0]->surface()->getNumVert() - 1);
      App->ShowMessage(estr);
      goto retry_bin;
  }

  if(pick_func)
  pick_func(pick_dat,(void*)m_Surfaces[0]->surface(),vertex);
VOUT("void SurfaceBox::SetBlendAlphas()")
}

void SurfaceBox::SetHome()
{
VIN("void SurfaceBox::SetHome()")
  for(int i=0;i<16;i++)
	m_home_matrix[i] = m_surfacemtx[i];

  m_home_center = center;
  m_home_eye    = m_eye;
VOUT("void SurfaceBox::SetHome()")
}

void SurfaceBox::Home()
{
VIN("void SurfaceBox::Home()")
  int i;
  for(i=0;i<16;i++)
	m_surfacemtx[i] = m_home_matrix[i];

  center = m_home_center;
  setCamera(m_home_eye.x(),m_home_eye.y(),m_home_eye.z());
  redraw();
VOUT("void SurfaceBox::Home()")
}

void SurfaceBox::set2D()
{
VIN("void SurfaceBox::set2D()")
  glViewport(0,0,m_width,m_height);
  glGetIntegerv(GL_VIEWPORT,m_viewport);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0,m_width,0,m_height,-0.1,0.1);
VOUT("void SurfaceBox::set2D()")
}

void SurfaceBox::set3D()
{
VIN("void SurfaceBox::set3D()")
  setCamera(m_home_eye.x(),m_home_eye.y(),m_home_eye.z());
VOUT("void SurfaceBox::set3D()")
}

void SurfaceBox::drawMessage()
{
VIN("void SurfaceBox::drawMessage()")
  if(m_message[0] == 0)
     return;

  set2D();
  m_GText.drawTextSmall(m_message,5,7,1.f,1.f,0.f);
VOUT("void SurfaceBox::drawMessage()")
}

void SurfaceBox::setMessage(const char *str)
{
VIN("void SurfaceBox::setMessage(const char *str)")
  if(!str || (str[0] == 0))
    m_message[0] = 0;
  else
    sprintf(m_message,"%s",str);
VOUT("void SurfaceBox::setMessage(const char *str)")
}

void SurfaceBox::setCamera(float x, float y, float z)
{
VIN("void SurfaceBox::setCamera(float x, float y, float z)")
  m_eye.set(x,y,z);

  if(!setCurrent())
	return;

  glViewport(0,0,m_width,m_height);
  glGetIntegerv(GL_VIEWPORT,m_viewport);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  gluPerspective(60.,(float)m_width/(float)m_height,NEAR_ZPLANE,FAR_ZPLANE);
  glTranslated(m_eye.x(),m_eye.y(),m_eye.z()); 
  glGetDoublev(GL_PROJECTION_MATRIX,m_projection);
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
VOUT("void SurfaceBox::setCamera(float x, float y, float z)")
}

void SurfaceBox::addRotation(float angle, float axx, float axy, float axz)
{
VIN("void SurfaceBox::addRotation(float angle, float axx, float axy, float axz)")
  glMatrixMode(GL_MODELVIEW);
//  glPushMatrix();
  glLoadIdentity();
  glRotatef(angle,axx,axy,axz);
  glMultMatrixd(m_surfacemtx);
  glGetDoublev(GL_MODELVIEW_MATRIX,m_surfacemtx);
//  glMatrixMode(GL_MODELVIEW);
//  glPopMatrix();
VOUT("void SurfaceBox::addRotation(float angle, float axx, float axy, float axz)")
}


void SurfaceBox::startClick(int x,int y, int m)
{
VIN("void SurfaceBox::startClick(int x,int y, int m)")
  drag_sx   = x; 
  drag_sy   = y; 
  start_sx  = x; 
  start_sy  = y; 
  drag_meta = m;

  if(mouse_func && !mouse_func(mouse_dat,x,y,m,0))
	return;
VOUT("void SurfaceBox::startClick(int x,int y, int m)")
}

void SurfaceBox::motion(int x, int y, int m)
{
VIN("void SurfaceBox::motion(int x, int y, int m)")
  setMessage(NULL);

  bool isflat;
  if(numSurfaces() > 0)
	isflat = m_Surfaces[0]->glsurface()->isFlat();
  else  isflat = false;	

  isflat = false;
  //else  return;	

  int dx  = x - drag_sx;
  int dy  = y - drag_sy;
  drag_sx = x;
  drag_sy = y;

  if(mouse_func && !mouse_func(mouse_dat,x,y,m,1))
	return;

  float fact;

  switch(m) {
  case 0: // rotate (Button1)
     if(!isflat) {
        if(dy) rotateScreenX( (float)dy * rotate_per_pixel);
        if(dx) rotateScreenY( (float)dx * rotate_per_pixel);
        redraw();
     }
     break;
  case 1: // zoom   (Button1 + shift)
     fact = (float)dx/(float)m_width * m_eye.norm(); 
     setCamera(m_eye.x(),m_eye.y(),m_eye.z() + fact);
     redraw();
     break;
  case 2: // translate (Button2 + shift)
     if(dx || dy)
	translate(dx,dy);
     break;
  default:
     break;
  }
VOUT("void SurfaceBox::motion(int x, int y, int m)")
}


void SurfaceBox::endClick(int x, int y, int m)
{
VIN("void SurfaceBox::endClick(int x, int y, int m)")
  if(mouse_func && !mouse_func(mouse_dat,x,y,m,2))
	return;
VOUT("void SurfaceBox::endClick(int x, int y, int m)")
}

void SurfaceBox::translate(int x, int y)
{
VIN("void SurfaceBox::translate(int x, int y)")
  float dx = (float)x * move_per_pixel;
  float dy = (float)-y * move_per_pixel;

  float xaxisx = m_projection[4*0+0] * dx;
  float xaxisy = m_projection[4*0+1] * dx;
  float xaxisz = m_projection[4*0+2] * dx;

  float yaxisx = m_projection[4*1+0] * dy;
  float yaxisy = m_projection[4*1+1] * dy;
  float yaxisz = m_projection[4*1+2] * dy;

  setCamera(m_eye.x()+xaxisx+yaxisx,
            m_eye.y()+xaxisy+yaxisy,
            m_eye.z()+xaxisz+yaxisz);
  redraw();
VOUT("void SurfaceBox::translate(int x, int y)")
}

void SurfaceBox::rotateScreenX(float deg)
{
VIN("void SurfaceBox::rotateScreenX(float deg)")
  if(!setCurrent())
	return;

  float axisx = m_projection[4*0+0];
  float axisy = m_projection[4*0+1];
  float axisz = m_projection[4*0+2];

  addRotation(deg,axisx,axisy,axisz);
VOUT("void SurfaceBox::rotateScreenX(float deg)")
}


void SurfaceBox::rotateScreenY(float deg)
{
VIN("void SurfaceBox::rotateScreenY(float deg)")
  if(!setCurrent())
	return;

  float axisx = m_projection[4*1+0];
  float axisy = m_projection[4*1+1];
  float axisz = m_projection[4*1+2];

  addRotation(deg,axisx,axisy,axisz);
VOUT("void SurfaceBox::rotateScreenY(float deg)")
}


void SurfaceBox::redraw()
{
VIN("void SurfaceBox::redraw()")
  if(!setCurrent())
	return;

  glClearColor(m_background.x(),m_background.y(),m_background.z(),1.f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_PROJECTION);
//  glPushMatrix();
  glLoadMatrixd(m_projection);

  render();

//  set2D();
  drawMessage();

//  glMatrixMode(GL_PROJECTION);
//  glPopMatrix();

  GLwDrawingAreaSwapBuffers(render_widget);
VOUT("void SurfaceBox::redraw()")
}


void SurfaceBox::select(int winx, int winy, int meta)
{
VIN("void SurfaceBox::select(int winx, int winy, int meta)")
// || (pick_style == PickNone) || !pick_func)

  if(!setCurrent())
	return;

  glSelectBuffer(BUFSIZE,m_selectBuffer);
  glRenderMode(GL_SELECT);
  glInitNames();
  glPushName(0);

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();

  glLoadIdentity();
  gluPickMatrix(winx,m_height-winy-1,5,5,m_viewport);
  glMultMatrixd(m_projection);

  renderSelect();

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  GLuint *ptr;
  u_int   nnames,z1,z2,minz,hits;
  int     i,srf,pol,minsrf,minpol;

  hits   = glRenderMode(GL_RENDER);
  minz   = (u_int)-1;
  minsrf = -1;
  minpol = -1;

  ptr = m_selectBuffer;
  for(i=0;i<hits;i++) {
	nnames = *ptr++; // will always be 2
        z1     = *ptr++; // min pseudo Z value
        z2     = *ptr++; // max pseudo Z value
        srf    = *ptr++; // surface #
        pol    = *ptr++; // polygon #
        if(z1 < minz) {
           minz   = z1;
           minsrf = srf;
           minpol = pol;
        }      
        if(z2 < minz) {
           minz   = z2;
           minsrf = srf;
           minpol = pol;
        }      
  }

  if((minsrf<0)||(minpol<0)||(minsrf>=m_numSurfaces))
	return;
  
  Surface *SS = m_Surfaces[minsrf]->surface();
  if(!SS) return;

  GLSurface *S = m_Surfaces[minsrf]->glsurface();

  if(!S || (minpol >= S->getNumPoly()))
	return;

  int   a,b,c,mini;
  float ax,ay,az;
  float bx,by,bz;
  float cx,cy,cz;
  float dx,dy,da,db,dc;
  char  pstr[1024];

  SWPickStyle p;
  //  bool doextract;
  int doextract;		// Lei 04/17/2003
				// 0 = NO; 2 = extract; 3 = extractKeep

  if (meta == 2) {		// Lei 04/17/2003
    p = PickPoint;
    doextract = 2;
  } else if (meta == 1) {
    p = PickPoint;
    doextract = 1;
  } else {
    p = pick_style;
    doextract = 0;
  }

  /*---<Lei 04/17/2003>---*
   | if(meta == 2) {      |
   |   p = PickPoint;     |
   |   doextract = true;  |
   | } else {             |
   |   p = pick_style;    |
   |   doextract = false; |
   | }                    |
   *----------------------*/

  switch(p) {
  case PickSurface:
        if(pick_func)
	   pick_func(pick_dat,(void*)S,0);
        break;
  case PickPoly:
        if(pick_func)
	   pick_func(pick_dat,(void*)S,minpol);
        break;
  case PickPoint:
        a = S->getFacets()[minpol][0];
        b = S->getFacets()[minpol][1];
        c = S->getFacets()[minpol][2];
        project(S->getVertices()[3*a+0],S->getVertices()[3*a+1],
                S->getVertices()[3*a+2],ax,ay,az);
        project(S->getVertices()[3*b+0],S->getVertices()[3*b+1],
                S->getVertices()[3*b+2],bx,by,bz);
        project(S->getVertices()[3*c+0],S->getVertices()[3*c+1],
                S->getVertices()[3*c+2],cx,cy,cz);

        dx = ax - winx;
        dy = ay - winy;
        da = sqrt(dx*dx+dy*dy);

        dx = bx - winx;
        dy = by - winy;
        db = sqrt(dx*dx+dy*dy);

        dx = cx - winx;
        dy = cy - winy;
        dc = sqrt(dx*dx+dy*dy);

        if(da < db) {
           if(da < dc) 
                mini = a;
           else 
                mini = c;
        } else {
           if(db < dc) 
                mini = b;
           else 
                mini = c;
        }

        if (doextract) {
          if(App->GetYesNo("Extract Surface?"))
            if(App->GetYesNo("ExtractKeep Surface?"))
              SS->extractKeep(mini);
            else
              SS->extract(mini);
          m_Surfaces[minsrf]->surfaceChanged();
	      redraw();
        } else 
	  if(pick_func) pick_func(pick_dat,(void*)SS,mini);

	/*-------------------<Lei 04/17/2003>-------------------*
         | if(doextract) {                                      |
         |    if(App->GetYesNo("Extract Surface?"))             |
         |       SS->extract(mini);                             |
         |       m_Surfaces[minsrf]->surfaceChanged();          |
         |       redraw();                                      |
         | } else                                               |
         |    if(pick_func) pick_func(pick_dat,(void*)SS,mini); |
	 *------------------------------------------------------*/

        break;
  case PickNone:
  default:
        break;
  }
VOUT("void SurfaceBox::select(int winx, int winy, int meta)")
}

bool SurfaceBox::setupForStyle(SbDrawStyle d)
{
VIN("bool SurfaceBox::setupForStyle(SbDrawStyle d)")
  bool col = false;
  switch(d) {
  case SbDrawStyleWireframe:
     glDisable(GL_BLEND);
     glDisable(GL_DEPTH_TEST);
     glLineWidth(1.f);
     glDisable(GL_CULL_FACE);
     glDisable(GL_COLOR_MATERIAL);
     glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
     glDisable(GL_LIGHTING);
     glColor3f(DEFCOLOR_R,DEFCOLOR_G,DEFCOLOR_B);
     glDisable(POLYGON_OFFSET);
     col = false;
     break;
  case SbDrawStyleHiddenLine:
     glDisable(GL_BLEND);
     glEnable(GL_DEPTH_TEST);
     glDepthFunc(GL_LESS);
     glLineWidth(1.f);
     glEnable(GL_CULL_FACE);
     glDisable(GL_COLOR_MATERIAL);
     glCullFace(GL_BACK);
     glPolygonMode(GL_FRONT,GL_LINE);
     glDisable(GL_LIGHTING);
     glColor3f(DEFCOLOR_R,DEFCOLOR_G,DEFCOLOR_B);
     glDisable(POLYGON_OFFSET);
     col = false;
     break;
  case SbDrawStyleSolid:
     glDisable(GL_BLEND);
     glEnable(GL_DEPTH_TEST);
     glDepthFunc(GL_LESS);
     glEnable(GL_CULL_FACE);
     glFrontFace(GL_CCW);
     glCullFace(GL_BACK);
     glEnable(POLYGON_OFFSET);
     glPolygonOffset(1.0,2);
     glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
//     glPolygonMode(GL_FRONT,GL_FILL);
     glEnable(GL_COLOR_MATERIAL);
     if(m_bLighting) {
       glEnable(GL_LIGHTING);
       glEnable(GL_LIGHT0);
     } else
       glDisable(GL_LIGHTING);
     glColor3f(DEFCOLOR_R,DEFCOLOR_G,DEFCOLOR_B);
     col = true;
     break;
  case SbDrawStyleBlend:
     glEnable(GL_BLEND);
     glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
     glEnable(GL_DEPTH_TEST);
     glDepthFunc(GL_LESS);
     glEnable(GL_CULL_FACE);
     glFrontFace(GL_CCW);
     glCullFace(GL_BACK);
     glEnable(POLYGON_OFFSET);
     glPolygonOffset(1.0,2);
     glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
// keith
//     glPolygonMode(GL_FRONT,GL_FILL);
     glEnable(GL_COLOR_MATERIAL);
     if(m_bLighting) {
       glEnable(GL_LIGHTING);
       glEnable(GL_LIGHT0);
     } else
       glDisable(GL_LIGHTING);
     glColor3f(DEFCOLOR_R,DEFCOLOR_G,DEFCOLOR_B);
     col = true;
     break;
  default:
     col = false;
  }

col = true;

VOUT("bool SurfaceBox::setupForStyle(SbDrawStyle d)")
return(col);
}


static void findRange(const float *dat, int nelem, float &fmin, float &fmax)
{
  if(!dat) { fmin = fmax = -1.f; return; }

  float v;
  fmin = fmax = *dat;
  for(int i=1;i<nelem;i++) {
       v = dat[i];
       if(v < fmin) fmin = v;
       if(v > fmax) fmax = v;
  }
}


static void findRange(const int *dat, int nelem, int &fmin, int &fmax)
{
  if(!dat) { fmin = fmax = -1; return; }

  int v;
  fmin = fmax = *dat;
  for(int i=1;i<nelem;i++) {
       v = dat[i];
       if(v < fmin) fmin = v;
       if(v > fmax) fmax = v;
  }
}



void SurfaceBox::superDump()
{
  int i,j;

  cout << "SurfaceBox: Super Dump" << endl;
  cout << "         W/H: " << m_width << " " << m_height << endl;
  cout << "   RenderWid: " << (unsigned long)render_widget << endl;
  cout << "     Context: " << (unsigned long)glx_context << endl;
flush(cout);
  if(m_bLighting) 
  cout << "    Lighting: ON" << endl;
  else
  cout << "    Lighting: OFF" << endl;
  cout << "         Eye: " << m_eye.x() << " " << m_eye.y() << " " << m_eye.z() << endl;
  cout << "      Center: " << center.x() << " " << center.y() << " " << center.z() << endl;
  cout << "  Primitives: " << (unsigned long)m_primitives << endl;
flush(cout);

  cout << "  SurfaceMtx: " << endl;
  for(i=0;i<4;i++) { 
    for(j=0;j<4;j++) 
     cout << m_surfacemtx[i*4+j] << " ";
    cout << endl;
  }
  cout << endl;
flush(cout);

  cout << "   Modelview: " << endl;
  for(i=0;i<4;i++) { 
    for(j=0;j<4;j++) 
     cout << m_modelview[i*4+j] << " ";
    cout << endl;
  }
  cout << endl;
flush(cout);


  cout << " Projection: " << endl;
  for(i=0;i<4;i++) { 
    for(j=0;j<4;j++) 
     cout << m_projection[i*4+j] << " ";
    cout << endl;
  }
  cout << endl;
flush(cout);

  cout << "   Viewport: " << m_viewport[0] << " " << m_viewport[1] << " " << m_viewport[2] <<" " << m_viewport[3] << endl;

flush(cout);
  if(m_bInitialized)
  cout << "       Init: TRUE" << endl;
  else
  cout << "       Init: FALSE!!!!!!" << endl;

  cout << "  N Surfaces: " << m_numSurfaces << endl;

  GLSurface *S;
  int np,nf,imin,imax;
  float fmin,fmax;
  for(i=0;i<m_numSurfaces;i++) {
     S = m_Surfaces[i]->glsurface();
     cout << "      Surface: " << i << endl;
     if(!S) {
     cout << "           IS NULL!!!" << endl;
     continue;
     }
flush(cout);
  
     np = S->getNumVert(); 
     nf = S->getNumPoly(); 
     cout << "                  NV/NF: " << np << " " << nf << endl;

     findRange(S->getVertices(),3*np,fmin,fmax);
     cout << "           Vertices Rng: " << fmin << " " << fmax << endl;
flush(cout);
     findRange(S->getNormals(),3*np,fmin,fmax);
     cout << "            Normals Rng: " << fmin << " " << fmax << endl;
flush(cout);
     findRange(S->getColors(),4*np,fmin,fmax);
     cout << "             Colors Rng: " << fmin << " " << fmax << endl;
flush(cout);
     findRange(S->getFacets().data(),3*nf,imin,imax);
     cout << "             Facets Rng: " << imin << " " << imax << endl;
flush(cout);
  }
cout << endl;
}

//
// Normal rendering, fastest
//
void SurfaceBox::render()
{
VIN("void SurfaceBox::render()")
#ifdef _TRACE_
  superDump();
#endif

  GLSurface *S;

  glMatrixMode(GL_MODELVIEW);
//  glPushMatrix();

  glLoadMatrixd(m_surfacemtx);
  glTranslated(-center.x(),-center.y(),-center.z());
  glGetDoublev(GL_MODELVIEW_MATRIX,m_modelview);

  bool doColor = setupForStyle(draw_style);
  bool haveC;

VTRACE("Render-1");

  for(int i=0;i<m_numSurfaces;i++) {
	if(!(S = m_Surfaces[i]->glsurface()))
		continue;

VTRACE("Render-2");

  glBegin(GL_TRIANGLES);
  for(int pp=0;pp<S->getNumPoly();pp++) {
#ifdef _TRACE_
if(pp % 10 == 9) { cout << "."; flush(cout); }
#endif
       int a = S->getFacets()[pp][0];
       int b = S->getFacets()[pp][1];
       int c = S->getFacets()[pp][2];

       haveC = (doColor && (S->getColors() != NULL));

       if(!haveC)
          glColor3f(DEFCOLOR_R,DEFCOLOR_G,DEFCOLOR_B);

       if(haveC)
         glColor4f(S->getColors()[4*a+0],S->getColors()[4*a+1],
                   S->getColors()[4*a+2],S->getColors()[4*a+3]);
       glNormal3f(S->getNormals()[3*a+0],S->getNormals()[3*a+1],
                  S->getNormals()[3*a+2]);
       glVertex3f(S->getVertices()[3*a+0],S->getVertices()[3*a+1],
                  S->getVertices()[3*a+2]);

       if(haveC)
         glColor4f(S->getColors()[4*b+0],S->getColors()[4*b+1],
                   S->getColors()[4*b+2],S->getColors()[4*b+3]);
       glNormal3f(S->getNormals()[3*b+0],S->getNormals()[3*b+1],
                  S->getNormals()[3*b+2]);
       glVertex3f(S->getVertices()[3*b+0],S->getVertices()[3*b+1],
                  S->getVertices()[3*b+2]);

       if(haveC)
         glColor4f(S->getColors()[4*c+0],S->getColors()[4*c+1],
                   S->getColors()[4*c+2],S->getColors()[4*c+3]);
       glNormal3f(S->getNormals()[3*c+0],S->getNormals()[3*c+1],
                  S->getNormals()[3*c+2]);
       glVertex3f(S->getVertices()[3*c+0],S->getVertices()[3*c+1],
                  S->getVertices()[3*c+2]);

  }
VTRACE("Render-3");
  glEnd();
VTRACE("Render-4");

/*
VTRACE("Render-2");
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3,GL_FLOAT,0,S->getVertices());

	if(doColor && (S->getColors())) {
	   glEnableClientState(GL_COLOR_ARRAY);
           glColorPointer(4,GL_FLOAT,0,S->getColors());
        } else {
	   glDisableClientState(GL_COLOR_ARRAY);
           glColor3f(DEFCOLOR_R,DEFCOLOR_G,DEFCOLOR_B);
        }

	if(m_bLighting && (S->getNormals())) {
	   glEnableClientState(GL_NORMAL_ARRAY);
	   glNormalPointer(GL_FLOAT,0,S->getNormals());
	} else
	   glDisableClientState(GL_NORMAL_ARRAY);

VTRACE("Render-3");
	glDrawElements(GL_TRIANGLES,3 * S->getNumPoly(),
		       GL_UNSIGNED_INT,(void*)S->getFacets().data());
VTRACE("Render-4");
*/


  }

VTRACE("Render-5");
  GLPrimitive *p;
  for(p=m_primitives;p;p=p->Next())
	p->render();

VTRACE("Render-6");
  // allow sub-classes to draw here (using valid GL matrices)
  if(redraw_func) 
     redraw_func(redraw_dat);

VTRACE("Render-7");
//  glMatrixMode(GL_MODELVIEW);
//  glPopMatrix();
VOUT("void SurfaceBox::render()")
}


//
// select mode rendering, slower. Name each poly
//
void SurfaceBox::renderSelect()
{
VIN("void SurfaceBox::renderSelect()")
  GLSurface *S;

  glMatrixMode(GL_MODELVIEW);
//  glPushMatrix();

  glLoadMatrixd(m_surfacemtx);
  glTranslated(-center.x(),-center.y(),-center.z());
  glGetDoublev(GL_MODELVIEW_MATRIX,m_modelview);

  (void)setupForStyle(draw_style);

  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);

  for(int i=0;i<m_numSurfaces;i++) {
	if(!(S = m_Surfaces[i]->glsurface()))
		continue;

	for(int p=0;p<S->getNumPoly();p++) {
	    glLoadName(i); // surface number
	    glPushName(p); // poly number
	    glBegin(GL_TRIANGLES);

	    int a = S->getFacets()[p][0];
	    int b = S->getFacets()[p][1];
	    int c = S->getFacets()[p][2];

	    glVertex3fv(&(S->getVertices()[3*a]));
	    glVertex3fv(&(S->getVertices()[3*b]));
	    glVertex3fv(&(S->getVertices()[3*c]));

            glEnd();
	    glPopName();
        }
  }

//  glMatrixMode(GL_MODELVIEW);
//  glPopMatrix();
VOUT("void SurfaceBox::renderSelect()")
}


void SurfaceBox::renderLine(int x1, int y1, int x2, int y2, int wid,
                            unsigned char r, unsigned char g, unsigned char b)
{
VIN("void SurfaceBox::renderLine(int x1, int y1, int x2, int y2, int wid,")
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  gluOrtho2D(0.0,m_width,0.0,m_height);
  glLineWidth(wid);

  glMatrixMode(GL_MODELVIEW);
//  glPushMatrix();
  glLoadIdentity();

  glColor3ub(r,g,b);
  glBegin(GL_LINES);
  glVertex2i(x1,m_height - y1 - 1);
  glVertex2i(x2,m_height - y2 - 1);
  glEnd();

//  glMatrixMode(GL_MODELVIEW);
//  glPopMatrix();

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
VOUT("void SurfaceBox::renderLine(int x1, int y1, int x2, int y2, int wid,")
}


SurfaceBox::SurfaceBox(Widget parent, LoadedSurface *lsurf,
		       Arg *args, Cardinal nargs)
{
VIN("SurfaceBox::SurfaceBox(Widget parent, LoadedSurface *lsurf,")
  int i,j;

  if(!m_bInitialized) {
     XtAppAddActions(App->AppContext(),sbactions,XtNumber(sbactions));
     m_bInitialized = true;
  }
  
  m_numSurfaces = 0;
  m_curSurface  = 0;
  redraw_func   = NULL;
  redraw_dat    = 0L;
  mouse_func    = NULL;
  mouse_dat     = 0L;
  pick_func     = NULL;
  pick_dat      = 0L;
  m_primitives  = NULL;
  m_bLighting   = true;
  glx_context   = NULL;
  m_alphaShell  = NULL;

  rotate_per_pixel = 1.f;
  move_per_pixel   = 1.f;

  m_background.set(0,0,0);

  pick_style    = PickNone;
  draw_style    = SbDrawStyleSolid;

  m_message[0] = 0;

  m_eye.set(0.,0.,0.);
  center.set(0.,0.,0.);

  for(i=0;i<4;i++) {
    m_viewport[i] = 0;
    for(j=0;j<4;j++)
       if(i==j) {
	  m_surfacemtx[4*i+j] = 1.;
	  m_modelview[4*i+j]  = 1.;
	  m_projection[4*i+j] = 1.;
       } else {
	  m_surfacemtx[4*i+j] = 0.;
	  m_modelview[4*i+j]  = 0.;
	  m_projection[4*i+j] = 0.;
       }
  }

  Arg arglist[30];
  Cardinal c;

  form = XmCreateDrawingArea(parent,"SurfaceBox",args,nargs);
  XtManageChild(form);
  XmChangeColor(form,SHELLCOLOR);

  c = 0;
  Set(XmNmarginWidth,0);
  Set(XmNmarginHeight,0);
  Set(XmNresizable,False);
  Set(XmNhorizontalSpacing,0);
  Set(XmNverticalSpacing,0);
  XtSetValues(form,arglist,c);

  char tstring[2000];
  sprintf(tstring,"<Configure>:configure(%ld)",(u_long)this);
  XtTranslations ftrans = XtParseTranslationTable(tstring);
  XtOverrideTranslations(form,ftrans);

  Dimension fw,fh;
  c = 0;
  Set(XmNwidth,&fw);
  Set(XmNheight,&fh);
  XtGetValues(form,arglist,c);
  m_width  = fw;
  m_height = fh;

  c = 0;
  Set(XmNwidth,fw);
  Set(XmNheight,fh);
  Set(XmNbackground, 0);
  Set(GLwNallocateBackground, TRUE);
  Set(GLwNinstallBackground, TRUE);
  //
  // These should be same as in KDApplication.C glXChooseVisual()
  //
  Set(GLwNredSize, 1);
  Set(GLwNgreenSize, 1);
  Set(GLwNblueSize, 1);
  Set(GLwNdepthSize,1);
  Set(GLwNalphaSize, 0);
  Set(GLwNrgba, True);
  Set(GLwNdoublebuffer, TRUE);
  Set(GLwNvisualInfo,App->glVisualInfo());
  Set(XmNtraversalOn, FALSE);
  Set(XmNtopAttachment,XmATTACH_FORM);
  Set(XmNleftAttachment,XmATTACH_FORM);
  Set(XmNrightAttachment,XmATTACH_FORM);
  Set(XmNbottomAttachment,XmATTACH_FORM);
  render_widget = XtCreateManagedWidget("SurfaceBox_render_widget",
                                        glwMDrawingAreaWidgetClass,
                                        form,arglist,c);
  XVisualInfo *glx_visinfo;
  XtVaGetValues(render_widget,GLwNvisualInfo,&glx_visinfo,NULL);

  if(!glx_visinfo) {
	cerr << "FATAL ERROR: could not get GL visual info" << endl;
        delete this;
        return;
  }

  glx_context = glXCreateContext(XtDisplay(render_widget),
                                 glx_visinfo,0,GL_TRUE);

  //
  // Btn1 - rotate
  // Shift + Btn1 Zoom
  // Shift+Btn2 - select
  // Shift + Btn3 - selectKeep (Lei 04/17/2003)
  //
  sprintf(tstring,
        "Shift<Btn1Down>:sbStartClick(%ld,1)\n\
         Shift<Btn1Motion>:sbMotion(%ld,1)\n\
         Shift<Btn1Up>:sbEndClick(%ld,1)\n\
         <Btn1Down>:sbStartClick(%ld,0)\n\
         <Btn1Motion>:sbMotion(%ld,0)\n\
         <Btn1Up>:sbEndClick(%ld,0)\n\
         Shift<Btn2Down>:sbStartClick(%ld,2)\n\
         Shift<Btn2Motion>:sbMotion(%ld,2)\n\
         Shift<Btn2Up>:sbEndClick(%ld,2)\n\
         <Btn2Down>:sbStartClick(%ld,3)\n\
         <Expose>:sbExpose(%ld)\n\
         <Configure>:sbResize(%ld)",
         (u_long)this,(u_long)this,(u_long)this,
         (u_long)this,(u_long)this,(u_long)this,
         (u_long)this, (u_long)this,(u_long)this,
         (u_long)this, (u_long)this,(u_long)this);

  ftrans = XtParseTranslationTable(tstring);
  XtOverrideTranslations(render_widget,ftrans);

  for(i=0;i<XtNumber(style_items);i++)
    if(style_items[i].type != KDSubMenu)
      style_items[i].cb_second_param_or_nitems = (void*)this;

  for(i=0;i<XtNumber(surf_items);i++)
    if(surf_items[i].type != KDSubMenu)
      surf_items[i].cb_second_param_or_nitems = (void*)this;

  pop = KDShell::AddPopup(render_widget,surf_items,XtNumber(surf_items));

  SetHome();

  if(lsurf) 
     addSurf(lsurf);

VOUT("SurfaceBox::SurfaceBox(Widget parent, LoadedSurface *lsurf,")
}


SurfaceBox::~SurfaceBox()
{
VIN("SurfaceBox::~SurfaceBox()")
  if(m_alphaShell) {
    m_alphaShell->Hide();
    delete m_alphaShell;
  }

  for(int i=0;i<m_numSurfaces;i++)
      if(m_Surfaces[i])
         m_Surfaces[i]->unreference();

  XtDestroyWidget(form);
VOUT("SurfaceBox::~SurfaceBox()")
}


void SurfaceBox::project(float wx, float wy, float wz,
                         float &sx, float &sy, float &sz)
{
VIN("void SurfaceBox::project(float wx, float wy, float wz,")
  double x,y,z;
  gluProject(wx,wy,wz,m_modelview,m_projection,m_viewport,&x,&y,&z);
  sx = (float)x;
  sy = (float)m_height - (float)y - 1.f;
  sz = (float)z;
VOUT("void SurfaceBox::project(float wx, float wy, float wz,")
}

void SurfaceBox::reverseProject(float sx, float sy, float sz,
                                float &wx, float &wy, float &wz)
{
VIN("void SurfaceBox::reverseProject(float sx, float sy, float sz,")
  double x,y,z;
  gluUnProject(sx,sy,sz,m_modelview,m_projection,m_viewport,&x,&y,&z);
  wx = (float)x;
  wy = (float)y;
  wz = (float)z;
VOUT("void SurfaceBox::reverseProject(float sx, float sy, float sz,")
}


void SurfaceBox::addSurf(LoadedSurface *lS)
{
VIN("void SurfaceBox::addSurf(LoadedSurface *lS)")
  if(!lS || !lS->glsurface())
     return;

  for(int i=0;i<m_numSurfaces;i++)
     if(m_Surfaces[i] == lS)
        return;

  if(m_numSurfaces >= (MAX_SURFACES_SB-1)) {
	cerr << "ERROR: AddSurface -- too many surfaces in view" << endl;
	return;
	}

  if(m_alphaShell) {
     m_alphaShell->Hide();
     delete m_alphaShell;
  }

  sprintf(&(m_names[m_numSurfaces][0]),"%s",App->FilenameTail(lS->Filename()));

  lS->reference();
  m_Surfaces[m_numSurfaces++] = lS;
  m_curSurface = m_numSurfaces - 1;

  calculateCenter();
  ViewAll();
VOUT("void SurfaceBox::addSurf(LoadedSurface *lS)")
}


void SurfaceBox::delSurf(LoadedSurface *S)
{
VIN("void SurfaceBox::delSurf(LoadedSurface *S)")
  int i,j;

  for(i=0;i<m_numSurfaces;i++)
	if(m_Surfaces[i] == S)
		break;

  if(i>=m_numSurfaces)
	return;

  S->unreference();

  for(j=i+1;j<m_numSurfaces;j++)
	m_Surfaces[j-1] = m_Surfaces[j];

  m_numSurfaces--;
  calculateCenter();
  redraw();
VOUT("void SurfaceBox::delSurf(LoadedSurface *S)")
}


