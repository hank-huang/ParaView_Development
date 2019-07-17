#include <stdio.h>
#include <sys/param.h>
#include <KDApplication.h>
#include <Xm/Form.h>
#include <TDL/Surface.h>
#include <SurfaceBox.h>
#include <ImagePane.h>
#include <LoadedSurface.h>
#include <LoadedVolume.h>
#include <MaskCreator.h>
#include <GL/gl.h>
#include <GL/glu.h>

DEFINECB(MaskCreator,resetSurfaceCB,resetSurface)
DEFINECB(MaskCreator,selectAllCB,selectAll)
DEFINECB(MaskCreator,undoDrawCB,undoDraw)
DEFINECB(MaskCreator,drawToggleCB,toggleDraw)
DEFINECB(MaskCreator,eraseToggleCB,toggleErase)
DEFINECB(MaskCreator,rotateToggleCB,toggleRotate)

static void redrawCB(u_long dat)
{
VIN("static void redrawCB(u_long dat)")
  if(dat) ((MaskCreator*)dat)->redraw();
VOUT("static void redrawCB(u_long dat)")
}

static bool mouseCB(u_long dat, int x, int y, int m, int tpe)
{
VIN("static bool mouseCB(u_long dat, int x, int y, int m, int tpe)")
  bool rv = false;
  if(dat) {
     MaskCreator *M = (MaskCreator*)dat;
     switch(tpe) {
	case 0: rv = M->mouseDown(x,y,m); break;
	case 1: rv = M->mouseMotion(x,y,m); break;
	case 2: rv = M->mouseUp(x,y,m); break;
     }
  }
VOUT("static bool mouseCB(u_long dat, int x, int y, int m, int tpe)")
  return(rv);
}

MaskCreator::MaskCreator(LoadedSurface &surf) :
    KDChildWindow(App->FilenameTail(surf.Filename()),200,200)
{
VIN("MaskCreator::MaskCreator(LoadedSurface &surf) :")
   m_surface    = &surf;
   m_surfaceBox = NULL;
   m_penWidth   = 5;
   m_bDrawing   = false;
   m_bErasing   = false;
   m_bRotating  = true;
 
   if(!m_surface) {
      delete this;
      return;
   }

   m_numVert = m_surface->surface()->getNumVert();

   if(m_surface->mask().getNelm() != m_numVert)
	m_surface->defaultMask();

   m_maskUndo.setDim(m_numVert);
   m_maskUndo = m_surface->mask();

   m_surface->reference();

   ATTACHMENTS_TBLR(0,100,0,75)
   m_surfaceBox = new SurfaceBox(this->Form(),m_surface,arglist,c);

   ATTACHMENTS_TBLR(0,100,77,100)
   Widget form = XmCreateForm(this->Form(),"MaskCreatorMenu",arglist,c);
   XtManageChild(form);
   XmChangeColor(form,LABELCOLOR);

   Widget w1,w2;
   int bh = App->BigFontHeight() + 5;

   ATTACHMENTS_TLRH(0,0,100,bh)
   w1 = App->Button(form,"Clear All Vertices",App->MedFontList(),
                    resetSurfaceCB,(XtPointer)this,arglist,c);

   ATTACHMENTS_WLRH(w1,0,100,bh)
   w1 = App->Button(form,"Select All Vertices",App->MedFontList(),
                    selectAllCB,(XtPointer)this,arglist,c);

   ATTACHMENTS_WLRH(w1,0,100,bh)
   w1 = App->LabelText(form,"Pen Width","20",10,
                       App->SmallFontList(),arglist,c);
   m_penWidthText = w1;

   ATTACHMENTS_WLRH(w1,0,100,bh)
   w1 = App->Button(form,"Undo Draw",App->MedFontList(),
                    undoDrawCB,(XtPointer)this,arglist,c);

   ATTACHMENTS_WLRH(w1,0,100,bh)
   w1 = App->Toggle(form,"Select Vertices",App->MedFontList(),m_bDrawing,
		    drawToggleCB,(XtPointer)this,arglist,c);
   m_drawToggle = w1;

   ATTACHMENTS_WLRH(w1,0,100,bh)
   w1 = App->Toggle(form,"Unselect Vertices",App->MedFontList(),m_bErasing,
		    eraseToggleCB,(XtPointer)this,arglist,c);
   m_eraseToggle = w1;

   ATTACHMENTS_WLRH(w1,0,100,bh)
   w1 = App->Toggle(form,"Rotate",App->MedFontList(),m_bRotating,
		    rotateToggleCB,(XtPointer)this,arglist,c);
   m_rotateToggle = w1;

   m_surfaceBox->SetRedrawFunction(redrawCB,(u_long)this);
   m_surfaceBox->SetMouseFunction(mouseCB,(u_long)this);

   Show();
VOUT("MaskCreator::MaskCreator(LoadedSurface &surf) :")
}

MaskCreator::~MaskCreator()
{
VIN("MaskCreator::~MaskCreator()")
  if(m_surface)
      m_surface->unreference();

  if(m_surfaceBox)
	delete m_surfaceBox;
VOUT("MaskCreator::~MaskCreator()")
}

void MaskCreator::redraw()
{
VIN("void MaskCreator::redraw()")
  int i;
  float x,y,z;

  // we know surfaceBox GL matrices are valid here
  glPushAttrib(GL_LIGHTING_BIT);
  glDisable(GL_LIGHTING);

  glPointSize(5.f);
  glColor3f(1.f,0.f,0.8f);
  glBegin(GL_POINTS);

  for(i=0;i<m_numVert;i++) 
    if(m_surface->mask()[i]) {
       x = m_surface->surface()->vertices()[i].x();
       y = m_surface->surface()->vertices()[i].y();
       z = m_surface->surface()->vertices()[i].z();
       glVertex3f(x,y,z);
    }
  glEnd();
  glLineWidth(1.f);
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  glPopAttrib();
VOUT("void MaskCreator::redraw()")
}


bool MaskCreator::mouseDown(int x, int y, int b)
{
VIN("bool MaskCreator::mouseDown(int x, int y, int b)")
  if(m_bRotating) return(true);

  m_maskUndo = m_surface->mask();

  m_penWidth = App->GetInt(m_penWidthText);
  if((m_penWidth < 1)||(m_penWidth > 100)) {
      m_penWidth = 20;
      App->SetTextString(m_penWidthText,"20");
      App->ShowMessage("Pen width must be between 1 and 100");
  }
  m_surfaceBox->setCurrent();

  m_mouseX = x;
  m_mouseY = y;
  m_width  = m_surfaceBox->getWidth();
  m_height = m_surfaceBox->getHeight();

VOUT("bool MaskCreator::mouseDown(int x, int y, int b)")
  return(false);
}

bool MaskCreator::mouseMotion(int x, int y, int b)
{
VIN("bool MaskCreator::mouseMotion(int x, int y, int b)")
  if(m_bRotating) return(true);

  glPushAttrib(GL_LIGHTING_BIT);
  glDisable(GL_LIGHTING);

  int buf;
  int wo2 = m_penWidth/2;

  glGetIntegerv(GL_DRAW_BUFFER,&buf);
  glDrawBuffer(GL_FRONT);

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  gluOrtho2D(0.0,m_width,0.0,m_height);

  glColor3f(1.f,0.f,0.f);
  glPolygonMode(GL_FRONT,GL_FILL);
/*
  glLineWidth(m_penWidth);
  glBegin(GL_LINES);
  glVertex2i(m_mouseX,m_height - m_mouseY - 1);
  glVertex2i(x,m_height - y - 1);
  glEnd();
*/
  y = m_height - y - 1;
  glRecti(x-wo2,y-wo2,x+wo2,y+wo2);

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  glDrawBuffer((GLenum)buf);

  glPopAttrib();

  m_mouseX = x;
  m_mouseY = y;

VOUT("bool MaskCreator::mouseMotion(int x, int y, int b)")
  return(false);
}

bool MaskCreator::mouseUp(int x, int y, int b)
{
VIN("bool MaskCreator::mouseUp(int x, int y, int b)")
  if(m_bRotating) return(true);

  updateMask();

  m_surfaceBox->redraw();

VOUT("bool MaskCreator::mouseUp(int x, int y, int b)")
  return(false);
}

void MaskCreator::resetSurface()
{
VIN("void MaskCreator::resetSurface()")
  if(App->GetYesNo("Reset surface to no vertices tagged. Are you sure?")) {
	m_surface->mask() = (u_char)0;
        m_maskUndo        = (u_char)0;
        m_surfaceBox->redraw();
  }
VOUT("void MaskCreator::resetSurface()")
}

void MaskCreator::selectAll()
{
VIN("void MaskCreator::selectAll()")
  m_maskUndo = m_surface->mask();
  m_surface->mask() = (u_char)255;
  m_surfaceBox->redraw();
VOUT("void MaskCreator::selectAll()")
}

void MaskCreator::undoDraw()
{
VIN("void MaskCreator::undoDraw()")
  m_surface->mask() = m_maskUndo;
  m_surfaceBox->redraw();
VOUT("void MaskCreator::undoDraw()")
}

void MaskCreator::updateMask()
{
VIN("void MaskCreator::updateMask()")
  int i,ix,iy,buf;
  float x,y,z,sx,sy,sz;
  float pixel[3];

  u_char val;

  if(m_bDrawing) val = 255U;
  else           val = 0U;

  glGetIntegerv(GL_READ_BUFFER,&buf);
  glReadBuffer(GL_FRONT);

  for(i=0;i<m_numVert;i++) {
      x = m_surface->surface()->vertices()[i].x();
      y = m_surface->surface()->vertices()[i].y();
      z = m_surface->surface()->vertices()[i].z();
      m_surfaceBox->project(x,y,z,sx,sy,sz);
      ix = (int)sx;
      iy = m_height - (int)sy - 1;
      glReadPixels(ix,iy,1,1,GL_RGB,GL_FLOAT,pixel);

      // Look for 'red' pixels that were drawn with mouse in mouseMotion()
      if((pixel[0] >0.95f)&&(pixel[1]<0.05f)&&(pixel[2]<0.05f)) {
	m_surface->mask()[i] = val;
      }
  }

  glReadBuffer((GLenum)buf);
VOUT("void MaskCreator::updateMask()")
}

void MaskCreator::toggleDraw()
{
VIN("void MaskCreator::toggleDraw()")
  m_bDrawing  = App->GetToggleValue(m_drawToggle);
  if(m_bDrawing) { m_bErasing = m_bRotating = false;       }
  else           { m_bErasing = false; m_bRotating = true; }

  App->SetToggleValue(m_eraseToggle,m_bErasing);
  App->SetToggleValue(m_rotateToggle,m_bRotating);
VOUT("void MaskCreator::toggleDraw()")
}


void MaskCreator::toggleErase()
{
VIN("void MaskCreator::toggleErase()")
  m_bErasing  = App->GetToggleValue(m_eraseToggle);
  if(m_bErasing) { m_bDrawing = m_bRotating = false;       }
  else           { m_bDrawing = false; m_bRotating = true; }

  App->SetToggleValue(m_drawToggle,m_bDrawing);
  App->SetToggleValue(m_rotateToggle,m_bRotating);
VOUT("void MaskCreator::toggleErase()")
}


void MaskCreator::toggleRotate()
{
VIN("void MaskCreator::toggleRotate()")
  m_bRotating  = App->GetToggleValue(m_rotateToggle);
  if(m_bRotating) { m_bErasing = m_bDrawing = false;       }
  else            { m_bErasing = false; m_bDrawing = true; }

  App->SetToggleValue(m_eraseToggle,m_bErasing);
  App->SetToggleValue(m_drawToggle,m_bDrawing);
VOUT("void MaskCreator::toggleRotate()")
}
