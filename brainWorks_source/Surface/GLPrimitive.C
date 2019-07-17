#include <iostream.h>
#include <GL/gl.h>
#include <TDL/Line.h>
#include <TDL/Point.h>
#include <KDApplication.h>
#include "GLPrimitive.h"


void GLPrimitive::clearList(GLPrimitive *list)
{
VIN("void GLPrimitive::clearList(GLPrimitive *list)")
  GLPrimitive *p = list;
  while(p) {
      list = list->_next;
      delete p;
      p = list;
  }
VOUT("void GLPrimitive::clearList(GLPrimitive *list)")
}

GLPrimitive *GLPrimitive::addToList(GLPrimitive *list, GLPrimitive *add)
{
VIN("GLPrimitive *GLPrimitive::addToList(GLPrimitive *list, GLPrimitive *add)")
if(!add) return(list);
add->_next = list;
VOUT("GLPrimitive *GLPrimitive::addToList(GLPrimitive *list, GLPrimitive *add)")
return(add);
}

GLPrimitive *GLPrimitive::removeFromList(GLPrimitive *list, GLPrimitive *rem)
{
VIN("GLPrimitive *GLPrimitive::removeFromList()")
if(!rem)  return(list);
if(!list) return(NULL);

GLPrimitive *first=list;
if(list == rem) {
	list = list->_next;
	delete rem;
	return(list);
	}

for(;list->_next;list=list->_next) 
	if(list->_next == rem) {
		list->_next = rem->_next;
		delete rem;
		//return(list);
		return(first);
		}
VOUT("GLPrimitive *GLPrimitive::removeFromList()")
return(list);
}

GLPrimitive::GLPrimitive(const Array1D<Point> &p)
{
VIN("GLPrimitive::GLPrimitive(const Array1D<Point> &p)")
_type      = GLLineStrip;
_vertices  = NULL;
_nvertices = 0;
_line_width = 1.0;
_next      = NULL;
_hook      = NULL;

if(p.getNelm() <= 0) return;

_nvertices = p.getNelm();
_vertices  = new GLfloat[3*_nvertices];
if(!_vertices) {
	cerr << "ERROR: GLPrimitive(Line) out of memory!" << endl;
	return;
	}
int i;
for(i=0;i<_nvertices;i++) {
	_vertices[3*i+0] = p[i].x();
	_vertices[3*i+1] = p[i].y();
	_vertices[3*i+2] = p[i].z();
	}

_colr = _colg = _colb = 1.0;
VOUT("GLPrimitive::GLPrimitive(const Array1D<Point> &p)")
}


GLPrimitive::GLPrimitive(Line *l)
{
VIN("GLPrimitive::GLPrimitive(Line *l)")
_type      = GLLine;
_vertices  = NULL;
_nvertices = 0;
_line_width = 1.0;
_next      = NULL;
_hook      = NULL;

if(l->numSegments() <= 0) return;
_nvertices = l->numSegments();

_vertices = new GLfloat[6*_nvertices];
if(!_vertices) {
	cerr << "ERROR: GLPrimitive(Line) out of memory!" << endl;
	return;
	}
int i;
Point a,b;
for(i=0;i<l->numSegments();i++) {
	l->getSegment(a,b,i);	
	_vertices[6*i+0] = a.x();
	_vertices[6*i+1] = a.y();
	_vertices[6*i+2] = a.z();
	_vertices[6*i+3] = b.x();
	_vertices[6*i+4] = b.y();
	_vertices[6*i+5] = b.z();
	}

_colr = _colg = _colb = 1.0;
VOUT("GLPrimitive::GLPrimitive(Line *l)")
}



GLPrimitive::~GLPrimitive()
{
VIN("GLPrimitive::~GLPrimitive()")
if((_type == GLLine)||(_type == GLLineStrip)) {
	if(_vertices) delete [] _vertices;
	}
else if(_type == GLSphere) {
	if(_quadric) gluDeleteQuadric(_quadric);
	}
VOUT("GLPrimitive::~GLPrimitive()")
}

void GLPrimitive::getExtents(float &x, float &y, float &z,
			float &cx, float &cy, float &cz)
{
VIN("void GLPrimitive::getExtents(float &x, float &y, float &z,")
int i;
float minx,miny,minz,maxx,maxy,maxz;
float xx,yy,zz;
float sumx,sumy,sumz,nsum;

switch(_type) {
	case GLLine:
	minx = miny = minz = 1E20;
	maxx = maxy = maxz = -1E20;
	sumx = sumy = sumz = nsum = 0.0;
	for(i=0;i<_nvertices;i++) {
		xx = _vertices[6*i+0];
		yy = _vertices[6*i+1];
		zz = _vertices[6*i+2];
		sumx += xx; sumy += yy; sumz += zz; nsum += 1.0;
		if(xx < minx) minx = xx;
		if(yy < miny) miny = yy;
		if(zz < minz) minz = zz;
		if(xx > maxx) maxx = xx;
		if(yy > maxy) maxy = yy;
		if(zz > maxz) maxz = zz;
		xx = _vertices[6*i+3];
		yy = _vertices[6*i+4];
		zz = _vertices[6*i+5];
		sumx += xx; sumy += yy; sumz += zz; nsum += 1.0;
		if(xx < minx) minx = xx;
		if(yy < miny) miny = yy;
		if(zz < minz) minz = zz;
		if(xx > maxx) maxx = xx;
		if(yy > maxy) maxy = yy;
		if(zz > maxz) maxz = zz;
		}
	x = maxx - minx;
	y = maxy - miny;
	z = maxz - minz;
	cx = sumx/nsum;
	cy = sumy/nsum;
	cz = sumz/nsum;
	break;
	case GLLineStrip:	
	minx = miny = minz = 1E20;
	maxx = maxy = maxz = -1E20;
	sumx = sumy = sumz = nsum = 0.0;
	for(i=0;i<_nvertices;i++) {
		xx = _vertices[3*i+0];
		yy = _vertices[3*i+1];
		zz = _vertices[3*i+2];
		sumx += xx; sumy += yy; sumz += zz; nsum += 1.0;
		if(xx < minx) minx = xx;
		if(yy < miny) miny = yy;
		if(zz < minz) minz = zz;
		if(xx > maxx) maxx = xx;
		if(yy > maxy) maxy = yy;
		if(zz > maxz) maxz = zz;
		}
	x = maxx - minx;
	y = maxy - miny;
	z = maxz - minz;
	cx = sumx/nsum;
	cy = sumy/nsum;
	cz = sumz/nsum;
	break;
	case GLSphere:
	minx = _centx - _radius;
	miny = _centy - _radius;
	minz = _centz - _radius;
	maxx = _centx + _radius;
	maxy = _centy + _radius;
	maxz = _centz + _radius;
	x = maxx - minx;
	y = maxy - miny;
	z = maxz - minz;
	cx = _centx;
	cy = _centy;
	cz = _centz;
	break;
	}
VOUT("void GLPrimitive::getExtents(float &x, float &y, float &z,")
}



void GLPrimitive::render()
{
VIN("void GLPrimitive::render()")
int i;
int islight = glIsEnabled(GL_LIGHTING);

switch(_type) {
	case GLLine:
	   glDisable(GL_LIGHTING);
	   glLineWidth(_line_width);
	   glBegin(GL_LINES);
	   for(i=0;i<_nvertices;i++) {
		glColor3f(_colr,_colg,_colb);
		glVertex3fv(&(_vertices[6*i]));
		glColor3f(_colr,_colg,_colb);
		glVertex3fv(&(_vertices[6*i+3]));
		}
	   glEnd();
	   break;
	case GLLineStrip:
	   glDisable(GL_LIGHTING);
	   glLineWidth(_line_width);
	   glBegin(GL_LINE_STRIP);
	   for(i=0;i<_nvertices;i++) {
		glColor3f(_colr,_colg,_colb);
		glVertex3fv(&(_vertices[3*i]));
		}
	   glEnd();
	   break;
	case GLSphere:
	   glEnable(GL_LIGHTING);
	   glMatrixMode(GL_MODELVIEW);
	   glPushMatrix();
	   glTranslatef(_centx,_centy,_centz);
	   glColor3f(_colr,_colg,_colb);
	   gluSphere(_quadric,_radius,10,10);
	   glPopMatrix();
	   break;
	}
if(islight)
  glEnable(GL_LIGHTING);

VOUT("void GLPrimitive::render()")
}
