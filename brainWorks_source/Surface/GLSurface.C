#include <iostream.h>
#include <X11/Intrinsic.h>
#include <TDL/Surface.h>
#include <KDApplication.h>
#include <GLSurface.h>

#define STRIP_HEAP_INCREMENT	50

//
// TriangleFace (internal)
//    linked list of all triangles
//    for generating triangle strips
//
typedef struct _TriangleFace {
        int    _used;
        int    _01,_12,_20,_adj;
	int    _next,_exit_edge;
        } TriangleFace;


void GLSurface::makeStripped()
{
VIN("void GLSurface::makeStripped()")
if(isstripped) return;
VOUT("void GLSurface::makeStripped()")
}

void GLSurface::setSurface(Surface *surface, int do_strips)
{
VIN("void GLSurface::setSurface(Surface *surface, int do_strips)")
isstripped = do_strips;

int i,j;

_facets.setDim(0,0);
_facets=(surface->facets());

if(!surface->hasUNorms()) 
    surface->genUnitNormals();

npoly = surface->getNumPoly();
selected.setDim(npoly);
if(npoly != 0) selected=0; 
 
if(surface->getNumVert() != nvert) {
	nvert = surface->getNumVert();
	if(vertices) delete [] vertices;
	if(normals)  delete [] normals;
	vertices = new GLfloat[3*nvert];
	normals  = new GLfloat[3*nvert];
	if(!vertices || !normals) {
		cerr << "ERROR: GLSurface::setSurface() out of memory!" << endl;
		cleanup();
		return;
		}
        if(_colors) delete [] _colors;
        _colors = NULL;

        if(nvert>0) 
          _colors = new GLfloat[4*nvert];

        for(i=0;i<nvert;i++) {
	   _colors[4*i+0] = _red;
	   _colors[4*i+1] = _green;
	   _colors[4*i+2] = _blue;
           _colors[4*i+3] = _alpha;
	   }
	}

for(i=0,j=0;i<nvert;i++) {
	vertices[j++] = surface->vertices()[i].x();
	vertices[j++] = surface->vertices()[i].y();
	vertices[j++] = surface->vertices()[i].z();
	}
normal_mode = GLSurface::PerPoint;

for(i=0,j=0;i<nvert;i++) {
	normals[j++] = surface->unormals()[i].x();
	normals[j++] = surface->unormals()[i].y();
	normals[j++] = surface->unormals()[i].z();
	}
 if(surface->hasCurvature()){
   _Curvature.setDim(0);
   _Curvature = surface->curvature();
 }
 else _Curvature.setDim(0);

cleanStrips();

if(surface) {
	_extentx = surface->findXMax() - surface->findXMin();
	_extenty = surface->findYMax() - surface->findYMin();
	_extentz = surface->findZMax() - surface->findZMin();

        Point xx;
        if(!surface->getSimpleCentroid(&xx)) {
           _centx   = xx.x();
           _centy   = xx.y();
           _centz   = xx.z();
           }
	else _centx = _centy = _centz = 0.0;
	}
else 	{
	_extentx = _extenty = _extentz = 100.0;
	_centx = _centy = _centz = 0.0;
	}

_max_extent = MAX3(_extentx,_extenty,_extentz);

if(do_strips) 
   generateStrips(surface);
VOUT("void GLSurface::setSurface(Surface *surface, int do_strips)")
}

void GLSurface::saveColors(ostream &fp)
{
VIN("void GLSurface::saveColors(ostream &fp)")
  int i,p = fp.precision(3);
  for(i=0;i<nvert;i++) {
     fp << _colors[4*i+0] << " " <<_colors[4*i+1] << " " 
        << _colors[4*i+2] << endl;
  }
  fp.precision(p);
VOUT("void GLSurface::saveColors(ostream &fp)")
}



bool GLSurface::loadColors(istream &fp)
{
VIN("bool GLSurface::loadColors(istream &fp)")
  int i;
  for(i=0;i<nvert;i++) {
     fp >> _colors[4*i+0] >> _colors[4*i+1] >> _colors[4*i+2];
     _colors[4*i+3] = _alpha;

     if(fp.fail()) break;
  }
VOUT("bool GLSurface::loadColors(istream &fp)")
  return(i>=nvert);
}



bool GLSurface::loadIntensity(istream &fp, float ff)
{
  float intens;
  int i;
  if(nvert <= 0) return(true);

  Array1D<double> func(nvert);
  double minfunc,maxfunc,rangefunc;
  float  a,b,c,g;
  float bgr,bgg,bgb;

  minfunc =  1E20;
  maxfunc = -1E20;
  for(i=0;i<nvert;i++) {
     fp >> func[i];
     if(fp.fail()) break;

     if(func[i] < minfunc) minfunc = func[i];
     if(func[i] > maxfunc) maxfunc = func[i];
  }
  if(i<nvert) {
     cerr << "ERROR: LoadAlpha : failed - could not load array" << endl;
     return(false);
  }

  rangefunc = maxfunc - minfunc;
  if(rangefunc == 0.) rangefunc = 1.;

  bgr = 0.6;
  bgg = 0.6;
  bgb = 0.6;

  for(i=0;i<nvert;i++) {
     intens = (float)((func[i]-minfunc)/rangefunc);
     a = _colors[4*i+0];
     b = _colors[4*i+1];
     c = _colors[4*i+2];

     //g = intens;
     g = ff*intens*intens*intens;
     if(g > 1.f) g = 1.f;
     if(g < 0.f) g = 0.f;
 
     a = g * a + (1.f - g) * bgr;
     b = g * b + (1.f - g) * bgg;
     c = g * c + (1.f - g) * bgb;

     _colors[4*i+0] = a;
     _colors[4*i+1] = b;
     _colors[4*i+2] = c;
  }

  return(true);
}




void GLSurface::getVertexColor(int idx, float &r, float &g, float &b)
{
VIN("void GLSurface::getVertexColor(int idx, float &r, float &g, float &b)")
  if((idx>=0)&&(idx<nvert)) {
	r = _colors[4*idx+0];
	g = _colors[4*idx+1];
	b = _colors[4*idx+2];
  } else {
    r = g = b = 0.f;
  }
VOUT("void GLSurface::getVertexColor(int idx, float &r, float &g, float &b)")
}

void GLSurface::setVertexColor(int idx, float r, float g, float b)
{
VIN("void GLSurface::setVertexColor(int idx, float r, float g, float b)")
  if((idx>=0)&&(idx<nvert)) {
	_colors[4*idx+0] = r;
	_colors[4*idx+1] = g;
	_colors[4*idx+2] = b;
	_colors[4*idx+3] = _alpha;
	}
VOUT("void GLSurface::setVertexColor(int idx, float r, float g, float b)")
}

void GLSurface::setAlpha(float alpha)
{
VIN("void GLSurface::setAlpha(float alpha)")
  _alpha = alpha;
  for(int i=0;i<nvert;i++)
     _colors[4*i+3] = _alpha;
VOUT("void GLSurface::setAlpha(float alpha)")
}

//////// added by Tong
void GLSurface::setColorSurface(float r, float g, float b, float alpha)
{
VIN("void GLSurface::setColorSurface(float r, float g, float b, float alpha)")
int i,j;
_red = r;
_green = g;
_blue = b;
_alpha = alpha;
 
for(i=0,j=0;i<nvert;i++) {
	_colors[j++] = r;
	_colors[j++] = g;
	_colors[j++] = b;
	_colors[j++] = _alpha;
	}
VOUT("void GLSurface::setColorSurface(float r, float g, float b, float alpha)")
}

void GLSurface::setColorSurface(GLSurface &S)
{
VIN("void GLSurface::setColorSurface(GLSurface &S)")
  if(nvert != S.nvert) return;

  for(int i=0;i<nvert*4;i++)
	_colors[i] = S._colors[i];
VOUT("void GLSurface::setColorSurface(GLSurface &S)")
}


void GLSurface::setColorSurface(Array1D<float> &r, Array1D<float> &g, Array1D<float> &b)
{
VIN("void GLSurface::setColorSurface(Array1D<float> &r, Array1D<float> &g, Array1D<float> &b)")
int np = nvert;

if((r.getNelm() != np)||(g.getNelm() != np)||(b.getNelm() != np)) {
      cerr << "WARNING: GLSurface::colorSurface() invalid array size" << endl;
      return;
      }

int i,j;
for(i=0,j=0;i<nvert;i++) {
	_colors[j++] = r[i];
	_colors[j++] = g[i];
	_colors[j++] = b[i];
        _colors[j++] = _alpha;
	}
VOUT("void GLSurface::setColorSurface(Array1D<float> &r, Array1D<float> &g, Array1D<float> &b)")
}
////////////////////////

void GLSurface::cleanup()
{
VIN("void GLSurface::cleanup()")
if(vertices) delete [] vertices;
vertices = NULL;
if(normals)  delete [] normals;
normals = NULL;

if(_colors) delete _colors;
_colors = NULL;

cleanStrips();

nvert = npoly = 0;
VOUT("void GLSurface::cleanup()")
}


void GLSurface::cleanStrips()
{
VIN("void GLSurface::cleanStrips()")
if(strip_heap) {
	int i,j;
	for(i=0;i<nstrips;i++)
	    if(strip_heap[i].vertices)
		delete [] strip_heap[i].vertices;
	delete [] strip_heap;
	}

strip_heap  = NULL;
nstrips     = 0;
nstrip_heap = 0;
VOUT("void GLSurface::cleanStrips()")
}


int GLSurface::newStripIndex()
{
VIN("int GLSurface::newStripIndex()")
if(nstrips >= nstrip_heap) {
	nstrip_heap += STRIP_HEAP_INCREMENT;
	TriangleStrip *new_strips = new TriangleStrip[nstrip_heap];
	int i;
	for(i=0;i<nstrips;i++) {
		new_strips[i].vertices  = strip_heap[i].vertices;
		new_strips[i].nvertices = strip_heap[i].nvertices;
		}
	if(strip_heap) delete [] strip_heap;
	strip_heap = new_strips;
	}
strip_heap[nstrips].vertices  = NULL;
strip_heap[nstrips].nvertices = 0;

VOUT("int GLSurface::newStripIndex()")
return(nstrips++);
}




//
// return the adjacent face with least
// adjacency count
//
int GLSurface::findNextFace(void *dat, int cur)
{
VIN("int GLSurface::findNextFace(void *dat, int cur)")
TriangleFace *f = (TriangleFace*)dat;
if(!f) return(-1);

int nvalid=0;
int f01,a01;
int f12,a12;
int f20,a20;
int rv;

f01 = f[cur]._01;
f12 = f[cur]._12;
f20 = f[cur]._20;

if((f01 >= 0)&&(!f[f01]._used))  {
	a01 = f[cur]._adj;
	nvalid++;
	}
else	a01 = 10000;
if((f12 >= 0)&&(!f[f12]._used))  {
	a12 = f[cur]._adj;
	nvalid++;
	}
else	a12 = 10000;
if((f20 >= 0)&&(!f[f20]._used))  {
	a20 = f[cur]._adj;
	nvalid++;
	}
else	a20 = 10000;

if(nvalid==0)
	return(-1);

if(a01 < a12) {
	if(a01 < a20)  {
		rv = f01;
		f[cur]._exit_edge = 1;
		}
	else	{
		if(a20 < a12) {
			rv = f20;
			f[cur]._exit_edge = 3;
			}
		else	{
			rv = f12;
			f[cur]._exit_edge = 2;
			}
		}
	}
else    {
	if(a12 < a20)  {
		rv = f12;
		f[cur]._exit_edge = 2;
		}
	else 	{
		if(a20 < a01) {
			rv = f20;
			f[cur]._exit_edge = 3;
			}
		else	{
			rv = f01;
			f[cur]._exit_edge = 1;
			}
		}
	}
return(rv);
VOUT("int GLSurface::findNextFace(void *dat, int cur)")
}

void GLSurface::generateStrips(Surface *surface)
{
VIN("void GLSurface::generateStrips(Surface *surface)")
cleanStrips();
if(!surface || (surface->getNumPoly() <= 0))
	return;

int i,j;
int a,b,c;
int aa,bb,cc;

npoly = surface->getNumPoly();

//
// generate linked list of all
// triangles in surface
//
TriangleFace *face_heap  = new TriangleFace[npoly];
for(i=0;i<npoly;i++) {
	face_heap[i]._01         = -1;
	face_heap[i]._12         = -1;
	face_heap[i]._20         = -1;
	face_heap[i]._adj        = 0;
	face_heap[i]._used       = 0;
	face_heap[i]._next       = -1;
	face_heap[i]._exit_edge  = -1;
	}
//
// identify adjacent triangles
//
int edgeno;
for(i=0;i<npoly;i++) {
	a = surface->facets()[i][0];
	b = surface->facets()[i][1];
	c = surface->facets()[i][2];
	//
	// find poly sharing edge a,b
	//
	if(face_heap[i]._01 < 0) {
		edgeno = 0;
		for(j=0;j<npoly;j++) {
			if(j==i) continue;
			aa = surface->facets()[j][0];
			bb = surface->facets()[j][1];
			cc = surface->facets()[j][2];

			if(((aa==a)&&(bb==b))||
			   ((aa==b)&&(bb==a))) 
				edgeno = 1;   // edge 01->01
			else if(((cc==a)&&(bb==b))||
			   ((cc==b)&&(bb==a))) 
				edgeno = 2;   // edge 01->12
			else if(((aa==a)&&(cc==b))||
			   ((aa==b)&&(cc==a))) 
				edgeno = 3;   // edge 01->20
			if(edgeno) break;
			}
		if(j<npoly) {
			face_heap[i]._01 = j;
			face_heap[i]._adj++;
			//
			// link it back to this one
			//
			switch(edgeno) {
				case 1: face_heap[j]._01 = i; break;
				case 2: face_heap[j]._12 = i; break;
				case 3: face_heap[j]._20 = i; break;
				}
			face_heap[j]._adj++;
			}
		}
	//
	// now for b,c
	//
	if(face_heap[i]._12 < 0) {
		edgeno = 0;
		for(j=0;j<npoly;j++) {
			if(j==i) continue;
			aa = surface->facets()[j][0];
			bb = surface->facets()[j][1];
			cc = surface->facets()[j][2];
			if(((aa==c)&&(bb==b))||
			   ((aa==b)&&(bb==c))) 
				edgeno = 1;   // edge 12->01
			else if(((cc==c)&&(bb==b))||
			   ((cc==b)&&(bb==c))) 
				edgeno = 2;   // edge 12->12
			else if(((aa==c)&&(cc==b))||
			   ((aa==b)&&(cc==c))) 
				edgeno = 3;   // edge 12->20
			if(edgeno) break;
			}
		if(j<npoly) {
			face_heap[i]._12 = j;
			face_heap[i]._adj++;
			//
			// link it back to this one
			//
			switch(edgeno) {
				case 1: face_heap[j]._01 = i; break;
				case 2: face_heap[j]._12 = i; break;
				case 3: face_heap[j]._20 = i; break;
				}
			face_heap[j]._adj++;
			}
		}
	//
	// now for c,a
	//
	if(face_heap[i]._20 < 0) {
		edgeno = 0;
		for(j=0;j<npoly;j++) {
			if(j==i) continue;
			aa = surface->facets()[j][0];
			bb = surface->facets()[j][1];
			cc = surface->facets()[j][2];
			if(((aa==a)&&(bb==c))||
			   ((aa==c)&&(bb==a))) 
				edgeno = 1;   // edge 20->01
			else if(((cc==a)&&(bb==c))||
			   ((cc==c)&&(bb==a))) 
				edgeno = 2;   // edge 20->12
			else if(((aa==a)&&(cc==c))||
			   ((aa==c)&&(cc==a))) 
				edgeno = 3;   // edge 20->20
			if(edgeno) break;
			}
		if(j<npoly) {
			face_heap[i]._20 = j;
			face_heap[i]._adj++;
			//
			// link it back to this one
			//
			switch(edgeno) {
				case 1: face_heap[j]._01 = i; break;
				case 2: face_heap[j]._12 = i; break;
				case 3: face_heap[j]._20 = i; break;
				}
			face_heap[j]._adj++;
			}
		}
	}

int xx,curi,nexti,p;
int vertcnt,npolcur;
int tria,trib,tric;
int need_swap,stripi;
int sharea,shareb,stack1;
int nsingleton,maxstrip,sumstrip;

nsingleton = 0;
maxstrip   = 0;
sumstrip   = 0;

gen_strips:

//
// find an un-used triangle
//
for(i=0;i<npoly;i++) 
	if(!face_heap[i]._used) 
		break;
//
// generate a strip from it
//
if(i<npoly) {
	npolcur = 1;
	curi    = i;
	//
	// Remove poly from useable list
	//
	face_heap[curi]._used    = 1;
	//
	// Decrease adjacency count in
	// neighboring triangles
	//
	if((xx=face_heap[curi]._01)>=0)
		face_heap[xx]._adj--;
	if((xx=face_heap[curi]._12)>=0)
		face_heap[xx]._adj--;
	if((xx=face_heap[curi]._20)>=0)
		face_heap[xx]._adj--;
	//
	// Keep finding adjacent polys
	//

	while((nexti = findNextFace(face_heap,curi))>=0) {

		face_heap[curi]._next = nexti;
		curi = nexti;
		face_heap[curi]._used = 1;
//
// NOTE: wait until after strip is created!
//
		//
		// Decrease adjacency count in
		// neighboring triangles
		//
		if((xx=face_heap[curi]._01)>=0)
			face_heap[xx]._adj--;
		if((xx=face_heap[curi]._12)>=0)
			face_heap[xx]._adj--;
		if((xx=face_heap[curi]._20)>=0)
			face_heap[xx]._adj--;
		npolcur++;
		}
	face_heap[curi]._next = -1;

	if(maxstrip < npolcur)
		maxstrip = npolcur;
	sumstrip += npolcur;

	stripi = newStripIndex();
	strip_heap[stripi].vertices  = new int[npolcur*2+2];
	if(strip_heap[stripi].vertices == NULL) {
		cerr << "ERROR: GLSurface::generateStrips() out of memory" <<
			endl;
		cleanup();
		return;
		}

	curi = i;
	vertcnt = 0;
	//
	// Get first poly (a,b,c)
	//
	a = surface->facets()[curi][0];
	b = surface->facets()[curi][1];
	c = surface->facets()[curi][2];
	if(npolcur == 1) { // singleton
		nsingleton++;
		strip_heap[stripi].vertices[0] = a;
		strip_heap[stripi].vertices[1] = b;
		strip_heap[stripi].vertices[2] = c;
		strip_heap[stripi].nvertices   = 3;
		goto gen_strips;
		}

	// Add first triangle:
	//    first point is point not
	//    coincident with next face
	//
	switch(face_heap[curi]._exit_edge) {
		case 1: sharea = a; shareb = b;
			strip_heap[stripi].vertices[vertcnt++] = c;
			break;
		case 2: sharea = b; shareb = c;
			strip_heap[stripi].vertices[vertcnt++] = a;
			break;
		case 3: sharea = c; shareb = a;
		 	strip_heap[stripi].vertices[vertcnt++] = b;
			break;
		}
	strip_heap[stripi].vertices[vertcnt++] = sharea;
	strip_heap[stripi].vertices[vertcnt++] = shareb;

	p = 1;
	while((curi=face_heap[curi]._next)>=0) {
		a = surface->facets()[curi][0];
		b = surface->facets()[curi][1];
		c = surface->facets()[curi][2];

		switch(face_heap[curi]._exit_edge) {
			case 1: sharea = a; shareb = b; break;
			case 2: sharea = b; shareb = c; break;
			case 3: sharea = c; shareb = a; break;
			}

		if(p == 1) {
		   tria = strip_heap[stripi].vertices[vertcnt-1];
		   trib = strip_heap[stripi].vertices[vertcnt-2];
		   }
		else {
		   tria = strip_heap[stripi].vertices[vertcnt-2];
		   trib = strip_heap[stripi].vertices[vertcnt-1];
		   }
		//
		// find not-already-added point in this poly
		//
		if((a != tria)&&(a != trib))
			tric = a;
		else if((b != tria)&&(b != trib))
			tric = b;
		else	tric = c;

		//
		// Make sure tria,trib,tric (the triangle
		// that will come next on the stack) has
		// correct exiting edge
		//
		if(p == 1) { // exit edge must be 3rd poly edge (ie c,a)
			if((sharea==tric)&&(shareb==tria))
				need_swap = 0;
			else	need_swap = 1;
			}
		else	{    // exit edge must be 2nd poly edge (ie b,c)
			if((sharea==trib)&&(shareb==tric))
				need_swap = 0;
			else	need_swap = 1;
			}

		if(need_swap) {
			stack1 = strip_heap[stripi].vertices[vertcnt-2];
			strip_heap[stripi].vertices[vertcnt++] = stack1;
			p = (p==1) ? 0 : 1;
			}
		strip_heap[stripi].vertices[vertcnt++] = tric;
		p = (p==1) ? 0 : 1;
		}

	strip_heap[stripi].nvertices = vertcnt;

	goto gen_strips;
	}

delete [] face_heap;

if(nstrips > 0) {
   float meanstrip;
   meanstrip = ((float)sumstrip)/((float)nstrips);
   cout << endl;
   cout << "Triangle stripping summary: " << endl;
   cout << "                    Number of strips: " << nstrips << endl;
   cout << "    Number of single triangle strips: " << nsingleton << endl;
   cout << "                Maximum strip length: " << maxstrip << 
						" triangles " << endl;
   cout << "                   Mean strip length: " << meanstrip << 
						" triangles " << endl;
   cout << endl;
   }
else {
   cout << " *** No triangle strips created ***" << endl;
   }
VOUT("void GLSurface::generateStrips(Surface *surface)")
}









