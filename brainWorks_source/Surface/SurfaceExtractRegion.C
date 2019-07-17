///////////////////////////////////////////////////////////////////////////
//
// File: SurfaceExtractRegion.C
//
// Author: Keith Doolittle
//
// Purpose:
//
///////////////////////////////////////////////////////////////////////////
#include <stdlib.h>
#include <iostream.h>
#include <KDApplication.h>
#include <KDShell.h>
#include <KDWorkingChild.h>
#include <Xm/Text.h>
#include <AppIcons.h>
#include <IconList.h>
#include <TDL/IDLdefines.h>
#include <TDL/SurfaceUtils.h>
#include <LoadedVolume.h>
#include <Main.h>
#include <SurfaceBox.h>
#include <DynArray.h>
#include <SurfaceExtractRegion.h>


class PEdge {
public:
   PEdge()                    { a = b = n = 0;       }
   PEdge(int x)               { a = b = n = x;       }
   PEdge(int x, int y, int z) { a = x; b = y; n = z; }

   inline const PEdge & operator = (const PEdge &A) {
      a = A.a;
      b = A.b;
      n = A.n;
      return(*this);
   }

   inline bool operator == (const PEdge &A) {
     return((a==A.a)&&(b==A.b));
   }

   inline bool operator > (const PEdge &A) { return(false); }
   inline bool operator < (const PEdge &A) { return(false); }

   int a,b,n;
};

#ifdef GNU_COMPILER
template class DynArray<PEdge>;
#endif

#ifdef RS6K
#pragma define(DynArray<PEdge>)
#endif

template<>
u_int DynArray<PEdge>::m_increment=50;


DEFINECB(SurfaceExtractRegion,FillHolesCB,autoFillHoles)

static void CalculateCallbackCB(Widget wid, XtPointer cld, XtPointer cad)
{
VIN("static void CalculateCallbackCB(Widget wid, XtPointer cld, XtPointer cad)")
  if(cld) ((SurfaceExtractRegion*)cld)->Callback(1);
VOUT("static void CalculateCallbackCB(Widget wid, XtPointer cld, XtPointer cad)")
}

static void DoneCallbackCB(Widget wid, XtPointer cld, XtPointer cad)
{
VIN("static void DoneCallbackCB(Widget wid, XtPointer cld, XtPointer cad)")
  if(cld) ((SurfaceExtractRegion*)cld)->Callback(-1);
VOUT("static void DoneCallbackCB(Widget wid, XtPointer cld, XtPointer cad)")
}

static void ResetCallbackCB(Widget wid, XtPointer cld, XtPointer cad)
{
VIN("static void ResetCallbackCB(Widget wid, XtPointer cld, XtPointer cad)")
  if(cld) ((SurfaceExtractRegion*)cld)->reset();
VOUT("static void ResetCallbackCB(Widget wid, XtPointer cld, XtPointer cad)")
}

static void ImageChangeCB(Icon *ic, void *dat)
{
VIN("static void ImageChangeCB(Icon *ic, void *dat)")
  if(dat) ((SurfaceExtractRegion*)dat)->imagesChanged();
VOUT("static void ImageChangeCB(Icon *ic, void *dat)")
}

static void SegValueAChangeCB(char *str, void *dat)
{
VIN("static void SegValueAChangeCB(char *str, void *dat)")
  int i;
  if(str && sscanf(str,"%d",&i) == 1)
	if(dat) ((SurfaceExtractRegion*)dat)->setSegValue((unsigned char)i);
VOUT("static void SegValueAChangeCB(char *str, void *dat)")
}

static void SegValueBChangeCB(char *str, void *dat)
{
VIN("static void SegValueBChangeCB(char *str, void *dat)")
  int i;
  if(str && sscanf(str,"%d",&i) == 1)
	if(dat) ((SurfaceExtractRegion*)dat)->setNotSegValue((unsigned char)i);
VOUT("static void SegValueBChangeCB(char *str, void *dat)")
}

static bool MouseCB(u_long dat, int x, int y, int btn, int fcn)
{
VIN("static bool MouseCB(u_long dat, int x, int y, int btn, int fcn)")
  if(dat) return(((SurfaceExtractRegion*)dat)->mouse(x,y,btn,fcn));
  else    return(true);
VOUT("static bool MouseCB(u_long dat, int x, int y, int btn, int fcn)")
}


bool SurfaceExtractRegion::mouse(int x, int y, int btn, int fcn)
{
VIN("bool SurfaceExtractRegion::mouse(int x, int y, int btn, int fcn)")

  // Fill hole on <Shift>Btn1 Down

  if((fcn==0)&&(btn==1))
      fillHole(x,y);

VOUT("bool SurfaceExtractRegion::mouse(int x, int y, int btn, int fcn)")
return(true);
}


SurfaceExtractRegion::~SurfaceExtractRegion()
{
VIN("SurfaceExtractRegion::~SurfaceExtractRegion()")
  for(int i=0;i<nsegNames;i++)
	delete [] segNames[i];

  if(imageList)     delete imageList;
  if(segList)       delete segList;
  if(segValueListA) delete segValueListA;
  if(segValueListB) delete segValueListB;
  if(surfBox)       delete surfBox;

  if(lsurf)        lsurf->unreference();
VOUT("SurfaceExtractRegion::~SurfaceExtractRegion()")
}


SurfaceExtractRegion::SurfaceExtractRegion(LoadedSurface *lS) :
	KDChildWindow("Surface Extract Region",
			4*App->Width()/5,4*App->Height()/5),
    SurfaceExtractRegionBase(lS->surface())
{
VIN("SurfaceExtractRegion::SurfaceExtractRegion(LoadedSurface *lS) :")
  imageList      = NULL;
  segList        = NULL;
  segValueListA  = NULL;
  segValueListB  = NULL;
  surfBox        = NULL;
  nsegNames      = 0;

  if(!lS || !lS->surface()) {
	delete this;
	return;
  }

  lsurf = lS;
  
  if(lsurf->referenceCount() > 0) {
	App->ShowMessage("Cannot perform this operation. Surface is in use");
        delete this;
        return;
  }

  lsurf->reference();

  int l,nv = surf->getNumVert();

  newLSurf = new LoadedSurface(newSurf);
  (void)new Icon(newLSurf);

  // initially, no vertices are removed
  delVert.setDim(nv);
  delVert = 0U;

  // define new vertex mapping (from orig surface to new surface)
  // originally identity

  old2NewIndex.setDim(nv);
  new2OldIndex.setDim(nv);

  for(l=0;l<nv;l++) {
	old2NewIndex[l] = l;
	new2OldIndex[l] = l;
  }


  // ----------------  GUI

  Widget btn;
  u_long flgs = KDSHELL_TITLEBAR;

  int nitems = 15;
  int w = 2*(25 * App->BigFontWidth());
  int h = nitems * (App->BigFontHeight() + 12);

  int x = App->Width()/2  - w/2;
  int y = App->Height()/2 - h/2;

  int ddy = 100/nitems;
  int yst = 0;
  int ynd = ddy;

  ATTACHMENTS_TBLR(yst,ynd,70,100)
  imageList = new IconList(this->Form(),"Image",IXImageIcon,
                                NULL,0,arglist,c);
  imageList->AddCallback(ImageChangeCB,(void*)this);
  yst = ynd;
  ynd += ddy;

  ATTACHMENTS_TBLR(yst,ynd,70,100)
  segList = new IconList(this->Form(),"Segmentation",IXImageIcon,
                                NULL,0,arglist,c);
  segList->AddCallback(ImageChangeCB,(void*)this);
  yst = ynd;
  ynd += ddy;

  ATTACHMENTS_TBLR(yst,ynd,70,85)
  minxw = App->LabelText(this->Form(),"Min X","0",10,
			App->SmallFontList(),arglist,c);
  ATTACHMENTS_TBLR(yst,ynd,85,100)
  maxxw = App->LabelText(this->Form(),"Max X","0",10,
			App->SmallFontList(),arglist,c);
  yst = ynd;
  ynd += ddy;

  ATTACHMENTS_TBLR(yst,ynd,70,85)
  minyw = App->LabelText(this->Form(),"Min Y","0",10,
			App->SmallFontList(),arglist,c);
  ATTACHMENTS_TBLR(yst,ynd,85,100)
  maxyw = App->LabelText(this->Form(),"Max Y","0",10,
			App->SmallFontList(),arglist,c);
  yst = ynd;
  ynd += ddy;

  ATTACHMENTS_TBLR(yst,ynd,70,85)
  minzw = App->LabelText(this->Form(),"Min Z","0",10,
			App->SmallFontList(),arglist,c);
  ATTACHMENTS_TBLR(yst,ynd,85,100)
  maxzw = App->LabelText(this->Form(),"Max Z","0",10,
			App->SmallFontList(),arglist,c);
  yst = ynd;
  ynd += ddy;


  ATTACHMENTS_TBLR(yst,ynd,70,100)
  btn = App->Label(this->Form(),"KEEP Vertices Where",
			App->MedFontList(),arglist,c);
  XmChangeColor(btn,SHELLCOLOR);
  yst = ynd;
  ynd += ddy;


  ATTACHMENTS_TBLR(yst,ynd,70,100)  
  wNonZeroSeg = App->Toggle(this->Form(),"Segmentation Non-Zero ",
		App->SmallFontList(),bNonZeroSeg,NULL,NULL,arglist,c);
  yst = ynd;
  ynd += ddy;

  ATTACHMENTS_TBLR(yst,ynd,70,100)  
  wNonZeroImage = App->Toggle(this->Form(),"Image Non-Zero ",
		App->SmallFontList(),bNonZeroImage,NULL,NULL,arglist,c);
  yst = ynd;
  ynd += ddy;

  ATTACHMENTS_TBLR(yst,ynd,70,90)  
  wSegEquals = App->Toggle(this->Form(),"Seg Value ==",
		App->SmallFontList(),bSegEquals,NULL,NULL,arglist,c);
  aX1 = 90;
  aX2 = 100;
  aY1 = yst;
  aY2 = ynd;
  yst = ynd;
  ynd += ddy;

  ATTACHMENTS_TBLR(yst,ynd,70,90)  
  wSegNotEquals = App->Toggle(this->Form(),"Seg Value !=",
		App->SmallFontList(),bSegNotEquals,NULL,NULL,arglist,c);
  bX1 = 90;
  bX2 = 100;
  bY1 = yst;
  bY2 = ynd;
  yst = ynd;
  ynd += ddy;


  ATTACHMENTS_TBLR(yst,ynd,70,100)
  wNeighborNotZero = App->Toggle(this->Form(),"Seg Neighbors Non-Zero",
		App->SmallFontList(),bNeighborNotZero,NULL,NULL,arglist,c);
  yst = ynd;
  ynd += ddy;

  ATTACHMENTS_TBLR(yst,ynd,70,100)
  wNormalNotZero = App->Toggle(this->Form(),"Normals Face Non-Zero in Seg",
		App->SmallFontList(),bNormalNotZero,NULL,NULL,arglist,c);
  yst = ynd;
  ynd += ddy;

  ATTACHMENTS_TBLR(yst,ynd,70,90)  
  wNormalFaces = App->Toggle(this->Form(),"Normals Face Seg Value >=",
		App->SmallFontList(),bNormalFaces,NULL,NULL,arglist,c);
  ATTACHMENTS_TBLR(yst,ynd,90,98)  
  wNormalFacesText = App->Text(this->Form(),"50",10,NULL,NULL,
		App->SmallFontList(),arglist,c);
  yst = ynd;
  ynd += ddy;


  ATTACHMENTS_TBLR(yst,ynd,70,100)
  wNormalSize = App->LabelText(this->Form(),"Length of Normals:",
		"1.0",10,App->SmallFontList(),arglist,c);
  yst = ynd;
  ynd += ddy;

  ATTACHMENTS_TBLR(yst,ynd,70,100)
  wFillHoles = App->ButtonText(this->Form(),"Fill Holes w/verts <= ",
                   "1000",10,70,App->SmallFontList(),App->SmallFontList(),
		   FillHolesCB,(XtPointer)this,arglist,c);
  yst = ynd;
  ynd += ddy;

  ATTACHMENTS_TBLR(yst,ynd,70,80)
  btn = App->Button(this->Form(),"Calculate",App->SmallFontList(),
				CalculateCallbackCB, (XtPointer)this,arglist,c); 
  ATTACHMENTS_TBLR(yst,ynd,80,90)
  btn = App->Button(this->Form(),"Done",App->SmallFontList(),DoneCallbackCB,
				(XtPointer)this,arglist,c); 
  XmChangeColor(btn,CANCELCOLOR);

  ATTACHMENTS_TBLR(yst,ynd,90,100)
  btn = App->Button(this->Form(),"RESET",App->SmallFontList(),ResetCallbackCB,
				(XtPointer)this,arglist,c); 


  ATTACHMENTS_TBLR(0,100,0,70)
  surfBox = new SurfaceBox(this->Form(),newLSurf,arglist,c);
  surfBox->SetPickStyle(PickNone);
  surfBox->SetMouseFunction(MouseCB,(u_long)this);

  this->Show();
  imagesChanged();

retry_extract:

  gendone = 0;
  while(!gendone) App->ProcessEvents();

  if(gendone > 0) {
     if(!mri || !seg) {
	App->ShowMessage("Please select valid Image and Segmentation volumes");
	goto retry_extract;
     }

     minx = App->GetInt(minxw);
     maxx = App->GetInt(maxxw);
     miny = App->GetInt(minyw);
     maxy = App->GetInt(maxyw);
     minz = App->GetInt(minzw);
     maxz = App->GetInt(maxzw);

     if((minx<0)||(minx>=maxx)||(maxx>=nx)) {
	App->ShowMessage("Invalid X coordinates for image region");
	goto retry_extract;
     }

     if((miny<0)||(miny>=maxy)||(maxy>=ny)) {
	App->ShowMessage("Invalid Y coordinates for image region");
	goto retry_extract;
     }

     if((minz<0)||(minz>=maxz)||(maxz>=nz)) {
	App->ShowMessage("Invalid Z coordinates for image region");
	goto retry_extract;
     }

     bNonZeroSeg      = App->GetToggleValue(wNonZeroSeg);
     bNonZeroImage    = App->GetToggleValue(wNonZeroImage);
     bSegEquals       = App->GetToggleValue(wSegEquals);
     bSegNotEquals    = App->GetToggleValue(wSegNotEquals);
     bNeighborNotZero = App->GetToggleValue(wNeighborNotZero);
     bNormalNotZero   = App->GetToggleValue(wNormalNotZero);
     bNormalFaces     = App->GetToggleValue(wNormalFaces);

     normalValueFaces = App->GetFloat(wNormalFacesText);

     seg->reference();
     mri->reference();
     calculate();
     seg->unreference();
     mri->unreference();

     goto retry_extract;
  }

// will get here only when 'Done' is pressed

this->Hide();
delete this;
VOUT("SurfaceExtractRegion::SurfaceExtractRegion(LoadedSurface *lS) :")
}

void SurfaceExtractRegion::reset()
{
VIN("void SurfaceExtractRegion::reset()")

  SurfaceExtractRegionBase::reset();
  
VOUT("void SurfaceExtractRegion::reset()")
}

void SurfaceExtractRegion::updateSurface()
{
VIN("void SurfaceExtractRegion::updateSurface()")

   GenericWorkProc("Updating Surface...");

   SurfaceExtractRegionBase::updateSurface();

   GenericWorkProc(NULL);

   newLSurf->surfaceChanged();
   surfBox->surfaceChanged();
VOUT("void SurfaceExtractRegion::updateSurface()")
}

// for qsort
static int compString(const void *s1, const void *s2)
{
VIN("static int compString(const void *s1, const void *s2)")
  if(!s1 || !s2) return(0);
  
  int i1,i2;
  sscanf((const char*)s1,"%d",&i1);
  sscanf((const char*)s2,"%d",&i2);
//cerr << "Comp( " << (const char*)s1 << ", " << (const char*)s2 << ")" << endl;
  if(i1 < i2) return(-1);
  else if(i1 == i2) return(0);
  else return(1);
VOUT("static int compString(const void *s1, const void *s2)")
}


void SurfaceExtractRegion::imagesChanged()
{
VIN("void SurfaceExtractRegion::imagesChanged()")
  if(!imageList || !segList)
	return;

  Icon *iIcon,*sIcon;
  iIcon = imageList->getSelectedIcon();
  sIcon = segList->getSelectedIcon();

  if(!iIcon || !sIcon || !(iIcon->Type() == IXImageIcon)||
     !(sIcon->Type() == IXImageIcon))
	return;

  LoadedVolume *oldseg = seg;

  mri = (LoadedVolume*)iIcon->Data();
  seg = (LoadedVolume*)sIcon->Data();

  int nx1 = mri->volume()->getSizeX();
  int ny1 = mri->volume()->getSizeY();
  int nz1 = mri->volume()->getSizeZ();
  int nx2 = seg->volume()->getSizeX();
  int ny2 = seg->volume()->getSizeY();
  int nz2 = seg->volume()->getSizeZ();

  if((nx1!=nx2)||(ny1!=ny2)||(nz1!=nz2)) {
	App->ShowMessage("Please select image/seg that have same dimensions");
	mri = seg = NULL;
	return;
  }

  nx   = nx1;
  ny   = ny1;
  nz   = nz1;

  char tstr[50];

  sprintf(tstr,"0");
  App->SetTextString(minxw,tstr);
  App->SetTextString(minyw,tstr);
  App->SetTextString(minzw,tstr);

  sprintf(tstr,"%d",nx-1);
  App->SetTextString(maxxw,tstr);

  sprintf(tstr,"%d",ny-1);
  App->SetTextString(maxyw,tstr);

  sprintf(tstr,"%d",nz-1);
  App->SetTextString(maxzw,tstr);

  // now define segmentation value list
  // iff seg has changed

  if(oldseg != seg) {

     int  i,j,k,isused[256];
     unsigned char b;
     char *segname;

     seg->volume()->convertToUnsignedChar();

     for(i=0;i<256;i++) isused[i] = 0;

     for(i=0;i<nsegNames;i++)
	delete [] segNames[i];
     nsegNames = 0;

     for(k=0;k<nz;k++)
       for(j=0;j<ny;j++)
         for(i=0;i<nx;i++) {
	   b = seg->volume()->u_char_data()[k][j][i];
           if(!isused[b] && (nsegNames < 255)) {
	      isused[b] = 1;
	      segname = new char[50];	
	      sprintf(segname,"%d",(int)b);
	      segNames[nsegNames++] = segname;
	      if(nsegNames == 1) {
		segValueEquals    = b; // once sort works, use first
		segValueNotEquals = b;
	      }
           }
         }

/*

QSORT not working properly

int ii;
cerr << "SegNames: Before :";
for(ii=0;ii<nsegNames;ii++)
cerr << segNames[ii] << " ";
cerr << endl;

     // sort the seg values in ascending order
     qsort((void*)segNames,nsegNames,sizeof(char*),compString);

cerr << "SegNames: After :";
for(ii=0;ii<nsegNames;ii++)
cerr << segNames[ii] << " ";
cerr << endl;
*/

     if(segValueListA) delete segValueListA;
     if(segValueListB) delete segValueListB;

     ATTACHMENTS_TBLR(aY1,aY2,aX1,aX2)
     segValueListA = new IconList(this->Form()," ",IXNoIcon,
                                segNames,nsegNames,arglist,c);
     segValueListA->AddExtraCallback(SegValueAChangeCB,(void*)this);

     ATTACHMENTS_TBLR(bY1,bY2,bX1,bX2)
     segValueListB = new IconList(this->Form()," ",IXNoIcon,
                                segNames,nsegNames,arglist,c);
     segValueListB->AddExtraCallback(SegValueBChangeCB,(void*)this);
  }
VOUT("void SurfaceExtractRegion::imagesChanged()")
}

bool SurfaceExtractRegion::fillPatch(int vert)
{
  int i,idx,n,nei;
  DynArray<int> unDel,didDel;

  int nmax = App->GetInt(wFillHoles);

  if(!surf->hasNbhd())
      surf->genNeighborhoods();

  didDel.clear();
  unDel.clear();

  unDel.add(vert);

  int ccnt=0;
  while(unDel.canPop()) {
    idx = unDel.pop();
    delVert[idx] = 0U;
    didDel.add(idx);
    for(i=0;i<surf->neighborhood()[idx].numNeighbors();i++) {
        nei = surf->neighborhood()[idx].getNeighbor(i);
        if(delVert[nei]) unDel.addUnique(nei);
    }
    if(++ccnt > nmax) {
       cerr << "** Patch Too Large.  Stopping fill" << endl;
       break;
    }
  }

  if(ccnt > nmax) {
     while(didDel.canPop())
       delVert[didDel.pop()] = 1U;
     return(false);
  } else
     return(true);
}

void SurfaceExtractRegion::fillPath(DynArray<int> &P,
                                    Array1D<u_char> &isEdge)
{
  int i,nf,a,b,c,targ;
  bool found;

  // find a deleted vertex that adjoins path

  found = false;
  nf    = surf->getNumPoly();

  for(i=0;i<nf && !found;i++) {
      a = surf->facets()[i][0];
      b = surf->facets()[i][1];
      c = surf->facets()[i][2];
      if(P.exists(old2NewIndex[a])) {
         if(delVert[b])      { targ = b; found = true; }
         else if(delVert[c]) { targ = c; found = true; }
      }
      if(P.exists(old2NewIndex[b])) {
         if(delVert[a])      { targ = a; found = true; }
         else if(delVert[c]) { targ = c; found = true; }
      }
      if(P.exists(old2NewIndex[c])) {
         if(delVert[b])      { targ = b; found = true; }
         else if(delVert[a]) { targ = a; found = true; }
      }
  }  

  if(found) {
     if(fillPatch(targ))
       for(i=0;i<P.size();i++)
         isEdge[P[i]] = 0U;
  }
}

void addPEdge(int a, int b, DynArray<PEdge> &D)
{
   int i,i1,i2;

   if(a < b) { i1 = a; i2 = b; }
   else      { i1 = b; i2 = a; }

   for(i=0;i<D.size();i++)
       if((D[i].a == i1)&&(D[i].b == i2)) {
           D[i].n++;
           return;
       }

   D.add(PEdge(i1,i2,1));
}

bool edgeExists(int a, int b, DynArray<PEdge> &D)
{
   int i,i1,i2;

   if(a < b) { i1 = a; i2 = b; }
   else      { i1 = b; i2 = a; }

   for(i=0;i<D.size();i++)
       if((D[i].a == i1)&&(D[i].b == i2))
           return(true);

   return(false);
}

int SurfaceExtractRegion::findPath(int idx,
                           Array1D<u_char> &isEdge,
                           Array1D<u_char> &vUsed,
                           DynArray<PEdge> &E,
                           DynArray<int> &P)
{
  int i,j,t,nei,nnei,np,targ;
  int startindx;

  startindx = idx;

  np = 0;
  P.clear();
  P.add(idx); np++;
  vUsed[idx] = 1U;

  bool done = false;
  targ = idx;

  while(!done) {
      nnei = newSurf->neighborhood()[targ].numNeighbors();
      t = -1;
      for(j=0;j<nnei;j++) {
          nei = newSurf->neighborhood()[targ].getNeighbor(j);
          if(isEdge[nei] && !vUsed[nei] && edgeExists(nei,targ,E)) {
             t = nei; 
             break;
          }
      }
      if(t >= 0) {
        vUsed[t] = 1;
        P.add(t); np++;
        targ = t;
      } else 
        done = true;
  }

  // make sure its closed by looking for start index

  nnei = newSurf->neighborhood()[targ].numNeighbors();
  t = -1;
  for(j=0;j<nnei;j++) {
      nei = newSurf->neighborhood()[targ].getNeighbor(j);
      if(isEdge[nei] && (nei == startindx)) {
          t = nei; 
          break;
      }
  }

  if(t >= 0) {
     cerr << ".";
     return(np);
  } else {
// this adds alot of time and only fills in a few more holes
//     for(i=0;i<np;i++)
//        vUsed[P[i]] = 0U;
     cerr << "X";
     return(-np);
  }
}

void SurfaceExtractRegion::autoFillHoles()
{
  GenericWorkProc("Filling Holes....");

  int i,a,b,c,nv,nf,cnt;

  nv = newSurf->getNumVert();
  nf = newSurf->getNumPoly();

  if(!newSurf->hasNbhd())
      newSurf->genNeighborhoods();

  Array1D<u_char> isEdge(nv);

  //------------------------------------------------------------

  // Tag edge verts

  isEdge = 0U;
  for(i=0;i<nf;i++) {
     a = newSurf->facets()[i][0];
     b = newSurf->facets()[i][1];
     c = newSurf->facets()[i][2];
     isEdge[a]++;
     isEdge[b]++;
     isEdge[c]++;
  }


  // remove polys attached by single vertex

  int na,nb,nc;
  for(i=0;i<nf;i++) {
     a = newSurf->facets()[i][0];
     b = newSurf->facets()[i][1];
     c = newSurf->facets()[i][2];
     na = isEdge[a];
     nb = isEdge[b];
     nc = isEdge[c];
     if((na==2)&&(nb==2)) delVert[new2OldIndex[c]] = 1U;
     if((nb==2)&&(nc==2)) delVert[new2OldIndex[a]] = 1U;
     if((na==2)&&(nc==2)) delVert[new2OldIndex[b]] = 1U;
  }

  updateSurface();

  //------------------------------------------------------------

  nv = newSurf->getNumVert();
  nf = newSurf->getNumPoly();

  if(!newSurf->hasNbhd())
      newSurf->genNeighborhoods();

  isEdge.setDim(nv);

  // Tag edge verts

  isEdge = 0U;
  for(i=0;i<nf;i++) {
     a = newSurf->facets()[i][0];
     b = newSurf->facets()[i][1];
     c = newSurf->facets()[i][2];
     isEdge[a]++;
     isEdge[b]++;
     isEdge[c]++;
  }

  int nedge = 0;
  for(i=0;i<nv;i++) 
     if(newSurf->neighborhood()[i].numNeighbors() != isEdge[i]) {
        isEdge[i] = 1;
        nedge++;
     } else
        isEdge[i] = 0;

  cerr << nedge << " vertex edge points found in new surface" << endl;



  DynArray<PEdge> fEdges;
  for(i=0;i<nf;i++) {
     a = newSurf->facets()[i][0];
     b = newSurf->facets()[i][1];
     c = newSurf->facets()[i][2];
     if(isEdge[a] && isEdge[b])
        addPEdge(a,b,fEdges);
     if(isEdge[b] && isEdge[c])
        addPEdge(b,c,fEdges);
     if(isEdge[a] && isEdge[c])
        addPEdge(a,c,fEdges);
  }

  // now invalidate any poly edges that have more than one reference
  // (which means it doesnt lie along an open boundary)

  DynArray<PEdge> P;
  for(i=0;i<fEdges.size();i++)
      if(fEdges[i].n == 1)
         P.add(fEdges[i]);

  // free old memory

  fEdges.clear();

  cerr << P.size() << " hole edges found" << endl;

  Array1D<u_char> vUsed(nv);
  DynArray<int> path,maxpath;
  int npath,nmaxpath;
  int npf,j,idx;

  vUsed    = 0U;
  npf      = 0;
  nmaxpath = 0;
  for(i=0;i<nv;i++) {
     if(!isEdge[i] || vUsed[i]) 
        continue; 
     npath = findPath(i,isEdge,vUsed,P,path);
     if(npath <= 0) continue;
     GenericWorkProc("Filling Holes....");
     npf++;
     if(npath > nmaxpath) {
        if(nmaxpath > 0) 
           fillPath(maxpath,isEdge);
        maxpath  = path;
        nmaxpath = npath;
     } else 
       fillPath(path,isEdge);
  }

  GenericWorkProc(NULL);

  updateSurface();
}


//
// look for previously deleted (delVert[i] = 1U) vertices
// and un-delete them if they are around the mouse click
//

void SurfaceExtractRegion::fillHole(int mx, int my)
{
VIN("void SurfaceExtractRegion::fillHole(int mx, int my)")
  int i,nv,ix,iy,mini;
  float x,y,z,xx,yy,zz,zz2;
  float minz;

  nv = surf->getNumVert();

  if(!surf->hasNorms())
    surf->genNormals();

  mini = -1;
  minz = 1E20;
  for(i=0;i<nv;i++) {
      if(!delVert[i]) continue;
      x = surf->vertices()[i].x();
      y = surf->vertices()[i].y();
      z = surf->vertices()[i].z();
      surfBox->project(x,y,z,xx,yy,zz);
      ix = (int)(xx+0.5);
      iy = (int)(yy+0.5);
      if((abs(ix-mx)<10)&&(abs(iy-my)<10)&&(zz < minz)) 
          { minz = zz; mini = i; }
  }
  if(mini >= 0) {
    if(fillPatch(mini))
       updateSurface();
  }

VOUT("void SurfaceExtractRegion::fillHole(int mx, int my)")
}


void SurfaceExtractRegion::calculate()
{
VIN("void SurfaceExtractRegion::calculate()")
   if(!surf || !mri || !seg) 
	return;

   GenericWorkProc("Calculating...");
   
   normalSize = App->GetFloat(wNormalSize);
   
   SurfaceExtractRegionBase::calculate(seg->volume(), mri->volume());

VOUT("void SurfaceExtractRegion::calculate()")
}



