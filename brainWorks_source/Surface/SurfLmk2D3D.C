///////////////////////////////////////////////////////////////////////////
//
// File: SurfLmk2D3D.C
//
// Author: Keith Doolittle
//
// Purpose:
//
///////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <string.h>
#include <sys/param.h>
#include <KDApplication.h>
#include <TDL/Surface.h>
#include <TDL/ByuSurface.h>
#include <TDL/SurfaceUtils.h>
#include <Xm/ScrolledW.h>
#include <Xm/RowColumn.h>
#include <Xm/Text.h>
#include <Xm/Form.h>
#include <Main.h>
#include <MainUtils.h>
#include <SurfaceBox.h>
#include <IconList.h>
#include <LoadedSurface.h>
#include <SurfaceViewer.h>
#include <GLPrimitive.h>
#include <GLSurface.h>
#include <SurfLmk2D3D.h>

#define LINE_WIDTH	5.0

#define INDICES_TITLE "Indices 1.0"

//
// default line color
//
#define DEFRED		0.0
#define DEFGREEN	0.0
#define DEFBLUE		1.0
//
// selected line color
//
#define SELRED		1.0
#define SELGREEN	0.0
#define SELBLUE		0.0

DEFINECB(SurfLmk2D3D,SetColorCB,setColor)
DEFINECBP(SurfLmk2D3D,SelectButtonCB,selectLandmark,wid)
DEFINECB(SurfLmk2D3D,DeleteCurrentCB,deleteCurrent)
DEFINECB(SurfLmk2D3D,LoadLandmarksCB,loadLandmarks)
DEFINECB(SurfLmk2D3D,LoadLandmarks2CB,loadLandmarks2)
DEFINECB(SurfLmk2D3D,SaveLandmarksCB,saveLandmarks)
DEFINECB(SurfLmk2D3D,SaveLandmarks2CB,saveLandmarks2)
DEFINECB(SurfLmk2D3D,LoadIndicesCB,loadIndices)
DEFINECB(SurfLmk2D3D,SaveIndicesCB,saveIndices)
DEFINECB(SurfLmk2D3D,ClearLandmarksCB,clearLandmarks)
DEFINECB(SurfLmk2D3D,ToggleChangeCB,toggleChange)
DEFINECB(SurfLmk2D3D,SetLmkSizeCB,setLmkSize)

static void PickCB(u_long dat, void *s, int indx)
{
  if(dat) ((SurfLmk2D3D*)dat)->pickPoint(indx);
}

void SurfLmk2D3D::toggleChange()
{
  m_bChange = (App->GetToggleValue(m_wToggle) == TRUE);
}

void SurfLmk2D3D::setLmkSize()
{
  float f;
  GLPrimitive *p;
  f = sphere_rad;
  if(GetFloatValue("Enter Lmk Size",0.00001,100,f)) {
     sphere_rad = f;
     for(p=lines;p;p=p->Next())
         p->setRadius(sphere_rad);
     for(p=lines2;p;p=p->Next())
         p->setRadius(sphere_rad);
     surf1->redraw();
     surf2->redraw();
  } 
}

void SurfLmk2D3D::loadIndices()
{
 int i,n,j;
 unsigned int idx;
 ifstream fp;
 char fname[1024];
 char lname[1024];
 char line[1024];

  if(!surface)
     return;

 retry_loadi:
  fname[0] = 0;
  if(!Browser->GetFilename("Load Landmark Indices From",fname,1024,"*.lmi"))
      return;

  fp.open(fname);
  if(!fp || fp.fail()) {
     sprintf(lname,"ERROR: cannot open '%s' for reading",fname);
     App->ShowMessage(lname);
     goto retry_loadi;
   }

   fp.getline(line,1023);
   if(fp.fail() || strncmp(line,INDICES_TITLE,strlen(INDICES_TITLE))) {
       sprintf(lname,"ERROR: '%s' is not a valid indices file",fname);
       App->ShowMessage(lname);
       return;
   }

   landmarks.clear();
   flandmarks.clear();
   indices.clear();
   while(!fp.eof()) {
       fp.getline(line,1023);
       if(sscanf(line,"%ld %n",&idx,&j) < 1)
          continue;
       if((idx >= surface->getNumVert())||(idx >= fsurface->getNumVert())) {
          fp.close();
          clearLandmarks();
          App->ShowMessage("ERROR: the landmark indices dont fit these surfaces");
          return;
       }
       sprintf(lname,"%s",&(line[j]));
       indices.add(idx);

       float x = surface->vertices()[idx].x();
       float y = surface->vertices()[idx].y();
       float z = surface->vertices()[idx].z();

       float x2 = fsurface->vertices()[idx].x();
       float y2 = fsurface->vertices()[idx].y();
       float z2 = fsurface->vertices()[idx].z();

       landmarks.add(x,y,z,lname,1);
       flandmarks.add(x2,y2,z2,lname,1);
       addLandmark(x,y,z,x2,y2,z2,lname);
   }
   fp.close();
}

void SurfLmk2D3D::saveIndices()
{
  ofstream fp;
  int i,n;
  const Point *P;
  char fname[1024];
  char lname[1024];
  Widget txt;

  n = indices.size();
  if(n <= 0) {
     App->ShowMessage("ERROR: No indices to save");
     return;
  }

 retry_savei:
  fname[0] = 0;
  if(!Browser->GetFilename("Save Landmark Indices As",fname,1024,"*.lmi"))
     return;

  fp.open(fname);
  if(!fp || fp.fail()) {
     sprintf(lname,"ERROR: cannot open '%s' for writing",fname);
     App->ShowMessage(lname);
     goto retry_savei;
  }
  fp << INDICES_TITLE << endl;

  for(i=0;i<n;i++) {
     txt = XtNameToWidget(XtParent(line_buttons[i]),"TextWidget");
     if(txt) App->GetString(txt,lname,1024);
     else sprintf(lname,"<no-name>");
     if(strlen(lname) < 1) sprintf(lname,"<no-name>");
     fp << indices[i] << " " << lname << endl;
  }
  fp.close();
}

void SurfLmk2D3D::loadLandmarks()
{
  clearLandmarks();

  if(!surface) 
      return;

  int i,j;
  const Point *P;
  char fname[1024];
  char loadStatus[1024];
  char lname[1024];

 retry_load:
  fname[0] = 0;
  if(Browser->GetFilename("Load 3D Landmarks From", fname, 1024))
  {
    if (!loadLandmarksFromFile(fname, 1, loadStatus))
    {
        App->ShowMessage(loadStatus);
        goto retry_load;
    }

     for (i = 0; i < landmarks.num(); i++)
     {
         if(!(P = landmarks.point(i)))
             continue;

         indices[i] = findLandmark(surface, P->x(), P->y(), P->z());
         if(indices[i] == (unsigned int)-1)
            continue;

        if (landmarks.name(i))
            sprintf(lname, "%s", landmarks.name(i));
        else
            sprintf(lname, "<no-name>");
         
         float x2 = fsurface->vertices()[indices[i]].x();
         float y2 = fsurface->vertices()[indices[i]].y();
         float z2 = fsurface->vertices()[indices[i]].z();

         addLandmark(
            P->x(), P->y(), P->z(),
            x2, y2, z2,
            lname);
     }
     surf1->redraw();
     surf2->redraw();
  }
}

void SurfLmk2D3D::loadLandmarks2()
{
  clearLandmarks();

  if(!fsurface) 
      return;

  int i,j;
  const Point *P;
  char fname[1024];
  char loadStatus[1024];
  char lname[1024];

 retry_load2:
  fname[0] = 0;
  if(Browser->GetFilename("Load 2D Landmarks From", fname, 1024))
  {
    if (!loadLandmarksFromFile(fname, 2, loadStatus))
    {
        App->ShowMessage(loadStatus);
        goto retry_load2;
    }

     for (i = 0; i < flandmarks.num(); i++)
     {
        if (!(P = flandmarks.point(i)))
            continue;
         
         indices[i] = findLandmark(fsurface, P->x(), P->y(), P->z());
         if(indices[i] == (unsigned int) -1)
            continue;

        if (flandmarks.name(i))
            sprintf(lname, "%s", flandmarks.name(i));
        else
            sprintf(lname, "<no-name>");

         float x2 = surface->vertices()[indices[i]].x();
         float y2 = surface->vertices()[indices[i]].y();
         float z2 = surface->vertices()[indices[i]].z();

         addLandmark(
            x2, y2, z2,
            P->x(), P->y(), P->z(),
            lname);
     }
     surf1->redraw();
     surf2->redraw();
  }
}



void SurfLmk2D3D::saveLandmarks()
{
VIN("void SurfLmk2D3D::saveLandmarks()")

  int i,n;
  const Point *P;
  char fname[1024];
  char lname[1024];
  Widget txt;

  n = landmarks.num();
  if(n <= 0) {
     App->ShowMessage("ERROR: No landmarks to save");
     return;
  }

 retry_save:
  fname[0] = 0;
  if(Browser->GetFilename("Save Landmarks As",fname,1024,"*.lmk"))  {

     // set the landmark names
     for(i=0;i<n;i++) {
        if(!(P=landmarks.point(i)))
           continue;

        txt = XtNameToWidget(XtParent(line_buttons[i]),"TextWidget");
        if(txt) App->GetString(txt,lname,1024);
        else sprintf(lname,"<no-name>");
        if(strlen(lname) < 1) sprintf(lname,"<no-name>");
        landmarks.set(i,P->x(),P->y(),P->z(),lname,1);
     }
     if (!saveLandmarksToFile(fname, 1, lname))
     {
        App->ShowMessage(lname);
        goto retry_save;
     }
  }
VOUT("void SurfLmk2D3D::saveLandmarks()")
}

void SurfLmk2D3D::saveLandmarks2()
{
  int i,n;
  const Point *P;
  char fname[1024];
  char lname[1024];
  Widget txt;

  n = flandmarks.num();
  if(n <= 0) {
     App->ShowMessage("ERROR: No landmarks to save");
     return;
  }

 retry_save2:
  fname[0] = 0;
  if(Browser->GetFilename("Save 2D Landmarks As",fname,1024,"*.lmk"))  {

     // set the landmark names
     for(i=0;i<n;i++) {
        if(!(P=flandmarks.point(i)))
           continue;

        txt = XtNameToWidget(XtParent(line_buttons[i]),"TextWidget");
        if(txt) App->GetString(txt,lname,1024);
        else sprintf(lname,"<no-name>");
        if(strlen(lname) < 1) sprintf(lname,"<no-name>");
        flandmarks.set(i,P->x(),P->y(),P->z(),lname,1);
     }

     if (!saveLandmarksToFile(fname, 2, lname))
     {
        App->ShowMessage(lname);
        goto retry_save2;
     }
  }
}

void SurfLmk2D3D::clearLandmarks()
{
VIN("void SurfLmk2D3D::clearLandmarks()")
  //
  // clear all GLPrimitive lines
  //
  while(lines)  lines  = GLPrimitive::removeFromList(lines,lines);
  while(lines2) lines2 = GLPrimitive::removeFromList(lines2,lines2);
  surf1->setPrimitives(NULL);
  surf2->setPrimitives(NULL);

  int i;
  for(i=0;i<nbuttons;i++) {
    XtUnmapWidget(XtParent(line_buttons[i]));
    XtDestroyWidget(XtParent(line_buttons[i]));
  }

  if(line_buttons) delete [] line_buttons; 
  line_buttons = NULL;

  nbuttons   = 0;
  current_button  = NULL;
  current_line    = NULL;
  current_line2   = NULL;

  SurfLmk2D3DBase::clearLandmarks();
  
  surf1->redraw();
  surf2->redraw();
VOUT("void SurfLmk2D3D::clearLandmarks()")
}

void SurfLmk2D3D::setColor()
{
  float r,g,b;
  if(!current_line) {
      App->ShowMessage("Please select a current landmark first");
      return;
  }

  current_line->getColor(r,g,b);
  if(GetColor(&r,&g,&b)) {
     current_line->setColor(r,g,b);
     current_line2->setColor(r,g,b);
     XmChangeColor(current_button,App->RGBToPixel(r,g,b));
     surf1->redraw();
     surf2->redraw();
     saver = r;
     saveg = g;
     saveb = b;
  }
}


GLPrimitive *SurfLmk2D3D::findLine(GLPrimitive *ls, Widget w)
{
  //
  // find current GLPrimitive sphere corresponding
  // to the current button Widget
  //
  GLPrimitive *s = NULL;
  GLPrimitive *p;
  for(p=ls;p;p=p->Next())
    if(p->getHook() == (void*)w) {
      s = p;
      break;
    }
  return(s);
}


void SurfLmk2D3D::selectLandmark(Widget w)
{
VIN("void SurfLmk2D3D::selectLandmark(Widget w)")
  if(!w) return;

  if(current_line)
     current_line->setColor(saver,saveg,saveb);

  if(current_line2)
     current_line2->setColor(saver,saveg,saveb);

  if(current_button)
     XmChangeColor(current_button,App->RGBToPixel(saver,saveg,saveb));

  current_button = w;

  current_line   = findLine(lines,w);
  current_line2  = findLine(lines2,w);

  if(current_line)
     current_line->getColor(saver,saveg,saveb);

  if(current_line)
     current_line->setColor(1.f,0.f,0.f);

  if(current_line2)
     current_line2->setColor(1.f,0.f,0.f);

  XmChangeColor(current_button,App->RGBToPixel(1.f,0.f,0.f));
  surf1->redraw();
  surf2->redraw();

VOUT("void SurfLmk2D3D::selectLandmark(Widget w)")
}


Widget SurfLmk2D3D::addButton(char *name)
{
VIN("Widget SurfLmk2D3D::addButton(char *name)")
  Dimension ww = App->BigFontWidth() * 10;
  Dimension bh = App->BigFontHeight() + 8;

  ATTACHMENTS_TLRH(0,0,100,bh)

  Widget bform = XmCreateForm(line_form,"line_form_button",arglist,c);
  Set(XmNresizeHeight,False);
  Set(XmNresizeWidth,False);
  XtManageChild(bform);
  XmChangeColor(bform,SHELLCOLOR);

  ATTACHMENTS_TBLR(0,100,0,15)
  Widget btn = App->Button(bform," ",App->SmallFontList(),
                           SelectButtonCB,(void*)this,arglist,c);

  ATTACHMENTS_TBLR(0,100,15,100)
  Widget txt = App->Text(bform,name,128,
			 NULL,NULL,App->SmallFontList(),arglist,c);
VOUT("Widget SurfLmk2D3D::addButton(char *name)")
  return(btn);
}


void SurfLmk2D3D::addLandmark(float x, float y, float z,
                              float x2, float y2, float z2,
                              char *name)
{
VIN("void SurfLmk2D3D::addLandmark()")

  if(!name) name = "<no-name>";

  //
  // create graphical line
  // and new button 
  //
  Widget  *newbuttons = new Widget[nbuttons+1];
  int i;
  for(i=0;i<nbuttons;i++) 
    newbuttons[i] = line_buttons[i];

  nbuttons++;
  if(line_buttons) delete [] line_buttons;
  line_buttons = newbuttons;
  line_buttons[i] = this->addButton(name);
  XmChangeColor(line_buttons[i],App->RGBToPixel(DEFRED,DEFGREEN,DEFBLUE));

  GLPrimitive *newline = new GLPrimitive(x,y,z,sphere_rad);
  newline->setHook((void*)line_buttons[i]);
  newline->setColor(DEFRED,DEFGREEN,DEFBLUE);
  newline->setLineWidth(LINE_WIDTH);

  GLPrimitive *newline2 = new GLPrimitive(x2,y2,z2,sphere_rad);
  newline2->setHook((void*)line_buttons[i]);
  newline2->setColor(DEFRED,DEFGREEN,DEFBLUE);
  newline2->setLineWidth(LINE_WIDTH);
  //
  // add new GL line to the list 
  //
  lines = GLPrimitive::addToList(lines,newline);
  surf1->setPrimitives(lines);

  lines2 = GLPrimitive::addToList(lines2,newline2);
  surf2->setPrimitives(lines2);

  selectLandmark(line_buttons[i]);
VOUT("void SurfLmk2D3D::addLandmark()")
}


void SurfLmk2D3D::deleteCurrent()
{
VIN("void SurfLmk2D3D::deleteCurrent()")
  if(!current_button && !current_line)
    return;

  if(current_button) {
    XtUnmapWidget(XtParent(current_button));
    XtDestroyWidget(XtParent(current_button));
  }

  lines = GLPrimitive::removeFromList(lines,current_line);
  surf1->setPrimitives(lines);

  lines2 = GLPrimitive::removeFromList(lines2,current_line2);
  surf2->setPrimitives(lines2);

  int i;
  for(i=0;i<nbuttons-1;i++)
    if(line_buttons[i] == current_button)
      break;

  if(i < nbuttons) {
     landmarks.remove(i);
     flandmarks.remove(i);
     indices.removeAtIndex(i);
  }

  for(;i<nbuttons-1;i++) {
    line_buttons[i] = line_buttons[i+1];
  }
 
  nbuttons--;

  current_button = NULL;
  current_line   = NULL;
  current_line2  = NULL;
  surf1->redraw();
  surf2->redraw();
VOUT("void SurfLmk2D3D::deleteCurrent()")
}



void SurfLmk2D3D::pickPoint(int indx)
{
VIN("void SurfLmk2D3D::pickPoint(int indx)")
 int i;
 float x,y,z,x2,y2,z2,r;
 char tstr[200];

 x  = surface->vertices()[indx].x();
 y  = surface->vertices()[indx].y();
 z  = surface->vertices()[indx].z();
 x2 = fsurface->vertices()[indx].x();
 y2 = fsurface->vertices()[indx].y();
 z2 = fsurface->vertices()[indx].z();

 sprintf(tstr,"%d : (%.2f, %.2f, %.2f) (%.2f, %.2f, %.2f)",indx,x,y,z,x2,y2,z2);
 App->SetLabelText(text_label,tstr,App->SmallFontTag());

 if(m_bChange) {
  for(i=0;i<nbuttons;i++)
    if(line_buttons[i] == current_button)
      break;

  if(i < nbuttons) {
     landmarks.set(i,x,y,z,landmarks.name(i),true);
     flandmarks.set(i,x2,y2,z2,landmarks.name(i),true);
     indices[i] = (unsigned int)indx;
     if(current_line)
        current_line->setCenter(x,y,z);
     if(current_line2)
        current_line2->setCenter(x2,y2,z2);
  }

 } else {
   addLandmark(x,y,z,x2,y2,z2);
   landmarks.add(x,y,z,"<no-name>",1);
   flandmarks.add(x2,y2,z2,"<no-name>",1);
   indices.add((unsigned int)indx);
 }

 surf1->redraw();
 surf2->redraw();
VOUT("void SurfLmk2D3D::pickPoint(int indx)")
}

//
// Same as SurfaceViewer(), only 1 view (surf1)
//
SurfLmk2D3D::SurfLmk2D3D(LoadedSurface *surf, LoadedSurface *fsurf) :
          KDChildWindow(App->FilenameTail(surf->Filename()),
                        App->Width(),App->Height()),
          SurfLmk2D3DBase(surf->surface(), fsurf->surface())
{
VIN("SurfLmk2D3D::SurfLmk2D3D(LoadedSurface *surf)")

  ATTACHMENTS_TBLR(0,92,0,40)
  surf1 = new SurfaceBox(this->Form(),surf,arglist,c);

  ATTACHMENTS_TBLR(0,92,40,80)
  surf2 = new SurfaceBox(this->Form(),fsurf,arglist,c);

  ATTACHMENTS_TBLR(92,100,0,80)
  text_label = App->Label(this->Form(),"",App->SmallFontList(),arglist,c);

  line_buttons    = NULL;
  lines           = NULL;
  lines2          = NULL;
  current_button  = NULL;
  current_line 	  = NULL;
  nbuttons        = 0;
  saver           = 1.;
  saveg           = 0.;
  saveb           = 0.;
  m_bChange       = false;

  sphere_rad = (surf->glsurface()->getMaxExtent())*0.005;

  int bh = App->MedFontHeight()+6;
    
  ATTACHMENTS_TLRH(0,80,100,bh)
  Widget btn = App->Button(this->Form(),"Load 3D LMK",
		   App->SmallFontList(),LoadLandmarksCB,(void*)this,
		   arglist,c);

  ATTACHMENTS_WLRH(btn,80,100,bh)
  btn = App->Button(this->Form(),"Load 2D LMK",
		   App->SmallFontList(),LoadLandmarks2CB,(void*)this,
		   arglist,c);


  ATTACHMENTS_WLRH(btn,80,100,bh)
  btn = App->Button(this->Form(),"Load Indices",
		   App->SmallFontList(),LoadIndicesCB,(void*)this,
		   arglist,c);

  ATTACHMENTS_WLRH(btn,80,100,bh)
  btn = App->Button(this->Form(),"Save 3D LMK",
		    App->SmallFontList(),SaveLandmarksCB,(void*)this,
		    arglist,c);

  ATTACHMENTS_WLRH(btn,80,100,bh)
  btn = App->Button(this->Form(),"Save 2D LMK",
		    App->SmallFontList(),SaveLandmarks2CB,(void*)this,
		    arglist,c);

  ATTACHMENTS_WLRH(btn,80,100,bh)
  btn = App->Button(this->Form(),"Save Indices",
		    App->SmallFontList(),SaveIndicesCB,(void*)this,
		    arglist,c);

  ATTACHMENTS_WLRH(btn,80,100,bh)
  btn = App->Button(this->Form(),"Clear",
		    App->SmallFontList(),ClearLandmarksCB,(void*)this,
		    arglist,c);


  ATTACHMENTS_WLRH(btn,80,100,bh)
  btn = App->Button(this->Form(),"Delete",
		    App->SmallFontList(),DeleteCurrentCB,(void*)this,
		    arglist,c);

  ATTACHMENTS_WLRH(btn,80,100,bh)
  btn = App->Button(this->Form(),"Set Color",
		    App->SmallFontList(),SetColorCB,(XtPointer)this,
		    arglist,c);

  ATTACHMENTS_WLRH(btn,80,100,bh)
  btn = App->Button(this->Form(),"Set Lmk Size",
		    App->SmallFontList(),SetLmkSizeCB,(XtPointer)this,
		    arglist,c);

  ATTACHMENTS_WLRH(btn,80,100,bh)
  btn = App->Toggle(this->Form(),"Change Current",App->SmallFontList(),
			m_bChange,ToggleChangeCB,(XtPointer)this,arglist,c);
  m_wToggle = btn;

  ATTACHMENTS_WLRB(btn,80,100,100)
  Set(XmNscrollBarDisplayPolicy,XmAUTOMATIC);
  Set(XmNscrollBarPlacement,XmBOTTOM_RIGHT);
  Set(XmNscrollingPolicy,XmAUTOMATIC);
  Set(XmNresizeWidth,False);
  Set(XmNspacing,0);
  btn = XmCreateScrolledWindow(this->Form(),"line_form_scrw",arglist,c);
  XtManageChild(btn);
  App->RecursiveColor(btn,SHELLCOLOR);

  Dimension ww,ww2;
  Widget    cw,sb;

  cw = XtNameToWidget(btn,"ClipWindow");
  sb = XtNameToWidget(btn,"VertScrollBar");

  c = 0;
  Set(XmNwidth,&ww);
  XtGetValues(cw,arglist,c);

  c = 0;
  Set(XmNwidth,&ww2);
  XtGetValues(sb,arglist,c);

  ww -= ww2;
  ww -= 5;

  if((ww < 100)||(ww>200)) 
      ww = 100;

  ATTACHMENTS_TBLR(0,100,10,100)
  Set(XmNwidth,ww);
  Set(XmNresizeWidth,False);
  Set(XmNorientation,XmVERTICAL);
  Set(XmNisHomogeneous,False);
  Set(XmNrubberPositioning,False);
  Set(XmNpacking,XmPACK_COLUMN);
  Set(XmNnumColumns,1);
  Set(XmNspacing,1);
  line_form = XmCreateRowColumn(btn, "line_form",arglist,c);
  XtManageChild(line_form);
  XmChangeColor(line_form,SHELLCOLOR);

  surf1->SetPickStyle(PickPoint,PickCB,(u_long)this);
  surf2->SetPickStyle(PickPoint,PickCB,(u_long)this);

  surf->surfaceChanged();
VOUT("SurfLmk2D3D::SurfLmk2D3D(LoadedSurface *surf)")
}

SurfLmk2D3D::~SurfLmk2D3D()
{
VIN("SurfLmk2D3D::~SurfLmk2D3D()")
  clearLandmarks();
  delete surf1;
  delete surf2;
VOUT("SurfLmk2D3D::~SurfLmk2D3D()")
}



