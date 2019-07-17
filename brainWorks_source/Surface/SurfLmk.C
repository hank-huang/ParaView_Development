///////////////////////////////////////////////////////////////////////////
//
// File: SurfLmk.C
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
#include <SurfLmk.h>

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

DEFINECB(SurfLmk,SetColorCB,setColor)
DEFINECBP(SurfLmk,SelectButtonCB,selectLandmark,wid)
DEFINECB(SurfLmk,DeleteCurrentCB,deleteCurrent)
DEFINECB(SurfLmk,LoadLandmarksCB,loadLandmarks)
DEFINECB(SurfLmk,SaveLandmarksCB,saveLandmarks)
DEFINECB(SurfLmk,LoadIndicesCB,loadIndices)
DEFINECB(SurfLmk,SaveIndicesCB,saveIndices)
DEFINECB(SurfLmk,ClearLandmarksCB,clearLandmarks)
DEFINECB(SurfLmk,ToggleChangeCB,toggleChange)
DEFINECB(SurfLmk,SetLmkSizeCB,setLmkSize)

static void PickCB(u_long dat, void *s, int indx)
{
  if(dat) ((SurfLmk*)dat)->pickPoint(indx);
}

void SurfLmk::toggleChange()
{
  m_bChange = (App->GetToggleValue(m_wToggle) == TRUE);
}

void SurfLmk::setLmkSize()
{
  float f;
  GLPrimitive *p;
  f = sphere_rad;
  if(GetFloatValue("Enter Lmk Size",0.00001,100,f)) {
     sphere_rad = f;
     for(p=lines;p;p=p->Next())
         p->setRadius(sphere_rad);
     surf1->redraw();
  }
}

void SurfLmk::loadIndices()
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
   indices.clear();
   while(!fp.eof()) {
       fp.getline(line,1023);
       if(sscanf(line,"%ld %n",&idx,&j) < 1)
          continue;
       if(idx >= surface->getNumVert()) {
          fp.close();
          clearLandmarks();
          App->ShowMessage("ERROR: the landmark indices dont fit this surface");
          return;
       }
       sprintf(lname,"%s",&(line[j]));
       indices.add(idx);
       float x = surface->vertices()[idx].x();
       float y = surface->vertices()[idx].y();
       float z = surface->vertices()[idx].z();
       landmarks.add(x,y,z,lname,1);
       addLandmark(x,y,z,lname);
   }
   fp.close();
}

void SurfLmk::saveIndices()
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


unsigned int SurfLmk::findLandmark(float x, float y, float z)
{
  unsigned int i,mini;
  double dx,dy,dz,d,mind;

  if(!surface) 
      return((unsigned int)-1);
   
  for(i=0;i<surface->getNumVert();i++) {
      dx = x - surface->vertices()[i].x();
      dy = y - surface->vertices()[i].y();
      dz = z - surface->vertices()[i].z();
      d  = dx*dx+dy*dy+dz*dz;
      if((i==0)||(d<mind)) {
         mind = d;
         mini = i;
      }
  }
  return(mini);
}

void SurfLmk::loadLandmarks()
{
VIN("void SurfLmk::loadLandmarks()")

  clearLandmarks();

  if(!surface) return;

  int i,j;
  const Point *P;
  char fname[1024];
  char lname[1024];

 retry_load:
  // fname[0] = 0;
  // know the surface name - use that for the lmk file name
  sprintf(fname, "%s.lmk", App->FilenameRemoveExtension(surfaceFileName));
  if(Browser->GetFilename("Load Landmarks From",fname,1024, fname))  {

     if(landmarks.load(fname) != ItXSuccess) {
        App->ShowMessage("ERROR: landmarks could not be loaded");
        goto retry_load;
     }

     indices.setSize(landmarks.num());

     for(i=0;i<landmarks.num();i++) {
         if(!(P=landmarks.point(i))) {
                indices[i] = (unsigned int)-1;
		continue;
         }

         if(landmarks.name(i))
            sprintf(lname,"%s",landmarks.name(i));
         else sprintf(lname,"<no-name>");
         addLandmark(P->x(),P->y(),P->z(),lname);
         indices[i] = findLandmark(P->x(),P->y(),P->z());
     }

     surf1->redraw();
  }

VOUT("void SurfLmk::loadLandmarks()")
}


void SurfLmk::saveLandmarks()
{
VIN("void SurfLmk::saveLandmarks()")

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

     if(landmarks.save(fname) != ItXSuccess) {
        App->ShowMessage("ERROR: landmarks did not save");
        goto retry_save;
     }
  }
VOUT("void SurfLmk::saveLandmarks()")
}


void SurfLmk::clearLandmarks()
{
VIN("void SurfLmk::clearLandmarks()")
  //
  // clear all GLPrimitive lines
  //
  while(lines) lines = GLPrimitive::removeFromList(lines,lines);
  surf1->setPrimitives(NULL);

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

  landmarks.clear();
  indices.clear();

  surf1->redraw();
VOUT("void SurfLmk::clearLandmarks()")
}

void SurfLmk::setColor()
{
  float r,g,b;
  if(!current_line) {
      App->ShowMessage("Please select a current landmark first");
      return;
  }

  current_line->getColor(r,g,b);
  if(GetColor(&r,&g,&b)) {
     current_line->setColor(r,g,b);
     XmChangeColor(current_button,App->RGBToPixel(r,g,b));
     surf1->redraw();
     saver = r;
     saveg = g;
     saveb = b;
  }
}


GLPrimitive *SurfLmk::findLine(Widget w)
{
  //
  // find current GLPrimitive sphere corresponding
  // to the current button Widget
  //
  GLPrimitive *s = NULL;
  GLPrimitive *p;
  for(p=lines;p;p=p->Next())
    if(p->getHook() == (void*)w) {
      s = p;
      break;
    }
  return(s);
}


void SurfLmk::selectLandmark(Widget w)
{
VIN("void SurfLmk::selectLandmark(Widget w)")
  if(!w) return;

  if(current_line)
     current_line->setColor(saver,saveg,saveb);

  if(current_button)
     XmChangeColor(current_button,App->RGBToPixel(saver,saveg,saveb));

  current_button = w;
  current_line   = findLine(w);

  if(current_line)
     current_line->getColor(saver,saveg,saveb);

  if(current_line)
     current_line->setColor(1.f,0.f,0.f);

  XmChangeColor(current_button,App->RGBToPixel(1.f,0.f,0.f));
  surf1->redraw();

VOUT("void SurfLmk::selectLandmark(Widget w)")
}


Widget SurfLmk::addButton(char *name)
{
VIN("Widget SurfLmk::addButton(char *name)")
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
VOUT("Widget SurfLmk::addButton(char *name)")
  return(btn);
}


void SurfLmk::addLandmark(float x, float y, float z,char *name)
{
VIN("void SurfLmk::addLandmark()")

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
  //
  // add new GL line to the list 
  //
  lines = GLPrimitive::addToList(lines,newline);
  surf1->setPrimitives(lines);

  selectLandmark(line_buttons[i]);
VOUT("void SurfLmk::addLandmark()")
}


void SurfLmk::deleteCurrent()
{
VIN("void SurfLmk::deleteCurrent()")
  if(!current_button && !current_line)
    return;

  if(current_button) {
    XtUnmapWidget(XtParent(current_button));
    XtDestroyWidget(XtParent(current_button));
  }

  lines = GLPrimitive::removeFromList(lines,current_line);
  surf1->setPrimitives(lines);

  int i;
  for(i=0;i<nbuttons-1;i++)
    if(line_buttons[i] == current_button)
      break;

  if(i < nbuttons) {
     landmarks.remove(i);
     indices.removeAtIndex(i);
  }

  for(;i<nbuttons-1;i++) {
    line_buttons[i] = line_buttons[i+1];
  }
 
  nbuttons--;

  current_button = NULL;
  current_line   = NULL;
  surf1->redraw();
VOUT("void SurfLmk::deleteCurrent()")
}



void SurfLmk::pickPoint(int indx)
{
VIN("void SurfLmk::pickPoint(int indx)")
 int i;
 float x,y,z,r;
 char tstr[200];
 x = surface->vertices()[indx].x();
 y = surface->vertices()[indx].y();
 z = surface->vertices()[indx].z();

 sprintf(tstr,"Vert %d : (%f, %f, %f)",indx,x,y,z);
 App->SetLabelText(text_label,tstr,App->SmallFontTag());

 if(m_bChange) {
  for(i=0;i<nbuttons;i++)
    if(line_buttons[i] == current_button)
      break;

  if(i < nbuttons) {
     landmarks.set(i,x,y,z,landmarks.name(i),true);
     indices[i] = (unsigned int)indx;
     if(current_line)
       current_line->setCenter(x,y,z);
  }

 } else {
   addLandmark(x,y,z);
   landmarks.add(x,y,z,"<no-name>",1);
   indices.add((unsigned int)indx);
 }

 surf1->redraw();
VOUT("void SurfLmk::pickPoint(int indx)")
}


//
// Same as SurfaceViewer(), only 1 view (surf1)
//
SurfLmk::SurfLmk(LoadedSurface *surf) :
          KDChildWindow(App->FilenameTail(surf->Filename()),
                        App->Width(),App->Height())
{
VIN("SurfLmk::SurfLmk(LoadedSurface *surf)")

  ATTACHMENTS_TBLR(0,92,0,80)
  surf1 = new SurfaceBox(this->Form(),surf,arglist,c);

  ATTACHMENTS_TBLR(92,100,0,80)
  text_label = App->Label(this->Form(),"",App->SmallFontList(),arglist,c);


  line_buttons    = NULL;
  lines           = NULL;
  current_button  = NULL;
  current_line 	  = NULL;
  nbuttons        = 0;
  saver           = 1.;
  saveg           = 0.;
  saveb           = 0.;
  strcpy(surfaceFileName, App->FilenameTail(surf->Filename()));
  surface         = surf->surface();
  m_bChange       = false;

  sphere_rad = (surf->glsurface()->getMaxExtent())*0.005;

  int bh = App->MedFontHeight()+6;

  ATTACHMENTS_TLRH(0,80,100,bh)
  Widget btn = App->Button(this->Form(),"Load LMK",
		   App->SmallFontList(),LoadLandmarksCB,(void*)this,
		   arglist,c);

  ATTACHMENTS_WLRH(btn,80,100,bh)
  btn = App->Button(this->Form(),"Load Indices",
		   App->SmallFontList(),LoadIndicesCB,(void*)this,
		   arglist,c);

  ATTACHMENTS_WLRH(btn,80,100,bh)
  btn = App->Button(this->Form(),"Save LMK",
		    App->SmallFontList(),SaveLandmarksCB,(void*)this,
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

  surf->surfaceChanged();
VOUT("SurfLmk::SurfLmk(LoadedSurface *surf)")
}


SurfLmk::~SurfLmk()
{
VIN("SurfLmk::~SurfLmk()")
  clearLandmarks();
  delete surf1;
VOUT("SurfLmk::~SurfLmk()")
}


