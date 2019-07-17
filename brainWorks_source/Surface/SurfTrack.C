///////////////////////////////////////////////////////////////////////////
//
// File: SurfTrack.C
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
#include <SurfTrack.h>
#include <GLPrimitive.h>
#include <GLSurface.h>
#define LINE_WIDTH	5.0

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

static char *type_names[] = {
  "Sulci",
  "Gyri",
  "Geodesics"
};

static char *pick_style_names[] = {
  "Endpoints",
  "Waypoints"
};

DEFINECB(SurfTrack,SetColorCB,setColor)
DEFINECB(SurfTrack,SetLineWidthCB,setLineWidth)
DEFINECB(SurfTrack,EndCurveCB,endCurve)
DEFINECB(SurfTrack,WriteVertexCB,writeVertex)

static void PickCB(u_long dat, void *s, int indx)
{
VIN("static void PickCB(u_long dat, void *s, int indx)")
  if(dat) ((SurfTrack*)dat)->pickPoint(indx);
VOUT("static void PickCB(u_long dat, void *s, int indx)")
}

 static void SelectButtonCB(Widget wid,XtPointer cld, XtPointer)
{
VIN(" static void SelectButtonCB(Widget wid,XtPointer cld, XtPointer)")
  if(cld) ((SurfTrack*)cld)->selectLine(wid);
VOUT(" static void SelectButtonCB(Widget wid,XtPointer cld, XtPointer)")
}

static void SaveIVLinesCB(Widget wid,XtPointer cld, XtPointer)
{
VIN("static void SaveIVLinesCB(Widget wid,XtPointer cld, XtPointer)")
  if(cld) ((SurfTrack*)cld)->saveIVLines();
VOUT("static void SaveIVLinesCB(Widget wid,XtPointer cld, XtPointer)")
}

static void SavePointCB (Widget wid, XtPointer cld, XtPointer)  
{
VIN("static void SavePointCB (Widget wid, XtPointer cld, XtPointer)  ")
  if (cld ) ((SurfTrack*)cld)->savePoint();
VOUT("static void SavePointCB (Widget wid, XtPointer cld, XtPointer)  ")
}

static void DeleteTrackCB(Widget wid,XtPointer cld, XtPointer)
{
VIN("static void DeleteTrackCB(Widget wid,XtPointer cld, XtPointer)")
  if(cld) ((SurfTrack*)cld)->userDeleteCurrent();
VOUT("static void DeleteTrackCB(Widget wid,XtPointer cld, XtPointer)")
}

static void LoadPointsCB(Widget,XtPointer cld, XtPointer)
{
VIN("static void LoadPointsCB(Widget,XtPointer cld, XtPointer)")
  if(cld) ((SurfTrack*)cld)->loadPoints();
VOUT("static void LoadPointsCB(Widget,XtPointer cld, XtPointer)")
}

static void LoadTracksCB(Widget,XtPointer cld, XtPointer)
{
VIN("static void LoadTracksCB(Widget,XtPointer cld, XtPointer)")
  if(cld) ((SurfTrack*)cld)->loadTracks();
VOUT("static void LoadTracksCB(Widget,XtPointer cld, XtPointer)")
}

static void SaveTracksCB(Widget,XtPointer cld, XtPointer)
{
VIN("static void SaveTracksCB(Widget,XtPointer cld, XtPointer)")
  if(cld) ((SurfTrack*)cld)->saveTracks();
VOUT("static void SaveTracksCB(Widget,XtPointer cld, XtPointer)")
}

static void ClearTracksCB(Widget,XtPointer cld, XtPointer)
{
VIN("static void ClearTracksCB(Widget,XtPointer cld, XtPointer)")
  if(cld) ((SurfTrack*)cld)->clearTracks();
VOUT("static void ClearTracksCB(Widget,XtPointer cld, XtPointer)")
}

static void TypeChangedCB(char *nm, void *dat)
{
VIN("static void TypeChangedCB(char *nm, void *dat)")
  if(dat) ((SurfTrack*)dat)->typeChanged(nm);
VOUT("static void TypeChangedCB(char *nm, void *dat)")
}

static void PickStyleChangedCB(char *nm, void *dat)
{
VIN("static void TypeChangedCB(char *nm, void *dat)")
  if(dat) ((SurfTrack*)dat)->change_pick_style(nm);
VOUT("static void TypeChangedCB(char *nm, void *dat)")
}

void SurfTrack::writeVertex()
{
VIN("void ImageViewer::writeVertex()")
  // write out the text string associated with the vertex number and position
  /* don't know why can't get this from text field of label widget????
  char vertexLabelTxt[1024];
  App->GetLabelText(m_vertexPointedTo, vertexLabelTxt);
  cout << m_vertexLabelTxt << endl;
  */

  if(surface)
  {
    cout << surface->Filename() << ":  ";
  }
  cout << m_vertexLabelTxt << endl;
VOUT("void ImageViewer::writeVertex()")
}

void SurfTrack::saveIVLines() 
{
VIN("void SurfTrack::saveIVLines() ")
  //changed to output coordinates- muge
  int i, j, np;
  if(nindices <= 0) {
    App->ShowMessage("No curves to save");
    return;
  }

  char fname[1024], fname2 [1030];
 retry_save_2:
  fname[0] = 0;

  if(Browser->GetFilename("Save Curve To",fname,1024))  {
    int iStrlen = strlen (fname);
    if (iStrlen > 6 && fname [iStrlen - 6] == '.')
      fname [iStrlen - 6] = 0;
    sprintf (fname2, "%s.coord", fname);

    ofstream fp(fname2,ios::out);
    if(!fp) {
      App->ShowMessage("ERROR: can't open file");
      goto retry_save_2;
    }

    for(i=0; i<nbuttons; i++)
      if (line_buttons[i]==current_button)  break;
    np = indices[i].getNelm();
    fp<<np<<endl;
    for(j=0;j<np;j++) {
      fp<<surface->vertices()[indices[i][j]]<<endl;
    }

    fp.close();

    sprintf (fname2, "%s.index", fname);
    fp.open (fname2, ios::out);
    if (fp)
      {
	//write the data
	fp << "1\n";
	fp << np << endl;
	for (j = 0; j < np; j++)
	  {
	    fp << indices [i] [j] << " ";
	    if (j % 10 == 9) fp << endl;
	  }
	fp << endl;

	fp.close ();
      }
  }
VOUT("void SurfTrack::saveIVLines() ")
}

void SurfTrack::savePoint() 
{
VIN("void SurfTrack::savePoint() ")
  Array1D<int> temp(1);
  if (point1!=-1) {
    temp[0]=point1;
    addIndices(temp);
  
    point1=-1;
    
    lines = GLPrimitive::removeFromList(lines,sphere);
    surf1->setPrimitives(lines);
    sphere = NULL;
  }
VOUT("void SurfTrack::savePoint() ")
}

//
// Parse a line for a number followed by a string (which can contain ws)
//
void SurfTrack::parseName(char *line, int &np, char *cname)
{
VIN("void SurfTrack::parseName(char *line, int &np, char *cname)")
  // FILE adds NL
  line[strlen(line)-1] = 0;

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
VOUT("void SurfTrack::parseName(char *line, int &np, char *cname)")
}

void SurfTrack::loadTracks()
{
VIN("void SurfTrack::loadTracks()")
  clearTracks();

  if(!surface) return;
  int nv = surface->getNumVert();

  char fname[1024];
 retry_load:
  fname[0] = 0;
  if(Browser->GetFilename("Load Curves From",fname,1024))  {
    FILE *fp = fopen(fname,"r");
    if(!fp) {
      App->ShowMessage("ERROR: can't open file");
      goto retry_load;
    }
    int nl,np,i,j;
    fscanf(fp,"%d",&nl);

    indices  = new Array1D<int>[nl];
    nindices = nl;

    char line[255];
    char cname[200];
    Array1D<Point> p;
    double xx,yy,zz,nrmx,nrmy,nrmz;

    for(i=0;i<nl;i++) {
      // this line ensures eating newlines
      while(!feof(fp) && (!fgets(line,255,fp) || (strlen(line)<3)));
      parseName(line,np,cname);
      if((np<=0)||(np>10000)) {
          fclose(fp);
          App->ShowMessage("ERROR: invalid curves file");
          return;
      }
      p.setDim(np);
      indices[i].setDim(np);
      for(j=0;j<np;j++) {
	int elm;
	fscanf(fp,"%d",&elm);

        if((elm<0)||(elm>=nv)) {
          fclose(fp);
          App->ShowMessage("ERROR: these tracks do not belong to this surface");
          clearTracks();
          return;   
        }

	indices[i][j] = elm;

	xx = surface->vertices()[elm].x();
	yy = surface->vertices()[elm].y();
	zz = surface->vertices()[elm].z();
	nrmx = surface->normals()[elm].x();
	nrmy = surface->normals()[elm].y();
	nrmz = surface->normals()[elm].z();

	p[j].set(xx+0.1*nrmx,yy+0.1*nrmy,zz+0.1*nrmz);	
      }
      addLine(p,cname);
    }

    // now load colors/lw (if they are there -- new version)
    float r,g,b,l;
    for(i=0;i<nbuttons;i++) {
        GLPrimitive *p = findLine(line_buttons[i]);
        if(p && (fscanf(fp,"%f %f %f %f",&r,&g,&b,&l) == 4)) {
           p->setColor(r,g,b);
           p->setLineWidth(l);
           XmChangeColor(line_buttons[i],App->RGBToPixel(r,g,b));
        }
    }

  fclose(fp);
  surf1->redraw();
  }


VOUT("void SurfTrack::loadTracks()")
}

void SurfTrack::loadPoints()
{
VIN("void SurfTrack::loadPoints()")
    clearTracks();
    
    if(!surface) return;
    int nv = surface->getNumVert();
    
    char fname[1024];
    
    ifstream indexFile;

 retry_load:
    fname[0] = 0;
    if(Browser->GetFilename("Load Track Indices From", fname, 1024))
    {
        indexFile.open(fname, ifstream::in);
        if(!indexFile) {
        App->ShowMessage("ERROR: can't open file");
        goto retry_load;
    }
    
    int indexCount = 0;
    double indexD;
    while (!indexFile.eof())
    {
        indexFile >> indexD;

        if (!indexFile.fail())
            indexCount++;
    }

    indexFile.close();
    
    if (indexCount < 2)
    {
        char errorMsg[1024];
        sprintf(errorMsg, "ERROR: less than two indices in file:  %s", fname);
        App->ShowMessage(errorMsg);
        return;
    }

    int indexStart, indexEnd;

    indexFile.clear();
    indexFile.open(fname, ifstream::in);

    indexFile >> indexD;
    indexStart = indexD;
    for (int i = 0; i < indexCount - 1; i++)
    {
      // these are integers, but in some files appear in exponential format
      indexFile >> indexD;
      indexEnd = indexD;

      if ((indexStart < 0) || (indexStart >= nv) ||
           (indexEnd < 0) || (indexEnd >= nv))
      {
          indexFile.close();
          App->ShowMessage("ERROR: these tracks do not belong to this surface");
          clearTracks();
          return;
      }

      // we have a start and an endpoint - find a track
      pickPoint(indexStart);
      pickPoint(indexEnd);

      indexStart = indexEnd;
    }

    indexFile.close();
    surf1->redraw();
  }

VOUT("void SurfTrack::loadPoints()")
}


void SurfTrack::saveTracks()
{
VIN("void SurfTrack::saveTracks()")
  int i;
  if(nindices <= 0) {
    App->ShowMessage("No curves to save");
    return;
  }

  char cname[200];
  char fname[1024];
 retry_save:
  fname[0] = 0;
  if(Browser->GetFilename("Save Curves To",fname,1024))  {
    FILE *fp = fopen(fname,"w");
    if(!fp) {
      App->ShowMessage("ERROR: cant open file");
      goto retry_save;
    }
    fprintf(fp,"%d\n",nindices);
    for(i=0;i<nindices;i++) {
      Widget txt = XtNameToWidget(XtParent(line_buttons[i]),"TextWidget");
      if(txt) App->GetString(txt,cname,200);
      else    sprintf(cname,"<no-name>"); // this shouldnt happen

      int np = indices[i].getNelm();	
      fprintf(fp,"%d %s\n",np,cname);
      for(int j=0;j<np;j++) {
	fprintf(fp,"%d ",indices[i][j]);
	if(j % 10 == 9) fprintf(fp,"\n");
      }
      fprintf(fp,"\n");
    }	

    // now save r,g,b,lw
    for(i=0;i<nbuttons;i++) {
        GLPrimitive *p = findLine(line_buttons[i]);
        float r,g,b,l;
        if(!p) {
           fprintf(fp,"1.0 0.0 0.0 1.0\n");
        } else {
           p->getColor(r,g,b);
           p->getLineWidth(l);
           fprintf(fp,"%.2f %.2f %.2f %.2f\n",r,g,b,l);
        }
    }
    fclose(fp);
  
    // add two more files.
    // create one long curve from all the curve pieces
    // write out an index file and coordinates file for the curve
    char fname2[1024];
    char* extensionPos = strrchr (fname, '.');
    if (extensionPos)
      extensionPos = '\0';

    sprintf (fname2, "%s.coord", fname);
    ofstream fpCoord(fname2,ios::out);
    if(!fpCoord)
    {
      char errorMsg[1024];
      sprintf(errorMsg, "ERROR: can't open coordinates file:  %d", fname2);
      App->ShowMessage(errorMsg);
      return;
    }

    sprintf (fname2, "%s.index", fname);
    ofstream fpIndex(fname2,ios::out);
    if(!fpIndex)
    {
      char errorMsg[1024];
      sprintf(errorMsg, "ERROR: can't open index file:  %d", fname2);
      App->ShowMessage(errorMsg);
      return;
    }

    // one curve
    fpIndex << "1\n";
    
    // need the index count
    int previousIndex = -1;
    int indexCount = 0;
    for (i = 0; i < nindices; i++)
    {
      int np = indices[i].getNelm();
      indexCount += np;
      // the indeces are numbered backwards in these curves
      // so to check if they are connected, see if
      // the last one of this curve is equal to the first one of the previous curve
      // since the points get thrown out, don't count them...
      if (indices[i][np-1] == previousIndex)
        indexCount--;
      previousIndex = indices[i][0];
    }   
    
    fpIndex << indexCount << " CompletePath" << endl;
    fpCoord << indexCount << endl;
    // use this for breaking up lines
    indexCount = 0;
    previousIndex = -1;
    for (i = 0; i < nindices; i++)
    {
      int np = indices[i].getNelm();
      for (int j = np - 1; j >= 0; j--)
      {
        // don't repeat indices
        if (previousIndex != indices[i][j])
        {
            fpCoord << surface->vertices()[indices[i][j]] << endl;
            fpIndex << indices[i][j] << " ";
            if(indexCount++ % 10 == 9)
                fpIndex << endl;
        }
        previousIndex = indices[i][j];
      }
    }

    fpIndex << endl;

    fpIndex.close ();
    fpCoord.close ();

  }
VOUT("void SurfTrack::saveTracks()")
}


void SurfTrack::clearTracks()
{
VIN("void SurfTrack::clearTracks()")
  point1 = -1;
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
  nindices = 0;
  if(indices) delete [] indices;
  indices  = NULL;
  current_button  = NULL;
  current_line    = NULL;

  surf1->redraw();
VOUT("void SurfTrack::clearTracks()")
}

void SurfTrack::setColor()
{
  float r,g,b;
  if(!current_line) {
      App->ShowMessage("Please select a current track first");
      return;
  }

  current_line->getColor(r,g,b);
  if(GetColor(&r,&g,&b)) {
     current_line->setColor(r,g,b);
     XmChangeColor(current_button,App->RGBToPixel(r,g,b));
     surf1->redraw();
  }
}

void SurfTrack::setLineWidth()
{
  float nl,l;
  char tstr[100];
  if(!current_line) {
      App->ShowMessage("Please select a current track first");
      return;
  }
  current_line->getLineWidth(l);

retry_lw:
  sprintf(tstr,"%.2f",l);
  if(App->GetText("New Line Width:",tstr,99)) {
       if((sscanf(tstr,"%f",&nl)!=1)||(nl<=0.)||(nl>50)) {
           App->ShowMessage("Please enter a line width > 0 and < 50");
           goto retry_lw;
       }
       current_line->setLineWidth(nl);
       surf1->redraw();
  }
}


void SurfTrack::endCurve()
{
  point1 = -1;
  mAddNew = true;
}

GLPrimitive *SurfTrack::findLine(Widget w)
{
  //
  // find current GLPrimitive line corresponding
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


void SurfTrack::selectLine(Widget w)
{
VIN("void SurfTrack::selectLine(Widget w)")
//  if(current_button) 
//    XmChangeColor(current_button,BLUE);

  current_button = w;
  current_line   = findLine(w);

//  if(current_button)
//    XmChangeColor(current_button,RED);

VOUT("void SurfTrack::selectLine(Widget w)")
}


void SurfTrack::typeChanged(char *new_type)
{
VIN("void SurfTrack::typeChanged(char *new_type)")
  if(!strcmp(new_type,"Sulci"))
    tracker.setTrackType(SurfaceTracker::Sulci);
  else if(!strcmp(new_type,"Gyri"))
    tracker.setTrackType(SurfaceTracker::Gyri);
  else if(!strcmp(new_type,"Geodesics"))
    tracker.setTrackType(SurfaceTracker::Geodesi);
VOUT("void SurfTrack::typeChanged(char *new_type)")
}

void SurfTrack::change_pick_style(char *new_type)
{
VIN("void SurfTrack::typeChanged(char *new_type)")
  // includes a clear of all existing tracks
  if(!strcasecmp(new_type,"Endpoints"))
  {
    if (mPickStyle != EndPoint)
    {
      XtSetSensitive(m_endCurveBtn,FALSE);
      clearTracks();
      mPickStyle = EndPoint;
    }
  }
  else if(!strcasecmp(new_type,"Waypoints"))
  {
    if (mPickStyle != WayPoint)
    {
      mAddNew = true;
      clearTracks();
      mPickStyle = WayPoint;
      XtSetSensitive(m_endCurveBtn,TRUE);
    }
  }
VOUT("void SurfTrack::typeChanged(char *new_type)")
}


Widget SurfTrack::addButton(char *name)
{
VIN("Widget SurfTrack::addButton(char *name)")
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
VOUT("Widget SurfTrack::addButton(char *name)")
  return(btn);
}


void SurfTrack::addIndices(Array1D<int> inds)
{
VIN("void SurfTrack::addIndices(Array1D<int> inds)")
  Array1D<int> *newindices = new Array1D<int>[nindices+1];
  int i;
  for(i=0;i<nindices;i++) 
    newindices[i] = indices[i];

  if(indices) delete [] indices;
  indices = newindices;
  indices[i] = inds;
  nindices++;
VOUT("void SurfTrack::addIndices(Array1D<int> inds)")
}


void SurfTrack::addToIndices(Array1D<int> inds, int index)
{
VIN("void SurfTrack::addIndices(Array1D<int> inds)")

  // shouldn't happen - but dont know what this silent failure will do
  if (index >= nindices) return;

  // adding in appended surf track mode (mPickType == WayPoint)
  Array1D<int> *newindices = new Array1D<int>[nindices];

  int i;
  for (i = 0; i < index; i++) 
    newindices[i] = indices[i];

  int newCount = inds.getNelm();
  int totCount = newCount + indices[index].getNelm();
  newindices[index].setDim(totCount);
  for(i = 0; i < newCount; i++)
    newindices[index][i] = inds[i];
  int ix = 0;
  for(i = newCount; i < totCount; i++)
    newindices[index][i] = indices[index][ix++];

  for (i = index + 1; i < nindices; i++) 
    newindices[i] = indices[i];
/*
  for(int i = 0; i < oldCount;i++)
    newindices[0][i] = indices[0][i];
  int ix = 0;
  for(int i = oldCount; i < newCount;i++)
    newindices[0][i] = inds[ix++];
*/
  delete [] indices;
  indices = newindices;
VOUT("void SurfTrack::addIndices(Array1D<int> inds)")
}

void SurfTrack::addLine(Array1D<Point> &p,char *name)
{
VIN("void SurfTrack::addLine(Array1D<Point> &p,char *name)")
  if(p.getNelm() <= 0) return;
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

  GLPrimitive *newline = new GLPrimitive(p);
  newline->setHook((void*)line_buttons[i]);
  newline->setColor(DEFRED,DEFGREEN,DEFBLUE);
  newline->setLineWidth(LINE_WIDTH);
  //
  // add new GL line to the list 
  //
  lines = GLPrimitive::addToList(lines,newline);
  surf1->setPrimitives(lines);

  selectLine(line_buttons[i]);
VOUT("void SurfTrack::addLine(Array1D<Point> &p,char *name)")
}

void SurfTrack::userDeleteCurrent()
{
VIN("void SurfTrack::userDeleteCurrent()")
  deleteCurrent();
  mAddNew = true;
  point1 = -1;
VOUT("void SurfTrack::userDeleteCurrent()")
}

void SurfTrack::deleteCurrent()
{
VIN("void SurfTrack::deleteCurrent()")
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
  for(;i<nbuttons-1;i++) {
    line_buttons[i] = line_buttons[i+1];
    indices[i].setDim(indices[i+1].getNelm());
    indices[i]      = indices[i+1];
  }
 
  indices[nbuttons-1].setDim(0);

  nbuttons--;
  nindices--;

  current_button = NULL;
  current_line   = NULL;
  surf1->redraw();
VOUT("void SurfTrack::deleteCurrent()")
}


void SurfTrack::getLineFromIndices(Array1D<Point>& path, int pathIx)
{
VIN("void SurfTrack::getLineFromIndices(Array1D<Point> path)")
  char estr[1024];
  Array1D<int>& reti = indices[pathIx];
  path.setDim(reti.getNelm());
  double xx,yy,zz,nrmx,nrmy,nrmz;
  for(int i=0;i<reti.getNelm();i++)
  {
    int elm = reti[i];
    xx = surface->vertices()[elm].x();
    yy = surface->vertices()[elm].y();
    zz = surface->vertices()[elm].z();
    path[i].set(xx,yy,zz);
  }
VOUT("void SurfTrack::getLineFromIndices(Array1D<Point> path)")
}


void SurfTrack::pickPoint(int indx)
{
VIN("void SurfTrack::pickPoint(int indx)")
  char estr[1024];
  float x,y,z,r;
  x = surface->vertices()[indx].x();
  y = surface->vertices()[indx].y();
  z = surface->vertices()[indx].z();
  sprintf(estr,"%d : (%.2f, %.2f, %.2f)",indx,x,y,z);
  if(point1 == -1) {
    //	r = 1.0;
    r = (surf1->getCurSurf()->glsurface()->getMaxExtent())*0.005;
    if(sphere) delete sphere;
    sphere = new GLPrimitive(x,y,z,r);
    sphere->setColor(1.0,0.0,0.0);
    lines = GLPrimitive::addToList(lines,sphere);
    surf1->setPrimitives(lines);
    point1 = indx;
  }
  else {
    point2 = indx;
    lines = GLPrimitive::removeFromList(lines,sphere);
    surf1->setPrimitives(lines);
    sphere = NULL;

    Array1D<int> reti;
    ItXECode rval = tracker.findPath(point1,point2,reti);
    if (mPickStyle == EndPoint)
      point1 = -1;
    else
      point1 = point2;

    if((rval != ItXSuccess) || (reti.getNelm() < 2)) {
      App->ShowMessage("ERROR: no path found");
      return;
    }

    if(!(surface->hasUNorms()))
      surface->genUnitNormals();

    Array1D<Point> p;

    // add to set of lines - if EndPoint
    if (mPickStyle == EndPoint || mAddNew)
    {
      this->addIndices(reti);
      this->getLineFromIndices(p, nindices - 1);
      this->addLine(p);
      mAddNew = false;
    }
    else  // add to current line
    {
      Array1D<int> currentIndices;
      if (indices) currentIndices = indices[nindices - 1];
      this->deleteCurrent();
      this->addIndices(currentIndices);
      this->addToIndices(reti, nindices - 1);
      this->getLineFromIndices(p, nindices - 1);
      this->addLine(p);
    }
  }
  sprintf(m_vertexLabelTxt, "Vertex:  %s", estr);
  App->SetLabelText(m_vertexPointedTo,m_vertexLabelTxt,App->SmallFontTag());
  surf1->redraw();
VOUT("void SurfTrack::pickPoint(int indx)")
}


//
// Same as SurfaceViewer(), only 1 view (surf1)
//
SurfTrack::SurfTrack(LoadedSurface *surf) :
          KDChildWindow(App->FilenameTail(surf->Filename()),
                        App->Width(),App->Height())
//		KDChildWindow("Surface",50,50)
{
VIN("SurfTrack::SurfTrack(LoadedSurface *surf)")
//  Hide();

  strcpy(m_vertexLabelTxt, "");
  c =0;
  Set(XmNleftAttachment,XmATTACH_FORM);
  Set(XmNrightAttachment,XmATTACH_FORM);
  Set(XmNtopAttachment,XmATTACH_FORM);
  Set(XmNbottomAttachment,XmATTACH_POSITION);
  Set(XmNbottomPosition,100);
  surf1 = new SurfaceBox(this->Form(),surf,arglist,c);
  nindices        = 0;
  indices         = NULL;
  line_buttons    = NULL;
  lines           = NULL;
  sphere          = NULL;
  current_button 	= NULL;
  current_line 	= NULL;
  nbuttons        = 0;
  point1 = point2 = -1;
  lsurface        = surf;
  surface         = surf->surface();
  mPickStyle      = EndPoint;
  mAddNew = true;

  c = 0;
  Set(XmNrightAttachment,XmATTACH_POSITION);
  Set(XmNrightPosition,80);
  XtSetValues(surf1->Main(),arglist,c);

  int bh = App->BigFontHeight()+6;

  c = 0;
  Set(XmNleftAttachment,XmATTACH_POSITION);
  Set(XmNleftPosition,80);
  Set(XmNtopAttachment,XmATTACH_FORM);
  Set(XmNrightAttachment,XmATTACH_FORM);
  Set(XmNbottomAttachment,XmATTACH_NONE);
  Set(XmNheight,bh);
  track_type = new IconList(this->Form(),"Track",IXNoIcon,
			    type_names,XtNumber(type_names),arglist,c);

  track_type->AddExtraCallback(TypeChangedCB,(void*)this);

  c = 0;
  Set(XmNleftAttachment,XmATTACH_POSITION);
  Set(XmNleftPosition,90);
  Set(XmNtopAttachment,XmATTACH_FORM);
  Set(XmNrightAttachment,XmATTACH_FORM);
  Set(XmNbottomAttachment,XmATTACH_NONE);
  Set(XmNheight,bh);
  pick_style_type = new IconList(this->Form(),"Pick Style",IXNoIcon,
                pick_style_names,XtNumber(pick_style_names),arglist,c);

  pick_style_type->AddExtraCallback(PickStyleChangedCB,(void*)this);

  c = 0;
  Set(XmNleftAttachment,XmATTACH_POSITION);
  Set(XmNleftPosition,80);
  Set(XmNtopAttachment,XmATTACH_WIDGET);
  Set(XmNtopWidget,track_type->Form());
  Set(XmNrightAttachment,XmATTACH_FORM);
  Set(XmNbottomAttachment,XmATTACH_NONE);
  Set(XmNheight,bh);

  Widget btn = App->Button(this->Form(),"Load Tracks",
			   App->MedFontList(),LoadTracksCB,(void*)this,
			   arglist,c);

  c = 0;
  Set(XmNleftAttachment,XmATTACH_POSITION);
  Set(XmNleftPosition,80);
  Set(XmNtopAttachment,XmATTACH_WIDGET);
  Set(XmNtopWidget,btn);
  Set(XmNrightAttachment,XmATTACH_FORM);
  Set(XmNbottomAttachment,XmATTACH_NONE);
  Set(XmNheight,bh);

  btn = App->Button(this->Form(),"Load Track EndPoints",
            App->MedFontList(),LoadPointsCB,(void*)this,
            arglist,c);

  c = 0;
  Set(XmNleftAttachment,XmATTACH_POSITION);
  Set(XmNleftPosition,80);
  Set(XmNtopAttachment,XmATTACH_WIDGET);
  Set(XmNtopWidget,btn);
  Set(XmNrightAttachment,XmATTACH_FORM);
  Set(XmNbottomAttachment,XmATTACH_NONE);
  Set(XmNheight,bh);

  btn = App->Button(this->Form(),"Save",
		    App->MedFontList(),SaveTracksCB,(void*)this,
		    arglist,c);

  c = 0;
  Set(XmNleftAttachment,XmATTACH_POSITION);
  Set(XmNleftPosition,80);
  Set(XmNtopAttachment,XmATTACH_WIDGET);
  Set(XmNtopWidget,btn);
  Set(XmNrightAttachment,XmATTACH_FORM);
  Set(XmNbottomAttachment,XmATTACH_NONE);
  Set(XmNheight,bh);

  btn = App->Button(this->Form(),"Clear",
		    App->MedFontList(),ClearTracksCB,(void*)this,
		    arglist,c);


  c = 0;
  Set(XmNleftAttachment,XmATTACH_POSITION);
  Set(XmNleftPosition,80);
  Set(XmNtopAttachment,XmATTACH_WIDGET);
  Set(XmNtopWidget,btn);
  Set(XmNrightAttachment,XmATTACH_FORM);
  Set(XmNbottomAttachment,XmATTACH_NONE);
  Set(XmNheight,bh);

  btn = App->Button(this->Form(),"Delete",
		    App->MedFontList(),DeleteTrackCB,(void*)this,
		    arglist,c);

  c = 0;
  Set(XmNleftAttachment,XmATTACH_POSITION);
  Set(XmNleftPosition,80);
  Set(XmNtopAttachment,XmATTACH_WIDGET);
  Set(XmNtopWidget,btn);
  Set(XmNrightAttachment,XmATTACH_FORM);
  Set(XmNbottomAttachment,XmATTACH_NONE);
  Set(XmNheight,bh);

  btn = App->Button(this->Form(),"Save One Path",
		    App->MedFontList(),SaveIVLinesCB,(void*)this,
		    arglist,c);

  c = 0;
  Set(XmNleftAttachment,XmATTACH_POSITION);
  Set(XmNleftPosition,80);
  Set(XmNtopAttachment,XmATTACH_WIDGET);
  Set(XmNtopWidget,btn);
  Set(XmNrightAttachment,XmATTACH_FORM);
  Set(XmNbottomAttachment,XmATTACH_NONE);
  Set(XmNheight,bh);

  btn = App->Button(this->Form(),"Save Point",
		    App->MedFontList(),SavePointCB,(XtPointer)this,
		    arglist,c);

  ATTACHMENTS_WLRH(btn,80,100,bh)
  btn = App->Button(this->Form(),"Set Color",
		    App->MedFontList(),SetColorCB,(XtPointer)this,
		    arglist,c);

  ATTACHMENTS_WLRH(btn,80,100,bh)
  btn = App->Button(this->Form(),"Set Line Width",
		    App->MedFontList(),SetLineWidthCB,(XtPointer)this,
		    arglist,c);

  ATTACHMENTS_WLRH(btn,80,100,bh)
  m_endCurveBtn = btn = App->Button(this->Form(),"End Curve",
            App->MedFontList(),EndCurveCB,(XtPointer)this,
            arglist,c);
  XtSetSensitive(m_endCurveBtn,FALSE);

  ATTACHMENTS_WLRH(btn,80,85,bh)
  Widget btn2 = m_writeVertexBtn =
    App->Button(this->Form(),"Stdout",App->SmallFontList(),
                    WriteVertexCB,(XtPointer)this,arglist,c);

  ATTACHMENTS_WLRH(btn,85,100,bh + 5)
  Set(XmNalignment,XmALIGNMENT_BEGINNING);
  strcpy(m_vertexLabelTxt, "");
  m_vertexPointedTo =
    App->Label(this->Form(),m_vertexLabelTxt,
               App->SmallFontList(),arglist,c);

  btn = btn2;
  ATTACHMENTS_WLRB(btn2,80,100,100)
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

  tracker.setSurface(surf->surface(),SurfaceTracker::Sulci);
  surf->surfaceChanged();
//  SetMinDimension(App->Width()/3,App->Height()/3);
//  SetSize(App->Width(),App->Height());
//  Show();
VOUT("SurfTrack::SurfTrack(LoadedSurface *surf)")
}


SurfTrack::~SurfTrack()
{
VIN("SurfTrack::~SurfTrack()")

  // need to delete this for surface unreference

  delete surf1;

  clearTracks();
 if(track_type) delete track_type;
VOUT("SurfTrack::~SurfTrack()")
}















