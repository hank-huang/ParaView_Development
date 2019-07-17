///////////////////////////////////////////////////////////////////////////
//
// File: SurfaceLoader.C
//
// Author: Keith Doolittle
//
// Purpose:
//
///////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <KDApplication.h>
#include <KDShell.h>
#include <X11/keysymdef.h>
#include <Xm/PushBG.h>
#include <Xm/CascadeBG.h>
#include <Xm/Form.h>
#include <Xm/Text.h>
#include <Xm/RowColumn.h>
#include <TDL/Surface.h>
#include <TDL/ByuSurface.h>
#include <TDL/SurfaceUtils.h>
#include <LoadedSurface.h>
#include <GLSurface.h>
#include <SurfaceLoader.h>
#include <FlatMap.h>
#include <FileBrowser.h>
#include <Icon.h>
#include <AppIcons.h>
#include <Main.h>

#define BACKCOLOR	WORKBACKGROUNDCOLOR

typedef struct {
	SurfaceLoaderType type;
	char	   *name;
	Widget widget;
	} surf_type;

static surf_type surface_types[] = {
	{ SLByuSurface, "Byu Surface            ",NULL },
	{ SLFlatMap,    "Flat Map               ",NULL }
	};

static void SurfaceLoaderFileSelectCB(Widget widget, XtPointer cld, XtPointer cad)
{
VIN("static void SurfaceLoaderFileSelectCB(Widget widget, XtPointer cld, XtPointer cad)")
SurfaceLoader *ild = (SurfaceLoader*)cld;
if(ild) ild->FileSelectCB();
VOUT("static void SurfaceLoaderFileSelectCB(Widget widget, XtPointer cld, XtPointer cad)")
}

static void SurfaceLoaderOkayCB(Widget widget, XtPointer cld, XtPointer cad)
{
VIN("static void SurfaceLoaderOkayCB(Widget widget, XtPointer cld, XtPointer cad)")
SurfaceLoader *ild = (SurfaceLoader*)cld;
if(ild) ild->OkayCB();
VOUT("static void SurfaceLoaderOkayCB(Widget widget, XtPointer cld, XtPointer cad)")
}

static void SurfaceLoaderCancelCB(Widget widget, XtPointer cld, XtPointer cad)
{
VIN("static void SurfaceLoaderCancelCB(Widget widget, XtPointer cld, XtPointer cad)")
SurfaceLoader *ild = (SurfaceLoader*)cld;
if(ild) ild->CancelCB();
VOUT("static void SurfaceLoaderCancelCB(Widget widget, XtPointer cld, XtPointer cad)")
}


void SurfaceLoader::FileSelectCB()
{
VIN("void SurfaceLoader::FileSelectCB()")

  int returnLineLength = MAXPATHLEN;
  if (multipleSelection)
    returnLineLength *= 4;
    
  char rstring[returnLineLength];
  rstring[0] = 0;

  char *pref;
  switch(getSurfaceType()) {
    case SLByuSurface: pref = "*.byu"; break;
    case SLFlatMap:    pref = "*.p";   break;
    default:           pref = "*.*";   break;
  }

  if(browser->GetFilename("Select Surface",rstring,returnLineLength,pref,NULL,multipleSelection)) {
	sprintf(filename,"%s",rstring);
	c = 0;
	Set(XmNvalue,filename);
	XtSetValues(file_text,arglist,c);
  }
VOUT("void SurfaceLoader::FileSelectCB()")
}

static void LoadWorkProc(char *str, void *dat)
{
VIN("static void LoadWorkProc(char *str, void *dat)")
  if(dat) ((SurfaceLoader*)dat)->WorkProc(str);
VOUT("static void LoadWorkProc(char *str, void *dat)")
}


void SurfaceLoader::WorkProc(char *str)
{
VIN("void SurfaceLoader::WorkProc(char *str)")
  if(str) shell->Working(str);
  else	shell->EndWorking();
VOUT("void SurfaceLoader::WorkProc(char *str)")
}


SurfaceLoader::SurfaceLoader(FileBrowser *brows)
{
VIN("SurfaceLoader::SurfaceLoader(FileBrowser *brows)")
  int w = App->BigFontWidth() * 30;
  int h = (App->BigFontHeight() + 15) * 4;
  int x = App->Width()/2 - w/2;
  int y = App->Height()/2 - h/2;
  u_long flgs =  KDSHELL_WORKICON | KDSHELL_INFOBAR | KDSHELL_TITLEBAR | KDSHELL_NOCLOSE;

  filename[0] = 0;
  multipleSelection = false;

  shell = new KDShell("Surface Load/Save",
		x,y,w,h,flgs,
		LOAD_IMAGE_ICON,
		MainWindow->Shell(),
		BACKCOLOR);

  c = 0;
  Set(XmNminWidth,w);
  Set(XmNminHeight,h);
  XtSetValues(shell->Shell(),arglist,c);

  c = 0;
  Set(XmNverticalSpacing,10);
  Set(XmNhorizontalSpacing,10);
  XtSetValues(shell->Main(),arglist,c);


  c = 0;
  Set(XmNleftAttachment,XmATTACH_FORM);
  Set(XmNtopAttachment,XmATTACH_FORM);
  Set(XmNwidth,w/3);
  Set(XmNresizeWidth,False);
  Widget filebutton = App->Button(shell->Main(),"File",
		App->MedFontList(),
		SurfaceLoaderFileSelectCB,(XtPointer)this,
		arglist,c);		

  c = 0;
  Set(XmNtopAttachment,XmATTACH_FORM);
  Set(XmNrightAttachment,XmATTACH_FORM);
  Set(XmNleftAttachment,XmATTACH_WIDGET);
  Set(XmNleftWidget,filebutton);
  Set(XmNborderWidth,0);
  file_text = App->Text(shell->Main(),filename,
		MAXPATHLEN,
		SurfaceLoaderOkayCB,(XtPointer)this,
		App->SmallFontList(),arglist,c);


  //
  // Create the 'Surface Type' menu
  //
  c = 0;
  Widget pulldown = XmCreatePulldownMenu(shell->Main(),"options",
		arglist,c);
  XmChangeColor(pulldown,MENUCOLOR);

  Widget btn;
  int nitems = XtNumber(surface_types);
  c = 0;
  Set(XmNfontList,App->MedFontList());
  Set(XmNleftAttachment,XmATTACH_FORM);
  Set(XmNrightAttachment,XmATTACH_FORM);
  for(int i=0;i<nitems;i++) {
	surface_types[i].widget = btn = XmCreatePushButtonGadget(pulldown,
			surface_types[i].name,
			arglist,c);
	//XtAddCallback(btn,XmNactivateCallback,TypeSelectCB,(XtPointer)this);
	XtManageChild(btn);
	}

  XmString label = XmStringCreate(" Surface Type                    ",App->MedFontTag());

  c = 0;
  Set(XmNleftAttachment,XmATTACH_FORM);
  Set(XmNrightAttachment,XmATTACH_FORM);
  Set(XmNtopAttachment,XmATTACH_WIDGET);
  Set(XmNtopWidget,filebutton);
  Set(XmNfontList,App->MedFontList());
  Set(XmNsubMenuId,pulldown);
  Set(XmNlabelString,label);
  surftype = XmCreateOptionMenu(shell->Main(),
			"surf_type_options",arglist,c);
  XtManageChild(surftype);
  XmChangeColor(surftype,MENUCOLOR);

  Widget lbl = XmOptionLabelGadget(surftype);
  Widget pdn = XmOptionButtonGadget(surftype);
  XmChangeColor(lbl,MENUCOLOR);
  XmChangeColor(pdn,MENUCOLOR);

  c = 0;
  Set(XmNleftAttachment,XmATTACH_FORM);
  Set(XmNrightAttachment,XmATTACH_SELF);
  Set(XmNfontList,App->MedFontList());
  XtSetValues(lbl,arglist,c);

  c = 0;
  Set(XmNrightAttachment,XmATTACH_FORM);
  Set(XmNleftAttachment,XmATTACH_SELF);
  XtSetValues(pdn,arglist,c);


  //--------------------------------------------------
  //
  // Create Load/Save and CANCEL buttons
  //
  //--------------------------------------------------
  Dimension bx,bw;
  bx = 50;
  bw = w/4;


  c = 0;
  Set(XmNx,bx);
  Set(XmNleftAttachment,XmATTACH_SELF);
  Set(XmNrightAttachment,XmATTACH_SELF);
  Set(XmNbottomAttachment,XmATTACH_FORM);
  Set(XmNtopAttachment,XmATTACH_WIDGET);
  Set(XmNtopWidget,surftype);
  Set(XmNwidth,bw);
  Set(XmNheight,App->MedFontHeight() + 5);
  Set(XmNresizeWidth,False);
  Set(XmNbottomOffset,10);
  okaybtn = App->Button(shell->Main(),"Save",
		App->MedFontList(),
		SurfaceLoaderOkayCB,(XtPointer)this,arglist,c);

  bx = w - 150;

  c = 0;
  Set(XmNx,bx);
  Set(XmNleftAttachment,XmATTACH_SELF);
  Set(XmNrightAttachment,XmATTACH_SELF);
  Set(XmNbottomAttachment,XmATTACH_FORM);
  Set(XmNtopAttachment,XmATTACH_WIDGET);
  Set(XmNtopWidget,surftype);
  Set(XmNwidth,bw);
  Set(XmNheight,App->MedFontHeight() + 5);
  Set(XmNresizeWidth,False);
  Set(XmNbottomOffset,10);
  Widget cancelbtn = App->Button(shell->Main(),"CANCEL",
		App->MedFontList(),
		SurfaceLoaderCancelCB,(XtPointer)this,arglist,c);
  XmChangeColor(cancelbtn,CANCELCOLOR);

  if(brows) {
	browser = brows;
	own_browser = 0;
	}
  else	{
	browser = new FileBrowser();
	own_browser = 1;
	}
VOUT("SurfaceLoader::SurfaceLoader(FileBrowser *brows)")
}


SurfaceLoader::~SurfaceLoader()
{
VIN("SurfaceLoader::~SurfaceLoader()")
  delete shell;
  if(own_browser) delete browser;
VOUT("SurfaceLoader::~SurfaceLoader()")
}


SurfaceLoaderType SurfaceLoader::getSurfaceType()
{
VIN("SurfaceLoaderType SurfaceLoader::getSurfaceType()")
  Widget hwid;
  int i,ni;
  SurfaceLoaderType ctype = SLByuSurface;

  c = 0;
  Set(XmNmenuHistory,&hwid);
  XtGetValues(surftype,arglist,c);

  ni = XtNumber(surface_types);
  for(i=0;i<ni;i++)
	if(surface_types[i].widget == hwid) 
		ctype = surface_types[i].type;

VOUT("SurfaceLoaderType SurfaceLoader::getSurfaceType()")
  return(ctype);
}


void SurfaceLoader::Save(LoadedSurface *loaded_surf)
{
VIN("void SurfaceLoader::Save(LoadedSurface *loaded_surf)")
  if(!loaded_surf)
	return;

  multipleSelection = false;
  filename[0] = 0;
  c = 0;
  Set(XmNvalue,"");
  XtSetValues(file_text,arglist,c);

  App->SetLabelText(okaybtn,"Save",App->MedFontTag());

  c = 0;
  Set(XmNtitle,"Save Surface");
  XtSetValues(shell->Shell(),arglist,c);

  shell->Show();

popup_ssaver:

  okay_cancel_clicked = 0;

  while(!okay_cancel_clicked) App->ProcessEvents();

  if(okay_cancel_clicked < 0)  {
	shell->Hide();
	return;
	}


  SurfaceLoaderType current_type = getSurfaceType();

  bool rval;

  char *rstring = XmTextGetString(file_text);
  sprintf(filename,"%s",rstring);
  XtFree(rstring);

  ByuSurface bsurf;

  switch(current_type) {
	case SLByuSurface:
		bsurf = *(loaded_surf->surface());
		bsurf.SetWorkProc(LoadWorkProc,(void*)this);
		rval = bsurf.Save(filename);
		bsurf.SetWorkProc(NULL,NULL);
		if(rval == False) {
			App->ShowMessage("ERROR: BYU surface did not save");
        		goto popup_ssaver;
			}
		break;
	case SLFlatMap:
		App->ShowMessage("SurfaceLoader::save flat map not coded yet");
        	goto popup_ssaver;
		break;
	default:
		App->ShowMessage("INTERNAL ERROR: invalid surface type");
        	goto popup_ssaver;
		break;
	}
  Icon *ic = Icon::getIconFromData((void*)loaded_surf);
  if(ic) ic->SetTitle(App->FilenameTail(filename));

  loaded_surf->surface()->SetFilename(filename);

  shell->EndWorking();
  shell->Hide();
VOUT("void SurfaceLoader::Save(LoadedSurface *loaded_surf)")
}


LoadedSurface **SurfaceLoader::Load(bool multiselect)
{
VIN("LoadedSurface *SurfaceLoader::Load()")
  multipleSelection = multiselect;
  filename[0] = 0;
  c = 0;
  Set(XmNvalue,filename);
  XtSetValues(file_text,arglist,c);

  App->SetLabelText(okaybtn,"Load",App->MedFontTag());

  c = 0;
  Set(XmNtitle,"Load Surface");
  XtSetValues(shell->Shell(),arglist,c);

  shell->Show();

popup_loader:

  okay_cancel_clicked = 0;

  while(!okay_cancel_clicked) App->ProcessEvents();

  if(okay_cancel_clicked < 0 || strlen(filename) == 0)  {
	shell->Hide();
	return((LoadedSurface**)NULL);
	}

  ByuSurface *bsurf;
  LoadedSurface **retsurf = NULL;
  SurfaceLoaderType current_type = getSurfaceType();

  if(current_type == SLByuSurface) {
  
    // determine how many file names are in the string
    // replace the ';' with a string terminator for easier processing
    int fileCount = 1;
    char *nextFile = filename;
    while (nextFile = index(nextFile, (int) ';'))
    {
        *nextFile = '\0';
        nextFile++;
        fileCount++;
    }
    
    retsurf = new LoadedSurface*[fileCount + 1];
    retsurf[fileCount] = NULL;
    nextFile = filename;
    char nextFileName[MAXPATHLEN];
    for (int fileIx = 0;
         fileIx < fileCount;
         fileIx++)
    {
        bsurf = new ByuSurface();
        bsurf->SetWorkProc(LoadWorkProc,(void*)this);
        
        if(!bsurf->Load(nextFile))
        {
            bsurf->SetWorkProc(NULL,NULL);
            delete bsurf;
            retsurf[fileIx] = NULL;
            sprintf(errstr,"ERROR: Could not load `%s'",nextFile);
            App->ShowMessage(errstr);
            break;
            // goto popup_loader;
        }

        bsurf->SetWorkProc(NULL,NULL);
        if((bsurf->getNumVert() <= 0)||(bsurf->getNumPoly() <= 0))
        {
            delete bsurf;
                sprintf(errstr,"ERROR: `%s' is not a valid surface",filename);
                App->ShowMessage(errstr);
            break;
            // goto popup_loader;
        }

        retsurf[fileIx] = new LoadedSurface(bsurf);
        retsurf[fileIx]->glsurface()->setFlat(checkFlat(bsurf));
        nextFile = &nextFile[strlen(nextFile) + 1];
    }
  } else if(current_type == SLFlatMap) {
    retsurf = new LoadedSurface*[2];
	retsurf[0] = LoadFlatMap(filename);
    retsurf[1] = NULL;
	if(!retsurf[0]) {
        	App->ShowMessage(errstr);
		goto popup_loader;
	}
	retsurf[0]->glsurface()->setFlat(true);
  }

  shell->EndWorking();
  shell->Hide();

VOUT("LoadedSurface *SurfaceLoader::Load()")
  return(retsurf);
}

bool SurfaceLoader::checkFlat(Surface *S)
{
VIN("bool SurfaceLoader::checkFlat(Surface *S)")
  int i,np;
  if(!S || ((np=S->getNumVert())<1))
     return(false);

  float mx,my,mz;
  float xx,xy,xz;
  float dx,dy,dz;
  float x,y,z;

  mx = my = mz =  1E20;
  xx = xy = xz = -1E20;
  for(i=0;i<np;i++) {
     x = S->vertices()[i].x();
     y = S->vertices()[i].y();
     z = S->vertices()[i].z();
     if(x<mx) mx = x; if(x>xx) xx = x;
     if(y<my) my = y; if(y>xy) xy = y;
     if(z<mz) mz = z; if(z>xz) xz = z;
  }
  dx = (xx - mx);
  dy = (xy - my);
  dz = (xz - mz);

VOUT("bool SurfaceLoader::checkFlat(Surface *S)")
  return((dx==0.)||(dy==0.)||(dz==0.));
}


#define FLATEXIT(a) { sprintf(errstr,"%s",a); goto flat_fail; }
#define GREY_RGB	0.5f,0.5f,0.5f
#define DK_GREY_RGB	0.2f,0.2f,0.2f
#define DK_GREEN_RGB	0.25f,0.54f,0.25f
#define RED_RGB		1.f,0.f,0.f
#define YELLOW_RGB	1.f,1.f,0.f
#define BLUE_RGB	0.f,0.f,1.f
#define PINK_RGB	1.f,0.f,1.f
#define CYAN_RGB	0.f,1.f,1.f
#define BR_GREEN_RGB	0.f,1.f,0.f
#define ORANGE_RGB	0.91f,0.54f,0.f
#define BROWN_RGB	0.63f,0.31f,0.f
#define PURPLE_RGB	0.69f,0.f,1.f
#define ORANGE_RAMP1_RGB 0.96f,0.37f,0.f
#define ORANGE_RAMP2_RGB 1.f,0.37f,0.11f
#define ORANGE_RAMP3_RGB 1.f,0.62,0.22f
#define ORANGE_RAMP4_RGB 1.f,0.67f,0.33f
#define ORANGE_RAMP5_RGB 1.f,0.73,0.44f
#define ORANGE_RAMP6_RGB 1.f,0.8,0.58f
#define ORANGE_RAMP7_RGB 1.f,0.85,0.69f



LoadedSurface *SurfaceLoader::LoadFlatMap(const char *filename)
{
VIN("LoadedSurface *SurfaceLoader::LoadFlatMap(const char *filename)")

    int i;
    
    FlatMap::flatvertex *vertices = NULL;
    ByuSurface *bsurf = FlatMap::LoadFlatMap(filename, &vertices, errstr);
    GLSurface  *gsurf = NULL;
    LoadedSurface *lsurf = NULL;
    
    if (!bsurf) FLATEXIT("ByuSurface Not Created")
    
    lsurf = new LoadedSurface(bsurf);
    if (!lsurf) FLATEXIT("out of memory")
    
    gsurf = lsurf->glsurface();

    if (gsurf)
        for (i = 0; i < bsurf->getNumVert(); i++)
        {
            switch(vertices[i].color)
            {
                case 1:   gsurf->setVertexColor(i,BLUE_RGB);         break;
                case 199: gsurf->setVertexColor(i,RED_RGB);          break;
                case 209: gsurf->setVertexColor(i,YELLOW_RGB);       break;
                case 211: gsurf->setVertexColor(i,ORANGE_RAMP1_RGB); break;
                case 212: gsurf->setVertexColor(i,ORANGE_RAMP2_RGB); break;
                case 213: gsurf->setVertexColor(i,ORANGE_RAMP3_RGB); break;
                case 214: gsurf->setVertexColor(i,ORANGE_RAMP4_RGB); break;
                case 215: gsurf->setVertexColor(i,ORANGE_RAMP5_RGB); break;
                case 216: gsurf->setVertexColor(i,ORANGE_RAMP6_RGB); break;
                case 217: gsurf->setVertexColor(i,ORANGE_RAMP7_RGB); break;
                case 218: gsurf->setVertexColor(i,BR_GREEN_RGB);     break;
                case 221: gsurf->setVertexColor(i,CYAN_RGB);         break;
                case 231: gsurf->setVertexColor(i,DK_GREY_RGB);      break;
                case 232: gsurf->setVertexColor(i,GREY_RGB);         break;
                case 233: gsurf->setVertexColor(i,DK_GREEN_RGB);     break;
                case 235: gsurf->setVertexColor(i,PINK_RGB);         break;
                case 236: gsurf->setVertexColor(i,ORANGE_RGB);       break;
                case 237: gsurf->setVertexColor(i,BROWN_RGB);        break;
                case 239: gsurf->setVertexColor(i,PURPLE_RGB);       break;
                default:
                        gsurf->setVertexColor(i,GREY_RGB); break;
            }
        }

  // Now, make sure flat map is oriented correctly for viewing

VOUT("LoadedSurface *SurfaceLoader::LoadFlatMap(const char *filename)")
  delete [] vertices;

  return (lsurf);

flat_fail:
    
    delete [] vertices;
    
    if (lsurf) delete lsurf;
    else if (bsurf) delete bsurf;

VOUT("LoadedSurface *SurfaceLoader::LoadFlatMap(const char *filename)")
  return((LoadedSurface*)NULL);
}


