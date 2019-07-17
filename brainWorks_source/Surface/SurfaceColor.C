///////////////////////////////////////////////////////////////////////////
//
// File: SurfaceColor.C
//
// Author: Keith Doolittle
//
// Purpose:
//
///////////////////////////////////////////////////////////////////////////
#include <iostream.h>
#include <KDApplication.h>
#include <KDShell.h>
#include <Xm/Form.h>
#include <Main.h>
#include <ColorLUT.h>
#include <HistLUT.h>
#include <SurfaceColor.h>
#include <IconList.h>
#include <SurfaceBox.h>
#include <GLSurface.h>


static void OKCB(Widget,XtPointer cld, XtPointer)
{
VIN("static void OKCB(Widget,XtPointer cld, XtPointer)")
if(cld) ((SurfaceColor*)cld)->setCB(1);
VOUT("static void OKCB(Widget,XtPointer cld, XtPointer)")
}

static void CancelCB(Widget,XtPointer cld, XtPointer)
{
VIN("static void CancelCB(Widget,XtPointer cld, XtPointer)")
if(cld) ((SurfaceColor*)cld)->setCB(-1);
VOUT("static void CancelCB(Widget,XtPointer cld, XtPointer)")
}

static void ApplyCB(Widget,XtPointer cld, XtPointer)
{
VIN("static void ApplyCB(Widget,XtPointer cld, XtPointer)")
if(cld) ((SurfaceColor*)cld)->setCB(2);
VOUT("static void ApplyCB(Widget,XtPointer cld, XtPointer)")
}

void SurfaceColor::applyFunction()
{
VIN("void SurfaceColor::applyFunction()")
char *ttl = "<unknown>";
if(mode == ColorNone)
	ttl = "<no function defined>";
else if(mode == ColorFunction)
	ttl = fcnname;

hist->setTitle(ttl);
//hist->setData(&values);
VOUT("void SurfaceColor::applyFunction()")
}

void SurfaceColor::colorSurface(LoadedSurface *ldSurface, ifstream &ifile)
{
	VIN("void ColorSurface(LoadedSurface *ldSurface, ifstream *ifile)")

	GLSurface *glsurf;
	Surface	  *surf;
	if (!ldSurface || (!(glsurf = ldSurface->glsurface())) ||
		(!(surf = ldSurface->surface())) || !ifile)
	{
		return;
	}

	int red,green,blue;
	int vertex;

	lut = new ColorLUT();
	lut->setRed();

	for (int index = 0; index < surf->getNumVert(); ++index)
	{
		ifile >> vertex >> red >> green >> blue;

		glsurf->setVertexColor(index, (float)red/255.0, (float)green/255.0, (float)blue/255.0);
	}

	VOUT("void ColorSurface(LoadedSurface *ldSurface, ifstream *ifile)")
}

void SurfaceColor::colorSurface(GLSurface *gsurf)
{
VIN("void SurfaceColor::colorSurface(GLSurface *gsurf)")
	glsurf = gsurf;
	if(!glsurf) {
		cerr << "ERROR: SurfaceColor::colorSurface() invalid surface" << endl;
		return;
	}

	nvert = glsurf->getNumVert();
	if((nvert <= 0)||(glsurf->getNumPoly() <= 0)) {
		cerr << "ERROR: SurfaceColor::colorSurface() invalid surface" << endl;
		return;
	}

	values.setDim(nvert);
	if(values.isEmpty()) {
		App->ShowMessage("ERROR: out of memory");
		return;
	}

	int w = App->BigFontWidth() * 20;
	int h = App->Height()/2;
	int x = App->Width()/2 - w/2;
	int y = App->Height()/2 - h/2;

	KDShell *shell = new KDShell("Color Surface",x,y,w,h,
								 KDSHELL_TITLEBAR,0,
								 MainWindow->Shell(),MESSAGECOLOR);

	int bh = App->BigFontHeight() + 8;

/*
c = 0;
Set(XmNtopAttachment,XmATTACH_FORM);
Set(XmNbottomAttachment,XmATTACH_NONE);
Set(XmNheight,bh);
Set(XmNleftAttachment,XmATTACH_FORM);
Set(XmNrightAttachment,XmATTACH_FORM);
fcn_list = new IconList(shell->Main(),"Function",IXNoIcon,
	fcn_names,XtNumber(fcn_names),arglist,c);
//fcn_list->AddExtraCallback(FcnChangedCB,(void*)this);
*/

	c = 0;
	Set(XmNtopAttachment,XmATTACH_NONE);
	Set(XmNheight,bh);
	Set(XmNbottomAttachment,XmATTACH_FORM);
	Set(XmNleftAttachment,XmATTACH_FORM);
	Set(XmNrightAttachment,XmATTACH_FORM);
	Widget frm = XmCreateForm(shell->Main(),"surfacecolor_frm",arglist,c);
	XtManageChild(frm);
	XmChangeColor(frm,SHELLCOLOR);

	c = 0;
	Set(XmNtopAttachment,XmATTACH_FORM);
	Set(XmNbottomAttachment,XmATTACH_FORM);
	Set(XmNleftAttachment,XmATTACH_POSITION);
	Set(XmNleftPosition,2);
	Set(XmNrightAttachment,XmATTACH_POSITION);
	Set(XmNrightPosition,31);
	Widget btn = App->Button(frm,"Apply",App->SmallFontList(),
							 ApplyCB,(void*)this, arglist,c);
	XmChangeColor(btn,BUTTON2COLOR);

	c = 0;
	Set(XmNtopAttachment,XmATTACH_FORM);
	Set(XmNbottomAttachment,XmATTACH_FORM);
	Set(XmNleftAttachment,XmATTACH_POSITION);
	Set(XmNleftPosition,33);
	Set(XmNrightAttachment,XmATTACH_POSITION);
	Set(XmNrightPosition,66);
	btn = App->Button(frm,"OK",App->SmallFontList(),
					  OKCB,(void*)this, arglist,c);
	c = 0;
	Set(XmNtopAttachment,XmATTACH_FORM);
	Set(XmNbottomAttachment,XmATTACH_FORM);
	Set(XmNleftAttachment,XmATTACH_POSITION);
	Set(XmNleftPosition,68);
	Set(XmNrightAttachment,XmATTACH_POSITION);
	Set(XmNrightPosition,98);
	btn = App->Button(frm,"CANCEL", App->SmallFontList(),
					  CancelCB,(void*)this, arglist,c);
	XmChangeColor(btn,CANCELCOLOR);

	lut = new ColorLUT();
        lut->setGreyscale();

	c = 0;
	Set(XmNtopAttachment,XmATTACH_WIDGET);
	Set(XmNbottomAttachment,XmATTACH_WIDGET);
	Set(XmNbottomWidget,frm);
	Set(XmNleftAttachment,XmATTACH_FORM);
	Set(XmNrightAttachment,XmATTACH_FORM);

	Array1D<float> red;
	Array1D<float> grn;
	Array1D<float> blu;

	shell->Show();

float dmin = 1E20;
float dmax = -1E20;
float drange;
//
// Set whatever is on surface now
//
	if(glsurf->hasCurvature()) {

		for(int k=0;k<nvert;k++) {
			values[k] = glsurf->curvature()[k];
                        if(values[k] < dmin) dmin = values[k];
                        if(values[k] > dmax) dmax = values[k];
                }

drange = dmax - dmin;
if(drange == 0.f) drange = 1.f;

		mode = ColorFunction;
                Set(XmNtopAttachment,XmATTACH_WIDGET);
//Set(XmNtopWidget,fcn_list->Form());
        Set(XmNbottomAttachment,XmATTACH_WIDGET);
        Set(XmNbottomWidget,frm);
        Set(XmNleftAttachment,XmATTACH_FORM);
        Set(XmNrightAttachment,XmATTACH_FORM);

	hist = new HistLUT(shell->Main(),NULL,lut,arglist,c,&values);
	hist->setTitle("---");
        hist->setForceRange(true,dmin,dmax);
        hist->setData(&values);
	}
	else	{
		values = 0.0;
		mode = ColorNone;
	}

 retry_coloring:

	cb = 0;
	while(!cb) App->ProcessEvents();
	if(cb > 0) {
		if((mode == ColorFunction)) {

			red.setDim(nvert);
			grn.setDim(nvert);
			blu.setDim(nvert);
			if(red.isEmpty()||grn.isEmpty()||blu.isEmpty()) {
				App->ShowMessage("ERROR: out of memory");
				goto endthis;
			}

			u_char *tbl = lut->gltable();

			for(int i=0;i<nvert;i++) {
		//		int   bini = hist->valueToBin(values[i]);
//int bini = (int)( (values[i] - dmin)/drange * 255.f );

int bini = (int)((float)hist->valueToBin(values[i])/(float)hist->numBins() * 255.f);
				if((bini<0)||(bini>255)) {
					red[i] = grn[i] = blu[i] = 0.0;
				}
       else	{
		   red[i] = ((float)tbl[3*bini+0])/255.0;
		   grn[i] = ((float)tbl[3*bini+1])/255.0;
		   blu[i] = ((float)tbl[3*bini+2])/255.0;
       }
     }
			glsurf->setColorSurface(red,grn,blu);

			App->ProcessEvents();
		}
		else if(mode == ColorNone) {
			glsurf->setColorSurface(0.7,0.7,0.9);
   }

   // apply (ok falls through)
		if(cb == 2) goto retry_coloring;
	}

 endthis:
	shell->Hide();
//delete fcn_list;
	delete shell;
VOUT("void SurfaceColor::colorSurface(GLSurface *gsurf)")
}





