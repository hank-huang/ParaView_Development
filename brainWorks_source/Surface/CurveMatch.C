///////////////////////////////////////////////////////////////////////////
//
// File: CurveMatch.C
//
// Author: Muge Bakircioglu
//
// Purpose:
//
///////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <iostream.h>
#include <sys/param.h>
#include <KDApplication.h>
#include <Main.h>
#include <SurfaceBox.h>
#include <IconList.h>
#include <GLPrimitive.h>
#include <CurveMatch.h>
#include <Xm/ScrolledW.h>
#include <Xm/RowColumn.h>
#include <Xm/Text.h>
#include <Xm/Form.h>
#include <TDL/Surface.h>
#include <TDL/ByuSurface.h>
#include <TDL/SurfaceUtils.h>



DEFINECB(CurveMatch,MatchCB,match)
DEFINECBP(CurveMatch,LoadCurve1CB,loadCurve,1)
DEFINECBP(CurveMatch,LoadCurve2CB,loadCurve,2)
DEFINECB(CurveMatch,SaveMatchCB,saveMatch)

// default line color
//
#define DEFRED          0.0
#define DEFGREEN        0.0  
#define DEFBLUE         1.0

#define SELRED          1.0
#define SELGREEN        0.0
#define SELBLUE         0.0

#define MATCHRED 1.0
#define MATCHGREEN 1.0
#define MATCHBLUE 1.0

#define LINE_WIDTH      5.0

void CurveMatch::loadCurve(int curveno) {
    if (_target.numSamples()!=0 && _template.numSamples()!=0) {
     _target.reset(); _template.reset();
     tempPoints.setDim(0); targPoints.setDim(0);
     _minFlag=-1;
     while(lines) lines=GLPrimitive::removeFromList(lines,lines);
     surfwin->setPrimitives(NULL);
     
    App->SetLabelText(_cost_wid,"0.0", App->SmallFontTag());
    App->SetLabelText(_curve1_name_wid,"<none>",
             App->SmallFontTag());
     App->SetLabelText(_curve2_name_wid, "<none>",
               App->SmallFontTag());
    }
    char fname[1024];
 retry_load:
   fname[0] = 0;
   GLPrimitive *newline;
   if(Browser->GetFilename("Select Curve File",fname,1023)) {
     // load curve here filename is 'fname'
     
     if(curveno == 1) { // curve 1
       App->SetLabelText(_curve1_name_wid,
             App->FilenameTail(fname),
             App->SmallFontTag());
     }
     else   {      // curve 2
       App->SetLabelText(_curve2_name_wid,
             App->FilenameTail(fname),
             App->SmallFontTag());
     }  
     
     ifstream fp(fname,ios::in);
     char readStatusMsg[1024];
     Array1D<Point> P;
     if(!readFile(curveno, fname, P, readStatusMsg)) {
       App->ShowMessage(readStatusMsg);
       goto retry_load;
     }
   
     setCurvesApart(P,curveno);
     newline= new GLPrimitive(P);
     newline->setLineWidth(LINE_WIDTH);
     if (curveno==1) 
       newline->setColor(DEFRED,DEFGREEN,DEFBLUE);
     else
       newline->setColor(SELRED,SELGREEN,SELBLUE);
     lines= GLPrimitive::addToList(lines,newline);
     surfwin->setPrimitives(lines);
     surfwin->surfaceChanged();
     surfwin->ViewAll();
     surfwin->redraw();
    
   }
 }
 
void CurveMatch::match()
{
VIN("void CurveMatch::match()")
  //use the original curves stored
  _target=origTarget;
  _template=origTemplate;
  newTarg.reset();
  int nsamples = App->GetInt(_nsamples_wid);
  int nLmks=App->GetInt(_nlandmarks_wid);
  int polyOrder=App->GetInt(_orderPoly_wid);
  double a_coeff=App->GetDouble(_arclengcoef_wid);
  double c_coeff=App->GetDouble(_curvcoef_wid);
  double t_coeff=App->GetDouble(_torscoef_wid);
  int nghbdsize=App->GetInt(_nghbdSize_wid);

  // Curve newTarg;
 
  Array1D<Point> storeMatch(2);
  int i;

  char statusMsg[1024];
  if (!checkMatchingParameters(nsamples, nLmks, polyOrder, nghbdsize, statusMsg))
  {
      App->ShowMessage(statusMsg);
      return;
  }

  _template.calculateMatchcost(_target,nsamples,nLmks,polyOrder,
			       a_coeff,c_coeff,t_coeff,nghbdsize,
			       newTarg,numNghbs,spacing,dist);
  char temp[100];
  sprintf(temp,"%10.3f", _template.matching(newTarg,numNghbs,spacing,dist));
  char temp1[100];
  strcpy(temp1,"Cost of current match ");
  strcat(temp1,temp);
  App->SetLabelText(_cost_wid,temp1, App->SmallFontTag());
  tempPoints.setDim(newTarg.numSamples());
  targPoints.setDim(newTarg.numSamples());

  _template.getMatch(newTarg,tempPoints,targPoints);

  while(lines) lines=GLPrimitive::removeFromList(lines,lines);
  surfwin->setPrimitives(NULL);

  Array1D<Point> P(nsamples);
  Array1D<Point> Q(nsamples);

  for (i=0; i<nsamples; ++i) {
    P[i]=_target.getPoint(i);
    Q[i]=_template.getPoint(i);
  }
  GLPrimitive *newline= new GLPrimitive(P);
  newline->setLineWidth(LINE_WIDTH);
  newline->setColor(DEFRED,DEFGREEN,DEFBLUE);
  lines= GLPrimitive::addToList(lines,newline);
  setCurvesApart(Q,2); 
  newline=new GLPrimitive(Q);
  newline->setLineWidth(LINE_WIDTH);
  newline->setColor(SELRED,SELGREEN,SELBLUE);
  lines= GLPrimitive::addToList(lines,newline);

  for (i=0; i<newTarg.numSamples(); ++i) {
    storeMatch[0]=tempPoints[i];
    storeMatch[1]=targPoints[i];
  
    switch(_minFlag) {
    case 0: storeMatch[0]+=Point(_offset,0.0,0.0); break;
    case 1: storeMatch[0]+=Point(0.0, _offset, 0.0); break;
    case 2: storeMatch[0]+=Point(0.0, 0.0, _offset); break;
    default: break;
    }
    newline = new GLPrimitive(storeMatch);
    lines = GLPrimitive::addToList(lines,newline);
  } 

  surfwin->setPrimitives(lines);
  surfwin->surfaceChanged();
  surfwin->ViewAll();
  surfwin->redraw();
VOUT("void CurveMatch::match()")
}
/**************************************************************/

void CurveMatch::saveMatch() 
{
VIN("void CurveMatch::saveMatch() ")
  char fname[1024];
retry_save:
  fname[0]=0;

  Array1D<Point> P(_template.numSamples());
  for (int i=0; i<_template.numSamples();   ++i) 
    P[i]=_template.getPoint(i);

   switch(_minFlag) {
   case 0: 
     //KWD P-=Point(_offset, 0.0,0.0); break;
     addToArray(P,Point(-_offset,0.0,0.0));
     break;
   case 1: 
     //KWD P-=Point(0.0, _offset, 0.0); break;
     addToArray(P,Point(0.0,-_offset,0.0));
     break;
   case 2: 
     //KWD P-=Point(0.0, 0.0, _offset); break;
     addToArray(P,Point(0.0,0.0,-_offset));
     break;
   default: break;
   }
	    
  if (tempPoints.getSize()!=0) {
    if (Browser->GetFilename("Save Match To", fname,1024)) {
      ofstream fp(fname,ios::out);
      if (!fp) {
	App->ShowMessage("ERROR: can't open file");
	goto retry_save;
      }
    
      for (int i=0; i<newTarg.numSamples(); ++i) 
	fp<<targPoints[i]<<"\t"<<tempPoints[i]<<endl;
      fp.close();
    }
  }
  else App->ShowMessage("ERROR: There is no match to save");
VOUT("void CurveMatch::saveMatch() ")
}
/*****************************************************************/

CurveMatch::CurveMatch() : KDChildWindow("Curve Matching",50,50)
{
VIN("CurveMatch::CurveMatch()")
lines = NULL;

c = 0;
Set(XmNfractionBase,100);
XtSetValues(this->Form(),arglist,c);
XmChangeColor(this->Form(),MODULECOLOR);

c = 0;
Set(XmNleftAttachment,XmATTACH_FORM);
Set(XmNtopAttachment,XmATTACH_FORM);
Set(XmNbottomAttachment,XmATTACH_FORM);
Set(XmNrightAttachment,XmATTACH_POSITION);
Set(XmNrightPosition,72);
surfwin = new SurfaceBox(this->Form(),NULL,arglist,c);

Dimension lh  = App->BigFontHeight() + 3;


c = 0;
Set(XmNleftAttachment,XmATTACH_POSITION);
Set(XmNleftPosition,72);
Set(XmNrightAttachment,XmATTACH_FORM);
Set(XmNtopAttachment,XmATTACH_FORM);
Set(XmNbottomAttachment,XmATTACH_NONE);
Set(XmNheight,lh);
Widget btn = App->Button(this->Form(),"Load Curve 1(target)",
		App->SmallFontList(),LoadCurve1CB,(XtPointer)this,
		arglist,c);
XmChangeColor(btn,RED);

c = 0;
Set(XmNleftAttachment,XmATTACH_POSITION);
Set(XmNleftPosition,72);
Set(XmNrightAttachment,XmATTACH_FORM);
Set(XmNtopAttachment,XmATTACH_WIDGET);
Set(XmNtopWidget,btn);
Set(XmNbottomAttachment,XmATTACH_NONE);
_curve1_name_wid = App->Label(this->Form(),"<none>",
		App->SmallFontList(),arglist,c);


c = 0;
Set(XmNleftAttachment,XmATTACH_POSITION);
Set(XmNleftPosition,72);
Set(XmNrightAttachment,XmATTACH_FORM);
Set(XmNtopAttachment,XmATTACH_WIDGET);
Set(XmNtopWidget,_curve1_name_wid);
Set(XmNbottomAttachment,XmATTACH_NONE);
Set(XmNheight,lh);
Widget  btn1 = App->Button(this->Form(),"Load Curve 2(template)",
		App->SmallFontList(),LoadCurve2CB,(XtPointer)this,
		arglist,c);

XmChangeColor(btn1,BLUE);
c = 0;
Set(XmNleftAttachment,XmATTACH_POSITION);
Set(XmNleftPosition,72);
Set(XmNrightAttachment,XmATTACH_FORM);
Set(XmNtopAttachment,XmATTACH_WIDGET);
Set(XmNtopWidget,btn1);
Set(XmNbottomAttachment,XmATTACH_NONE);
_curve2_name_wid = App->Label(this->Form(),"<none>",
		App->SmallFontList(),arglist,c);



c = 0;
Set(XmNleftAttachment,XmATTACH_POSITION);
Set(XmNleftPosition,72);
Set(XmNrightAttachment,XmATTACH_FORM);
Set(XmNtopAttachment,XmATTACH_WIDGET);
Set(XmNtopWidget,_curve2_name_wid);
Set(XmNbottomAttachment,XmATTACH_NONE);
Set(XmNheight,lh);
_nsamples_wid = App->LabelText(this->Form(),"Number of Samples","100",10,
		App->MedFontList(),arglist,c);

c = 0;

Set(XmNleftAttachment,XmATTACH_POSITION);
Set(XmNleftPosition,72);
Set(XmNrightAttachment,XmATTACH_FORM);
Set(XmNtopAttachment,XmATTACH_WIDGET);
Set(XmNtopWidget,_nsamples_wid);
Set(XmNbottomAttachment,XmATTACH_NONE);
Set(XmNheight,lh);
_nlandmarks_wid = App->LabelText(this->Form(),"Number of Landmarks","12",10,
		App->MedFontList(),arglist,c);

c = 0;
Set(XmNleftAttachment,XmATTACH_POSITION);
Set(XmNleftPosition,72);
Set(XmNrightAttachment,XmATTACH_FORM);
Set(XmNtopAttachment,XmATTACH_WIDGET);
Set(XmNtopWidget,_nlandmarks_wid);
Set(XmNbottomAttachment,XmATTACH_NONE);
Set(XmNheight,lh);
_orderPoly_wid = App->LabelText(this->Form(),"Order of Poly","8",10,
		App->MedFontList(),arglist,c);

c = 0;
Set(XmNleftAttachment,XmATTACH_POSITION);
Set(XmNleftPosition,72);
Set(XmNrightAttachment,XmATTACH_FORM);
Set(XmNtopAttachment,XmATTACH_WIDGET);
Set(XmNtopWidget,_orderPoly_wid);
Set(XmNbottomAttachment,XmATTACH_NONE);
Set(XmNheight,lh);
_nghbdSize_wid = App->LabelText(this->Form(),"Nghbd Size (% #Samp)","18",10,App->MedFontList(),arglist,c);



c = 0;
Set(XmNleftAttachment,XmATTACH_POSITION);
Set(XmNleftPosition,72);
Set(XmNrightAttachment,XmATTACH_FORM);
Set(XmNtopAttachment,XmATTACH_WIDGET); 
Set(XmNtopWidget,_nghbdSize_wid);
Set(XmNbottomAttachment,XmATTACH_NONE);
Set(XmNheight,lh);
_arclengcoef_wid = App->LabelText(this->Form(),"Arclength Coeff","1",10,  App->MedFontList(),arglist,c);

c = 0;
Set(XmNleftAttachment,XmATTACH_POSITION);
Set(XmNleftPosition,72);
Set(XmNrightAttachment,XmATTACH_FORM);
Set(XmNtopAttachment,XmATTACH_WIDGET);
Set(XmNtopWidget,_arclengcoef_wid);
Set(XmNbottomAttachment,XmATTACH_NONE);
Set(XmNheight,lh);
_curvcoef_wid = App->LabelText(this->Form(),"Curvature Coeff","1",10,  App->MedFontList(),arglist,c);

c = 0;
Set(XmNleftAttachment,XmATTACH_POSITION);
Set(XmNleftPosition,72);
Set(XmNrightAttachment,XmATTACH_FORM);
Set(XmNtopAttachment,XmATTACH_WIDGET);   
Set(XmNtopWidget,_curvcoef_wid);
Set(XmNbottomAttachment,XmATTACH_NONE);
Set(XmNheight,lh);
_torscoef_wid = App->LabelText(this->Form(),"Torsion Coeff","1",10, 
App->MedFontList(),arglist,c);

c = 0;
Set(XmNleftAttachment,XmATTACH_POSITION);
Set(XmNleftPosition,72);
Set(XmNrightAttachment,XmATTACH_FORM);
Set(XmNtopAttachment,XmATTACH_NONE);
Set(XmNheight,lh);
Set(XmNbottomAttachment,XmATTACH_FORM);
Widget btn3 = App->Button(this->Form(),"Match",App->MedFontList(),
		MatchCB,(XtPointer)this,arglist,c);

c = 0;
Set(XmNleftAttachment,XmATTACH_POSITION);
Set(XmNleftPosition,72);
Set(XmNrightAttachment,XmATTACH_FORM);
Set(XmNtopAttachment,XmATTACH_NONE);
Set(XmNbottomAttachment,XmATTACH_WIDGET);
Set(XmNbottomWidget,btn3)
Set(XmNheight,lh);
Widget btn2 = App->Button(this->Form(),"Save matches",
		App->MedFontList(),SaveMatchCB,(XtPointer)this,
		arglist,c);


c = 0;
Set(XmNleftAttachment,XmATTACH_POSITION);
Set(XmNleftPosition,72);
Set(XmNrightAttachment,XmATTACH_FORM);
Set(XmNbottomAttachment,XmATTACH_WIDGET);
Set(XmNbottomWidget,btn2);
Set(XmNtopAttachment,XmATTACH_NONE);
_cost_wid =  App->Label(this->Form(),"Cost of current match 0.0",App->MedFontList(),arglist,c);

SetSize(App->Width(),App->Height());
Show();
VOUT("CurveMatch::CurveMatch()")
}

CurveMatch::~CurveMatch()
{
VIN("CurveMatch::~CurveMatch()")
while(lines) lines = GLPrimitive::removeFromList(lines,lines);
delete surfwin;
VOUT("CurveMatch::~CurveMatch()")
}

