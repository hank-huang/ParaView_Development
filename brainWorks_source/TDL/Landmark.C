#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <TDL/Point.h>
#include <TDL/Landmark.h>

#define DALLOC	10

#define NEW_ANALYZE_HEADER_STR	"#   X     Y     Z       Val"
#define WARP_HEADER_FORMAT_STR	"Landmarks-%s"
#define WARP_HEADER_STR  "Landmarks-1.0"

Landmark::Landmark(int nx, int ny, int nz, double dx, double dy, double dz)
{
sprintf(filename,"");
sprintf(warpVersion,"");
points   = NULL;
names    = NULL;
valid    = NULL;
n        = 0;
nalloced = 0;

//
// check argument ranges for validity
//
if(dx <= 0.0) dx = 1.0;
if(dy <= 0.0) dy = 1.0;
if(dz <= 0.0) dz = 1.0;
if(nx < 1)    nx = 1;
if(ny < 1)    ny = 1;
if(nz < 1)    nz = 1;

pixeldx  = dx;
pixeldy  = dy;
pixeldz  = dz;
imagenx  = nx;
imageny  = ny;
imagenz  = nz;
}


void Landmark::clear()
{
for(int i = 0; i < n; i++)
    delete names[i];
delete [] points;
delete [] names;
delete [] valid;

points   = NULL;
names    = NULL;
valid    = NULL;
n        = 0;
nalloced = 0;
}


ItXECode Landmark::remove(int i)
{
if((i<0)||(i >= n)) 
	return(ItXError);

for(;i<n-1;i++) {
	memcpy((void*)&(points[i]),
		(const void*)&(points[i+1]),
		sizeof(Point));
	names[i] = names[i+1];
	valid[i] = valid[i+1];
	}
n--;
return(ItXSuccess);
}


Landmark::~Landmark()
{
clear();
}


/* KWD
int Landmark::getLandmarkIndex(char *search_name, int case_matters)
{
int i,j;

for(j=0;j<n;j++) {
	if(case_matters) 
		i = strcmp(search_name,names[j]);
	else    i = strcasecmp(search_name,names[j]);
	if(i == 0) return(j);
	}
return(-1);
}


int Landmark::getLandmarkIndex(double x, double y, double z, double r)
{
int i;
double dx,dy,dz,tr;

for(i=0;i<n;i++) {
	if(r == 0.0) {
		if((points[i].x() == x)&&(points[i].y() == y)&&
		   (points[i].z() == z))
			return(i);
		}
	else	{
		dx = points[i].x() - x;
		dy = points[i].y() - y;
		dz = points[i].z() - z;
		tr = sqrt(dx*dx + dy*dy + dz*dz);
		if(tr <= r) return(i);
		}
	}
return(-1);
}
*/


int Landmark::add(double x, double y, double z, const char *nm, int isvalid)
{
int onalloced = nalloced;
if(nalloced <= n) {
	nalloced += DALLOC;
	Point *npoints = new Point[nalloced];
	if(points && (n>0)) {
	    	memcpy((void*)npoints,(const void*)points,
			sizeof(Point)*n);
	    	delete [] points;
		}
	points = npoints;


	int *nvalid = new int[nalloced];
	if(valid && (n>0)) {
	    	memcpy((void*)nvalid,(const void*)valid,
			sizeof(int)*n);
	    	delete [] valid;
		}
	valid = nvalid;

	char **nnames = new char*[nalloced];
	if(names && (onalloced>0)) {
	    	memcpy((void*)nnames,(const void*)names,
			sizeof(char*)*onalloced);
	    	delete [] names;
		}

//keith check this here
	names = nnames;
	for(int i=onalloced;i<nalloced;i++) {
		names[i] = 0;
		}
	}

points[n].set(x,y,z);
valid[n] = isvalid;
int l = strlen(nm);
if(l <= 0) {
	names[n] = new char[3];
 	sprintf(names[n],"<>");
	}
else	{
    // truncate the name to our max len
    char tmpLandmarkName[MAXLANDMARKNAMELEN];
    strncpy(tmpLandmarkName, nm, MAXLANDMARKNAMELEN);
    tmpLandmarkName[MAXLANDMARKNAMELEN] = '\0';
    
    // remove carriage return and anything after...
    if (char* carriageReturn = strchr(tmpLandmarkName, '\n'))
        *carriageReturn = '\0';
    
    // allocate and copy...
	names[n] = new char[strlen(tmpLandmarkName)+1];
 	sprintf(names[n],"%s",tmpLandmarkName);
	}

n++;
return(n-1);
}


#define MAXLOADEDLINE 512

ItXECode Landmark::loadAnalyze(FILE *fp)
{
ItXECode rv = loadStandard(fp);
int i;
double fx,fy,fz;
char   sstring[MAXLOADEDLINE];
//
// Re-order
//
if(rv == ItXSuccess)
  for(i=0;i<n;i++) {
	sprintf(sstring,"Landmark_%03d",i+1);
	fx = points[i].x();
	fy = points[i].y();
	fz = points[i].z();
	ItXECode nrv = this->set(i,fz,fx,fy,sstring,1);
	if(nrv != ItXSuccess) {
		rv = nrv;
		break;
		}
	}
return(rv);
}



ItXECode Landmark::loadAnalyzeAVW(FILE *fp)
{
int i,x,y,z;
char sstring[MAXLOADEDLINE],nstring[50];
double fx,fy,fz,ftmp;

i = 0;
while(!feof(fp) && !ferror(fp)) {
	if(fgets(sstring,MAXLOADEDLINE - 1,fp)==NULL)
		break;
    sstring[MAXLOADEDLINE] = '\0';
	if(sscanf(sstring,"%d %d %d %lf\n",&x,&y,&z,&ftmp) == 4) {
		sprintf(nstring,"Landmark_%03d",++i);
		fx = (double)x * pixeldx;
		fy = (double)y * pixeldy;
		fz = (double)z * pixeldz;
		(void)this->add(fx,fy,fz,nstring,1);
		}
	}
if(i > 0) return(ItXSuccess);
else	  return(ItXError);
}


ItXECode Landmark::loadWarp(FILE *fp)
{
int nmax;

//
// Max # Landmarks
//
if((fscanf(fp,"%d\n",&nmax) != 1)||(nmax <= 0))
	return(ItXError);

int j;
ItXECode rv = ItXSuccess;

for(j=0;j<nmax && (rv == ItXSuccess) && !feof(fp) && !ferror(fp);j++) 
	rv = this->readWarp(fp);

if(j > 0) return(ItXSuccess);
else	  return(ItXError);
}



ItXECode Landmark::load(char *fname, int imgnx, int imgny, int imgnz,
			double pixdx, double pixdy, double pixdz)
{
FILE *fp;
double dx,dy,dz;
int   nx,ny,nz;
char sstring[MAXLOADEDLINE];
ItXECode rv;

if((fp=fopen(fname,"r"))==(FILE*)NULL)
	return(ItXInvalidFile);

sprintf(filename,fname);

pixeldx = pixdx;
pixeldy = pixdy;
pixeldz = pixdz;
imagenx = imgnx;
imageny = imgny;
imagenz = imgnz;

// get first line of file
fgets(sstring,MAXLOADEDLINE - 1,fp);
sstring[MAXLOADEDLINE] = '\0';

if(sscanf(sstring,"# BrainWorks Landmarks : Image Size (%d,%d,%d) Pixel Size (%lf,%lf,%lf)\n", &nx,&ny,&nz,&dx,&dy,&dz) == 6) {
	pixeldx = dx;
	pixeldy = dy;
	pixeldz = dz;
	imagenx = nx;
	imageny = ny;
	imagenz = nz;
	rv = loadBrainWorks(fp);
	fclose(fp);
	sprintf(filename,fname);
	return(rv);
	}


int    xx,yy,zz;
double fx,fy,fz;

//
// check for Analyze/Warp in sstring here
//

//if(!strncmp(sstring,WARP_HEADER_STR,
//        strlen(WARP_HEADER_STR))) {

if(sscanf(sstring, WARP_HEADER_FORMAT_STR, warpVersion) == 1) {
     fprintf(stderr,"Assuming '%s' are Warp landmarks\n",fname);
     rv = loadWarp(fp);
     }
else if(!strncmp(sstring,NEW_ANALYZE_HEADER_STR,
		strlen(NEW_ANALYZE_HEADER_STR))) {
     fprintf(stderr,"Assuming '%s' are Analyze AVW landmarks\n",fname);
     rv = loadAnalyzeAVW(fp);
     }
else if(sscanf(sstring,"%d %d %d %lf\n",&xx,&yy,&zz,&fx)==4) {
     rewind(fp);
     fprintf(stderr,"Assuming '%s' are Analyze AVW landmarks\n",fname);
     rv = loadAnalyzeAVW(fp);
     }
else if(sscanf(sstring,"%lf %lf %lf\n",&fx,&fy,&fz) == 3) {
     fprintf(stderr,"Assuming '%s' are Analyze landmarks\n",fname);
     rewind(fp);
     rv = loadAnalyze(fp);
     }
else {
     fprintf(stderr,"\nERROR: unknown landmark format in `%s'\n",fname);
     rv = ItXError;
     }

fclose(fp);

sprintf(filename,fname);
return(rv);
}

ItXECode Landmark::loadBrainWorks(FILE *fp)
{
int i = 0;
ItXECode rv = ItXSuccess;
while(rv == ItXSuccess) {
	rv = this->readBrainWorks(fp);
	i++;
	}
fclose(fp);
if(i > 0) return(ItXSuccess);
else	  return(ItXError);
}


ItXECode Landmark::loadStandard(FILE *fp)
{
int i = 0;
ItXECode rv = ItXSuccess;
while(rv == ItXSuccess) {
	rv = this->readStandard(fp);
	i++;
	}
fclose(fp);
if(i > 0) return(ItXSuccess);
else	  return(ItXError);
}



ItXECode Landmark::save(char *fname)
{
FILE *fp;
int i;
ItXECode rv;

if((fp=fopen(fname,"w"))==(FILE*)NULL)
	return(ItXError);

if (strlen(warpVersion))
  fprintf(fp,"%s\n",WARP_HEADER_STR);
else
  fprintf(fp,"%s\n",WARP_HEADER_STR);
fprintf(fp,"%d\n",n);
for(i=0;i<n;i++)  
   if((rv=this->saveWarp(fp,i)) != ItXSuccess) {
		fclose(fp);
		return(rv);
		}

fclose(fp);
sprintf(filename,fname);
return(ItXSuccess);
}




void Landmark::fillArray(Array2D<double> &d,
	double dx, double dy, double dz)
{
int i,cnt = 0;
for(i=0;i<n;i++) 
	if(valid[i]) cnt++;

if((dx==0.0)||(dy==0.0)||(dz==0.0)) {
	cerr << "ERROR: Landmark::fillMatrix() invalid pixel dimensions" << endl;
	exit(1);
	}

d.setDim(cnt,3);
for(i=0,cnt=0;i<n;i++) 
  if(valid[i]) {
	d[cnt][0] = points[i].x()/dx;
	d[cnt][1] = points[i].y()/dy;
	d[cnt++][2] = points[i].z()/dz;
	}
}

void Landmark::fillMatrix(Matrix<double> &d, 
		double dx, double dy, double dz)
{
int i,cnt;

cnt = 0;
for(i=0;i<n;i++) 
	if(valid[i]) cnt++;

if((dx==0.0)||(dy==0.0)||(dz==0.0)) {
   cerr << "ERROR: Landmark::fillMatrix() invalid pixel dimensions" << endl;
   exit(1);
}

d.setDim(cnt,3);
for(i=0,cnt=0;i<n;i++) 
   if(valid[i]) {
	d[cnt][0]   = points[i].x()/dx;
	d[cnt][1]   = points[i].y()/dy;
	d[cnt++][2] = points[i].z()/dz;
	}
}





ItXECode Landmark::set(int i, double xpt, double ypt, 
			double zpt, const char *nm, int isvalid)
{
if((i<0)||(i>=n))
	return(ItXError);
points[i].set(xpt,ypt,zpt);
valid[i] = isvalid;
if(nm) sprintf(names[i],"%s",nm);
else   names[i][0] = 0;
return(ItXSuccess);
}


ItXECode Landmark::saveStandard(FILE *fp, int i)
{
double fx = ((double)imagenx - points[i].x()/pixeldx - 1.0);
double fy = ((double)imageny - points[i].y()/pixeldy - 1.0);
double fz = ((double)imagenz - points[i].z()/pixeldz - 1.0);
int rv = fprintf(fp,"%lf %lf %lf %s\n",fx,fy,fz,names[i]);
if(rv > 0) return(ItXSuccess);
else	   return(ItXError);
}


ItXECode Landmark::saveWarp(FILE *fp, int i)
{
if((pixeldx <= 0.0)||(pixeldy <= 0.0)||(pixeldz <= 0.0)) {
	cerr << "ERROR: Landmark::save : invalid pixel dimensions" << endl;
	return(ItXError);
	}

double fx = points[i].x() / pixeldx;
double fy = points[i].y() / pixeldy;
double fz = points[i].z() / pixeldz;

int rv = fprintf(fp,"%c%s%c\n%lf %lf %lf %d 1\n",
		'"',names[i],'"',fx,fy,fz,valid[i]);
if(rv > 0) return(ItXSuccess);
else	   return(ItXError);
}


ItXECode Landmark::readWarp(FILE *fp)
{
double fx,fy,fz;
char *sptr,sstring[MAXLOADEDLINE];
int da,db;

if(!fp) return(ItXInvalidFile);
//
// read landmark name
//
if(fgets(sstring,MAXLOADEDLINE - 1,fp)==NULL)
	return(ItXError);
sstring[MAXLOADEDLINE] = '\0';
if(strlen(sstring) <= 0) 
	return(ItXError);

sptr = sstring;
if(*sptr == '"') sptr++;
if(sptr[strlen(sptr)-2] == '"')
	sptr[strlen(sptr)-2] = 0;

//
// NULL landmark names should be caught here
// ending loading process
//
if(strlen(sptr) <= 0)
	return(ItXError);

if(fscanf(fp,"%lf %lf %lf %d %d\n",&fx,&fy,&fz,&da,&db) != 5)
	return(ItXError);

double x = fx * pixeldx;
double y = fy * pixeldy;
double z = fz * pixeldz;

int nindx = this->add(x,y,z,sptr,da);
if(nindx >= 0)
	return(ItXSuccess);
else	return(ItXError);
}



ItXECode Landmark::readBrainWorks(FILE *fp)
{
double fxp,fyp,fzp;
char  line[MAXLOADEDLINE],*nm;
int   n;

if(!fp) return(ItXInvalidFile);

redo_read:
if(feof(fp) || ferror(fp))
	return(ItXInvalidFile);

(void)fgets(line,MAXLOADEDLINE - 1,fp);
line[MAXLOADEDLINE] = '\0';
if(strlen(line) <= 1) 
	return(ItXError);

line[strlen(line)-1] = 0;
if(sscanf(line,"%lf %lf %lf %n",&fxp,&fyp,&fzp,&n) != 3) {
	fprintf(stderr,"ERROR: line `%s' for landmarks is invalid\n",line);
	goto redo_read;
	}
if(strlen(line) < n)
     nm = "<>";
else nm = &(line[n]);

double x = ((double)imagenx - 1.0) * pixeldx - fxp;
double y = ((double)imageny - 1.0) * pixeldy - fyp;
double z = ((double)imagenz - 1.0) * pixeldz - fzp;

int nindx = this->add(x,y,z,nm,1);

if(nindx >= 0) return(ItXSuccess);
else	       return(ItXError);
}

ItXECode Landmark::readStandard(FILE *fp)
{
double fxp,fyp,fzp;
char  line[MAXLOADEDLINE],*nm;
int   n;

if(!fp) return(ItXInvalidFile);

redo_read:
if(feof(fp) || ferror(fp))
	return(ItXInvalidFile);

(void)fgets(line,MAXLOADEDLINE - 1,fp);
line[MAXLOADEDLINE] = '\0';

if(strlen(line) <= 1) 
	return(ItXError);

line[strlen(line)-1] = 0;

if(sscanf(line,"%lf %lf %lf %n",&fxp,&fyp,&fzp,&n) != 3) {
	fprintf(stderr,"ERROR: line `%s' for landmarks is invalid\n",line);
	goto redo_read;
	}
if(strlen(line) < n)
     nm = "<>";
else nm = &(line[n]);

double x = ((double)imagenx - fxp - 1.0) * pixeldx;
double y = ((double)imageny - fyp - 1.0) * pixeldy;
double z = ((double)imagenz - fzp - 1.0) * pixeldz;

int nindx = this->add(x,y,z,nm,1);

if(nindx >= 0) return(ItXSuccess);
else	       return(ItXError);
}

