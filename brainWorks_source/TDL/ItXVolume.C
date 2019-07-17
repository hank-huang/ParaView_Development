#include <iostream.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <TDL/IDLdefines.h>
#include <TDL/ItXVolume.h>
#include <TDL/ItXVolumeUtils.h>
#include <TDL/analyze_hdr.h>

#define ErrorExit(str) {cerr << str << endl; exit(1);}
//
// class constants
//
int ItXVolume::minZoom = 1;
int ItXVolume::maxZoom = 8;

const char *ItXVolume::dataTypeName(ItXDataType tpe)
{
switch(tpe) {
	case ItXVolume::UnsignedChar:
		return("Unsigned Char");
	case ItXVolume::UnsignedShort:
		return("Unsigned Short");
	case ItXVolume::Short:
		return("Short");
	case ItXVolume::Float:
		return("Float");
	case ItXVolume::Int:
		return("Int");
	case ItXVolume::RGB:
		return("RGB");
	}
return("<UNKNOWN>");
}


//--------------------------------------------------------------
// 
// setDefaults (protected)
//
//	set basic default parameters
//
//--------------------------------------------------------------
void ItXVolume::setDefaults()
{
nx               = 0;
ny               = 0;
nz               = 0;
data_type        = ItXVolume::UnsignedChar;
pixel_dimensionx = 1.0;
pixel_dimensiony = 1.0;
pixel_dimensionz = 1.0;
offset_x         = 0.0;
offset_y         = 0.0;
offset_z         = 0.0;

u_char_dat.setDim(0,0,0);
u_short_dat.setDim(0,0,0);
short_dat.setDim(0,0,0);
float_dat.setDim(0,0,0);
int_dat.setDim(0,0,0);
}

//--------------------------------------------------------------
//
// Basic constructor, zero size array (public)
//
//--------------------------------------------------------------
ItXVolume::ItXVolume()
{
sprintf(filename,"image");

setDefaults();

data_type   = UnsignedChar;
}

ItXVolume::ItXVolume(ItXDataType tpe)
{
sprintf(filename,"image");

setDefaults();

data_type   = tpe;
}


//--------------------------------------------------------------
//
//--------------------------------------------------------------
ItXVolume::ItXVolume(const char *file, ItXVolumeType tpe)
{
setDefaults();
if(this->load(file,tpe) != ItXSuccess)
   setDefaults();
}

ItXECode ItXVolume::loadDataRGB(Array3D<u_char> &dat,istream &file,bool doswap)
{
int x,y,z,xoff,yoff;
int nnx = dat.getXsize()/3;
int nny = dat.getYsize();
int nnz = dat.getZsize();
Array2D<u_char> pln(3,nnx*nny);

for(z=0;z<nnz;z++) {
   file.read((char*)pln.data(),pln.getSize());
   for(y=0;y<nny;y++) {
     yoff = y*nnx;
     for(x=0;x<nnx;x++) {
       xoff = x*3;
       dat[z][y][xoff+0] = pln[0][yoff+x];
       dat[z][y][xoff+1] = pln[1][yoff+x];
       dat[z][y][xoff+2] = pln[2][yoff+x];
     }
   }
}
if(!file) return(ItXError);
else      return(ItXSuccess);
}


ItXECode ItXVolume::saveDataRGB(Array3D<u_char> &dat, ostream &file)
{
int x,y,z,xoff,yoff;
int nnx = dat.getXsize()/3;
int nny = dat.getYsize();
int nnz = dat.getZsize();
Array2D<u_char> pln(3,nnx*nny);

for(z=0;z<nnz;z++) {
   for(y=0;y<nny;y++) {
     yoff = y*nnx;
     for(x=0;x<nnx;x++) {
       xoff = x*3;
       pln[0][yoff+x] = dat[z][y][xoff+0];
       pln[1][yoff+x] = dat[z][y][xoff+1];
       pln[2][yoff+x] = dat[z][y][xoff+2];
     }
   }
   file.write((char*)pln.data(),pln.getSize());
}
if(!file) return(ItXError);
else      return(ItXSuccess);
}




ItXECode ItXVolume::loadData(Array3D<u_char> &dat, istream &file, bool doswap)
{
file.read((char*)dat.data(),dat.getSize());
if(!file) return(ItXError);
// no byte-swapping needed for char volumes
return(ItXSuccess);
}

ItXECode ItXVolume::loadData(Array3D<u_short> &dat, istream &file, bool doswap)
{
file.read((char*)dat.data(),dat.getSize());
if(!file) return(ItXError);
if(doswap) {
   for(int k=0;k<dat.getZsize();k++)
     for(int j=0;j<dat.getYsize();j++)
       for(int i=0;i<dat.getXsize();i++)
	dat[k][j][i] = (unsigned short)swap_short((short)dat[k][j][i]);
}
return(ItXSuccess);
}


ItXECode ItXVolume::loadData(Array3D<short> &dat, istream &file, bool doswap)
{
file.read((char*)dat.data(),dat.getSize());
if(!file) return(ItXError);
if(doswap) {
   for(int k=0;k<dat.getZsize();k++)
     for(int j=0;j<dat.getYsize();j++)
       for(int i=0;i<dat.getXsize();i++)
	dat[k][j][i] = (short)swap_short((short)dat[k][j][i]);
}
return(ItXSuccess);
}


ItXECode ItXVolume::loadData(Array3D<float> &dat, istream &file, bool doswap)
{
file.read((char*)dat.data(),dat.getSize());
if(!file) return(ItXError);
if(doswap) {
   for(int k=0;k<dat.getZsize();k++)
     for(int j=0;j<dat.getYsize();j++)
       for(int i=0;i<dat.getXsize();i++)
	dat[k][j][i] = swap_float(dat[k][j][i]);
}
return(ItXSuccess);
}

ItXECode ItXVolume::loadData(Array3D<int> &dat, istream &file, bool doswap)
{
file.read((char*)dat.data(),dat.getSize());
if(!file) return(ItXError);
if(doswap) {
   for(int k=0;k<dat.getZsize();k++)
     for(int j=0;j<dat.getYsize();j++)
       for(int i=0;i<dat.getXsize();i++)
	dat[k][j][i] = swap_int(dat[k][j][i]);
}
return(ItXSuccess);
}

ItXECode ItXVolume::saveData(Array3D<u_char> &dat, ostream &file)
{
file.write((char*)dat.data(),dat.getSize());
if(!file) return(ItXError);
return(ItXSuccess);
}


ItXECode ItXVolume::saveData(Array3D<u_short> &dat, ostream &file)
{
file.write((char*)dat.data(),dat.getSize());
if(!file) return(ItXError);
return(ItXSuccess);
}

ItXECode ItXVolume::saveData(Array3D<short> &dat, ostream &file)
{
file.write((char*)dat.data(),dat.getSize());
if(!file) return(ItXError);
return(ItXSuccess);
}

ItXECode ItXVolume::saveData(Array3D<float> &dat, ostream &file)
{
file.write((char*)dat.data(),dat.getSize());
if(!file) return(ItXError);
return(ItXSuccess);
}

ItXECode ItXVolume::saveData(Array3D<int> &dat, ostream &file)
{
file.write((char*)dat.data(),dat.getSize());
if(!file) return(ItXError);
return(ItXSuccess);
}

void ItXVolume::flipDataUnsignedChar()
{
u_char *dat = u_char_dat.data();
u_char *ndat = new u_char[nx*ny*nz];

u_char *nptr = ndat + nx*ny*nz-1;

int i,j,k;
for(k=0;k<nz;k++)
  for(j=0;j<ny;j++)
    for(i=0;i<nx;i++,dat++,nptr--)
        *nptr = *dat;

memcpy((void*)u_char_dat.data(),(const void*)ndat,nx*ny*nz);
delete [] ndat;
}

void ItXVolume::flipDataRGB()
{
u_char *dat = u_char_dat.data();
u_char *ndat = new u_char[3*nx*ny*nz];

u_char *nptr = ndat + 3*nx*ny*nz-1;

int i,j,k;
for(k=0;k<nz;k++)
  for(j=0;j<ny;j++)
    for(i=0;i<nx;i++,dat+=3,nptr-=3) {
        *(nptr+2) = *(dat+0);
        *(nptr+1) = *(dat+1);
        *(nptr+0) = *(dat+2);
    }

memcpy((void*)u_char_dat.data(),(const void*)ndat,3*nx*ny*nz);
delete [] ndat;
}

void ItXVolume::flipDataFloat()
{
float *dat = float_dat.data();
float *ndat = new float[nx*ny*nz];

float *nptr = ndat + nx*ny*nz-1;

int i,j,k;
for(k=0;k<nz;k++)
  for(j=0;j<ny;j++)
    for(i=0;i<nx;i++,dat++,nptr--)
        *nptr = *dat;

memcpy((void*)float_dat.data(),(const void*)ndat,nx*ny*nz*sizeof(float));
delete [] ndat;
}


void ItXVolume::flipDataUnsignedShort()
{
u_short *dat  = u_short_dat.data();
u_short *ndat = new u_short[nx*ny*nz];
u_short *nptr = ndat + nx*ny*nz-1;

int i,j,k;
for(k=0;k<nz;k++)
  for(j=0;j<ny;j++)
    for(i=0;i<nx;i++,dat++,nptr--)
        *nptr = *dat;

memcpy((void*)u_short_dat.data(),(const void*)ndat,nx*ny*nz*sizeof(u_short));
delete [] ndat;
}


void ItXVolume::flipDataShort()
{
short *dat  = short_dat.data();
short *ndat = new short[nx*ny*nz];
short *nptr = ndat + nx*ny*nz-1;

int i,j,k;
for(k=0;k<nz;k++)
  for(j=0;j<ny;j++)
    for(i=0;i<nx;i++,dat++,nptr--)
        *nptr = *dat;

memcpy((void*)short_dat.data(),(const void*)ndat,nx*ny*nz*sizeof(short));
delete [] ndat;
}

void ItXVolume::flipDataInt()
{
int *dat  = int_dat.data();
int *ndat = new int[nx*ny*nz];
int *nptr = ndat + nx*ny*nz-1;

int i,j,k;
for(k=0;k<nz;k++)
  for(j=0;j<ny;j++)
    for(i=0;i<nx;i++,dat++,nptr--)
        *nptr = *dat;

memcpy((void*)int_dat.data(),(const void*)ndat,nx*ny*nz*sizeof(int));
delete [] ndat;
}


void ItXVolume::flipData()
{
  switch(data_type) {
	case ItXVolume::UnsignedChar:
		flipDataUnsignedChar();
		break;
	case ItXVolume::UnsignedShort:
		flipDataUnsignedShort();
		break;
	case ItXVolume::Short:
		flipDataShort();
		break;
	case ItXVolume::Float:
		flipDataFloat();
		break;
	case ItXVolume::Int:
		flipDataInt();
		break;
	case ItXVolume::RGB:
		flipDataRGB();
		break;
	default:
		cerr << "WARNING: cannot re-orient data" << endl;
		cerr << "         Unhandled data type" << endl;
		break;
	}
}

void ItXVolume::setFilename(const char *fname)
{
if(fname) sprintf(filename,"%s",fname);
}


//--------------------------------------------------------------
//
//--------------------------------------------------------------
ItXECode ItXVolume::load(const char *file, ItXVolumeType tpe)
{
ItXECode rval;
sprintf(filename,"%s",file);
switch(tpe) {
	case ItXVolume::RawVolume :
	  rval = ItXInvalidDataType;
	  break;
	case ItXVolume::AnalyzeVolume : 
	  if((rval=loadAnalyze()) == ItXSuccess)
	  	flipData();
	  else setDefaults();
	  break;
	default: 
	  rval = ItXInvalidImageType;
	  break;
	}

return(rval);
}


//----------------------------------------------------
//
// Constructor for loading a raw volume (IXRawVolume)
//
//----------------------------------------------------

ItXVolume::ItXVolume(const char *file, int xdim, int ydim, int zdim,
            ItXDataType tpe)
{
setDefaults();
setSize(xdim,ydim,zdim,tpe);

sprintf(filename,"%s",file);

if(loadRaw() != ItXSuccess) {
	cerr << "ERROR: ItXVolume ctor load of " << file 
	     << " failed" << endl;
	}
}


ItXECode ItXVolume::load(const char *file, int xdim, int ydim, int zdim,
			 float pixx, float pixy, float pixz,
            		 ItXDataType tpe)
{
data_type   = tpe;
nx          = xdim;
ny          = ydim;
nz          = zdim;
pixel_dimensionx = pixx;
pixel_dimensiony = pixy;
pixel_dimensionz = pixz;

sprintf(filename,"%s",file);

ItXECode rval;
if((rval=loadRaw()) != ItXSuccess) {
	cerr << "ERROR: ItXVolume::ItXVolume() load of " << file 
	     << " failed" << endl;
	}
return(rval);
}


ItXECode ItXVolume::loadRaw()
{
ItXECode rval;

ifstream inFile(filename);

if(!inFile) {
	cerr << "ERROR: cant open " << filename << endl;
	return(ItXInvalidFile);
	}

switch(data_type) {
	case ItXVolume::UnsignedChar  : 
		u_short_dat.setDim(0,0,0);  
		short_dat.setDim(0,0,0);  
		float_dat.setDim(0,0,0);  
		int_dat.setDim(0,0,0);  
		u_char_dat.setDim(nz,ny,nx);
		if(u_char_dat.isEmpty()) 
			rval = ItXNoMemory;
		else rval = loadData(u_char_dat,inFile);
		break;
	case ItXVolume::UnsignedShort : 
		u_char_dat.setDim(0,0,0); 
		short_dat.setDim(0,0,0); 
		float_dat.setDim(0,0,0);  
		int_dat.setDim(0,0,0);  
		u_short_dat.setDim(nz,ny,nx); 
		if(u_short_dat.isEmpty()) 
			rval = ItXNoMemory;
		else rval = loadData(u_short_dat,inFile);
		break;
	case ItXVolume::Short : 
		u_char_dat.setDim(0,0,0); 
		u_short_dat.setDim(0,0,0); 
		float_dat.setDim(0,0,0);  
		int_dat.setDim(0,0,0);  
		short_dat.setDim(nz,ny,nx); 
		if(short_dat.isEmpty()) 
		     rval = ItXNoMemory;
		else rval = loadData(short_dat,inFile);
		break;
	case ItXVolume::Float : 
		u_char_dat.setDim(0,0,0); 
		u_short_dat.setDim(0,0,0);  
		short_dat.setDim(0,0,0);  
		int_dat.setDim(0,0,0);  
		float_dat.setDim(nz,ny,nx); 
		if(float_dat.isEmpty()) 
			rval = ItXNoMemory;
		else rval = loadData(float_dat,inFile);
		break;
	case ItXVolume::RGB: 
		u_short_dat.setDim(0,0,0);  
		short_dat.setDim(0,0,0);  
		float_dat.setDim(0,0,0);  
		int_dat.setDim(0,0,0);  
		u_char_dat.setDim(nz,ny,nx*3);
		if(u_char_dat.isEmpty()) 
			rval = ItXNoMemory;
		else rval = loadDataRGB(u_char_dat,inFile);
		break;
	case ItXVolume::Int : 
		u_char_dat.setDim(0,0,0); 
		u_short_dat.setDim(0,0,0);  
		short_dat.setDim(0,0,0);  
		float_dat.setDim(0,0,0);  
		int_dat.setDim(nz,ny,nx); 
		if(int_dat.isEmpty()) 
			rval = ItXNoMemory;
		else rval = loadData(int_dat,inFile);
		break;
	default: 
		rval = ItXInvalidDataType;
	}
return(rval);
}



ItXECode ItXVolume::saveRaw(const char *fname)
{
ItXECode rval;

ofstream outFile(fname);

if(!outFile) {
	cerr << "ERROR: cant open " << fname << endl;
	return(ItXInvalidFile);
	}

switch(data_type) {
	case ItXVolume::UnsignedChar: 
		rval = saveData(u_char_dat,outFile);
		break;
	case ItXVolume::UnsignedShort: 
		rval = saveData(u_short_dat,outFile);
		break;
	case ItXVolume::Short: 
		rval = saveData(short_dat,outFile);
		break;
	case ItXVolume::Float: 
		rval = saveData(float_dat,outFile);
		break;
	case ItXVolume::Int: 
		rval = saveData(int_dat,outFile);
		break;
	case ItXVolume::RGB: 
		rval = saveDataRGB(u_char_dat,outFile);
		break;
	default: 
		rval = ItXInvalidDataType;
	}
return(rval);
}


//----------------------------------------------------
//
// Constructor for creating a blank volume
//
//----------------------------------------------------
ItXVolume::ItXVolume(int xdim, int ydim, int zdim,
            ItXDataType tpe)
{
sprintf(filename,"image");

setDefaults();
setSize(xdim,ydim,zdim,tpe);
}


void ItXVolume::setSize(int xdim,int ydim, int zdim, ItXDataType tpe)
{
data_type = tpe;
nx        = xdim;
ny        = ydim;
nz        = zdim;

switch(data_type) {
	case ItXVolume::UnsignedChar  : 
		u_short_dat.setDim(0,0,0);  
		short_dat.setDim(0,0,0);  
		float_dat.setDim(0,0,0);  
		int_dat.setDim(0,0,0);  
		u_char_dat.setDim(nz,ny,nx);  
		u_char_dat = 0;  
		break;
	case ItXVolume::UnsignedShort : 
		u_char_dat.setDim(0,0,0);
		short_dat.setDim(0,0,0);
		float_dat.setDim(0,0,0);  
		int_dat.setDim(0,0,0);  
		u_short_dat.setDim(nz,ny,nx); 
		u_short_dat = 0; 
		break;
	case ItXVolume::Short : 
		u_char_dat.setDim(0,0,0);
		u_short_dat.setDim(0,0,0);
		float_dat.setDim(0,0,0);  
		int_dat.setDim(0,0,0);  
		short_dat.setDim(nz,ny,nx); 
		short_dat = 0; 
		break;
	case ItXVolume::Float : 
		u_char_dat.setDim(0,0,0);
		u_short_dat.setDim(0,0,0);
		short_dat.setDim(0,0,0);
		int_dat.setDim(0,0,0);
		float_dat.setDim(nz,ny,nx); 
		float_dat = 0; 
		break;
	case ItXVolume::RGB  : 
		u_short_dat.setDim(0,0,0);  
		short_dat.setDim(0,0,0);  
		float_dat.setDim(0,0,0);  
		int_dat.setDim(0,0,0);  
		u_char_dat.setDim(nz,ny,nx*3);  
		u_char_dat = 0;  
		break;
	case ItXVolume::Int : 
		u_char_dat.setDim(0,0,0);
		u_short_dat.setDim(0,0,0);
		short_dat.setDim(0,0,0);
		float_dat.setDim(0,0,0);
		int_dat.setDim(nz,ny,nx); 
		int_dat = 0; 
		break;
	default: 
		cerr << "ERROR: ItXVolume::setSize() : defaulting to u_char" 
		     << endl;
		u_short_dat.setDim(0,0,0);  
		float_dat.setDim(0,0,0);  
		int_dat.setDim(0,0,0);  
		u_char_dat.setDim(nz,ny,nx);  
		u_char_dat = 0;  
		data_type = ItXVolume::UnsignedChar;
		break;
	}
}


ItXECode ItXVolume::save(const char *fname, ItXVolumeType itpe,
		   ItXDataType tpe)
{
if(tpe == ItXVolume::UnknownDataType) 
	tpe = data_type;

if(fname) sprintf(filename,fname);

ItXECode rval,convval;

if(data_type == tpe) {
	flipData();
	switch(itpe) {
		case ItXVolume::AnalyzeVolume:
			rval = this->saveAnalyze(filename);
		break;
		case ItXVolume::RawVolume:
			rval = this->saveRaw(filename);
			break;
		default: rval = ItXInvalidDataType;
		}
	flipData();
	}
else    {
	ItXVolume simg = *this;
	switch(tpe) {
		case ItXVolume::UnsignedChar:
			convval = simg.convertToUnsignedChar();
			break;
		case ItXVolume::UnsignedShort:
			convval = simg.convertToUnsignedShort();
			break;
		case ItXVolume::Short:
			convval = simg.convertToShort();
			break;
		case ItXVolume::Float:
			convval = simg.convertToFloat();
			break;
		default: 
			convval = ItXInvalidDataType;
			rval    = ItXInvalidDataType;
			break;
		}

	if(convval != ItXSuccess)
		rval = convval;
	else {
	     simg.flipData();
	     switch(itpe) {
		case ItXVolume::AnalyzeVolume:
			rval = simg.saveAnalyze(filename);
			break;
		case ItXVolume::RawVolume:
			rval = simg.saveRaw(filename);
			break;
		default: 
			rval = ItXInvalidDataType;
		}
	     }
	}
return(rval);
}


//=====================================================================
//
// Analyze format specific routines
//
//=====================================================================

void ItXVolume::constructAnalyzeNames(
                     const char *filename,
                           char *ret_prefix,
                           char *ret_hdrname,
                           char *ret_imgname)
{
int l = strlen(filename);

sprintf(ret_prefix,filename);

if((l>4)&&((!strncmp(&(ret_prefix[l-4]),".img",3))||
   (!strncmp(&(ret_prefix[l-4]),".hdr",3))) )
   ret_prefix[l-4] = 0;

sprintf(ret_hdrname,"%s.hdr",ret_prefix);
sprintf(ret_imgname,"%s.img",ret_prefix);
}

// make header values all BIG_ENDIAN for saving

bool ItXVolume::swapAnalyzeHeader(void *hdr)
{
int   i;
short s;
struct dsr *h = (struct dsr*)hdr;

// byte swap only if needed

s = h->dime.datatype;
if((s>=0)&&(s<DT_ALL))
    return(false);

h->hk.sizeof_hdr     = swap_int(h->hk.sizeof_hdr);
h->hk.extents        = swap_int(h->hk.extents);
h->hk.session_error  = swap_short(h->hk.session_error);

for(i=0;i<8;i++)
   h->dime.dim[i]    = swap_short(h->dime.dim[i]);
h->dime.unused1      = swap_short(h->dime.unused1);
h->dime.datatype     = swap_short(h->dime.datatype);
h->dime.bitpix       = swap_short(h->dime.bitpix);
h->dime.dim_un0      = swap_short(h->dime.dim_un0);

for(i=0;i<8;i++)
   h->dime.pixdim[i] = swap_float(h->dime.pixdim[i]);
h->dime.vox_offset   = swap_float(h->dime.vox_offset);
h->dime.roi_scale    = swap_float(h->dime.roi_scale);
h->dime.funused1     = swap_float(h->dime.funused1);
h->dime.funused2     = swap_float(h->dime.funused2);
h->dime.cal_max      = swap_float(h->dime.cal_max);
h->dime.cal_min      = swap_float(h->dime.cal_min);
h->dime.compressed   = swap_int(h->dime.compressed);
h->dime.verified     = swap_int(h->dime.verified);
h->dime.glmax        = swap_int(h->dime.glmax);
h->dime.glmin        = swap_int(h->dime.glmin);

h->hist.views        = swap_int(h->hist.views);
h->hist.vols_added   = swap_int(h->hist.vols_added);
h->hist.start_field  = swap_int(h->hist.start_field);
h->hist.field_skip   = swap_int(h->hist.field_skip);
h->hist.omax         = swap_int(h->hist.omax);
h->hist.omin         = swap_int(h->hist.omin);
h->hist.smax         = swap_int(h->hist.smax);
h->hist.smin         = swap_int(h->hist.smin);

return(true);
}

//--------------------------------------------------------------
//
// Read .hdr file for Analyze and set current volume properties
//
//--------------------------------------------------------------
ItXECode ItXVolume::readAnalyzeHeader(const char *hdrname, bool &need_swap)
{
struct dsr hdr;

ifstream inFile(hdrname);

if(!inFile) {
        cerr << "ERROR: cant read " << hdrname << endl;
        return(ItXInvalidFile);
        }

inFile.read((char*)&hdr,sizeof(struct dsr));

if(inFile.fail())
        return(ItXShortRead);

need_swap = swapAnalyzeHeader((void*)&hdr);

nx = hdr.dime.dim[1];
ny = hdr.dime.dim[2];
nz = hdr.dime.dim[3];

pixel_dimensionx = hdr.dime.pixdim[1];
pixel_dimensiony = hdr.dime.pixdim[2];
pixel_dimensionz = hdr.dime.pixdim[3];

if(pixel_dimensionx == 0.f)
   pixel_dimensionx = 1.f;

if(pixel_dimensiony == 0.f)
   pixel_dimensiony = 1.f;

if(pixel_dimensionz == 0.f)
   pixel_dimensionz = 1.f;
   
switch(hdr.hist.orient)
{
    case 0:
        planarOrientation = ItXTransverseUnflipped;
        break;
    case 1:
        planarOrientation = ItXCoronalUnflipped;
        break;
    case 2:
        planarOrientation = ItXSagittalUnflipped;
        break;
    case 3:
        planarOrientation = ItXTransverseFlipped;
        break;
    case 4:
        planarOrientation = ItXCoronalFlipped;
        break;
    case 5:
        planarOrientation = ItXSagittalFlipped;
        break;
    default:
        planarOrientation = ItxUnknownOrientationType;
        break;
}
ItXECode rval = ItXSuccess;

switch(hdr.dime.datatype) {
        case DT_UNSIGNED_CHAR:
                data_type = ItXVolume::UnsignedChar;
                break;
        case DT_SIGNED_SHORT:
                data_type = ItXVolume::Short;
                break;
	case DT_FLOAT:
                data_type = ItXVolume::Float;
                break;
        case DT_RGB:
                data_type = ItXVolume::RGB;
                break;
	case DT_SIGNED_INT:
                data_type = ItXVolume::Int;
                break;
	case DT_UNKNOWN:
		switch(hdr.dime.bitpix) {
			case 8:  data_type = ItXVolume::UnsignedChar; break;
			case 16: data_type = ItXVolume::UnsignedShort; break;
			case 24: data_type = ItXVolume::RGB; break;
			case 32: data_type = ItXVolume::Float; break;
			default: rval = ItXInvalidImageType; break;
			}
		break;
        default:
                rval = ItXInvalidImageType;
                break;
        }

return(rval);
}


//--------------------------------------------------------------
//
// Save .hdr file for Analyze using current volume properties
//
//--------------------------------------------------------------
ItXECode  ItXVolume::saveAnalyzeHeader(
                        const char *hdrname)
{
struct    dsr hdr;
float fmax,fmin;

getMinMax(fmin,fmax);

ofstream outFile(hdrname);

if(!outFile) {
        cerr << "ERROR: cant save to " << hdrname << endl;
        return(ItXInvalidFile);
        }

memset((void*)&hdr,0,sizeof(struct dsr));

memset(&hdr,0, sizeof(struct dsr));
for(int i=0;i<8;i++)
    hdr.dime.pixdim[i]=0.0;

hdr.dime.vox_offset = 0.0;
hdr.dime.roi_scale  = 1.0;
hdr.dime.funused1   = 0.0;
hdr.dime.funused2   = 0.0;
hdr.dime.cal_max    = 0.0;
hdr.dime.cal_min    = 0.0;

hdr.hk.regular    = 'r';
hdr.hk.sizeof_hdr = sizeof(struct dsr);

/* all Analyze images are taken as 4 dimensional */
hdr.dime.dim[0] = 4;  
hdr.dime.dim[1] = nx;
hdr.dime.dim[2] = ny;
hdr.dime.dim[3] = nz;
hdr.dime.dim[4] = 1;

hdr.dime.pixdim[1] = pixel_dimensionx;
hdr.dime.pixdim[2] = pixel_dimensiony;
hdr.dime.pixdim[3] = pixel_dimensionz;

hdr.dime.vox_offset = 0.0;
hdr.hist.orient     = (char) planarOrientation; 
strcpy(hdr.dime.vox_units," ");
strcpy(hdr.dime.cal_units," ");
hdr.dime.cal_max = 0.0;
hdr.dime.cal_min = 0.0;

switch(data_type){
        case ItXVolume::UnsignedChar:
                hdr.dime.datatype = DT_UNSIGNED_CHAR;
                hdr.dime.bitpix   = 8;
		hdr.dime.glmax    = 255;
		hdr.dime.glmin    = 0;
                break;
        case ItXVolume::UnsignedShort:
                // Analyze has no u_short
                hdr.dime.datatype = DT_SIGNED_SHORT;
                hdr.dime.bitpix   = 16;
		hdr.dime.glmax    = USHRT_MAX;
		hdr.dime.glmin    = 0;
                break;
        case ItXVolume::Short:
                hdr.dime.datatype = DT_SIGNED_SHORT;
                hdr.dime.bitpix   = 16;
		hdr.dime.glmax    = SHRT_MAX;
		hdr.dime.glmin    = SHRT_MIN;
                break;
        case ItXVolume::Float:
                hdr.dime.datatype = DT_FLOAT;
                hdr.dime.bitpix   = 32;
		hdr.dime.glmax    = INT_MAX;
		hdr.dime.glmin    = INT_MIN;
                break;
        case ItXVolume::RGB:
                hdr.dime.datatype = DT_RGB;
                hdr.dime.bitpix   = 24;
		hdr.dime.glmax    = 255;
		hdr.dime.glmin    = 0;
                break;
        case ItXVolume::Int:
                hdr.dime.datatype = DT_SIGNED_INT;
                hdr.dime.bitpix   = 32;
		hdr.dime.glmax    = INT_MAX;
		hdr.dime.glmin    = INT_MIN;
                break;
        default:
                return(ItXInvalidDataType);
        }

outFile.write((const char*)&hdr,sizeof(struct dsr));

if(outFile.fail())
        return(ItXShortWrite);
else    return(ItXSuccess);
}


//--------------------------------------------------------------
//
// loadAnalyze (protected)
//
//      Load the volume from an Analyze format file
//
//  returns ItXSucess on sucess
//
//--------------------------------------------------------------
ItXECode ItXVolume::loadAnalyze()
{
char imgname[MAXPATHLEN];
char hdrname[MAXPATHLEN];
char fprefix[MAXPATHLEN];
ItXECode rval;
bool doswap;

constructAnalyzeNames(filename,fprefix,hdrname,imgname);

sprintf(filename,"%s",fprefix);

if((rval=readAnalyzeHeader(hdrname,doswap)) != ItXSuccess)
        return(rval);

ifstream inFile(imgname);

if(!inFile) {
	cerr << "ERROR: cant open " << imgname << endl;
	return(ItXInvalidFile);
	}

switch(data_type) {
        case ItXVolume::UnsignedChar:
                u_short_dat.setDim(0,0,0);
                short_dat.setDim(0,0,0);
                float_dat.setDim(0,0,0);
                int_dat.setDim(0,0,0);
                u_char_dat.setDim(nz,ny,nx);
                if(u_char_dat.isEmpty()) {
		     setDefaults();
                     rval = ItXNoMemory;
                   }
		else rval = loadData(u_char_dat,inFile,doswap);
                break;

        case ItXVolume::UnsignedShort:
                u_char_dat.setDim(0,0,0);
                short_dat.setDim(0,0,0);
                float_dat.setDim(0,0,0);
                int_dat.setDim(0,0,0);
                u_short_dat.setDim(nz,ny,nx);
                if(u_short_dat.isEmpty()) {
		     setDefaults();
                     rval  = ItXNoMemory;
                     }
		else rval = loadData(u_short_dat,inFile,doswap);
                break;

        case ItXVolume::Short:
                u_char_dat.setDim(0,0,0);
                u_short_dat.setDim(0,0,0);
                float_dat.setDim(0,0,0);
                int_dat.setDim(0,0,0);
                short_dat.setDim(nz,ny,nx);
                if(short_dat.isEmpty()) {
		     setDefaults();
                     rval  = ItXNoMemory;
                     }
		else rval = loadData(short_dat,inFile,doswap);
                break;

        case ItXVolume::Float:
                u_char_dat.setDim(0,0,0);
                u_short_dat.setDim(0,0,0);
                short_dat.setDim(0,0,0);
                int_dat.setDim(0,0,0);
                float_dat.setDim(nz,ny,nx);
                if(float_dat.isEmpty()) {
		     setDefaults();
                     rval  = ItXNoMemory;
                     }
		else rval = loadData(float_dat,inFile,doswap);
                break;

        case ItXVolume::RGB:
                u_short_dat.setDim(0,0,0);
                short_dat.setDim(0,0,0);
                float_dat.setDim(0,0,0);
                int_dat.setDim(0,0,0);
                u_char_dat.setDim(nz,ny,nx*3);
                if(u_char_dat.isEmpty()) {
		     setDefaults();
                     rval = ItXNoMemory;
                   }
		else rval = loadDataRGB(u_char_dat,inFile,doswap);
                break;
        case ItXVolume::Int:
                u_char_dat.setDim(0,0,0);
                u_short_dat.setDim(0,0,0);
                short_dat.setDim(0,0,0);
                float_dat.setDim(0,0,0);
                int_dat.setDim(nz,ny,nx);
                if(int_dat.isEmpty()) {
		     setDefaults();
                     rval  = ItXNoMemory;
                     }
		else rval = loadData(int_dat,inFile,doswap);
                break;

        default:
                rval      = ItXInvalidDataType;
		setDefaults();
                break;
        }
return(rval);
}


//--------------------------------------------------------------
//
// saveAnalyze (protected)
//
//      Save the volume as an Analyze format volume
//      with both a header and an image file
//
//  returns ItXSucess on sucess
//
//--------------------------------------------------------------
ItXECode ItXVolume::saveAnalyze(const char *fname)
{
char      imgname[MAXPATHLEN];
char      hdrname[MAXPATHLEN];
char      fprefix[MAXPATHLEN];

constructAnalyzeNames(fname,fprefix,hdrname,imgname);

sprintf(filename,"%s",fprefix);

ItXECode rval = saveAnalyzeHeader(hdrname);
if(rval != ItXSuccess)
        return(rval);

ofstream outFile(imgname);

if(!outFile) {
	cerr << "ERROR: cant open " << imgname << endl;
	return(ItXInvalidFile);
	}

return(saveRaw(imgname));

}



//=====================================================================
//
// Datatype conversion routines
//
//=====================================================================
ItXECode ItXVolume::convertToUnsignedChar()
{
u_short smin,*sptr;
short   xmin,*xptr;
float   fmin,fmax,*fptr;
u_char  *cptr;
float   frange;
int     i,j,k;

if(data_type != ItXVolume::UnsignedChar) {
	u_char_dat.setDim(nz,ny,nx);
	if(u_char_dat.isEmpty()) {
		cerr << "ERROR: ItXVolume::convertToUnsignedChar()" << endl;
		cerr << "       Out of memory!" << endl;
		setDefaults();
		return(ItXNoMemory);
		}
	}

switch(data_type) {
	case ItXVolume::UnsignedChar:
		break;
	case ItXVolume::UnsignedShort:
		this->getMinMax(fmin,fmax);
		smin = (u_short)fmin;
		frange = fmax - fmin;
		if(frange == 0.0) frange = 1.0;
		sptr = u_short_dat.data();
		cptr = u_char_dat.data();
		for(k=0;k<nz;k++)
		  for(j=0;j<ny;j++)
		   for(i=0;i<nx;i++,sptr++,cptr++) 
			*cptr = (u_char)(((float)(*sptr-smin))/frange*255.0 + 0.5);
		u_short_dat.setDim(0,0,0);
		break;
	case ItXVolume::Short:
		this->getMinMax(fmin,fmax);
		xmin = (short)fmin;
		frange = fmax - fmin;
		if(frange == 0.0) frange = 1.0;
		xptr = short_dat.data();
		cptr = u_char_dat.data();
		for(k=0;k<nz;k++)
		  for(j=0;j<ny;j++)
		   for(i=0;i<nx;i++,xptr++,cptr++) 
			*cptr = (u_char)(((float)(*xptr-xmin))/frange*255.0 + 0.5);
		short_dat.setDim(0,0,0);
		break;
	case ItXVolume::Float:
		this->getMinMax(fmin,fmax);
		frange = fmax - fmin;
		if(frange == 0.0) frange = 1.0;
		fptr = float_dat.data();
		cptr = u_char_dat.data();
		for(k=0;k<nz;k++)
		  for(j=0;j<ny;j++)
		   for(i=0;i<nx;i++,fptr++,cptr++) 
			*cptr = (u_char)((*fptr-fmin)/frange*255.0 + 0.5);
		float_dat.setDim(0,0,0);
		break;
	default: 
		cerr << "WARNING: ItXVolume::convertToUnsignedChar()" << endl;
		cerr << "         Invalid datatype. " << endl;
		u_char_dat = 0;
		break;
	}
data_type = ItXVolume::UnsignedChar;
return(ItXSuccess);
}



ItXECode ItXVolume::convertToUnsignedShort()
{
u_short *sptr;
short   *xptr;
u_char  *cptr;
float    fmin,fmax,frange,*fptr;
int      i,j,k;

if(data_type != ItXVolume::UnsignedShort) {
	u_short_dat.setDim(nz,ny,nx);
	if(u_short_dat.isEmpty()) {
		cerr << "ERROR: ItXVolume::convertToUnsignedShort()" << endl;
		cerr << "       Out of memory!" << endl;
		setDefaults();
		return(ItXNoMemory);
		}
	}

switch(data_type) {
	case ItXVolume::UnsignedShort:
		break;
	case ItXVolume::UnsignedChar:
		sptr = u_short_dat.data();
		cptr = u_char_dat.data();
		for(k=0;k<nz;k++)
		  for(j=0;j<ny;j++)
		   for(i=0;i<nx;i++,cptr++,sptr++) 
			*sptr = (u_short) *cptr;
		u_char_dat.setDim(0,0,0);
		break;
	case ItXVolume::Short:
		this->getMinMax(fmin,fmax);
		frange = fmax - fmin;
		if(frange == 0.0) frange = 1.0;
		xptr = short_dat.data();
		sptr = u_short_dat.data();
		for(k=0;k<nz;k++)
		  for(j=0;j<ny;j++)
		   for(i=0;i<nx;i++,xptr++,sptr++) 
			*sptr = (u_short)((*xptr-fmin)/frange * 
				          (float)USHRT_MAX + 0.5);
		short_dat.setDim(0,0,0);
		break;
	case ItXVolume::Float:
		this->getMinMax(fmin,fmax);
		frange = fmax - fmin;
		if(frange == 0.0) frange = 1.0;
		fptr = float_dat.data();
		sptr = u_short_dat.data();
		for(k=0;k<nz;k++)
		  for(j=0;j<ny;j++)
		   for(i=0;i<nx;i++,fptr++,sptr++) 
			*sptr = (u_short)((*fptr-fmin)/frange * 
				          (float)USHRT_MAX + 0.5);
		float_dat.setDim(0,0,0);
		break;
	default: 
		cerr << "WARNING: ItXVolume::convertToUnsignedShort()" << endl;
		cerr << "         Invalid datatype. " << endl;
		u_short_dat = 0;
		break;
	}
data_type = ItXVolume::UnsignedShort;
return(ItXSuccess);
}


ItXECode ItXVolume::convertToShort()
{
u_short *sptr;
short   *xptr;
u_char  *cptr;
float    fmin,fmax,frange,*fptr;
int      i,j,k;

if(data_type != ItXVolume::Short) {
	short_dat.setDim(nz,ny,nx);
	if(short_dat.isEmpty()) {
		cerr << "ERROR: ItXVolume::convertToUnsignedShort()" << endl;
		cerr << "       Out of memory!" << endl;
		setDefaults();
		return(ItXNoMemory);
		}
	}

switch(data_type) {
	case ItXVolume::Short:
		break;
	case ItXVolume::UnsignedChar:
		xptr = short_dat.data();
		cptr = u_char_dat.data();
		for(k=0;k<nz;k++)
		  for(j=0;j<ny;j++)
		   for(i=0;i<nx;i++,cptr++,xptr++) 
			*xptr = (short)*cptr;
		u_char_dat.setDim(0,0,0);
		break;
	case ItXVolume::UnsignedShort:
		this->getMinMax(fmin,fmax);
		frange = fmax - fmin;
		if(frange == 0.0) frange = 1.0;
		xptr = short_dat.data();
		sptr = u_short_dat.data();
		for(k=0;k<nz;k++)
		  for(j=0;j<ny;j++)
		   for(i=0;i<nx;i++,xptr++,sptr++) 
			*xptr = (short)((*sptr-fmin)/frange * 
				          (float)SHRT_MAX + 0.5);
		u_short_dat.setDim(0,0,0);
		break;
	case ItXVolume::Float:
		this->getMinMax(fmin,fmax);
		frange = fmax - fmin;
		if(frange == 0.0) frange = 1.0;
		fptr = float_dat.data();
		xptr = short_dat.data();
		for(k=0;k<nz;k++)
		  for(j=0;j<ny;j++)
		   for(i=0;i<nx;i++,fptr++,xptr++) 
			*xptr = (short)((*fptr-fmin)/frange * 
				          (float)SHRT_MAX + 0.5);
		float_dat.setDim(0,0,0);
		break;
	default: 
		cerr << "WARNING: ItXVolume::convertToShort()" << endl;
		cerr << "         Invalid datatype. " << endl;
		short_dat = 0;
		break;
	}
data_type = ItXVolume::Short;
return(ItXSuccess);
}


ItXECode ItXVolume::convertToFloat()
{
u_short *sptr;
short   *xptr;
u_char  *cptr;
float   *fptr;
int      i,j,k;

if(data_type != ItXVolume::Float) {
	float_dat.setDim(nz,ny,nx);
	if(float_dat.isEmpty()) {
		cerr << "ERROR: ItXVolume::convertToFloat()" << endl;
		cerr << "       Out of memory!" << endl;
		setDefaults();
		return(ItXNoMemory);
		}
	}

switch(data_type) {
	case ItXVolume::Float:
		break;
	case ItXVolume::UnsignedChar:
		fptr = float_dat.data();
		cptr = u_char_dat.data();
		for(k=0;k<nz;k++)
		  for(j=0;j<ny;j++)
		   for(i=0;i<nx;i++,cptr++,fptr++) 
			*fptr = (float)*cptr;
		u_char_dat.setDim(0,0,0);
		break;
	case ItXVolume::UnsignedShort:
		fptr = float_dat.data();
		sptr = u_short_dat.data();
		for(k=0;k<nz;k++)
		  for(j=0;j<ny;j++)
		   for(i=0;i<nx;i++,sptr++,fptr++) 
			*fptr = (float)*sptr;
		u_short_dat.setDim(0,0,0);
		break;
	case ItXVolume::Short:
		fptr = float_dat.data();
		xptr = short_dat.data();
		for(k=0;k<nz;k++)
		  for(j=0;j<ny;j++)
		   for(i=0;i<nx;i++,xptr++,fptr++) 
			*fptr = (float)*xptr;
		short_dat.setDim(0,0,0);
		break;
	default: 
		cerr << "WARNING: ItXVolume::convertToFloat()" << endl;
		cerr << "         Invalid datatype. " << endl;
		float_dat = 0;
		break;
	}
data_type = ItXVolume::Float;
return(ItXSuccess);
}

ItXVolume & ItXVolume::operator= (ItXVolume const &I)
{
data_type   = I.data_type;
nx          = I.nx;
ny          = I.ny;
nz          = I.nz;

switch(data_type) {
	case ItXVolume::UnsignedChar:
		u_char_dat.setDim(nz,ny,nx);
		if(u_char_dat.isEmpty()) {
			cout << "FATAL ERROR: out of memory" << endl;
			exit(1);
			}
		u_char_dat  = I.u_char_dat;
		u_short_dat.setDim(0,0,0);
		short_dat.setDim(0,0,0);
		float_dat.setDim(0,0,0);
		int_dat.setDim(0,0,0);
		break;
	case ItXVolume::UnsignedShort:
		u_short_dat.setDim(nz,ny,nx);
		if(u_short_dat.isEmpty()) {
			cout << "FATAL ERROR: out of memory" << endl;
			exit(1);
			}
		u_short_dat = I.u_short_dat;
		u_char_dat.setDim(0,0,0);
		short_dat.setDim(0,0,0);
		float_dat.setDim(0,0,0);
		int_dat.setDim(0,0,0);
		break;
	case ItXVolume::Short:
		short_dat.setDim(nz,ny,nx);
		if(short_dat.isEmpty()) {
			cout << "FATAL ERROR: out of memory" << endl;
			exit(1);
			}
		short_dat = I.short_dat;
		u_char_dat.setDim(0,0,0);
		u_short_dat.setDim(0,0,0);
		float_dat.setDim(0,0,0);
		int_dat.setDim(0,0,0);
		break;
	case ItXVolume::Float:
		float_dat.setDim(nz,ny,nx);
		if(float_dat.isEmpty()) {
			cout << "FATAL ERROR: out of memory" << endl;
			exit(1);
			}
		float_dat = I.float_dat;
		u_char_dat.setDim(0,0,0);
		u_short_dat.setDim(0,0,0);
		short_dat.setDim(0,0,0);
		int_dat.setDim(0,0,0);
		break;
	case ItXVolume::RGB:
		u_char_dat.setDim(nz,ny,3*nx);
		if(u_char_dat.isEmpty()) {
			cout << "FATAL ERROR: out of memory" << endl;
			exit(1);
			}
		u_char_dat  = I.u_char_dat;
		u_short_dat.setDim(0,0,0);
		short_dat.setDim(0,0,0);
		float_dat.setDim(0,0,0);
		int_dat.setDim(0,0,0);
		break;
	case ItXVolume::Int:
		int_dat.setDim(nz,ny,nx);
		if(int_dat.isEmpty()) {
			cout << "FATAL ERROR: out of memory" << endl;
			exit(1);
			}
		float_dat = I.float_dat;
		u_char_dat.setDim(0,0,0);
		u_short_dat.setDim(0,0,0);
		short_dat.setDim(0,0,0);
		float_dat.setDim(0,0,0);
		break;
	default:
		cerr << "ERROR: ItXVolume::operator= invalid datatype" << endl;
		break;
	}

pixel_dimensionx = I.pixel_dimensionx;
pixel_dimensiony = I.pixel_dimensiony;
pixel_dimensionz = I.pixel_dimensionz;

offset_x      = I.offset_x;
offset_y      = I.offset_y;
offset_z      = I.offset_z;

center_x      = I.center_x;
center_y      = I.center_y;
center_z      = I.center_z;

sprintf(filename,"%s",I.filename);

return *this;
}

ItXVolume & ItXVolume::operator= (float f)
{
switch(data_type) {
	case ItXVolume::UnsignedChar:
		u_char_dat = (unsigned char)f;
		break;
	case ItXVolume::UnsignedShort:
		u_short_dat = (unsigned short)f;
		break;
	case ItXVolume::Short:
		short_dat = (short)f;
		break;
	case ItXVolume::Float:
		float_dat = f;
		break;
	case ItXVolume::RGB:
		u_char_dat = (unsigned char)f;
		break;
	case ItXVolume::Int:
		int_dat = (int)f;
		break;
	default:
		cerr << "ERROR: ItXVolume::operator= invalid datatype" << endl;
		break;
	}
return *this;
}


float ItXVolume::pointValue(int x, int y, int z, int rgb)
{
//
// ADL handles out of bounds
//
float rval;
switch(data_type) {
	case ItXVolume::UnsignedChar:
		rval = (float)u_char_dat[z][y][x];
		break;
	case ItXVolume::UnsignedShort:
		rval = (float)u_short_dat[z][y][x];
		break;
	case ItXVolume::Short:
		rval = (float)short_dat[z][y][x];
		break;
	case ItXVolume::Float:
		rval = float_dat[z][y][x];
		break;
	case ItXVolume::RGB:
		rval = (float)u_char_dat[z][y][3*x+rgb];
		break;
	case ItXVolume::Int:
		rval = (float)int_dat[z][y][x];
		break;
	default:
		rval = 0.0;
	}
return(rval);
}


void ItXVolume::setPointValue(int x, int y, int z, float f, int rgb)
{
//
// ADL handles out of bounds
//
switch(data_type) {
	case ItXVolume::UnsignedChar:
		u_char_dat[z][y][x] = (u_char)f;
		break;
	case ItXVolume::UnsignedShort:
		u_short_dat[z][y][x] = (u_short)f;
		break;
	case ItXVolume::Short:
		short_dat[z][y][x] = (short)f;
		break;
	case ItXVolume::Float:
		float_dat[z][y][x] = f;
		break;
	case ItXVolume::RGB:
		u_char_dat[z][y][3*x+rgb] = (u_char)f;
		break;
	case ItXVolume::Int:
		int_dat[z][y][x] = (int)f;
		break;
	}
}



void ItXVolume::getMinMax(float &minret, float &maxret)
{
u_char  cval,cmin,cmax,*cptr;
u_short sval,smin,smax,*sptr;
short   xval,xmin,xmax,*xptr;
float   fval,fmin,fmax,*fptr;
int     ival,imin,imax,*iptr;
int     i,sz;

sz = nx*ny*nz;
switch(data_type) {
	case ItXVolume::UnsignedChar:
		cptr = u_char_dat.data();
		cmin = cmax = *cptr;
		for(i=0;i<sz;i++,cptr++) {
			cval = *cptr;
			if(cmin > cval) cmin = cval;
			if(cmax < cval) cmax = cval;
			}
		minret = (float)cmin;
		maxret = (float)cmax;
		break;
	case ItXVolume::UnsignedShort:
		sptr = u_short_dat.data();
		smin = smax = *sptr;
		for(i=0;i<sz;i++,sptr++) {
			sval = *sptr;
			if(smin > sval) smin = sval;
			if(smax < sval) smax = sval;
			}
		minret = (float)smin;
		maxret = (float)smax;
		break;
	case ItXVolume::Short:
		xptr = short_dat.data();
		xmin = xmax = *xptr;
		for(i=0;i<sz;i++,xptr++) {
			xval = *xptr;
			if(xmin > xval) xmin = xval;
			if(xmax < xval) xmax = xval;
			}
		minret = (float)xmin;
		maxret = (float)xmax;
		break;
	case ItXVolume::Float:
		fptr = float_dat.data();
		fmin = fmax = *fptr;
		for(i=0;i<sz;i++,fptr++) {
			fval = *fptr;
			if(fmin > fval) fmin = fval;
			if(fmax < fval) fmax = fval;
			}
		minret = fmin;
		maxret = fmax;
		break;
	case ItXVolume::RGB:
		cptr = u_char_dat.data();
		cmin = cmax = *cptr;
		for(i=0;i<sz*3;i++,cptr++) {
			cval = *cptr;
			if(cmin > cval) cmin = cval;
			if(cmax < cval) cmax = cval;
			}
		minret = (float)cmin;
		maxret = (float)cmax;
		break;
	case ItXVolume::Int:
		iptr = int_dat.data();
		imin = imax = *iptr;
		for(i=0;i<sz;i++,iptr++) {
			ival = *iptr;
			if(imin > ival) imin = ival;
			if(imax < ival) imax = ival;
			}
		minret = (float)imin;
		maxret = (float)imax;
		break;
	default:
		break;
	}
}


