///////////////////////////////////////////////////////////////////////////
//
// File: ItXVolumeUtils.C
//
// Author: Keith Doolittle
//
// Purpose:
//
///////////////////////////////////////////////////////////////////////////



#include <iostream.h>
#include <string.h>
#include <ADL/Array3D.h>
#include <TDL/IDLdefines.h>
#include <TDL/ItXVolume.h>
#include <TDL/ItXVolumeUtils.h>
#include <TDL/Resample.h>
#include <TDL/Convolution.h>

#define InBounds(aa,bb,cc) ((aa>=0)&&(aa<nx)&&(bb>=0)&&(bb<ny)&&(cc>=0)&&(cc<nz))


static float blur_mask[27] = {
        1,1,1,
        1,1,1,
        1,1,1,

        1,1,1,
        1,1,1,
        1,1,1,

        1,1,1,
        1,1,1,
        1,1,1,
        };


static float sharp_mask[27] = {
        -1,-1,-1,
        -1,-1,-1,
        -1,-1,-1,

        -1,-1,-1,
        -1,26,-1,
        -1,-1,-1,

        -1,-1,-1,
        -1,-1,-1,
        -1,-1,-1,
        };

ItXECode ItXVolumeUtils::resample(ItXVolume &vol, 
		int newnx, int newny, int newnz)
{
if((newnx <= 0)||(newny <= 0)||(newnz <= 0))
        return(ItXInvalidDimensions);

switch(vol.data_type) {
        case ItXVolume::UnsignedChar :
                return(resampleUChar(vol,newnx,newny,newnz));
        case ItXVolume::UnsignedShort :
                return(resampleUShort(vol,newnx,newny,newnz));
        case ItXVolume::Short :
                return(resampleShort(vol,newnx,newny,newnz));
        case ItXVolume::Float :
                return(resampleFloat(vol,newnx,newny,newnz));
        case ItXVolume::Int :
                return(resampleInt(vol,newnx,newny,newnz));
	default:
		break;
        }

return(ItXInvalidDataType);
}


ItXECode ItXVolumeUtils::resample(ItXVolume &vol, 
		float scalex, float scaley, float scalez)
{
int newnx = (int)(scalex * (float)vol.nx);
int newny = (int)(scaley * (float)vol.ny);
int newnz = (int)(scalez * (float)vol.nz);

return(resample(vol,newnx,newny,newnz));
}

ItXECode ItXVolumeUtils::convolute(ItXVolume &vol,char *type,float para)
{
switch(vol.data_type) {
        case ItXVolume::UnsignedChar :
                return(convoluteUChar(vol,type,para));
        case ItXVolume::UnsignedShort :
                return(convoluteUShort(vol,type,para));
        case ItXVolume::Short :
                return(convoluteShort(vol,type,para));
        case ItXVolume::Float :
                return(convoluteFloat(vol,type,para));
        case ItXVolume::Int :
                return(convoluteInt(vol,type,para));
	default:
		break;
        }

return(ItXInvalidDataType);
}

//--------------------------------------------------------------
//
// resize  (public)
//
//       Resize volume to new dimensions.  If 'keep_data'
//       is not 0, then old data will be cropped/padded
//       into new volume
//
//  Returns: ItXSucess on sucess
//
//--------------------------------------------------------------
ItXECode ItXVolumeUtils::resize(ItXVolume &vol ,
		int newnx, int newny, int newnz, 
		int keep_data)
{
ItXECode rval;
switch(vol.data_type) {
   case ItXVolume::UnsignedChar:
        if(keep_data) rval = resizeUChar(vol,newnx,newny,newnz);
        else {
             vol.u_char_dat.setDim(newnz,newny,newnx);
             if(vol.u_char_dat.isEmpty()) {
                        vol.setDefaults();
                        rval = ItXNoMemory;
                        }
             else {
                  vol.nx = newnx;
                  vol.ny = newny;
                  vol.nz = newnz;
                  rval = ItXSuccess;
                  }
             }
        break;
   case ItXVolume::UnsignedShort:
        if(keep_data) rval = resizeUShort(vol,newnx,newny,newnz);
        else {
             vol.u_short_dat.setDim(newnz,newny,newnx);
             if(vol.u_short_dat.isEmpty()) {
                        vol.setDefaults();
                        rval = ItXNoMemory;
                        }
             else {
                  vol.nx = newnx;
                  vol.ny = newny;
                  vol.nz = newnz;
                  rval = ItXSuccess;
                  }
             }
        break;

   case ItXVolume::Short:
        if(keep_data) rval = resizeShort(vol,newnx,newny,newnz);
        else {
             vol.short_dat.setDim(newnz,newny,newnx);
             if(vol.short_dat.isEmpty()) {
                        vol.setDefaults();
                        rval = ItXNoMemory;
                        }
             else {
                  vol.nx = newnx;
                  vol.ny = newny;
                  vol.nz = newnz;
                  rval = ItXSuccess;
                  }
             }
        break;

   case ItXVolume::Float:
        if(keep_data) rval = resizeFloat(vol,newnx,newny,newnz);
        else {
             vol.float_dat.setDim(newnz,newny,newnx);
             if(vol.float_dat.isEmpty()) {
                        vol.setDefaults();
                        rval = ItXNoMemory;
                        }
             else {
                  vol.nx = newnx;
                  vol.ny = newny;
                  vol.nz = newnz;
                  rval = ItXSuccess;
                  }
             }
	break;
   case ItXVolume::Int:
        if(keep_data) rval = resizeInt(vol,newnx,newny,newnz);
        else {
             vol.int_dat.setDim(newnz,newny,newnx);
             if(vol.int_dat.isEmpty()) {
                        vol.setDefaults();
                        rval = ItXNoMemory;
                        }
             else {
                  vol.nx = newnx;
                  vol.ny = newny;
                  vol.nz = newnz;
                  rval = ItXSuccess;
                  }
             }
	break;
   default:
        rval = ItXInvalidDataType;
        break;
   }
return(rval);
}


ItXECode ItXVolumeUtils::blur(ItXVolume &vol)
{
ItXECode rval;
switch(vol.data_type) {
        case ItXVolume::UnsignedChar :
                rval = applyMaskUChar(vol,blur_mask,3,3,3);
                break;
        case ItXVolume::UnsignedShort :
                rval = applyMaskUShort(vol,blur_mask,3,3,3);
                break;
        case ItXVolume::Short :
                rval = applyMaskShort(vol,blur_mask,3,3,3);
                break;
        case ItXVolume::Float :
                rval = applyMaskFloat(vol,blur_mask,3,3,3);
                break;
        case ItXVolume::Int :
                rval = applyMaskInt(vol,blur_mask,3,3,3);
                break;
        default:
                rval = ItXInvalidDataType;
                break;
        }
return(rval);
}



ItXECode ItXVolumeUtils::sharpen(ItXVolume &vol)
{
ItXECode rval;
switch(vol.data_type) {
        case ItXVolume::UnsignedChar :
                rval = applyMaskUChar(vol,sharp_mask,3,3,3);
                break;
        case ItXVolume::UnsignedShort :
                rval = applyMaskUShort(vol,sharp_mask,3,3,3);
                break;
        case ItXVolume::Short :
                rval = applyMaskShort(vol,sharp_mask,3,3,3);
                break;
        case ItXVolume::Float :
                rval = applyMaskFloat(vol,sharp_mask,3,3,3);
                break;
        case ItXVolume::Int :
                rval = applyMaskInt(vol,sharp_mask,3,3,3);
                break;
        default:
                rval = ItXInvalidDataType;
                break;
        }
return(rval);
}


ItXECode ItXVolumeUtils::applyMask(ItXVolume &vol, float *msk, int mx,int my, int mz)
{
ItXECode rval;
switch(vol.data_type) {
        case ItXVolume::UnsignedChar :
                rval = applyMaskUChar(vol,msk,mx,my,mz);
                break;
        case ItXVolume::UnsignedShort :
                rval = applyMaskUShort(vol,msk,mx,my,mz);
                break;
        case ItXVolume::Short :
                rval = applyMaskShort(vol,msk,mx,my,mz);
                break;
        case ItXVolume::Float :
                rval = applyMaskFloat(vol,msk,mx,my,mz);
                break;
        case ItXVolume::Int :
                rval = applyMaskInt(vol,msk,mx,my,mz);
                break;
        default:
                rval = ItXInvalidDataType;
                break;
        }
return(rval);
}



ItXECode ItXVolumeUtils::linearScale(ItXVolume &vol,float fmin, float fmax)
{
ItXECode rval;
switch(vol.data_type) {
        case ItXVolume::UnsignedChar:
                if((fmin < 0.0)||(fmax > 255.0)||(fmin >= fmax)) {
                   cerr << "ERROR: ItXVolumeUtils::linearScale values out of range" 
			<< endl;
                   rval = ItXInvalidParms;
                   }
                else {
                   linearScaleUChar(vol,(u_char)fmin,(u_char)fmax);
                   rval = ItXSuccess;
                   }
                break;
        case ItXVolume::UnsignedShort:
                if((fmin < 0.0)||(fmax > (float)(USHRT_MAX))||(fmin >= fmax)) {
                   cerr << "ERROR: ItXVolumeUtils::linearScale values out of range"
 			<< endl;
                   rval = ItXInvalidParms;
                   }
                else {
                   linearScaleUShort(vol,(u_short)fmin,(u_short)fmax);
                   rval = ItXSuccess;
                   }
                break;

        case ItXVolume::Short:
                if((fmin < SHRT_MIN)||(fmax > SHRT_MAX)||(fmin >= fmax)) {
                   cerr << "ERROR: ItXVolumeUtils::linearScale values out of range"
 			<< endl;
                   rval = ItXInvalidParms;
                   }
                else {
                   linearScaleShort(vol,(short)fmin,(short)fmax);
                   rval = ItXSuccess;
                   }
                break;
        case ItXVolume::Float:
                if(fmin >= fmax) {
                   cerr << "ERROR: ItXVolumeUtils::linearScale values out of range"
 			<< endl;
                   rval = ItXInvalidParms;
                   }
                else {
                   linearScaleFloat(vol,fmin,fmax);
                   rval = ItXSuccess;
                   }
                break;
        case ItXVolume::Int:
                if(fmin >= fmax) {
                   cerr << "ERROR: ItXVolumeUtils::linearScale values out of range"
 			<< endl;
                   rval = ItXInvalidParms;
                   }
                else {
                   linearScaleInt(vol,(int)fmin,(int)fmax);
                   rval = ItXSuccess;
                   }
                break;
        default:
                rval = ItXInvalidDataType;
                break;
        }
return(rval);
}



ItXECode ItXVolumeUtils::resampleUChar(ItXVolume &vol, int newnx, int newny, int newnz)
{
float xfact = ((float)newnx)/((float)vol.nx);
float yfact = ((float)newny)/((float)vol.ny);
float zfact = ((float)newnz)/((float)vol.nz);

Array3D<u_char> newdat(newnz,newny,newnx);

if(newdat.isEmpty())
        return(ItXNoMemory);

u_char *ptr = newdat.data();

int   work_k = (int)( (float)newnz / 10.0 );

////////////////////////////////////////////////////
// modified by G. Tong

 vol.ShowWorking("Resizing...");

   Resample<unsigned char>::resample(Resample<unsigned char>::TRI,vol.u_char_dat.data()
				     ,vol.nx,vol.ny,vol.nz,ptr,newnx,newny,newnz);
 
////////////////////////////////////////////////////////

vol.ShowWorking((char*)NULL);

vol.u_char_dat.setDim(newnz,newny,newnx);
//
// See if resize was successful
//
if(vol.u_char_dat.isEmpty()) {
        vol.setDefaults();
        return(ItXNoMemory);
        }

vol.nx = newnx;
vol.ny = newny;
vol.nz = newnz;

vol.u_char_dat = newdat;

//
// Set new pixel sizes (offset stays same)
//
vol.pixel_dimensionx /= xfact;
vol.pixel_dimensiony /= yfact;
vol.pixel_dimensionz /= zfact;

return(ItXSuccess);
}


ItXECode ItXVolumeUtils::convoluteUChar(ItXVolume &vol, char *type,float para)
{
float xsize = (float)vol.nx;
float ysize = (float)vol.ny;
float zsize = (float)vol.nz;

Array3D<u_char> newdat(zsize,ysize,xsize);

if(newdat.isEmpty())
        return(ItXNoMemory);

u_char *ptr = newdat.data();

//int   work_k = (int)( (float)newnz / 10.0 );


 vol.ShowWorking("Convolving...");

 if(!strcmp(type,"GAUSS"))
 Convolution<unsigned char>::convolution(Convolution<unsigned char>::GAUSS,para,
					 vol.u_char_dat.data(),vol.nx,vol.ny,vol.nz,ptr);
 if(!strcmp(type,"CUBIC"))
 Convolution<unsigned char>::convolution(Convolution<unsigned char>::CUBIC,para,
					 vol.u_char_dat.data(),vol.nx,vol.ny,vol.nz,ptr);
 if(!strcmp(type,"LINEAR"))
 Convolution<unsigned char>::convolution(Convolution<unsigned char>::LINEAR,para,
					 vol.u_char_dat.data(),vol.nx,vol.ny,vol.nz,ptr);

 vol.ShowWorking("Done");

vol.u_char_dat.setDim(zsize,ysize,xsize);
//
// See if resize was successful
//
if(vol.u_char_dat.isEmpty()) {
        vol.setDefaults();
        return(ItXNoMemory);
        }

vol.nx = (int)xsize;
vol.ny = (int)ysize;
vol.nz = (int)zsize;

vol.u_char_dat = newdat;

return(ItXSuccess);
}






ItXECode ItXVolumeUtils::convoluteUShort(ItXVolume &vol, char *type,float para)
{
float xsize = (float)vol.nx;
float ysize = (float)vol.ny;
float zsize = (float)vol.nz;

Array3D<u_short> newdat(zsize,ysize,xsize);

if(newdat.isEmpty())
        return(ItXNoMemory);

u_short *ptr = newdat.data();

//int   work_k = (int)( (float)zsize / 10.0 );

 vol.ShowWorking("Convolving...");

 if(!strcmp(type,"GAUSS"))
   Convolution<unsigned short>::convolution( Convolution<unsigned short>::GAUSS,para,
					     vol.u_short_dat.data(),vol.nx,vol.ny,vol.nz,ptr);
 if(!strcmp(type,"CUBIC"))
   Convolution<unsigned short>::convolution( Convolution<unsigned short>::CUBIC,para,
					     vol.u_short_dat.data(),vol.nx,vol.ny,vol.nz,ptr);
 if(!strcmp(type,"LINEAR"))
   Convolution<unsigned short>::convolution( Convolution<unsigned short>::LINEAR,para,
					     vol.u_short_dat.data(),vol.nx,vol.ny,vol.nz,ptr);

vol.ShowWorking("Done");

vol.u_short_dat.setDim(zsize,ysize,xsize);

//
// See if resize was successful
//
if(vol.u_short_dat.isEmpty()) {
        vol.setDefaults();
        return(ItXNoMemory);
        }

vol.nx = (int)xsize;
vol.ny = (int)ysize;
vol.nz = (int)zsize;

vol.u_short_dat = newdat;

return(ItXSuccess);
}


ItXECode ItXVolumeUtils::convoluteShort(ItXVolume &vol, char *type,float para)
{
float xsize = (float)vol.nx;
float ysize = (float)vol.ny;
float zsize = (float)vol.nz;

Array3D<short> newdat(zsize,ysize,xsize);

if(newdat.isEmpty())
        return(ItXNoMemory);

 short *ptr = newdat.data();

 vol.ShowWorking("Convolving...");

 if(!strcmp(type,"GAUSS"))
   Convolution<short>::convolution( Convolution<short>::GAUSS,para,
					     vol.short_dat.data(),vol.nx,vol.ny,vol.nz,ptr);
 if(!strcmp(type,"CUBIC"))
   Convolution<short>::convolution( Convolution<short>::CUBIC,para,
					     vol.short_dat.data(),vol.nx,vol.ny,vol.nz,ptr);
 if(!strcmp(type,"LINEAR"))
   Convolution<short>::convolution( Convolution<short>::LINEAR,para,
					     vol.short_dat.data(),vol.nx,vol.ny,vol.nz,ptr);

vol.ShowWorking("Done");

vol.short_dat.setDim(zsize,ysize,xsize);

//
// See if resize was successful
//
if(vol.short_dat.isEmpty()) {
        vol.setDefaults();
        return(ItXNoMemory);
        }

vol.nx = (int)xsize;
vol.ny = (int)ysize;
vol.nz = (int)zsize;

vol.short_dat = newdat;

return(ItXSuccess);
}




ItXECode ItXVolumeUtils::convoluteFloat(ItXVolume &vol, char *type,float para)
{
float xsize = (float)vol.nx;
float ysize = (float)vol.ny;
float zsize = (float)vol.nz;

Array3D<float> newdat(zsize,ysize,xsize);

if(newdat.isEmpty())
        return(ItXNoMemory);

float *ptr = newdat.data();

//int   work_k = (int)( (float)zsize / 10.0 );

 vol.ShowWorking("Convolving...");

if(!strcmp(type,"GAUSS"))
  Convolution<float>::convolution( Convolution<float>::GAUSS,para,vol.float_dat.data()
				     ,vol.nx,vol.ny,vol.nz,ptr);
if(!strcmp(type,"CUBIC"))
  Convolution<float>::convolution( Convolution<float>::CUBIC,para,vol.float_dat.data()
				     ,vol.nx,vol.ny,vol.nz,ptr);
if(!strcmp(type,"LINEAR"))
  Convolution<float>::convolution( Convolution<float>::LINEAR,para,vol.float_dat.data()
				     ,vol.nx,vol.ny,vol.nz,ptr);

 vol.ShowWorking("Done");

vol.float_dat.setDim(zsize,ysize,xsize);

//
// See if resize was successful
//
if(vol.float_dat.isEmpty()) {
        vol.setDefaults();
        return(ItXNoMemory);
        }

vol.nx = (int)xsize;
vol.ny = (int)ysize;
vol.nz = (int)zsize;

vol.float_dat = newdat;

return(ItXSuccess);
}


ItXECode ItXVolumeUtils::convoluteInt(ItXVolume &vol, char *type,float para)
{
float xsize = (float)vol.nx;
float ysize = (float)vol.ny;
float zsize = (float)vol.nz;

Array3D<int> newdat(zsize,ysize,xsize);

if(newdat.isEmpty())
        return(ItXNoMemory);

int *ptr = newdat.data();
//int   work_k = (int)( (float)zsize / 10.0 );

 vol.ShowWorking("Convolving...");

if(!strcmp(type,"GAUSS"))
  Convolution<int>::convolution( Convolution<int>::GAUSS,para,vol.int_dat.data() ,vol.nx,vol.ny,vol.nz,ptr);
if(!strcmp(type,"CUBIC"))
  Convolution<int>::convolution( Convolution<int>::CUBIC,para,vol.int_dat.data() ,vol.nx,vol.ny,vol.nz,ptr);
if(!strcmp(type,"LINEAR"))
  Convolution<int>::convolution( Convolution<int>::LINEAR,para,vol.int_dat.data() ,vol.nx,vol.ny,vol.nz,ptr);

 vol.ShowWorking("Done");

vol.int_dat.setDim(zsize,ysize,xsize);

//
// See if resize was successful
//
if(vol.int_dat.isEmpty()) {
        vol.setDefaults();
        return(ItXNoMemory);
        }

vol.nx = (int)xsize;
vol.ny = (int)ysize;
vol.nz = (int)zsize;

vol.int_dat = newdat;

return(ItXSuccess);
}




ItXECode ItXVolumeUtils::resampleUShort(ItXVolume &vol, int newnx, int newny, int newnz)
{
float xfact = ((float)newnx)/((float)vol.nx);
float yfact = ((float)newny)/((float)vol.ny);
float zfact = ((float)newnz)/((float)vol.nz);

Array3D<u_short> newdat(newnz,newny,newnx);

if(newdat.isEmpty())
        return(ItXNoMemory);

u_short *ptr = newdat.data();

int   work_k = (int)( (float)newnz / 10.0 );

////////////////////////////////////////////////////
// modified by G. Tong

 vol.ShowWorking("Resizing...");

 Resample<unsigned short>::resample(Resample<unsigned short>::TRI,vol.u_short_dat.data()
				    ,vol.nx,vol.ny,vol.nz,ptr,newnx,newny,newnz);

/////////////////////////////////////////////////////

vol.ShowWorking((char*)NULL);

vol.u_short_dat.setDim(newnz,newny,newnx);

//
// See if resize was successful
//
if(vol.u_short_dat.isEmpty()) {
        vol.setDefaults();
        return(ItXNoMemory);
        }

vol.nx = newnx;
vol.ny = newny;
vol.nz = newnz;

vol.u_short_dat = newdat;
//
// Set new pixel sizes (offset stays same)
//
vol.pixel_dimensionx /= xfact;
vol.pixel_dimensiony /= yfact;
vol.pixel_dimensionz /= zfact;

return(ItXSuccess);
}


ItXECode ItXVolumeUtils::resampleShort(ItXVolume &vol, int newnx, int newny, int newnz)
{
float xfact = ((float)newnx)/((float)vol.nx);
float yfact = ((float)newny)/((float)vol.ny);
float zfact = ((float)newnz)/((float)vol.nz);

Array3D<short> newdat(newnz,newny,newnx);

if(newdat.isEmpty())
        return(ItXNoMemory);

short *ptr = newdat.data();

int   work_k = (int)( (float)newnz / 10.0 );

////////////////////////////////////////////////////
// modified by G. Tong

 vol.ShowWorking("Resizing...");

 Resample<short>::resample(Resample<short>::TRI,vol.short_dat.data(),
                           vol.nx,vol.ny,vol.nz,ptr,newnx,newny,newnz);

/////////////////////////////////////////////////////

vol.ShowWorking((char*)NULL);

vol.short_dat.setDim(newnz,newny,newnx);

//
// See if resize was successful
//
if(vol.short_dat.isEmpty()) {
        vol.setDefaults();
        return(ItXNoMemory);
        }

vol.nx = newnx;
vol.ny = newny;
vol.nz = newnz;

vol.short_dat = newdat;
//
// Set new pixel sizes (offset stays same)
//
vol.pixel_dimensionx /= xfact;
vol.pixel_dimensiony /= yfact;
vol.pixel_dimensionz /= zfact;

return(ItXSuccess);
}




ItXECode ItXVolumeUtils::resampleFloat(ItXVolume &vol, int newnx, int newny, int newnz)
{
float xfact = ((float)newnx)/((float)vol.nx);
float yfact = ((float)newny)/((float)vol.ny);
float zfact = ((float)newnz)/((float)vol.nz);

Array3D<float> newdat(newnz,newny,newnx);

if(newdat.isEmpty())
        return(ItXNoMemory);

float *ptr = newdat.data();

int   work_k = (int)( (float)newnz / 10.0 );

////////////////////////////////////////////////////// 
// modified by G. Tong

 vol.ShowWorking("Resizing...");

 Resample<float>::resample(Resample<float>::TRI,vol.float_dat.data()
			   ,vol.nx,vol.ny,vol.nz,ptr,newnx,newny,newnz);

/////////////////////////////////////////////////////////
 
vol.ShowWorking((char*)NULL);

vol.float_dat.setDim(newnz,newny,newnx);

//
// See if resize was successful
//
if(vol.float_dat.isEmpty()) {
        vol.setDefaults();
        return(ItXNoMemory);
        }

vol.nx = newnx;
vol.ny = newny;
vol.nz = newnz;

vol.float_dat = newdat;
//
// Set new pixel sizes (offset stays same)
//
vol.pixel_dimensionx /= xfact;
vol.pixel_dimensiony /= yfact;
vol.pixel_dimensionz /= zfact;

return(ItXSuccess);
}


ItXECode ItXVolumeUtils::resampleInt(ItXVolume &vol, int newnx, int newny, int newnz)
{
float xfact = ((float)newnx)/((float)vol.nx);
float yfact = ((float)newny)/((float)vol.ny);
float zfact = ((float)newnz)/((float)vol.nz);

Array3D<int> newdat(newnz,newny,newnx);

if(newdat.isEmpty())
        return(ItXNoMemory);

int *ptr = newdat.data();

int   work_k = (int)( (float)newnz / 10.0 );

////////////////////////////////////////////////////// 
// modified by G. Tong

 vol.ShowWorking("Resizing...");

 Resample<int>::resample(Resample<int>::TRI,vol.int_dat.data(),
			vol.nx,vol.ny,vol.nz,ptr,newnx,newny,newnz);

/////////////////////////////////////////////////////////
 
vol.ShowWorking((char*)NULL);

vol.int_dat.setDim(newnz,newny,newnx);

//
// See if resize was successful
//
if(vol.int_dat.isEmpty()) {
        vol.setDefaults();
        return(ItXNoMemory);
        }

vol.nx = newnx;
vol.ny = newny;
vol.nz = newnz;

vol.int_dat = newdat;
//
// Set new pixel sizes (offset stays same)
//
vol.pixel_dimensionx /= xfact;
vol.pixel_dimensiony /= yfact;
vol.pixel_dimensionz /= zfact;

return(ItXSuccess);
}




void ItXVolumeUtils::linearScaleUShort(ItXVolume &vol, u_short bmin, u_short bmax)
{
float frange = (bmax - bmin);
if(frange == 0.0) frange = 1.0;

u_short dmin,dval;
float drange,fmin,fmax;

vol.getMinMax(fmin,fmax);

dmin = (u_short)fmin;

int i,j,k;

drange = (fmax - fmin);
if(drange == 0.0) drange = 1.0;
for(k=0;k<vol.nz;k++)
  for(j=0;j<vol.ny;j++)
    for(i=0;i<vol.nx;i++) {
        dval = vol.u_short_dat[k][j][i];
        vol.u_short_dat[k][j][i] = 
		(u_short)((float)(dval - dmin)/drange * frange)+bmin;
        }
}



void ItXVolumeUtils::linearScaleShort(ItXVolume &vol, short bmin, short bmax)
{
float frange = (bmax - bmin);
if(frange == 0.0) frange = 1.0;

short dmin,dval;
float drange,fmin,fmax;

vol.getMinMax(fmin,fmax);

dmin = (short)fmin;

int i,j,k;

drange = (fmax - fmin);
if(drange == 0.0) drange = 1.0;
for(k=0;k<vol.nz;k++)
  for(j=0;j<vol.ny;j++)
    for(i=0;i<vol.nx;i++) {
        dval = vol.short_dat[k][j][i];
        vol.short_dat[k][j][i] = 
		(short)((float)(dval - dmin)/drange * frange)+bmin;
        }
}




void ItXVolumeUtils::linearScaleUChar(ItXVolume &vol, u_char bmin, u_char bmax)
{
float frange = (bmax - bmin);
if(frange == 0.0) frange = 1.0;

u_short dmin,dval;
float fmin,fmax,drange;

vol.getMinMax(fmin,fmax);
dmin = (u_char)fmin;

int i,j,k;

drange = (fmax - fmin);
if(drange == 0.0) drange = 1.0;
for(k=0;k<vol.nz;k++)
  for(j=0;j<vol.ny;j++)
    for(i=0;i<vol.nx;i++) {
        dval = vol.u_char_dat[k][j][i];
        vol.u_char_dat[k][j][i] = 
		(u_char)((float)(dval - dmin)/drange * frange)+bmin;
        }
}




void ItXVolumeUtils::linearScaleFloat(ItXVolume &vol, float bmin, float bmax)
{
float frange = (bmax - bmin);
if(frange == 0.0) frange = 1.0;

float dmin,dmax,dval;
float drange;

vol.getMinMax(dmin,dmax);

int i,j,k;

drange = (dmax - dmin);
if(drange == 0.0) drange = 1.0;
for(k=0;k<vol.nz;k++)
  for(j=0;j<vol.ny;j++)
    for(i=0;i<vol.nx;i++) {
        dval = vol.float_dat[k][j][i];
        vol.float_dat[k][j][i] = ((float)(dval - dmin)/drange * frange)+bmin;
        }
}


void ItXVolumeUtils::linearScaleInt(ItXVolume &vol, int bmin, int bmax)
{
double frange = (bmax - bmin);
if(frange == 0.0) frange = 1.0;

double dmin,dmax,dval;
double drange;
float  fmin,fmax;

vol.getMinMax(fmin,fmax);
dmin = fmin;
dmax = fmax;

int i,j,k;

drange = (dmax - dmin);
if(drange == 0.0) drange = 1.0;
for(k=0;k<vol.nz;k++)
  for(j=0;j<vol.ny;j++)
    for(i=0;i<vol.nx;i++) {
        dval = (double)vol.int_dat[k][j][i];
        vol.int_dat[k][j][i] = (int)((dval - dmin)/drange * frange)+bmin;
        }
}



ItXECode ItXVolumeUtils::applyMaskUChar(ItXVolume &vol, float *mask, 
					int mx, int my, int mz)
{
Array3D<float>  fdat(vol.nz,vol.ny,vol.nx);

if(fdat.isEmpty())
        return(ItXNoMemory);

float dmin,dmax,drange;
vol.getMinMax(dmin,dmax);
drange = dmax - dmin;
if(drange == 0.f) drange = 1.f;

int mpl = my*mz;
int dx  = mx/2;
int dy  = my/2;
int dz  = mz/2;

float mval;
float fmin = 1E20;
float fmax = -1E20;
float frange,fval;

int i,j,k,ii,jj,kk;
int xx,yy,zz;
int mskx,msky,mskz;

char work_string[100];
float work_perc = 100.0/((float)vol.nz);

for(k=0;k<vol.nz;k++) {
  sprintf(work_string,"%d%% Complete",(int)((float)k * work_perc));
  vol.ShowWorking(work_string);
  for(j=0;j<vol.ny;j++)
    for(i=0;i<vol.nx;i++) {

        float msum = 0.0;
        int   mcnt = 0;
        for(kk=-dz;kk<=dz;kk++) {
          zz = k + kk;
          if((zz < 0)||(zz >= vol.nz)) continue;
          mskz = kk + dz;

          for(jj=-dy;jj<=dy;jj++) {
            yy = j + jj;
            if((yy < 0)||(yy >= vol.ny)) continue;
            msky = jj + dy;

            for(ii=-dx;ii<=dx;ii++) {
              xx = i + ii;
              if((xx < 0)||(xx >= vol.nx)) continue;
              mskx = ii + dx;

              mval = *(mask + mskx + msky*my + mskz*mpl);
              if(mval) {
                msum += (mval * (float)vol.u_char_dat[zz][yy][xx]);
                mcnt++;
              }
              }
            }
          }
        if(mcnt) fdat[k][j][i] = msum/((float)mcnt);
        else     fdat[k][j][i] = 0;

        if(fdat[k][j][i] < fmin)      fmin = fdat[k][j][i];
        else if(fdat[k][j][i] > fmax) fmax = fdat[k][j][i];
        }
  }

frange = fmax - fmin;
if(frange == 0.0) frange = 1.0;
for(k=0;k<vol.nz;k++)
  for(j=0;j<vol.ny;j++)
    for(i=0;i<vol.nx;i++) {
        fval = (float)fdat[k][j][i];
//        vol.u_char_dat[k][j][i] = (u_char)((fval-fmin)/frange * 255.0 + 0.5);
//        vol.u_char_dat[k][j][i] = (u_char)(dmin + (fval-fmin)/frange * drange + 0.5);
        if(fval < 0.f) fval = 0.f;
        else if(fval > 255.f) fval = 255.f;
        vol.u_char_dat[k][j][i] = (u_char)(fval + 0.5);
        }

vol.ShowWorking((char*)NULL);

return(ItXSuccess);
}



ItXECode ItXVolumeUtils::applyMaskUShort(ItXVolume &vol, float *mask, 
			int mx, int my, int mz)
{
Array3D<float>  fdat(vol.nz,vol.ny,vol.nx);

if(fdat.isEmpty())
        return(ItXNoMemory);

float dmin,dmax,drange;
vol.getMinMax(dmin,dmax);
drange = dmax - dmin;
if(drange == 0.f) drange = 1.f;


int mpl = my * mz;
int dx = mx/2;
int dy = my/2;
int dz = mz/2;

float mval;
float fmin = 1E20;
float fmax = -1E20;
float fval,frange;

int i,j,k,ii,jj,kk;
int xx,yy,zz;
int mskx,msky,mskz;

char work_string[100];
float work_perc = 100.0/((float)vol.nz);

for(k=0;k<vol.nz;k++) {
  sprintf(work_string,"%d%% Complete",(int)((float)k * work_perc));
  vol.ShowWorking(work_string);
  for(j=0;j<vol.ny;j++)
    for(i=0;i<vol.nx;i++) {
        float msum = 0.0;
        int   mcnt = 0;
        for(kk=-dz;kk<=dz;kk++) {
          zz = k + kk;
          if((zz < 0)||(zz >= vol.nz)) continue;
          mskz = kk + dz;

          for(jj=-dy;jj<=dy;jj++) {
            yy = j + jj;
            if((yy < 0)||(yy >= vol.ny)) continue;
            msky = jj + dy;

            for(ii=-dx;ii<=dx;ii++) {
              xx = i + ii;
              if((xx < 0)||(xx >= vol.nx)) continue;
              mskx = ii + dx;

              mval = *(mask + mskx + msky*my + mskz*mpl);
              if(mval) {
                msum += (mval * (float)vol.u_short_dat[zz][yy][xx]);
                mcnt++;
              }
              }
            }
          }
        if(mcnt) fdat[k][j][i] = msum/((float)mcnt);
        else     fdat[k][j][i] = 0.0;

        if(fdat[k][j][i] < fmin)      fmin = fdat[k][j][i];
        else if(fdat[k][j][i] > fmax) fmax = fdat[k][j][i];
        }
  }

frange = fmax - fmin;
if(frange == 0.0) frange = 1.0;
for(k=0;k<vol.nz;k++)
  for(j=0;j<vol.ny;j++)
    for(i=0;i<vol.nx;i++) {
        fval = (float)fdat[k][j][i];
//        vol.u_short_dat[k][j][i] = (u_short)((fval-fmin)/frange * USHRT_MAX + 0.5);
//        vol.u_short_dat[k][j][i] = (u_short)(dmin + (fval-fmin)/frange * drange + 0.5);
        if(fval < 0.f) fval = 0.f;
        else if(fval > USHRT_MAX) fval = USHRT_MAX;
        vol.u_short_dat[k][j][i] = (u_short)(fval + 0.5);
        }

vol.ShowWorking((char*)NULL);

return(ItXSuccess);
}


ItXECode ItXVolumeUtils::applyMaskShort(ItXVolume &vol, float *mask, 
			int mx, int my, int mz)
{
Array3D<float>  fdat(vol.nz,vol.ny,vol.nx);

if(fdat.isEmpty())
        return(ItXNoMemory);

float dmin,dmax,drange;
vol.getMinMax(dmin,dmax);
drange = dmax - dmin;
if(drange == 0.f) drange = 1.f;

int mpl = my * mz;
int dx = mx/2;
int dy = my/2;
int dz = mz/2;

float mval;
float fmin = 1E20;
float fmax = -1E20;
float fval,frange;

int i,j,k,ii,jj,kk;
int xx,yy,zz;
int mskx,msky,mskz;

char work_string[100];
float work_perc = 100.0/((float)vol.nz);

for(k=0;k<vol.nz;k++) {
  sprintf(work_string,"%d%% Complete",(int)((float)k * work_perc));
  vol.ShowWorking(work_string);
  for(j=0;j<vol.ny;j++)
    for(i=0;i<vol.nx;i++) {
        float msum = 0.0;
        int   mcnt = 0;
        for(kk=-dz;kk<=dz;kk++) {
          zz = k + kk;
          if((zz < 0)||(zz >= vol.nz)) continue;
          mskz = kk + dz;

          for(jj=-dy;jj<=dy;jj++) {
            yy = j + jj;
            if((yy < 0)||(yy >= vol.ny)) continue;
            msky = jj + dy;

            for(ii=-dx;ii<=dx;ii++) {
              xx = i + ii;
              if((xx < 0)||(xx >= vol.nx)) continue;
              mskx = ii + dx;

              mval = *(mask + mskx + msky*my + mskz*mpl);
              if(mval) {
                msum += (mval * (float)vol.short_dat[zz][yy][xx]);
                mcnt++;
              }
              }
            }
          }
        if(mcnt) fdat[k][j][i] = msum/((float)mcnt);
        else     fdat[k][j][i] = 0.0;

        if(fdat[k][j][i] < fmin)      fmin = fdat[k][j][i];
        else if(fdat[k][j][i] > fmax) fmax = fdat[k][j][i];
        }
  }

frange = fmax - fmin;
if(frange == 0.0) frange = 1.0;
for(k=0;k<vol.nz;k++)
  for(j=0;j<vol.ny;j++)
    for(i=0;i<vol.nx;i++) {
        fval = (float)fdat[k][j][i];
//        vol.short_dat[k][j][i] = (short)((fval-fmin)/frange * SHRT_MAX + 0.5);
//        vol.short_dat[k][j][i] = (short)(dmin + (fval-fmin)/frange * drange + 0.5);

        if(fval < SHRT_MIN) fval = SHRT_MIN;
        else if(fval > SHRT_MAX) fval = SHRT_MAX;
        vol.short_dat[k][j][i] = (short)(fval + 0.5);
        }

vol.ShowWorking((char*)NULL);

return(ItXSuccess);
}


ItXECode ItXVolumeUtils::applyMaskFloat(ItXVolume &vol, float *mask, 
			int mx, int my, int mz)
{
Array3D<float>  fdat(vol.nz,vol.ny,vol.nx);

if(fdat.isEmpty())
        return(ItXNoMemory);

int mpl = my * mz;
int dx = mx/2;
int dy = my/2;
int dz = mz/2;

float mval;
float fmin = 1E20;
float fmax = -1E20;
float fval,frange;

int i,j,k,ii,jj,kk;
int xx,yy,zz;
int mskx,msky,mskz;

char work_string[100];
float work_perc = 100.0/((float)vol.nz);

for(k=0;k<vol.nz;k++) {
  sprintf(work_string,"%d%% Complete",(int)((float)k * work_perc));
  vol.ShowWorking(work_string);
  for(j=0;j<vol.ny;j++)
    for(i=0;i<vol.nx;i++) {
        float msum = 0.0;
        int   mcnt = 0;
        for(kk=-dz;kk<=dz;kk++) {
          zz = k + kk;
          if((zz < 0)||(zz >= vol.nz)) continue;
          mskz = kk + dz;

          for(jj=-dy;jj<=dy;jj++) {
            yy = j + jj;
            if((yy < 0)||(yy >= vol.ny)) continue;
            msky = jj + dy;

            for(ii=-dx;ii<=dx;ii++) {
              xx = i + ii;
              if((xx < 0)||(xx >= vol.nx)) continue;
              mskx = ii + dx;

              mval = *(mask + mskx + msky*my + mskz*mpl);
              if(mval) {
                msum += (mval * vol.float_dat[zz][yy][xx]);
                mcnt++;
              }
              }
            }
          }

        if(mcnt) fdat[k][j][i] = msum/((float)mcnt);
        else     fdat[k][j][i] = 0.0;

        if(fdat[k][j][i] < fmin)      fmin = fdat[k][j][i];
        else if(fdat[k][j][i] > fmax) fmax = fdat[k][j][i];
        }
  }

vol.float_dat = fdat;

vol.ShowWorking((char*)NULL);

return(ItXSuccess);
}


ItXECode ItXVolumeUtils::applyMaskInt(ItXVolume &vol, float *mask, 
			int mx, int my, int mz)
{
Array3D<int>  fdat(vol.nz,vol.ny,vol.nx);

if(fdat.isEmpty())
        return(ItXNoMemory);

int mpl = my * mz;
int dx = mx/2;
int dy = my/2;
int dz = mz/2;

double mval;
double fmin = 1E20;
double fmax = -1E20;
double fval,frange;

int i,j,k,ii,jj,kk;
int xx,yy,zz;
int mskx,msky,mskz;

char work_string[100];
float work_perc = 100.0/((float)vol.nz);

for(k=0;k<vol.nz;k++) {
  sprintf(work_string,"%d%% Complete",(int)((float)k * work_perc));
  vol.ShowWorking(work_string);
  for(j=0;j<vol.ny;j++)
    for(i=0;i<vol.nx;i++) {
        double msum = 0.0;
        int   mcnt = 0;
        for(kk=-dz;kk<=dz;kk++) {
          zz = k + kk;
          if((zz < 0)||(zz >= vol.nz)) continue;
          mskz = kk + dz;

          for(jj=-dy;jj<=dy;jj++) {
            yy = j + jj;
            if((yy < 0)||(yy >= vol.ny)) continue;
            msky = jj + dy;

            for(ii=-dx;ii<=dx;ii++) {
              xx = i + ii;
              if((xx < 0)||(xx >= vol.nx)) continue;
              mskx = ii + dx;

              mval = *(mask + mskx + msky*my + mskz*mpl);
              if(mval) {
                msum += (mval * vol.float_dat[zz][yy][xx]);
                mcnt++;
              }
              }
            }
          }

        if(mcnt) fdat[k][j][i] = (int)(msum/((double)mcnt)+0.5);
        else     fdat[k][j][i] = 0;

        if(fdat[k][j][i] < fmin)      fmin = fdat[k][j][i];
        else if(fdat[k][j][i] > fmax) fmax = fdat[k][j][i];
        }
  }

vol.int_dat = fdat;

vol.ShowWorking((char*)NULL);

return(ItXSuccess);
}



//--------------------------------------------------------------
//
// resizeUChar  (protected)
//
//       Crop or pad existing u_char volume to new dimensions.
//       Data is centered at new volume center
//
//  Returns: ItXSucess on sucess
//
//--------------------------------------------------------------

ItXECode ItXVolumeUtils::resizeUChar(ItXVolume &vol, int newnx, int newny, int newnz)
{
Array3D<u_char> ndat(newnz,newny,newnx);

if(ndat.isEmpty())
        return(ItXNoMemory);

int dx = (newnx - vol.nx)/2;
int dy = (newny - vol.ny)/2;
int dz = (newnz - vol.nz)/2;
int xx,yy,zz;

for(int k=0;k<newnz;k++) {
  zz = k - dz;
  for(int j=0;j<newny;j++) {
    yy = j - dy;
    for(int i=0;i<newnx;i++) {
        xx = i - dx;
        if((xx<0)||(xx>=vol.nx)||(yy<0)||(yy>=vol.ny)||(zz<0)||(zz>=vol.nz))
                ndat[k][j][i] = 0;
        else    ndat[k][j][i] = vol.u_char_dat[zz][yy][xx];
        }
    }
  }

vol.nx = newnx;
vol.ny = newny;
vol.nz = newnz;

vol.u_char_dat.setDim(vol.nz,vol.ny,vol.nx);

if(vol.u_char_dat.isEmpty()) {
        vol.setDefaults();
        return(ItXNoMemory);
        }
vol.u_char_dat = ndat;
return(ItXSuccess);
}






//--------------------------------------------------------------
//
// resizeUShort  (protected)
//
//       Crop or pad existing u_short volume to new dimensions.
//       Data is centered at new volume center
//
//  Returns: ItXSucess on sucess
//
//--------------------------------------------------------------
ItXECode ItXVolumeUtils::resizeUShort(ItXVolume &vol, int newnx, int newny, int newnz)
{
Array3D<u_short> ndat(newnz,newny,newnz);

if(ndat.isEmpty())
        return(ItXNoMemory);

int dx = (newnx - vol.nx)/2;
int dy = (newny - vol.ny)/2;
int dz = (newnz - vol.nz)/2;
int xx,yy,zz;

for(int k=0;k<newnz;k++) {
  zz = k - dz;
  for(int j=0;j<newny;j++) {
    yy = j - dy;
    for(int i=0;i<newnx;i++) {
        xx = i - dx;
        if((xx<0)||(xx>=vol.nx)||(yy<0)||(yy>=vol.ny)||(zz<0)||(zz>=vol.nz))
                ndat[k][j][i] = 0;
        else    ndat[k][j][i] = vol.u_short_dat[zz][yy][xx];
        }
    }
  }

vol.nx = newnx;
vol.ny = newny;
vol.nz = newnz;

vol.u_short_dat.setDim(vol.nz,vol.ny,vol.nx);

if(vol.u_short_dat.isEmpty()) {
        vol.setDefaults();
        return(ItXNoMemory);
        }
vol.u_short_dat = ndat;
return(ItXSuccess);
}


//--------------------------------------------------------------
//
// resizeShort  (protected)
//
//       Crop or pad existing short volume to new dimensions.
//       Data is centered at new volume center
//
//  Returns: ItXSucess on sucess
//
//--------------------------------------------------------------
ItXECode ItXVolumeUtils::resizeShort(ItXVolume &vol, int newnx, int newny, int newnz)
{
Array3D<short> ndat(newnz,newny,newnz);

if(ndat.isEmpty())
        return(ItXNoMemory);

int dx = (newnx - vol.nx)/2;
int dy = (newny - vol.ny)/2;
int dz = (newnz - vol.nz)/2;
int xx,yy,zz;

for(int k=0;k<newnz;k++) {
  zz = k - dz;
  for(int j=0;j<newny;j++) {
    yy = j - dy;
    for(int i=0;i<newnx;i++) {
        xx = i - dx;
        if((xx<0)||(xx>=vol.nx)||(yy<0)||(yy>=vol.ny)||(zz<0)||(zz>=vol.nz))
                ndat[k][j][i] = 0;
        else    ndat[k][j][i] = vol.short_dat[zz][yy][xx];
        }
    }
  }

vol.nx = newnx;
vol.ny = newny;
vol.nz = newnz;

vol.short_dat.setDim(vol.nz,vol.ny,vol.nx);

if(vol.u_short_dat.isEmpty()) {
        vol.setDefaults();
        return(ItXNoMemory);
        }
vol.short_dat = ndat;
return(ItXSuccess);
}


//--------------------------------------------------------------
//
// resizeFloat  (protected)
//
//       Crop or pad existing float volume to new dimensions.
//       Data is centered at new volume center
//
//  Returns: ItXSucess on sucess
//
//--------------------------------------------------------------

ItXECode ItXVolumeUtils::resizeFloat(ItXVolume &vol, int newnx, int newny, int newnz)
{
Array3D<float> ndat(newnz,newny,newnz);

if(ndat.isEmpty())
        return(ItXNoMemory);

int dx = (newnx - vol.nx)/2;
int dy = (newny - vol.ny)/2;
int dz = (newnz - vol.nz)/2;
int xx,yy,zz;

for(int k=0;k<newnz;k++) {
  zz = k - dz;
  for(int j=0;j<newny;j++) {
    yy = j - dy;
    for(int i=0;i<newnx;i++) {
        xx = i - dx;
        if((xx<0)||(xx>=vol.nx)||(yy<0)||(yy>=vol.ny)||(zz<0)||(zz>=vol.nz))
                ndat[k][j][i] = 0;
        else    ndat[k][j][i] = vol.float_dat[zz][yy][xx];
        }
    }
  }

vol.nx = newnx;
vol.ny = newny;
vol.nz = newnz;

vol.float_dat.setDim(vol.nz,vol.ny,vol.nx);

if(vol.float_dat.isEmpty()) {
        vol.setDefaults();
        return(ItXNoMemory);
        }
vol.float_dat = ndat;
return(ItXSuccess);
}

ItXECode ItXVolumeUtils::resizeInt(ItXVolume &vol, int newnx, int newny, int newnz)
{
Array3D<int> ndat(newnz,newny,newnz);

if(ndat.isEmpty())
        return(ItXNoMemory);

int dx = (newnx - vol.nx)/2;
int dy = (newny - vol.ny)/2;
int dz = (newnz - vol.nz)/2;
int xx,yy,zz;

for(int k=0;k<newnz;k++) {
  zz = k - dz;
  for(int j=0;j<newny;j++) {
    yy = j - dy;
    for(int i=0;i<newnx;i++) {
        xx = i - dx;
        if((xx<0)||(xx>=vol.nx)||(yy<0)||(yy>=vol.ny)||(zz<0)||(zz>=vol.nz))
                ndat[k][j][i] = 0;
        else    ndat[k][j][i] = vol.int_dat[zz][yy][xx];
        }
    }
  }

vol.nx = newnx;
vol.ny = newny;
vol.nz = newnz;

vol.int_dat.setDim(vol.nz,vol.ny,vol.nx);

if(vol.int_dat.isEmpty()) {
        vol.setDefaults();
        return(ItXNoMemory);
        }
vol.int_dat = ndat;
return(ItXSuccess);
}





const char *ItXVolumeUtils::filenameTail(const char *fil)
{
int l,i;

if(!fil || ((l=strlen(fil)) <= 2))
        return(fil);

for(i=l-2;i>0;i--)
  if(fil[i] == '/')
        return(&(fil[i+1]));

return(fil);
}

void ItXVolumeUtils::trilinearInterp(const ItXVolume &img,
				float x, float y, float z,
                                float &v1, float &v2, float &v3)
{
float rval;
switch(img.dataType()) {
	case ItXVolume::UnsignedChar:
		 v1 = (float)ItXVolumeUtils::trilinearInterp(
				img.u_char_dat,x,y,z);
		 break;
	case ItXVolume::UnsignedShort:
		 v1 = (float)ItXVolumeUtils::trilinearInterp(
				img.u_short_dat,x,y,z);
		 break;
	case ItXVolume::Short:
		 v1 = (float)ItXVolumeUtils::trilinearInterp(
				img.short_dat,x,y,z);
		 break;
	case ItXVolume::Float:
		 v1 = ItXVolumeUtils::trilinearInterp(
				img.float_dat,x,y,z);
		 break;
	case ItXVolume::Int:
		 v1 = ItXVolumeUtils::trilinearInterp(
				img.int_dat,x,y,z);
		 break;
	case ItXVolume::RGB:
		 ItXVolumeUtils::trilinearInterp(
				img.u_char_dat,x,y,z,v1,v2,v3);
		 break;
	default: v1 = 0.0;
	}
}


u_char ItXVolumeUtils::trilinearInterp(const Array3D<u_char> &dat, 
			         float x, float y, float z)
{
int x1,y1,z1,x2,y2,z2;
float dx,dy,dz,tx,ty,tz;
float l1,l2,l3,l4,l5,l6,l7,l8;
u_char rval;

l1 = l2 = l3 = l4 = 0.0;
l5 = l6 = l7 = l8 = 0.0;

x1 = (int)floor(x);
y1 = (int)floor(y);
z1 = (int)floor(z);
dx = x - (float)x1;
dy = y - (float)y1;
dz = z - (float)z1;
tx = 1.0 - dx;
ty = 1.0 - dy;
tz = 1.0 - dz;

x2 = x1 + 1;
y2 = y1 + 1;
z2 = z1 + 1;

int nx = dat.getXsize();
int ny = dat.getYsize();
int nz = dat.getZsize();

const u_char *const *const *image_data = dat.address();

if(InBounds(x1,y1,z1)) l1 = image_data[z1][y1][x1];
if(InBounds(x2,y1,z1)) l2 = image_data[z1][y1][x2];
if(InBounds(x2,y2,z1)) l3 = image_data[z1][y2][x2];
if(InBounds(x1,y2,z1)) l4 = image_data[z1][y2][x1];
if(InBounds(x1,y1,z2)) l5 = image_data[z2][y1][x1];
if(InBounds(x2,y1,z2)) l6 = image_data[z2][y1][x2];
if(InBounds(x2,y2,z2)) l7 = image_data[z2][y2][x2];
if(InBounds(x1,y2,z2)) l8 = image_data[z2][y2][x1];

rval = (u_char)(l1*tx*ty*tz + l2*dx*ty*tz +
                l3*dx*dy*tz + l4*tx*dy*tz +
                l5*tx*ty*dz + l6*dx*ty*dz +
                l7*dx*dy*dz + l8*tx*dy*dz + 0.5);
return(rval);
}

u_short ItXVolumeUtils::trilinearInterp(const Array3D<u_short> &dat,
				       float x, float y, float z)
{
int x1,y1,z1,x2,y2,z2;
float dx,dy,dz,tx,ty,tz;
float l1,l2,l3,l4,l5,l6,l7,l8;
u_short rval;

l1 = l2 = l3 = l4 = 0.0;
l5 = l6 = l7 = l8 = 0.0;

x1 = (int)floor(x);
y1 = (int)floor(y);
z1 = (int)floor(z);
dx = x - (float)x1;
dy = y - (float)y1;
dz = z - (float)z1;
tx = 1.0 - dx;
ty = 1.0 - dy;
tz = 1.0 - dz;

x2 = x1 + 1;
y2 = y1 + 1;
z2 = z1 + 1;

int nx = dat.getXsize();
int ny = dat.getYsize();
int nz = dat.getZsize();

const u_short *const *const *image_data = dat.address();

if(InBounds(x1,y1,z1)) l1 = image_data[z1][y1][x1];
if(InBounds(x2,y1,z1)) l2 = image_data[z1][y1][x2];
if(InBounds(x2,y2,z1)) l3 = image_data[z1][y2][x2];
if(InBounds(x1,y2,z1)) l4 = image_data[z1][y2][x1];
if(InBounds(x1,y1,z2)) l5 = image_data[z2][y1][x1];
if(InBounds(x2,y1,z2)) l6 = image_data[z2][y1][x2];
if(InBounds(x2,y2,z2)) l7 = image_data[z2][y2][x2];
if(InBounds(x1,y2,z2)) l8 = image_data[z2][y2][x1];

rval = (u_short)(l1*tx*ty*tz + l2*dx*ty*tz +
                l3*dx*dy*tz + l4*tx*dy*tz +
                l5*tx*ty*dz + l6*dx*ty*dz +
                l7*dx*dy*dz + l8*tx*dy*dz + 0.5);
return(rval);
}


short ItXVolumeUtils::trilinearInterp(const Array3D<short> &dat,
				       float x, float y, float z)
{
int x1,y1,z1,x2,y2,z2;
float dx,dy,dz,tx,ty,tz;
float l1,l2,l3,l4,l5,l6,l7,l8;
short rval;

l1 = l2 = l3 = l4 = 0.0;
l5 = l6 = l7 = l8 = 0.0;

x1 = (int)floor(x);
y1 = (int)floor(y);
z1 = (int)floor(z);
dx = x - (float)x1;
dy = y - (float)y1;
dz = z - (float)z1;
tx = 1.0 - dx;
ty = 1.0 - dy;
tz = 1.0 - dz;

x2 = x1 + 1;
y2 = y1 + 1;
z2 = z1 + 1;

int nx = dat.getXsize();
int ny = dat.getYsize();
int nz = dat.getZsize();

const short *const *const *image_data = dat.address();

if(InBounds(x1,y1,z1)) l1 = image_data[z1][y1][x1];
if(InBounds(x2,y1,z1)) l2 = image_data[z1][y1][x2];
if(InBounds(x2,y2,z1)) l3 = image_data[z1][y2][x2];
if(InBounds(x1,y2,z1)) l4 = image_data[z1][y2][x1];
if(InBounds(x1,y1,z2)) l5 = image_data[z2][y1][x1];
if(InBounds(x2,y1,z2)) l6 = image_data[z2][y1][x2];
if(InBounds(x2,y2,z2)) l7 = image_data[z2][y2][x2];
if(InBounds(x1,y2,z2)) l8 = image_data[z2][y2][x1];

rval = (short)(l1*tx*ty*tz + l2*dx*ty*tz +
               l3*dx*dy*tz + l4*tx*dy*tz +
               l5*tx*ty*dz + l6*dx*ty*dz +
               l7*dx*dy*dz + l8*tx*dy*dz + 0.5);
return(rval);
}






float ItXVolumeUtils::trilinearInterp(const Array3D<float> &dat,
				       float x, float y, float z)
{
int x1,y1,z1,x2,y2,z2;
float dx,dy,dz,tx,ty,tz;
float l1,l2,l3,l4,l5,l6,l7,l8;
float rval;

l1 = l2 = l3 = l4 = 0.0;
l5 = l6 = l7 = l8 = 0.0;

x1 = (int)floor(x);
y1 = (int)floor(y);
z1 = (int)floor(z);
dx = x - (float)x1;
dy = y - (float)y1;
dz = z - (float)z1;
tx = 1.0 - dx;
ty = 1.0 - dy;
tz = 1.0 - dz;

x2 = x1 + 1;
y2 = y1 + 1;
z2 = z1 + 1;

int nx = dat.getXsize();
int ny = dat.getYsize();
int nz = dat.getZsize();

const float *const *const *image_data = dat.address();

if(InBounds(x1,y1,z1)) l1 = image_data[z1][y1][x1];
if(InBounds(x2,y1,z1)) l2 = image_data[z1][y1][x2];
if(InBounds(x2,y2,z1)) l3 = image_data[z1][y2][x2];
if(InBounds(x1,y2,z1)) l4 = image_data[z1][y2][x1];
if(InBounds(x1,y1,z2)) l5 = image_data[z2][y1][x1];
if(InBounds(x2,y1,z2)) l6 = image_data[z2][y1][x2];
if(InBounds(x2,y2,z2)) l7 = image_data[z2][y2][x2];
if(InBounds(x1,y2,z2)) l8 = image_data[z2][y2][x1];

rval = (float)(l1*tx*ty*tz + l2*dx*ty*tz +
               l3*dx*dy*tz + l4*tx*dy*tz +
               l5*tx*ty*dz + l6*dx*ty*dz +
               l7*dx*dy*dz + l8*tx*dy*dz + 0.5);
return(rval);
}

float ItXVolumeUtils::trilinearInterp(const Array3D<double> &dat,
				       float x, float y, float z)
{
int x1,y1,z1,x2,y2,z2;
double dx,dy,dz,tx,ty,tz;
double l1,l2,l3,l4,l5,l6,l7,l8;
double rval;

l1 = l2 = l3 = l4 = 0.0;
l5 = l6 = l7 = l8 = 0.0;

x1 = (int)floor(x);
y1 = (int)floor(y);
z1 = (int)floor(z);
dx = (double)x - (double)x1;
dy = (double)y - (double)y1;
dz = (double)z - (double)z1;
tx = 1.0 - dx;
ty = 1.0 - dy;
tz = 1.0 - dz;

x2 = x1 + 1;
y2 = y1 + 1;
z2 = z1 + 1;

int nx = dat.getXsize();
int ny = dat.getYsize();
int nz = dat.getZsize();

const double *const *const *image_data = dat.address();

if(InBounds(x1,y1,z1)) l1 = image_data[z1][y1][x1];
if(InBounds(x2,y1,z1)) l2 = image_data[z1][y1][x2];
if(InBounds(x2,y2,z1)) l3 = image_data[z1][y2][x2];
if(InBounds(x1,y2,z1)) l4 = image_data[z1][y2][x1];
if(InBounds(x1,y1,z2)) l5 = image_data[z2][y1][x1];
if(InBounds(x2,y1,z2)) l6 = image_data[z2][y1][x2];
if(InBounds(x2,y2,z2)) l7 = image_data[z2][y2][x2];
if(InBounds(x1,y2,z2)) l8 = image_data[z2][y2][x1];

rval = (float)(l1*tx*ty*tz + l2*dx*ty*tz +
               l3*dx*dy*tz + l4*tx*dy*tz +
               l5*tx*ty*dz + l6*dx*ty*dz +
               l7*dx*dy*dz + l8*tx*dy*dz + 0.5);
return(rval);
}

int ItXVolumeUtils::trilinearInterp(const Array3D<int> &dat,
				       float x, float y, float z)
{
int x1,y1,z1,x2,y2,z2;
double dx,dy,dz,tx,ty,tz;
int l1,l2,l3,l4,l5,l6,l7,l8;
int rval;

l1 = l2 = l3 = l4 = 0;
l5 = l6 = l7 = l8 = 0;

x1 = (int)floor(x);
y1 = (int)floor(y);
z1 = (int)floor(z);
dx = x - (float)x1;
dy = y - (float)y1;
dz = z - (float)z1;
tx = 1.0 - dx;
ty = 1.0 - dy;
tz = 1.0 - dz;

x2 = x1 + 1;
y2 = y1 + 1;
z2 = z1 + 1;

int nx = dat.getXsize();
int ny = dat.getYsize();
int nz = dat.getZsize();

const int *const *const *image_data = dat.address();

if(InBounds(x1,y1,z1)) l1 = image_data[z1][y1][x1];
if(InBounds(x2,y1,z1)) l2 = image_data[z1][y1][x2];
if(InBounds(x2,y2,z1)) l3 = image_data[z1][y2][x2];
if(InBounds(x1,y2,z1)) l4 = image_data[z1][y2][x1];
if(InBounds(x1,y1,z2)) l5 = image_data[z2][y1][x1];
if(InBounds(x2,y1,z2)) l6 = image_data[z2][y1][x2];
if(InBounds(x2,y2,z2)) l7 = image_data[z2][y2][x2];
if(InBounds(x1,y2,z2)) l8 = image_data[z2][y2][x1];

rval = (int)(l1*tx*ty*tz + l2*dx*ty*tz +
             l3*dx*dy*tz + l4*tx*dy*tz +
             l5*tx*ty*dz + l6*dx*ty*dz +
             l7*dx*dy*dz + l8*tx*dy*dz + 0.5);
return(rval);
}


void ItXVolumeUtils::trilinearInterp(const Array3D<u_char> &dat, 
			         float x, float y, float z,
                                 float &r, float &g, float &b)
{
int x1,y1,z1,x2,y2,z2;
float dx,dy,dz,tx,ty,tz;
float l1,l2,l3,l4,l5,l6,l7,l8;

x1 = (int)floor(x);
y1 = (int)floor(y);
z1 = (int)floor(z);
dx = x - (float)x1;
dy = y - (float)y1;
dz = z - (float)z1;
tx = 1.0 - dx;
ty = 1.0 - dy;
tz = 1.0 - dz;

x2 = x1 + 1;
y2 = y1 + 1;
z2 = z1 + 1;

int nx = dat.getXsize();
int ny = dat.getYsize();
int nz = dat.getZsize();

l1 = l2 = l3 = l4 = 0.0;
l5 = l6 = l7 = l8 = 0.0;
if(InBounds(x1,y1,z1)) l1 = dat[z1][y1][3*x1+0];
if(InBounds(x2,y1,z1)) l2 = dat[z1][y1][3*x2+0];
if(InBounds(x2,y2,z1)) l3 = dat[z1][y2][3*x2+0];
if(InBounds(x1,y2,z1)) l4 = dat[z1][y2][3*x1+0];
if(InBounds(x1,y1,z2)) l5 = dat[z2][y1][3*x1+0];
if(InBounds(x2,y1,z2)) l6 = dat[z2][y1][3*x2+0];
if(InBounds(x2,y2,z2)) l7 = dat[z2][y2][3*x2+0];
if(InBounds(x1,y2,z2)) l8 = dat[z2][y2][3*x1+0];

r = (l1*tx*ty*tz + l2*dx*ty*tz +
     l3*dx*dy*tz + l4*tx*dy*tz +
     l5*tx*ty*dz + l6*dx*ty*dz +
     l7*dx*dy*dz + l8*tx*dy*dz + 0.5);


l1 = l2 = l3 = l4 = 0.0;
l5 = l6 = l7 = l8 = 0.0;
if(InBounds(x1,y1,z1)) l1 = dat[z1][y1][3*x1+1];
if(InBounds(x2,y1,z1)) l2 = dat[z1][y1][3*x2+1];
if(InBounds(x2,y2,z1)) l3 = dat[z1][y2][3*x2+1];
if(InBounds(x1,y2,z1)) l4 = dat[z1][y2][3*x1+1];
if(InBounds(x1,y1,z2)) l5 = dat[z2][y1][3*x1+1];
if(InBounds(x2,y1,z2)) l6 = dat[z2][y1][3*x2+1];
if(InBounds(x2,y2,z2)) l7 = dat[z2][y2][3*x2+1];
if(InBounds(x1,y2,z2)) l8 = dat[z2][y2][3*x1+1];

g = (l1*tx*ty*tz + l2*dx*ty*tz +
     l3*dx*dy*tz + l4*tx*dy*tz +
     l5*tx*ty*dz + l6*dx*ty*dz +
     l7*dx*dy*dz + l8*tx*dy*dz + 0.5);

l1 = l2 = l3 = l4 = 0.0;
l5 = l6 = l7 = l8 = 0.0;
if(InBounds(x1,y1,z1)) l1 = dat[z1][y1][3*x1+2];
if(InBounds(x2,y1,z1)) l2 = dat[z1][y1][3*x2+2];
if(InBounds(x2,y2,z1)) l3 = dat[z1][y2][3*x2+2];
if(InBounds(x1,y2,z1)) l4 = dat[z1][y2][3*x1+2];
if(InBounds(x1,y1,z2)) l5 = dat[z2][y1][3*x1+2];
if(InBounds(x2,y1,z2)) l6 = dat[z2][y1][3*x2+2];
if(InBounds(x2,y2,z2)) l7 = dat[z2][y2][3*x2+2];
if(InBounds(x1,y2,z2)) l8 = dat[z2][y2][3*x1+2];

b = (l1*tx*ty*tz + l2*dx*ty*tz +
     l3*dx*dy*tz + l4*tx*dy*tz +
     l5*tx*ty*dz + l6*dx*ty*dz +
     l7*dx*dy*dz + l8*tx*dy*dz + 0.5);
}

void ItXVolumeUtils::keepOneValue(ItXVolume &img, float val)
{
if(img.dataType() == ItXVolume::UnsignedChar) {
	u_char cval = (u_char)val;
	for(int k=0;k<img.getSizeZ();k++)
	  for(int j=0;j<img.getSizeY();j++)
	    for(int i=0;i<img.getSizeX();i++)
		if(img.u_char_dat[k][j][i] != cval)
			img.u_char_dat[k][j][i] = 0;
	}
else if(img.dataType() == ItXVolume::UnsignedShort) {
	u_short cval = (u_short)val;
	for(int k=0;k<img.getSizeZ();k++)
	  for(int j=0;j<img.getSizeY();j++)
	    for(int i=0;i<img.getSizeX();i++)
		if(img.u_short_dat[k][j][i] != cval)
			img.u_short_dat[k][j][i] = 0;
	}
else if(img.dataType() == ItXVolume::Short) {
	short xval = (short)val;
	for(int k=0;k<img.getSizeZ();k++)
	  for(int j=0;j<img.getSizeY();j++)
	    for(int i=0;i<img.getSizeX();i++)
		if(img.short_dat[k][j][i] != xval)
			img.short_dat[k][j][i] = 0;
	}
else if(img.dataType() == ItXVolume::Float) {
	for(int k=0;k<img.getSizeZ();k++)
	  for(int j=0;j<img.getSizeY();j++)
	    for(int i=0;i<img.getSizeX();i++)
		if(img.float_dat[k][j][i] != val)
			img.float_dat[k][j][i] = 0;
	}
}


ItXECode ItXVolumeUtils::extract(ItXVolume &in, ItXVolume &out,
				 int sx, int sy, int sz)
{
if(in.data_type != out.data_type) {
	cerr << "ERROR: ItXVolumeUtils::extract needs same data type" << endl;
	return(ItXInvalidDataType);
	}

if((sx < 0)||(sy < 0)||(sz < 0))
	return(ItXInvalidParms);

int do_fill = 0;
int ex = sx + out.getSizeX() - 1;
if(ex >= in.getSizeX())
	ex = in.getSizeX() - 1, 
	do_fill = 1;

int ey = sy + out.getSizeY() - 1;
if(ey >= in.getSizeY())
	ey = in.getSizeY() - 1, 
	do_fill = 1;

int ez = sz + out.getSizeZ() - 1;
if(ez >= in.getSizeZ())
	ez = in.getSizeZ() - 1, 
	do_fill = 1;

int x,y,z,xx,yy,zz;
switch(in.data_type) {
	case ItXVolume::UnsignedChar:
		if(do_fill) out.u_char_dat = (u_char)0;
		for(z=sz,zz=0;z<=ez;z++,zz++)
  		  for(y=sy,yy=0;y<=ey;y++,yy++)
    		    for(x=sx,xx=0;x<=ex;x++,xx++)
			out.u_char_dat[zz][yy][xx] = in.u_char_dat[z][y][x];
		break;
	case ItXVolume::UnsignedShort:
		if(do_fill) out.u_short_dat = (u_short)0;
		for(z=sz,zz=0;z<=ez;z++,zz++)
  		  for(y=sy,yy=0;y<=ey;y++,yy++)
    		    for(x=sx,xx=0;x<=ex;x++,xx++)
			out.u_short_dat[zz][yy][xx] = in.u_short_dat[z][y][x];
		break;
	case ItXVolume::Short:
		if(do_fill) out.short_dat = (u_short)0;
		for(z=sz,zz=0;z<=ez;z++,zz++)
  		  for(y=sy,yy=0;y<=ey;y++,yy++)
    		    for(x=sx,xx=0;x<=ex;x++,xx++)
			out.short_dat[zz][yy][xx] = in.short_dat[z][y][x];
		break;
	case ItXVolume::Float:
		if(do_fill) out.float_dat = (float)0.0;
		for(z=sz,zz=0;z<=ez;z++,zz++)
  		  for(y=sy,yy=0;y<=ey;y++,yy++)
    		    for(x=sx,xx=0;x<=ex;x++,xx++)
			out.float_dat[zz][yy][xx] = in.float_dat[z][y][x];
		break;
	case ItXVolume::Int:
		if(do_fill) out.int_dat = 0;
		for(z=sz,zz=0;z<=ez;z++,zz++)
  		  for(y=sy,yy=0;y<=ey;y++,yy++)
    		    for(x=sx,xx=0;x<=ex;x++,xx++)
			out.int_dat[zz][yy][xx] = in.int_dat[z][y][x];
		break;
	default:
		return(ItXInvalidDataType);
	}
return(ItXSuccess);
}


void ItXVolumeUtils::clear(ItXVolume &in, float val)
{
if(in.dataType() == ItXVolume::UnsignedChar)
	in.u_char_dat = (u_char)val;
else if(in.dataType() == ItXVolume::UnsignedShort)
	in.u_short_dat = (u_short)val;
else if(in.dataType() == ItXVolume::Short)
	in.short_dat = (short)val;
else if(in.dataType() == ItXVolume::Float)
	in.float_dat = val;
else if(in.dataType() == ItXVolume::Int)
	in.int_dat = val;

}

//function for global thresholding,returns 0/255 image- muge 
ItXECode ItXVolumeUtils::globalThreshUChar(ItXVolume &vol, int thresh, ItXVolume &result) 
{
 
	if (vol.dataType() != ItXVolume:: UnsignedChar)
		return (ItXInvalidDataType); 
	for (int k=0; k<vol.nz; ++k)
		for (int j=0; j<vol.ny; ++j)
			for (int i=0; i<vol.nx; ++i) {
				if ((int) vol.u_char_dat[k][j][i]<thresh)
					result.u_char_dat[k][j][i] = (u_char) 0;
				else result.u_char_dat[k][j][i] = (u_char) 255;
	
			}
	return(ItXSuccess);
}

ItXECode ItXVolumeUtils::invertVolume(ItXVolume &vol) {

	if (vol.dataType() != ItXVolume:: UnsignedChar)
			return (ItXInvalidDataType); 
	for (int k=0; k<vol.nz; ++k)
		for (int j=0; j<vol.ny; ++j)
			for (int i=0; i<vol.nx; ++i)
				vol.u_char_dat[k][j][i] =(u_char) (255 - (int) vol.u_char_dat[k][j][i]);
					
	return(ItXSuccess);
}
	
//function to apply a greyscale erosion
//uses a cube structuring element of size 2dx+1, 2dy+1, 2dz+1
 
ItXECode ItXVolumeUtils::erosionUChar(ItXVolume &vol, int size ){
  
  if(vol.dataType() != ItXVolume:: UnsignedChar)
		return (ItXInvalidDataType); 

  ItXVolume result(vol.getSizeX(), vol.getSizeY(),  vol.getSizeZ(),
		 ItXVolume::UnsignedChar);

  int i,j,k,ii,jj,kk;
  int xx,yy,zz;
  int mskx,msky,mskz;

  char work_string[100];
  float work_perc = 100.0/ (size);

  //set the boundary to 0
  // Assumes the data is 0;
 

  for (int n=0; n<size; n++) {
	  sprintf(work_string,"%d%% Complete", (int)(work_perc * n ));
	  vol.ShowWorking(work_string);
	  for(k=1; k<vol.nz-1; k++) {
  
		  for(j=1; j<vol.ny-1; j++) {
			  for(i=1; i<vol.nx-1; i++) { 
				  result.u_char_dat[k][j][i] = MIN3(vol.u_char_dat[k][j][i],
													vol.u_char_dat[k][j][i-1],
													vol.u_char_dat[k][j][i+1]);
			  }
		  }
	  }

	  for(k=1; k<vol.nz-1; k++) {
		  for(j=1; j<vol.ny-1; j++) {
			  for(i=1; i<vol.nx-1; i++) { 
				  vol.u_char_dat[k][j][i] = MIN3(result.u_char_dat[k][j][i],
												 result.u_char_dat[k][j-1][i],
												 result.u_char_dat[k][j+1][i]);
			  }
		  }
	  }

	  for(k=1; k<vol.nz-1; k++) {
		  for(j=1; j<vol.ny-1; j++) {
			  for(i=1; i<vol.nx-1; i++) { 
				  result.u_char_dat[k][j][i] = MIN3(vol.u_char_dat[k][j][i],
													vol.u_char_dat[k-1][j][i],
													vol.u_char_dat[k+1][j][i]);
			  }
		  }
	  }
    
	  vol=result;
  }
  return(ItXSuccess);
}

//slow erosion
ItXECode ItXVolumeUtils::erosionSlowUChar(ItXVolume &vol, int dx, int dy, int dz, 
				  ItXVolume &result){

  if(vol.dataType() != ItXVolume:: UnsignedChar)
		return (ItXInvalidDataType); 

  int i,j,k,ii,jj,kk;
  int xx,yy,zz;
  int mskx,msky,mskz;

  char work_string[100];
  float work_perc = 100.0/((float)vol.nz);

  for(k=0;k<vol.nz;k++) {
    sprintf(work_string,"%d%% Complete",(int)((float)k * work_perc));
    vol.ShowWorking(work_string);
    for(j=0;j<vol.ny;j++) {
      for(i=0;i<vol.nx;i++) { 
	float cur_min=(float) vol.u_char_dat[k][j][i];
	for (kk=-dz; kk<=dz; kk++) {
	  zz= k + kk;
	  if((zz < 0)||(zz >= vol.nz)) continue;
	  for(jj=-dy;jj<=dy;jj++) {
            yy = j + jj;
            if((yy < 0)||(yy >= vol.ny))   continue;  
	    for(ii=-dx;ii<=dx;ii++) {
              xx = i + ii;
              if((xx < 0) || (xx >= vol.nx))  continue; 
	      if ((float) vol.u_char_dat[zz][yy][xx] < cur_min) {
	      	cur_min = (float) vol.u_char_dat[zz][yy][xx];
	      }
	    } //end ii
	  } // end jj
	} //end kk
      result.u_char_dat[k][j][i]=(u_char) cur_min;
    }  //end i
  } //end j
} //end k
return(ItXSuccess);
}




//function to apply grayscale dilation
//uses cubic structuring element of size 2dx+1,2dy+1,2dz+1
ItXECode ItXVolumeUtils::dilationSlowUChar(ItXVolume &vol, int dx, int dy, int dz, 
				  ItXVolume &result){

  if(vol.dataType() != ItXVolume:: UnsignedChar)
		return (ItXInvalidDataType); 

  int i,j,k,ii,jj,kk;
  int xx,yy,zz;
  int mskx,msky,mskz;

  char work_string[100];
  float work_perc = 100.0/((float)vol.nz);

  for(k=0;k<vol.nz;k++) {
    sprintf(work_string,"%d%% Complete",(int)((float)k * work_perc));
    vol.ShowWorking(work_string);
    for(j=0;j<vol.ny;j++) {
      for(i=0;i<vol.nx;i++) { 
	float cur_max=(float) vol.u_char_dat[k][j][i];
	for (kk=-dz; kk<=dz; kk++) {
	  zz= k + kk;
	  if((zz < 0)||(zz >= vol.nz)) continue;
	  for(jj=-dy;jj<=dy;jj++) {
            yy = j + jj;
            if((yy < 0)||(yy >= vol.ny))   continue;  
	    for(ii=-dx;ii<=dx;ii++) {
              xx = i + ii;
              if((xx < 0) || (xx >= vol.nx))  continue; 
	      if ((float) vol.u_char_dat[zz][yy][xx] > cur_max) {
	      	cur_max = (float) vol.u_char_dat[zz][yy][xx];
	      }
	    } //end ii
	  } // end jj
	} //end kk
      result.u_char_dat[k][j][i]=(u_char) cur_max;
    }  //end i
  } //end j
} //end k
return(ItXSuccess);
}


ItXECode ItXVolumeUtils::openingUChar(ItXVolume &vol, int dx) 
{
  if(vol.dataType() != ItXVolume:: UnsignedChar)
		return (ItXInvalidDataType); 
 
  ItXVolumeUtils::erosionUChar(vol,dx);
  ItXVolumeUtils::dilationUChar(vol,dx);
  return (ItXSuccess);

}


ItXECode ItXVolumeUtils::RegionGrow(ItXVolume &vol, ItXVolume &result, 
				    int cx, int cy, int cz, int toleranceL,int toleranceU)
{
  if(vol.dataType() != ItXVolume:: UnsignedChar)
		return (ItXInvalidDataType); 
				
  int height = vol.ny; 
  int width  = vol.nx;  
  int depth  = vol.nz;
	
  int minSeed = (int) vol.u_char_dat[cz][cy][cx] -  toleranceL;
  int maxSeed = (int) vol.u_char_dat[cz][cy][cx] + toleranceU;
  if (minSeed < 0 )  {
	cout<<" Min seed <0.  Cannnot region grow "<<endl;
	return(ItXError);
}

  if (maxSeed >255)  { 
	cout<<" Max seed >255. Set to 255" <<endl;
	maxSeed = 255;
}
  result.u_char_dat=(u_char) 0;

  struct shpt {
	  short x;
	  short y;
	  short z;
  };
  shpt * stack = new shpt[height*width*depth/30]; //create Point Stack
  
  int tempx, tempy, tempz;
  int stacksize=1;
  
  stack[0].x=cx;
  stack[0].y=cy;
  stack[0].z=cz;	
 
  int maxstacksize = height*width*depth/30;
  
 
  while(stacksize != 0){
	  if(stacksize > maxstacksize - 7)  {
      shpt * temp = stack;
	  stack = new shpt[2*maxstacksize];
      for(int i=0; i <maxstacksize; i++){
		  stack[i] = temp[i];
      }
      delete [] temp;
      maxstacksize = 2*maxstacksize;
      //cerr << maxstacksize << endl;
    }
    
    //if the node on the stack has the right color, mark it as part of image, take it off the stack, and then add the pixels around it to the stack if they aren't there already		
	   
    tempx = stack[stacksize-1].x;
    tempy = stack[stacksize-1].y;
    tempz = stack[stacksize-1].z;
    int color = (int) vol.u_char_dat[tempz][tempy][tempx];
    stacksize--;
	
    if((float) result.u_char_dat[tempz][tempy][tempx] != 255 
       && color >= minSeed && color <=maxSeed) {
		result.u_char_dat[tempz][tempy][tempx] = (u_char) 255;
		if (tempx !=0 && 
			(float) result.u_char_dat[tempz][tempy][tempx-1] != 255 ) {
			stack[stacksize].x = tempx-1;
			stack[stacksize].y = tempy;
			stack[stacksize].z = tempz;
			stacksize++;
		}
		if (tempy!= 0 
			&& (float) result.u_char_dat[tempz][tempy-1][tempx] !=255 ) {
			stack[stacksize].x = tempx;
			stack[stacksize].y = tempy-1;
			stack[stacksize].z = tempz;
			stacksize++;
		}
      if (tempz !=0 && 
		  (float) result.u_char_dat[tempz-1][tempy][tempx] !=255 ) {
		  stack[stacksize].x=tempx;
		  stack[stacksize].y=tempy;
		  stack[stacksize].z=tempz-1;
		  stacksize++;
      }
      if (tempx !=width-1 && 
		  (float) result.u_char_dat[tempz][tempy][tempx+1] !=255 ) {
		  stack[stacksize].x = tempx+1;
		  stack[stacksize].y = tempy;
		  stack[stacksize].z = tempz;
		  stacksize++;
      }
      if (tempy != height-1 &&
		  (float) result.u_char_dat[tempz][tempy+1][tempx] !=255 ) {
		  stack[stacksize].x = tempx;
		  stack[stacksize].y = tempy+1;
		  stack[stacksize].z = tempz;
		  stacksize++;
      }
      if ( tempz != depth-1 &&
		   (float) result.u_char_dat[tempz+1][tempy][tempx] !=255 ) {
		  stack[stacksize].x = tempx;
		  stack[stacksize].y = tempy;
		  stack[stacksize].z = tempz+1;
		  stacksize++;
      }
	  
    } //end if right color
    
  } //end while 
	
  delete [] stack;

 
}

ItXECode ItXVolumeUtils::removeSkull(ItXVolume &vol, ItXVolume &result, 
				     int cx, int cy,  int cz, 
				     int toleranceL, int toleranceU) {

  if(vol.dataType() != ItXVolume:: UnsignedChar)
		return (ItXInvalidDataType); 

  ItXVolume temp1(vol.nx, vol.ny, vol.nz,
                   ItXVolume::UnsignedChar) ;

  temp1=vol;
  ItXVolumeUtils::erosionUChar(temp1,4);
  ItXVolumeUtils::RegionGrow(temp1,result,cx,cy,cz, toleranceL,toleranceU);

  ItXVolumeUtils::dilationUChar(result,5);
  ItXVolumeUtils::invertVolume(result);

  cx=vol.nx-2;
  cy=vol.ny-2;
  cz=0;

  ItXVolume slice_in(vol.nx, vol.ny, 1, ItXVolume::UnsignedChar);
  ItXVolume slice_out(vol.nx, vol.ny, 1,ItXVolume::UnsignedChar);

  int i,j,k;
  for(k=0; k<vol.nz; ++k) {
	  for(j=0; j<vol.ny; ++j)
		  for(i=0; i<vol.nx; ++i)
			  slice_in.u_char_dat[0][j][i] = result.u_char_dat[k][j][i];
	  
	  ItXVolumeUtils::RegionGrow(slice_in,slice_out,cx,cy,cz,0,0);
	  for(j=0; j<vol.ny; ++j)
		  for(i=0; i<vol.nx; ++i)  
			temp1.u_char_dat[k][j][i] = slice_out.u_char_dat[0][j][i];
	  
  }

  
  for (k=0; k<vol.nz; ++k)
	  for(j=0; j<vol.ny; ++j)
		  for(i=0; i<vol.nx; ++i) {
			  if ((int) temp1.u_char_dat[k][j][i] == 0)
				  result.u_char_dat[k][j][i] = vol.u_char_dat[k][j][i];
			  else result.u_char_dat[k][j][i] = (u_char) 0;
		  }
  
  return(ItXSuccess);
}

ItXECode ItXVolumeUtils::dilationUChar(ItXVolume &vol, int size ){
  
  if(vol.dataType() != ItXVolume:: UnsignedChar)
		return (ItXInvalidDataType); 

  ItXVolume result(vol.getSizeX(), vol.getSizeY(),  vol.getSizeZ(),
		 ItXVolume::UnsignedChar);

  int i,j,k,ii,jj,kk;
  int xx,yy,zz;
  int mskx,msky,mskz;

  char work_string[100];
  float work_perc = 100.0/ (size);


  // Assumes the data is 0;


  for (int n=0; n<size; n++) {
	  sprintf(work_string,"%d%% Complete", (int)(work_perc * n ));
	  vol.ShowWorking(work_string);
	  for(k=1; k<vol.nz-1; k++) {
		  for(j=1; j<vol.ny-1; j++) {
			  for(i=1; i<vol.nx-1; i++) { 
				  result.u_char_dat[k][j][i] = MAX3(vol.u_char_dat[k][j][i],
													vol.u_char_dat[k][j][i-1],
													vol.u_char_dat[k][j][i+1]);
			  }
		  }
	  }

	  for(k=1; k<vol.nz-1; k++) {
		  for(j=1; j<vol.ny-1; j++) {
			  for(i=1; i<vol.nx-1; i++) { 
				  vol.u_char_dat[k][j][i] = MAX3(result.u_char_dat[k][j][i],
												 result.u_char_dat[k][j-1][i],
												 result.u_char_dat[k][j+1][i]);
			  }
		  }
	  }
	  
	  for(k=1; k<vol.nz-1; k++) {
		  for(j=1; j<vol.ny-1; j++) {
			  for(i=1; i<vol.nx-1; i++) { 
				  result.u_char_dat[k][j][i] = MAX3(vol.u_char_dat[k][j][i],
													vol.u_char_dat[k-1][j][i],
													vol.u_char_dat[k+1][j][i]);
			  }
		  }
	  }
    
	  vol=result;
  }
  return(ItXSuccess);
}

ItXECode ItXVolumeUtils::AND(ItXVolume &vol, ItXVolume &mask) {
	
  if(vol.dataType() != ItXVolume:: UnsignedChar)
		return (ItXInvalidDataType); 
	
	if (vol.getSizeZ() != mask.getSizeZ() || vol.getSizeY() != mask.getSizeY() ||
		vol.getSizeX() != vol.getSizeX() )
		return(ItXError);

	int i,j,k;
	for (k=0; k<vol.getSizeZ(); ++k)
		for (j=0; j<vol.getSizeY(); ++j)
			for ( i=0; i<vol.getSizeX(); ++i) {
				if ((int) mask.u_char_dat[k][j][i] !=255)
					vol.u_char_dat[k][j][i]=(u_char) 0;
			}

	return (ItXSuccess);
}

//function to pad a volume
ItXECode ItXVolumeUtils::pad(ItXVolume &vol,ItXVolume &result,int px1,int py1,int pz1,
			       int px2,int py2,int pz2){

int vx,vy,vz;
int nx,ny,nz;

vx = vol.getSizeX();
vy = vol.getSizeY();
vz = vol.getSizeZ();

nx = vx+px1+px2;
ny = vy+py1+py2;
nz = vz+pz1+pz2;

result.setSize(nx,ny,nz,vol.dataType());
ItXVolumeUtils::clear(result,0.f);

int i,j,k;
int i1,j1,k1;

for (i=0,i1=px1;i<vx;i++,i1++)
 for (j=0,j1=py1;j<vy;j++,j1++)	 
   for (k=0,k1=pz1;k<vz;k++,k1++)
     result.setPointValue(i1,j1,k1,vol.pointValue(i,j,k));

 return(ItXSuccess);
}



ItXECode ItXVolumeUtils::OR(ItXVolume &vol,ItXVolume &mask) {
                        
                                
  if(vol.dataType() != ItXVolume:: UnsignedChar)
		return (ItXInvalidDataType); 

  if(mask.dataType() != ItXVolume:: UnsignedChar)
		return (ItXInvalidDataType); 

   if (vol.getSizeZ() != mask.getSizeZ() ||
        vol.getSizeY() !=mask.getSizeY()   ||
        vol.getSizeX() != vol.getSizeX() )
                return(ItXError);

        int i,j,k;
        for (k=0; k<vol.getSizeZ(); ++k)
                for (j=0; j<vol.getSizeY(); ++j)
                        for ( i=0; i<vol.getSizeX(); ++i) {
                                if ((int) mask.u_char_dat[k][j][i] ==255)

                                vol.u_char_dat[k][j][i]=(u_char)255;
                        }
                
        return (ItXSuccess);
        
}
