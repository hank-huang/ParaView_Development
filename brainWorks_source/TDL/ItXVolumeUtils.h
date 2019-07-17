#ifndef _ImageUtils_h_
#define _ImageUtils_h_
#include <OS.h>
#include <TDL/ItXVolume.h>
//////////////////////////////////////////////////////////////////////////
//
// File: ItXVolumeUtils.h
//
// Author: Keith Doolittle
/////////////////////////////////////////////////////////////////
class ItXVolumeUtils {
public:
        //---------------------------------------------------------
        // Apply generic mask to volume
        //---------------------------------------------------------
	static ItXECode applyMask(ItXVolume &vol, float *msk, 
			int mx,int my, int mz);

        //---------------------------------------------------------
        // Blur data
        //---------------------------------------------------------
        static ItXECode blur(ItXVolume&);

        //---------------------------------------------------------
        // Sharpen data
        //---------------------------------------------------------
        static ItXECode sharpen(ItXVolume&);

	static void zeroMemory(char *ptr, int l) {
                int   i;
                char *cptr = ptr;
		for(i=0;i<l;i++,cptr++) 
			*cptr = 0U; 
	}

	inline static void zeroMemory(u_char *ptr, int l) {
		zeroMemory((char*)ptr,l);
	}

        //---------------------------------------------------------
        // Linearly rescale volume values
        //---------------------------------------------------------
        static ItXECode linearScale(ItXVolume&,float bmin, float bmax);

        //---------------------------------------------------------
        // Array resize
        // If keep_old_data != 0 then original data is padded
        // or cropped evenly from center of volume
        //---------------------------------------------------------
        static ItXECode resize(ItXVolume&, int new_nx, int new_ny, int new_nz,
                        int keep_old_data = 0);

	static ItXECode pad(ItXVolume &vol,ItXVolume &result,int px1,int py1,int pz1,
			    int px2,int py2,int pz2);

        //---------------------------------------------------------
        // resample
        //      shrink/stretch image to new dimensions
        //---------------------------------------------------------
        static ItXECode resample(ItXVolume&, int newnx, int newny, int newnz);
        static ItXECode resample(ItXVolume&, float scalex, float scaley, float scalez);


       //---------------------------------------------------------
        // Convole image with some function
        //---------------------------------------------------------
	static ItXECode convolute(ItXVolume&,char *type, float para);
		
        //---------------------------------------------------------
        // Point by point voxel scaling
        //---------------------------------------------------------
        static void scale(ItXVolume&,float s);

	static const char *filenameTail(const char *fname);

        //
        // Trilinear Interpolation:
        //
        // if img->dataType() == RGB fr,fg,fb are interpolated
        //                       r,g,b values, else only 'fr' is valid
        //
	static void trilinearInterp(const ItXVolume &img,
				     float x, float y, float z,
                                     float &fr, float &fg, float &fb);
	
	static void keepOneValue(ItXVolume &img, float val);

	static ItXECode extract(ItXVolume &in, ItXVolume &out,
				int startx, int starty, int startz);

	static void clear(ItXVolume &in, float value = 0.0);
	
	static ItXECode globalThreshUChar(ItXVolume &vol, int thresh, 
				      ItXVolume &result);
	static ItXECode invertVolume(ItXVolume &vol);
	static ItXECode erosionUChar(ItXVolume &vol, int size);
	static ItXECode dilationUChar(ItXVolume &vol, int size);
	static ItXECode erosionSlowUChar(ItXVolume &vol, int dx, 
                                   int dy, int dz,  ItXVolume &result);
	static ItXECode dilationSlowUChar(ItXVolume &vol, int dx, 
                                   int dy, int dz, 
					  ItXVolume &result);
	
	static ItXECode openingUChar(ItXVolume &vol, int dx);

	static ItXECode RegionGrow(ItXVolume &vol, ItXVolume &result, 
				   int cx,int cy, int cz, 
                                   int toleranceL,int toleranceU);
	
	static ItXECode removeSkull(ItXVolume &vol, ItXVolume &result,
				    int cx, int cy, int cz, 
				    int  toleranceL, int toleranceU);
	static ItXECode AND(ItXVolume &vol, ItXVolume &mask);
	static ItXECode OR(ItXVolume &vol, ItXVolume &mask);

	static u_char  trilinearInterp(const Array3D<u_char> &dat,
				       float x, float y, float z);

	static u_short trilinearInterp(const Array3D<u_short> &dat,
				float x, float y, float z);

	static short   trilinearInterp(const Array3D<short> &dat,
				float x, float y, float z);

	static float   trilinearInterp(const Array3D<float> &dat,
				float x, float y, float z);

	static float   trilinearInterp(const Array3D<double> &dat,
				float x, float y, float z);

	static int     trilinearInterp(const Array3D<int> &dat,
				float x, float y, float z);

	static void    trilinearInterp(const Array3D<u_char> &dat,
				float x, float y, float z,
                                float &fr, float &fg, float &fb);

private:

	static ItXECode resampleUChar(ItXVolume&,int, int, int);
        static ItXECode resampleUShort(ItXVolume&,int, int, int);
        static ItXECode resampleShort(ItXVolume&,int, int, int);
        static ItXECode resampleFloat(ItXVolume&,int, int, int);
        static ItXECode resampleInt(ItXVolume&,int, int, int);

	static ItXECode convoluteUChar(ItXVolume&,char *,float);
        static ItXECode convoluteUShort(ItXVolume&,char *,float);
        static ItXECode convoluteShort(ItXVolume&,char *,float);
        static ItXECode convoluteFloat(ItXVolume&,char *,float);
        static ItXECode convoluteInt(ItXVolume&,char *,float);
		
        static void     linearScaleUShort(ItXVolume&, u_short, u_short);
        static void     linearScaleShort(ItXVolume&, short, short);
        static void     linearScaleUChar(ItXVolume&, u_char, u_char);
        static void     linearScaleFloat(ItXVolume&, float, float);
        static void     linearScaleInt(ItXVolume&, int, int);

        static ItXECode applyMaskUChar(ItXVolume&, float*, int, int, int);
        static ItXECode applyMaskUShort(ItXVolume&, float*, int, int, int);
        static ItXECode applyMaskShort(ItXVolume&, float*, int, int, int);
        static ItXECode applyMaskFloat(ItXVolume&, float*, int, int, int);
        static ItXECode applyMaskInt(ItXVolume&, float*, int, int, int);

        static ItXECode resizeUChar(ItXVolume&, int, int, int);
        static ItXECode resizeUShort(ItXVolume&, int, int, int);
        static ItXECode resizeShort(ItXVolume&, int, int, int);
        static ItXECode resizeFloat(ItXVolume&, int, int, int);
        static ItXECode resizeInt(ItXVolume&, int, int, int);
};

#endif





