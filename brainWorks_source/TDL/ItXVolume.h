#ifndef _ItXVolume_h_
#define _ItXVolume_h_
///////////////////////////////////////////////////////////////////////////
//
// File: ItXVolume.h
//
// Author: Keith Doolittle
//
// Purpose: ItXVolume Class Definition
//  This class is a base class for representation of image volumes.  
//  It includes routines for reading/writing volumes of different
//  types (ie Analyze), different data types, and contains basic 
//  scaling/resizing routines.
///////////////////////////////////////////////////////////////////////////
#include <OS.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/param.h>
#include <ADL/Array3D.h>
#include <ADL/Array2D.h>
#include <TDL/IDLdefines.h>
#include <TDL/ItXWorkClass.h>
#include <TDL/Point.h>

class ItXVolume : public ItXWorkClass {
 public:
  friend class ItXVolumeUtils;
  friend class ItXVolumeDev;
  friend class XFormApply;
  friend class LoadedVolume;

  //---------------------------------------------------------
  // The type of file to load/save the volume as
  // (in addition to a possible '.itx' header)
  //---------------------------------------------------------
  enum ItXVolumeType   { RawVolume, AnalyzeVolume };

  //---------------------------------------------------------
  // The data type
  //---------------------------------------------------------
  enum ItXDataType    { UnsignedChar, UnsignedShort, Float, Short, RGB,
                        Int, UnknownDataType };

  enum ItXDataOrientation{  ItXTransverseUnflipped,
                            ItXCoronalUnflipped,
                            ItXSagittalUnflipped,
                            ItXTransverseFlipped,
                            ItXCoronalFlipped,
                            ItXSagittalFlipped,
                            ItxUnknownOrientationType };

  enum ItXInterpolation {
    Bilinear, NearestNeighbor
  };

  static int minZoom; 
  static int maxZoom; 

  //=========================================================
  //============ Constructors ===============================
  //=========================================================

  //---------------------------------------------------------
  // Basic constructor, zero size array
  //---------------------------------------------------------
  ItXVolume();
  ItXVolume(ItXDataType tpe);

  //---------------------------------------------------------
  // Constructor where image is loaded from file
  // (dimensions and data type are derived from header(s))
  //---------------------------------------------------------
  ItXVolume(const char *file, ItXVolumeType tpe);

  //---------------------------------------------------------
  // Constructor for loading a raw volume (IXRawVolume)
  //---------------------------------------------------------
  ItXVolume(const char *file, int nx, int ny, int nz, 
	    ItXDataType tpe);

  //---------------------------------------------------------
  // Constructor for creating a blank volume
  //---------------------------------------------------------
  ItXVolume(int xdim, int ydim, int zdim, 
	    ItXDataType tpe);

  void setSize(int xdim, int ydim, int zdim, ItXDataType tpe);

  //=========================================================
  //============ Load/Save routines =========================
  //=========================================================
  ItXECode load(const char *filename, ItXVolumeType tpe);
  ItXECode load(const char *filename, int nx, int ny, int nz, 
		float pixel_dx, float pixel_dy, float pixel_dz,	
		ItXDataType tpe);

  //---------------------------------------------------------
  // Save method.  Uses type for writing appropriate
  // headers (with the possibility of additional '.itx'
  // header for extra params, history, etc).
  //---------------------------------------------------------
  ItXECode save(const char *filename, ItXVolumeType itpe,
		ItXDataType tpe = ItXVolume::UnknownDataType);

  // assignment operator
  ItXVolume & operator= (ItXVolume const &I);

  ItXVolume & operator= (float constval);
  //---------------------------------------------------------
  // The following operators will work in float to avoid
  // under/overflow, then linearly scaled (if needed) back to
  // original datatype
  //---------------------------------------------------------

  // addition operator
  ItXVolume & operator+= (ItXVolume const &I);

  // subtraction operator
  ItXVolume & operator-= (ItXVolume const &I);

  // multiplication operator
  ItXVolume & operator*= (ItXVolume const &I); 
  ItXVolume & operator*= (float const f);
	
  // division operator
  ItXVolume & operator/= (ItXVolume const &I); 
  ItXVolume & operator/= (float const f);


  //---------------------------------------------------------
  // Image data (one will have zero dimension)
  //---------------------------------------------------------
  const Array3D<u_char>  &u_char_data()  const { return(u_char_dat);  }
  const Array3D<u_short> &u_short_data() const { return(u_short_dat); }
  const Array3D<short>   &short_data()   const { return(short_dat);   }
  const Array3D<float>   &float_data()   const { return(float_dat);   }
  const Array3D<int>     &int_data()     const { return(int_dat);     }

  Array3D<u_char>  &u_char_data()         { return(u_char_dat);  }
  Array3D<u_short> &u_short_data()        { return(u_short_dat); }
  Array3D<short>   &short_data()          { return(short_dat);   }
  Array3D<float>   &float_data()          { return(float_dat);   }
  Array3D<int>     &int_data()            { return(int_dat);     }

  void setFilename(const char*);

  //---------------------------------------------------------
  // Dimensions of image data
  //---------------------------------------------------------
  inline int getSizeX() const	{ return(nx); }
  inline int getSizeY() const	{ return(ny); }
  inline int getSizeZ() const	{ return(nz); }

  inline const char *Filename()   { return(filename); }
  static const char *dataTypeName(ItXDataType);

  //---------------------------------------------------------
  // Sub-volume offsets in millimeters
  //---------------------------------------------------------
  inline float getOffsetX()	{ return(offset_x);      }
  inline float getOffsetY()	{ return(offset_y);      }
  inline float getOffsetZ()	{ return(offset_z);      }

  void setOffsetX(float f) 	{ offset_x = f; }
  void setOffsetY(float f) 	{ offset_y = f; }
  void setOffsetZ(float f) 	{ offset_z = f; }

  //---------------------------------------------------------
  // Pixel size
  //---------------------------------------------------------
  inline float getPixelDimensionX()       { return(pixel_dimensionx); }
  inline float getPixelDimensionY()       { return(pixel_dimensiony); }
  inline float getPixelDimensionZ()       { return(pixel_dimensionz); }

  void setPixelDimensionX(float f)	
    { if(f > 0.0) pixel_dimensionx = f; 
    else cerr << "ItXVolume::setPixelDimensionX() invalid value" 
	      << endl;
    }

  void setPixelDimensionY(float f)	
    { if(f > 0.0) pixel_dimensiony = f; 
    else cerr << "ItXVolume::setPixelDimensionY() invalid value" 
	      << endl;
    }

  void setPixelDimensionZ(float f)	
    { if(f > 0.0) pixel_dimensionz = f; 
    else cerr << "ItXVolume::setPixelDimensionZ() invalid value" 
	      << endl;
    }

  //
  // get voxel value at coordinates (x,y,z)
  //
  void  setPointValue(int x, int y, int z, float f, int rgb=0);
  float pointValue(int x, int y, int z, int rgb=0);

  void  getMinMax(float &minret, float &maxret);

  //---------------------------------------------------------
  // Data type 
  //---------------------------------------------------------
  inline ItXDataType     dataType()    const { return(data_type); }

  ItXDataOrientation getPlanarOrientation() const { return planarOrientation; }
  
  //---------------------------------------------------------
  // Datatype conversion routines
  //---------------------------------------------------------
  ItXECode convertToUnsignedChar();
  ItXECode convertToUnsignedShort();
  ItXECode convertToShort();
  ItXECode convertToFloat();

  //---------------------------------------------------------
  // Get/set center point routines
  //---------------------------------------------------------
  void setCenter(float  x, float  y, float  z);
  void getCenter(float &x, float &y, float &z);

  //---------------------------------------------------------

 protected:
  Array3D<u_char>  u_char_dat;	
  Array3D<u_short> u_short_dat;	
  Array3D<short>   short_dat;	
  Array3D<float>   float_dat;
  Array3D<int>     int_dat;

  ItXDataType    data_type;

  //---------------------------------------------------------
  // Volume data.  Data type determines valid
  // array, others will have zero dimension
  //---------------------------------------------------------
  char filename[MAXPATHLEN];

  //---------------------------------------------------------
  // Convenience access to volume dimensions
  // instead of using Array3D 
  //---------------------------------------------------------
  int   nx,ny,nz;

  float pixel_dimensionx;
  float pixel_dimensiony;
  float pixel_dimensionz;

  float offset_x;
  float offset_y;
  float offset_z;

  ItXDataOrientation planarOrientation;
  
  //---------------------------------------------------------
  // Current center point of volume
  //---------------------------------------------------------
  float center_x;
  float center_y;
  float center_z;
  
  void     setDefaults();
  bool	   swapAnalyzeHeader(void*);
  ItXECode loadData(Array3D<u_char>&,istream&,bool = false);
  ItXECode loadDataRGB(Array3D<u_char>&,istream&,bool = false);
  ItXECode loadData(Array3D<u_short>&,istream&,bool = false);
  ItXECode loadData(Array3D<short>&,istream&,bool = false);
  ItXECode loadData(Array3D<float>&,istream&,bool = false);
  ItXECode loadData(Array3D<int>&,istream&,bool = false);
  ItXECode saveData(Array3D<u_char>&,ostream&);
  ItXECode saveData(Array3D<u_short>&,ostream&);
  ItXECode saveData(Array3D<short>&,ostream&);
  ItXECode saveData(Array3D<float>&,ostream&);
  ItXECode saveData(Array3D<int>&,ostream&);
  ItXECode saveDataRGB(Array3D<u_char>&,ostream&);
  ItXECode loadRaw();
  ItXECode saveRaw(const char *file);
  ItXECode resampleUChar(int newnx, int newny, int newnz);
  ItXECode resampleUShort(int newnx, int newny, int newnz);
  void     linearScaleUShort(u_short bmin, u_short bmax);
  void     linearScaleUChar(u_char bmin, u_char bmax);
  ItXECode applyMaskUChar(float *mask, int mx, int my, int mz);
  ItXECode applyMaskUShort(float *mask, int mx, int my, int mz);
  ItXECode resizeUChar(int new_nx, int new_ny, int new_nz);
  ItXECode resizeUShort(int new_nx, int new_ny, int new_nz);
  void	   flipData();
  void     flipDataUnsignedChar();
  void     flipDataUnsignedShort();
  void     flipDataShort();
  void     flipDataFloat();
  void     flipDataInt();
  void     flipDataRGB();

  //
  // Analyze file format routines
  //
  ItXECode loadAnalyze();
  ItXECode saveAnalyze(const char *file);
  void constructAnalyzeNames(const char *filename,
			     char *prefix,
			     char *ret_hdrname,
			     char *ret_imgname);
  ItXECode  readAnalyzeHeader(const char *hdrname, bool &need_swap);
  ItXECode  saveAnalyzeHeader(const char *hdrname);
};


#endif





