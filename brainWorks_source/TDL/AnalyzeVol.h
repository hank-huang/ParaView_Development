#ifndef __ANALYZEVOL_H__
#define __ANALYZEVOL_H__
#include <OS.h>
#include <string.h>
#include <strings.h>

#include <ADL/Message.h>
#include <TDL/IDLdefines.h>

/*
 * 
 * (c) Copyright, 1986-1995
 * Biomedical Imaging Resource
 * Mayo Foundation
 *
 * dbh.h
 *
 *
 * database sub-definitions
 */

/* Acceptable values for hdr.dime.datatype */
#define DT_NONE				0
#define DT_UNKNOWN			0
#define DT_BINARY			1
#define DT_UNSIGNED_CHAR		2
#define DT_SIGNED_SHORT			4
#define DT_SIGNED_INT			8
#define DT_FLOAT			16
#define DT_COMPLEX			32
#define DT_DOUBLE			64
#define DT_RGB				128
#define DT_ALL				255

struct hdr_key {          		/*      header_key       */
  int sizeof_hdr;              		/* 0 + 4     */
  char data_type[10];          		/* 4 + 10    */
  char db_name[18];            		/* 14 + 18   */
  int extents;                 		/* 32 + 4    */
  short int session_error;     		/* 36 + 2    */
  char regular;                		/* 38 + 1    */
  char hkey_un0;               		/* 39 + 1    */
};           				/* total=40  */
	
struct img_dimension {     		/*      image_dimension  */
  short int dim[8];           		/* 0 + 16    */
  char vox_units[4];			/* 16 + 4    */
  char cal_units[8];			/* 20 + 4    */
  short int unused1;			/* 24 + 2    */
  short int datatype;			/* 30 + 2    */
  short int bitpix;                   	/* 32 + 2    */
  short int dim_un0;                    /* 34 + 2    */
  float pixdim[8];                      /* 36 + 32   */
  /* 
     pixdim[] specifies the voxel dimensions:
     pixdim[1] - voxel width
     pixdim[2] - voxel height
     pixdim[3] - interslice distance
     ..etc
     */
  float vox_offset;                    	/* 68 + 4    */
  float roi_scale;                      /* 72 + 4    */
  float funused1;                      	/* 76 + 4    */
  float funused2;                      	/* 80 + 4    */
  float cal_max;                    	/* 84 + 4    */
  float cal_min;                     	/* 88 + 4    */
  int compressed;                     	/* 92 + 4    */
  int verified;                     	/* 96 + 4    */
  int glmax, glmin;                 	/* 100 + 8   */
};          				/* total=108 */

struct data_hist {         		/*      data_history     */
  char descrip[80];                	/* 0 + 80    */
  char aux_file[24];               	/* 80 + 24   */
  char orient;                     	/* 104 + 1   */
  char originator[10];             	/* 105 + 10  */
  char generated[10];              	/* 115 + 10  */
  char scannum[10];                	/* 125 + 10  */
  char patient_id[10];             	/* 135 + 10  */
  char exp_date[10];               	/* 145 + 10  */
  char exp_time[10];               	/* 155 + 10  */
  char hist_un0[3];                	/* 165 + 3   */
  int views;                       	/* 168 + 4   */
  int vols_added;                  	/* 172 + 4   */
  int start_field;                 	/* 176 + 4   */
  int field_skip;                  	/* 180 + 4   */
  int omax,omin;                   	/* 184 + 8   */
  int smax,smin;                   	/* 192 + 8   */
};                     			/* total=200 */

struct adsr {                		/*      dsr              */
  struct hdr_key hk;          		/* 0 + 40    */
  struct img_dimension dime;   		/* 40 + 108  */
  struct data_hist hist;      		/* 148 + 200 */
};                     			/* total=348 */


class AnalyzeVol {
private:
  struct adsr _hdr;
  unsigned short ***_p;

public:
  AnalyzeVol();	        // Ctor
  ~AnalyzeVol();	// Dtor

  int load(const char *, Message *message=NULL);
  int loadHDR(const char *);
  int loadIMG(const char *, Message *message=NULL);

  int save(const char *);
  int saveHDR(const char *);
  int saveIMG(const char *);

  int setDim(int zsize, int ysize, int xsize);
  void printHeader(void);
  
  inline unsigned short ***address(void) { return _p; }
  inline unsigned short *data(int z=0, int y=0, int x=0) { 
    return &_p[z][y][x];
  }

  // accessor functions
  inline int getXsize(void) { return _hdr.dime.dim[1]; }
  inline int getYsize(void) { return _hdr.dime.dim[2]; }
  inline int getZsize(void) { return _hdr.dime.dim[3]; }
  inline int getNelm(void) { return getXsize()*getYsize()*getZsize(); } 
  inline float getVoxelXScale(void) { return _hdr.dime.pixdim[1]; }
  inline float getVoxelYScale(void) { return _hdr.dime.pixdim[2]; }
  inline float getVoxelZScale(void) { return _hdr.dime.pixdim[3]; }
  inline char *getPatientName(void) { return _hdr.hk.db_name; }
  inline int isEmpty(void) { return !getNelm(); }
  inline int getGlmax() { return _hdr.dime.glmax; }
  inline int getGlmin() { return _hdr.dime.glmin; }
  inline void clear() { bzero(data(), 2*getNelm()); }
  inline int getSmin() { return _hdr.hist.smin; }
  inline int getSmax() { return _hdr.hist.smax; }
  inline int getOmin() { return _hdr.hist.omin; }
  inline int getOmax() { return _hdr.hist.omax; }

  // mutator functions
  inline void setXscale(float xscale) { _hdr.dime.pixdim[1] = xscale; }
  inline void setYscale(float yscale) { _hdr.dime.pixdim[2] = yscale; } 
  inline void setZscale(float zscale) { _hdr.dime.pixdim[3] = zscale; } 
  inline void setPatientName(char *s) { strncpy(_hdr.hk.db_name, s, 18); }
  inline void setGlmax(int i) { _hdr.dime.glmax = i; }
  inline void setGlmin(int i) { _hdr.dime.glmin = i; }
  inline void  setSmin(int i) { _hdr.hist.smin = i; }
  inline void  setSmax(int i) { _hdr.hist.smax = i; }
  inline void  setOmin(int i) { _hdr.hist.omin = i; }
  inline void  setOmax(int i) { _hdr.hist.omax = i; }
};

#endif
