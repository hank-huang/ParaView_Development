#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <iostream.h>
#include <TDL/IDLdefines.h>
#include <TDL/AnalyzeVol.h>

// I don't want subroutines to exit from the program [csb:19970729.0923MST]
#if 1
#define throwError(STRING) { cerr << STRING << endl; }
#else
void throwError(char *str) {
  cout << str << endl;
  exit(1);
}
#endif

AnalyzeVol::AnalyzeVol() {
  _p = 0;
  bzero(&_hdr, sizeof(_hdr));
}

AnalyzeVol::~AnalyzeVol() {
  setDim(0, 0, 0);
}

int
AnalyzeVol::setDim(int zsize, int ysize, int xsize) {
  if ( _p ) {
    delete [] _p[0][0];
    delete [] _p[0];
    delete [] _p;
    _p = 0;
    _hdr.dime.dim[1] = 0;
    _hdr.dime.dim[2] = 0;
    _hdr.dime.dim[3] = 0;
  }

  if (xsize < 0 || ysize < 0 || zsize < 0)
    throwError("AnalyzeVol::setDim negative dimensions attempted.");

  if (xsize && ysize && zsize) {        // Check for null array.
    unsigned short **ptr1;
    unsigned short  *ptr2;
    _hdr.dime.dim[1] = xsize;
    _hdr.dime.dim[2] = ysize;
    _hdr.dime.dim[3] = zsize;
    _p   = new unsigned short** [zsize];
    ptr1 = new unsigned short*  [zsize*ysize];
    ptr2 = new unsigned short   [zsize*ysize*xsize];
    if(!_p || !ptr1 || !ptr2)
      throwError("AnalyzeVol::setDim New failed.");
    
    ptr1[0] = ptr2;
    int length = zsize*ysize;
    for(int i=1; i<length; i++)
      ptr1[i] = ptr1[i-1] + xsize;

    _p[0] = ptr1;
    for(int z=1; z<zsize; z++)
      _p[z] = _p[z-1] + ysize;
  }
  return 1;
}

int AnalyzeVol::load(const char *img_pathname, Message *message) {
  int compressed = (img_pathname[strlen(img_pathname)-1] == 'Z');
  char 	hdr_pathname[1024];
  strcpy(hdr_pathname, img_pathname);
  if (compressed) {// take off trailing .Z
    //    hdr_pathname[strlen(img_pathname)-2] = 0;
    hdr_pathname[strlen(hdr_pathname)-3] = 'r';
    hdr_pathname[strlen(hdr_pathname)-4] = 'd';
    hdr_pathname[strlen(hdr_pathname)-5] = 'h';
  } else {
    hdr_pathname[strlen(hdr_pathname)-1] = 'r';
    hdr_pathname[strlen(hdr_pathname)-2] = 'd';
    hdr_pathname[strlen(hdr_pathname)-3] = 'h';
  }

  if (loadHDR(hdr_pathname) == 0)
    return 0;
  if (loadIMG(img_pathname, message) == 0)
    return 0;
  //	printHeader();	
  return 1;
}

int AnalyzeVol::loadHDR(const char *hdr_pathname) {
  FILE 	*fp;
  
  int compressed = (hdr_pathname[strlen(hdr_pathname)-1] == 'Z');

  if (compressed) {
    char cmd[256];
    strcpy(cmd, "zcat ");
    strncat(cmd, hdr_pathname, sizeof(cmd)-strlen(cmd));
    if((fp=popen(cmd, "r")) == NULL) {
      char errorMessage[256];
      sprintf(errorMessage, "unable to open pipe: %s\n", cmd);
      throwError(errorMessage);
      return 0;
    }
  } else 
    if((fp=fopen(hdr_pathname, "r"))==0) {
      char errorMessage[256];
      sprintf(errorMessage, "unable to read: %s\n", hdr_pathname);
      throwError(errorMessage);
      return 0;
    }
  
  
  fread(&_hdr, sizeof(struct adsr), 1, fp);
  
  //do some byte swapping within the _hdr struct
  //take care of the hdr_key struct first.
  _hdr.hk.sizeof_hdr =big_endian_int( _hdr.hk.sizeof_hdr);
  _hdr.hk.extents = big_endian_int(_hdr.hk.extents);
  _hdr.hk.session_error = big_endian_short(_hdr.hk.session_error);
  
  //take care of the img_dimension struct second
  for(int k=0;k<8;k++){
    _hdr.dime.dim[k]=big_endian_short(_hdr.dime.dim[k]);
    _hdr.dime.pixdim[k]=big_endian_float(_hdr.dime.pixdim[k]);
  }
  _hdr.dime.unused1=big_endian_short(_hdr.dime.unused1);
  _hdr.dime.datatype=big_endian_short(_hdr.dime.datatype);
  _hdr.dime.bitpix=big_endian_short(_hdr.dime.bitpix);
  _hdr.dime.dim_un0=big_endian_short(_hdr.dime.dim_un0);
  _hdr.dime.vox_offset=big_endian_float(_hdr.dime.vox_offset);
  _hdr.dime.roi_scale=big_endian_float(_hdr.dime.roi_scale);
  _hdr.dime.funused1=big_endian_float(_hdr.dime.funused1);
  _hdr.dime.funused2=big_endian_float(_hdr.dime.funused2);
  _hdr.dime.cal_max=big_endian_float(_hdr.dime.cal_max);
  _hdr.dime.cal_min=big_endian_float(_hdr.dime.cal_min);
  _hdr.dime.compressed=big_endian_int(_hdr.dime.compressed);

  _hdr.dime.verified=big_endian_int(_hdr.dime.verified);
  _hdr.dime.glmax=big_endian_int(_hdr.dime.glmax);
  _hdr.dime.glmin=big_endian_int(_hdr.dime.glmin);
  //take care of byte swapping the data_hist struct
  _hdr.hist.views = big_endian_int(_hdr.hist.views);
  _hdr.hist.vols_added = big_endian_int(_hdr.hist.vols_added);
  _hdr.hist.start_field= big_endian_int(_hdr.hist.start_field);
  _hdr.hist.field_skip = big_endian_int(_hdr.hist.field_skip);
  _hdr.hist.omax = big_endian_int(_hdr.hist.omax);
  _hdr.hist.omin = big_endian_int(_hdr.hist.omin);
  _hdr.hist.smax = big_endian_int(_hdr.hist.smax);
  _hdr.hist.smin = big_endian_int(_hdr.hist.smin);

if (compressed)
    pclose(fp);
  else
    fclose(fp);
  return 1;
}

int AnalyzeVol::loadIMG(const char *img_pathname, Message *message) {
  int		i;
  FILE	*fp;
  
  if (message)
    message->showMessageWindow("Loading");
  if (message)
    message->percentComplete(0.0);
  
  setDim(getZsize(), getYsize(), getXsize());

  int compressed = (img_pathname[strlen(img_pathname)-1] == 'Z');

  if (compressed) {
    char cmd[256];
    strcpy(cmd, "zcat ");
    strncat(cmd, img_pathname, sizeof(cmd)-strlen(cmd));
    if((fp=popen(cmd, "r")) == NULL) {
      char errorMessage[256];
      sprintf(errorMessage, "can't open pipe: %s\n", cmd);
      throwError(errorMessage);
      //printf;
      return 0;
    }
  } else 
    if((fp=fopen(img_pathname, "r"))==0) {
      char errorMessage[256];
      sprintf(errorMessage, "can't open: %s\n", img_pathname);
      throwError(errorMessage);
      //printf("can't open: %s\n", img_pathname);
      return 0;
    }
  
  unsigned short* ptr; //JAVA needs 16 bits
  if (_hdr.dime.bitpix == 8) {
    int             n;
    int             length = getXsize() * getYsize(); 

    // changed by kem -- march 5, 1997 from char to unsigned char
    unsigned char*           buffer = new unsigned char[length];
    
    for (i=0; i<getZsize(); i++) {
      fread(buffer, 1*length, 1, fp);
      ptr = &_p[i][0][0];
      for (n=0; n<length; n++) {
	ptr[n] = buffer[n]; //NO need to SWAP!
      }
      if (message)
	message->percentComplete((float)i/(float)(getZsize()-1));
    }
    delete [] buffer;
  } else if (_hdr.dime.bitpix == 16) {
    int n;
    int length = getXsize() * getYsize();
    unsigned short *buffer = new unsigned short[length];

    for (i=0; i<getZsize(); i++) {
      fread(buffer, 2, length, fp);
      ptr = &_p[i][0][0];
      for (n=0; n<length; n++)
  	ptr[n] = big_endian_short(buffer[n]);
      //ptr = &_p[i][0][0];
      //fread(ptr, 2, length, fp);
      //fread(&_p[i][0][0], getXsize()*getYsize()*2, 1, fp);
      if (message)
	message->percentComplete((float)i/(float)(getZsize()-1));
    }
    delete [] buffer;
  } else {
// I don't want this subroutine to exit the program [csb:19970729.0923MST]
#if 1
    throwError("error: can't read pixel data.");
#else
    printf(".\n");
    exit(1);
#endif
  }

  if (compressed)
    pclose(fp);
  else
    fclose(fp);
  
  return 1;
}

int AnalyzeVol::save(const char *img_pathname) {
  char 	hdr_pathname[1024];

  sprintf(hdr_pathname, "%s", img_pathname);
  hdr_pathname[strlen(hdr_pathname)-1] = 'r';
  hdr_pathname[strlen(hdr_pathname)-2] = 'd';
  hdr_pathname[strlen(hdr_pathname)-3] = 'h';
  
  if (saveHDR(hdr_pathname) == 0)
    return 0;
  if (saveIMG(img_pathname) == 0)
    return 0;
  
  return 1;
}

int AnalyzeVol::saveHDR(const char *hdr_pathname) {
  FILE *fp;
  
  _hdr.dime.datatype = DT_SIGNED_SHORT;
  _hdr.dime.bitpix = 16;
  
  _hdr.dime.dim[0] = 4;
  _hdr.dime.dim[4] = 1;
  _hdr.hk.regular = 'r';
  _hdr.hk.sizeof_hdr = sizeof(struct adsr);
  
  _hdr.dime.roi_scale = 1.0;
  _hdr.dime.glmin = 0;
  //_hdr.dime.glmax = 0;

  _hdr.dime.vox_offset  = 0.0;
  _hdr.dime.roi_scale   = 1.0;
  _hdr.dime.funused1    = 0.0;
  _hdr.dime.funused2    = 0.0;
  _hdr.dime.cal_max     = 0.0;
  _hdr.dime.cal_min     = 0.0;
  
  strcpy(_hdr.dime.vox_units, " ");
  strcpy(_hdr.dime.cal_units, " ");
  
  /*   Planar Orientation;    */
  /*   Movie flag OFF: 0 = transverse, 1 = coronal, 2 = sagittal
       Movie flag ON:  3 = transverse, 4 = coronal, 5 = sagittal  */
  
  _hdr.hist.orient     = 0;
  
  if((fp=fopen(hdr_pathname, "w"))==0) {
    char errorMessage[256];
    sprintf(errorMessage, "unable to save: %s\n", hdr_pathname);
    throwError(errorMessage);
    return 0;
  }
  //swapping time again
  adsr temp = _hdr;

  //take care of the hdr_key struct first.
  temp.hk.sizeof_hdr =big_endian_int( _hdr.hk.sizeof_hdr);
  temp.hk.extents = big_endian_int(_hdr.hk.extents);
  temp.hk.session_error = big_endian_short(_hdr.hk.session_error);
  
  //take care of the img_dimension struct second
  for(int k=0;k<8;k++){
    temp.dime.dim[k]=big_endian_short(_hdr.dime.dim[k]);
    temp.dime.pixdim[k]=big_endian_float(_hdr.dime.pixdim[k]);
  }
  temp.dime.unused1=big_endian_short(_hdr.dime.unused1);
  temp.dime.datatype=big_endian_short(_hdr.dime.datatype);
  temp.dime.bitpix=big_endian_short(_hdr.dime.bitpix);
  temp.dime.dim_un0=big_endian_short(_hdr.dime.dim_un0);
  temp.dime.vox_offset=big_endian_float(_hdr.dime.vox_offset);
  temp.dime.roi_scale=big_endian_float(_hdr.dime.roi_scale);
  temp.dime.funused1=big_endian_float(_hdr.dime.funused1);
  temp.dime.funused2=big_endian_float(_hdr.dime.funused2);
  temp.dime.cal_max=big_endian_float(_hdr.dime.cal_max);
  temp.dime.cal_min=big_endian_float(_hdr.dime.cal_min);
  temp.dime.compressed=big_endian_int(_hdr.dime.compressed);

  temp.dime.verified=big_endian_int(_hdr.dime.verified);
  temp.dime.glmax=big_endian_int(_hdr.dime.glmax);
  temp.dime.glmin=big_endian_int(_hdr.dime.glmin);
  //take care of byte swapping the data_hist struct
  temp.hist.views = big_endian_int(_hdr.hist.views);
  temp.hist.vols_added = big_endian_int(_hdr.hist.vols_added);
  temp.hist.start_field= big_endian_int(_hdr.hist.start_field);
  temp.hist.field_skip = big_endian_int(_hdr.hist.field_skip);
  temp.hist.omax = big_endian_int(_hdr.hist.omax);
  temp.hist.omin = big_endian_int(_hdr.hist.omin);
  temp.hist.smax = big_endian_int(_hdr.hist.smax);
  temp.hist.smin = big_endian_int(_hdr.hist.smin);

  

  fwrite(&temp, sizeof(struct adsr), 1, fp);
  fclose(fp);
  
  return 1;
}

int AnalyzeVol::saveIMG(const char *img_pathname) {
  FILE *fp;	

  if((fp=fopen(img_pathname, "w"))==0) {
    char errorMessage[256];
    sprintf(errorMessage, "unable to save: %s\n", img_pathname);
    throwError(errorMessage);
    return 0;
  }

  //  fwrite(data(), 2*getNelm(), 1, fp); changed this for byteswapping
  unsigned short swapshort;
  for(int g=0; g<getNelm();g++)
	{
	swapshort = big_endian_short(*(data()+g));
	fwrite(&swapshort, sizeof(swapshort),1,fp);
	}

fclose(fp);
  
  return 1;
}

void AnalyzeVol::printHeader(void) {
  // header key
  printf("data_type: %s\n", _hdr.hk.data_type); 
  printf("db_name: %s\n", _hdr.hk.db_name); 
  
  // image dimensions 
  printf("bitpix: %3d\n", _hdr.dime.bitpix);
  printf("xsize: %3d\n", _hdr.dime.dim[1]); 
  printf("ysize: %3d\n", _hdr.dime.dim[2]); 
  printf("zsize: %3d\n", _hdr.dime.dim[3]); 
  printf("xscale: %.6f\n", _hdr.dime.pixdim[1]); 
  printf("yscale: %.6f\n", _hdr.dime.pixdim[2]); 
  printf("zscale: %.6f\n", _hdr.dime.pixdim[3]); 
  printf("vox_units: %s\n", _hdr.dime.vox_units);
  printf("cal_units: %s\n", _hdr.dime.cal_units);
  printf("datatype: %d\n", _hdr.dime.datatype);
  printf("vox_offset: %f\n", _hdr.dime.vox_offset);
  printf("roi_scale: %f\n", _hdr.dime.roi_scale);
  printf("cal_min: %f\n", _hdr.dime.cal_min);
  printf("cal_max: %f\n", _hdr.dime.cal_max);
  printf("compressed: %d\n", _hdr.dime.compressed);
  printf("verified: %d\n", _hdr.dime.verified);
  
  // data history 
  printf("descrip: %s\n", _hdr.hist.descrip); 
  printf("orient: %c\n", _hdr.hist.orient); 
  printf("originator: %s\n", _hdr.hist.originator); 
  printf("scannum: %s\n", _hdr.hist.scannum); 
  printf("patient_id: %s\n", _hdr.hist.patient_id); 
  printf("exp_date: %s\n", _hdr.hist.exp_date); 
  printf("exp_time: %s\n", _hdr.hist.exp_time); 
  printf("hist_un0: %s\n", _hdr.hist.hist_un0); 
  printf("views: %d\n", _hdr.hist.views); 
  printf("vols_added: %d\n", _hdr.hist.vols_added); 
  printf("start_field: %d\n", _hdr.hist.start_field); 
  printf("field_skip: %d\n", _hdr.hist.field_skip); 
  printf("omin: %d\n", _hdr.hist.omin); 
  printf("omax: %d\n", _hdr.hist.omax); 
  printf("smin: %d\n", _hdr.hist.smin); 
  printf("smax: %d\n", _hdr.hist.smax); 
}










