#ifndef _Landmark_h_
#define _Landmark_h_

#include <OS.h>
#include <sys/param.h>
#include <ADL/Array2D.h>
#include <ADL/Matrix.h>
#include <TDL/IDLdefines.h>
#include <TDL/Point.h>

#define MAXLANDMARKNAMELEN	100
#define MAXWARPVERSIONLEN  10

class Landmark;

class Landmark {
protected:
	int	 n,nalloced;
	char    filename[MAXPATHLEN];
	double  pixeldx,pixeldy,pixeldz;
	int     imagenx,imageny,imagenz;

	Point	*points;
	char   **names;
	int     *valid;

	ItXECode loadBrainWorks(FILE *fp);
	ItXECode loadStandard(FILE *fp);
	ItXECode loadAnalyze(FILE *fp);
	ItXECode loadWarp(FILE *fp);
	ItXECode loadAnalyzeAVW(FILE *fp);

    char     warpVersion[MAXWARPVERSIONLEN];
	ItXECode saveStandard(FILE *fp, int i);
	ItXECode saveWarp(FILE *fp, int i);
	ItXECode readWarp(FILE *fp);
	ItXECode readBrainWorks(FILE *fp);
	ItXECode readStandard(FILE *fp);
public:
	Landmark(int nx=1, int ny=1, int nz=1, 
		 double dx=1.0, double dy=1.0, double dz=1.0);
	~Landmark();

	ItXECode load(char *fname, int nx=1, int ny=1, int nz=1,
			double pixdx=1.0, double pixdy=1.0, double pixdz=1.0);
	ItXECode save(char *fname);

	//
	// fill Matix/Array with point values, scaling back to voxels
	// with pixel dimensions if needed
	//
	void fillArray(Array2D<double> &d,
		double pixelx = 1.0, double pixely = 1.0, double pixelz = 1.0);
	void fillMatrix(Matrix<double> &d,
		double pixelx = 1.0, double pixely = 1.0, double pixelz = 1.0);

	inline int    num() 		 { return(n); }
	inline int    imageDimensionX()  { return(imagenx); }
	inline int    imageDimensionY()  { return(imageny); }
	inline int    imageDimensionZ()  { return(imagenz); }
	inline double pixelDimensionX()  { return(pixeldx); }
	inline double pixelDimensionY()  { return(pixeldy); }
	inline double pixelDimensionZ()  { return(pixeldz); }


	void setDimensions(int x, int y, int z) {
		imagenx = x;	
		imageny = y;	
		imagenz = z;	
		}
	//

	void setResolution(float x, float y, float z) {
		pixeldx = x;
		pixeldy = y;
		pixeldz = z;
		}

	inline const char *Filename()		{ return(filename); }
	inline const Point *point(int i)	
		{ if((i>=0)&&(i<n)) return(&(points[i])); 
		  else return((Point*)NULL); }

	inline const char  *name(int i)	
		{ if((i>=0)&&(i<n)) return(names[i]);  else return(NULL); }

	inline int  isValid(int i)	
		{ if((i>=0)&&(i<n)) return(valid[i]);  else return(0); }

	void  setValid(int i, int isvalid)	
		{ if((i>=0)&&(i<n)) valid[i] = isvalid; }


	void      clear();
	ItXECode  remove(int i);
	//
	// returns new landmark index, -1 if error
	//
	int       add(double x, double y, double z, const char *nm, int isvalid);

	ItXECode  set(int i, double x, double y, double z, const char *nm, int isvalid);

	//
	// returns landmark index
	//
/*
	int       getLandmarkIndex(double x,double y, double z, 
				double searchrad = 0.0);
	int       getLandmarkIndex(char *name, int case_matters = 0);
*/

	};

#endif
