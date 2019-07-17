//-----------------------------------------------------------------------------------//
//  Author and Copyright Owner of this Code: Mirza Faisal Beg                        //
//  Ph.D. Candidate in Biomedical Engineering, The Johns Hopkins University          //
//  School of Medicine. mfbeg@bme.jhu.edu                                            //
//  Code written towards thesis: Computing Metrics on Medical Images via Metrics     //
//  on Flows of Diffeomorphisms using Image Intensity and Landmark data.             //
//-----------------------------------------------------------------------------------//


#ifndef __MFARRAY3D_H__
#define __MFARRAY3D_H__

////////////////////////FORTRAN STYLE INDEXING CLASS//////////////////////
//
// Purpose: A generic class to create a 3D data array
//          and provide access functions for data
//
// Author: M. Faisal Beg, JHU
//
//	   Stand-along usage of this class (without Analyze image functions) is default
//         To turn on Analyze Load/save functions, use flag -DANALYZE_IMAGE during compilation
///////////////////////////////////////////////////////////////////////////

#include <OS.h>
#include "MFByteSwap.h"

template <class T> class MFArray3D {

public:

	MFArray3D();

	MFArray3D(int Nz,int Ny,int Nx);

	MFArray3D(int Nzl,int Nzh, int Nyl, int Nyh, int Nxl, int Nxh);

	MFArray3D( int Nzmin, int Nzbegin, int Nzend, int Nzmax,
		 int Nymin, int Nybegin, int Nyend, int Nymax, 
		 int Nxmin, int Nxbegin, int Nxend, int Nxmax );


	MFArray3D(const MFArray3D<T> &A);         // Copy constructor.
					//called when Cat A = B;
					//or Cat A(B);

	~MFArray3D();

	 void deleteMemory(); 		//delete mem adn cleanup object



        MFArray3D<T>& operator = (const MFArray3D<T> &A); //Copy Operator
				//if inner-outer indices same, copies them all
				//if only inner same, copies only the inner-data.

        MFArray3D<T>& operator = (T k); //assign the entire array a constant value

        void operator *= (T k); //Multiply all elements by k
        void operator += (T k); //Add all elements by k







	void setDim(int Nz, int Ny, int Nx);

	void setIndices(int Zmin, int Zbegin, int Zend, int Zmax,
			int Ymin, int Ybegin, int Yend, int Ymax, int Xmin, int Xbegin, int Xend, int Xmax);






	inline T* data(int z, int y, int x);

	inline T const * data(int z, int y, int x) const;

        inline T *** address();

        inline T const * const * const * address() const;

	inline T ** operator [] (int i);

	inline T const * const * operator [] (int i) const;
	

	inline int isEmpty() const;	

	inline int isInnerSizeSame(const MFArray3D<T> & A) const; 
	inline int isInnerOuterSizeSame(const MFArray3D<T> & A) const; 

	inline int isInnerIndicesSame(const MFArray3D<T> & A) const; 
	inline int isInnerOuterIndicesSame(const MFArray3D<T> & A) const; 


	inline int getNelm()  const;
	inline int getNelm_Inner()  const;
	inline int getXmin() const;
	inline int getYmin() const;
	inline int getZmin() const;

	inline int getXbegin() const;
	inline int getYbegin() const;
	inline int getZbegin() const;

	inline int getXmax() const;
	inline int getYmax() const;
	inline int getZmax() const;

	inline int getXend() const;
	inline int getYend() const;
	inline int getZend() const;

	inline int getNxTotal() const; //returns #elements between indices min and max in X direction
        inline int getNxInner() const;//returns #elements between indices begin and end in X direction
        inline int getNyTotal() const;
        inline int getNyInner() const;
        inline int getNzTotal() const;
        inline int getNzInner() const;

	inline T getMaxValue() const;
	inline T getMaxAbsValue() const;
	inline T getMaxValue(int &zpos, int &ypos, int &xpos) const; //get max value and its position

        inline T getMinValue() const;
        inline T getMinAbsValue() const;
	inline T getMinValue(int &zpos, int &ypos, int &xpos) const; //get min value and its position

        inline double getMeanValue() const;
        inline double getStdDeviation() const;

	inline T getTrilinear(double z,double y,double x) const;	//return value of array at point(z,y,x)
	inline T getNearestNeighbour(double z,double y,double x) const;	//return value of array at point(z,y,x)

	void getHistogram(int *histogram, int num_bins, int max_value = 255, int min_value = 0) const; 
					  //return a pointer to array of num_bins integers storing the count of pixels
					  //that take on each color value in [0,255]. Assuming pixels in range [0,255].
					  //Assuming histogram array to have num_bins elements.

        int show() const;
        int show(char *ch) const;
        int printXslice(int k) const;   //show the Sagittal slice at x =k
        int printYslice(int k) const;   //show the Coronal Slice at y = k
        int printZslice(int k) const;   //show the Transverse Slice at z = k




	int save(ofstream &fout, streampos byteoffset = 0) const;
	int save(const char * file, streampos byteoffset = 0) const;
	int saveAsText(char * file, streampos byteoffset = 0) const;
        int saveAsVTK(char * file, streampos byteoffset = 0) const;

	int load(ifstream &fin, streampos byteoffset = 0);
	int load(const char *file, streampos byteoffset = 0);


protected:

		T ***_p; // Data Pointer.

		//this is where memory is allocated and pointers set to
		// point to data etc.

	        void setup(int Xmin, int Xbegin, int Xend, int Xmax,
                        int Ymin, int Ybegin, int Yend, int Ymax,
                        int Zmin, int Zbegin, int Zend, int Zmax);

		void  copyInner(const MFArray3D<T> & A);
		void  copyInnerOuter(const MFArray3D<T> & A);



private:

	////////////////////////////////////////////////
	//
	//   -----------Indices Layout--------------------
	//     <-PAD-><--DATA...DATA...DATA--><-PAD->
	//    _xmin..._xbegin............._xend..._xmax
        //            <------_nx------------>
	//    <-------------------_nxTotal----------->
	//   
	//   _xend =  _xbegin  +  _nx  - 1;
        //   _xmax  = _xmin    +  _nxTotal -  1;
	//
	/////////////////////////////////////////////////

 
	int _nxTotal;	// Number of elements in x direction.
	int _nyTotal;	// Number of elements in y direction. 
	int _nzTotal;	// Number of elements in z direction. 
	int _xmin;    // Array starts in x direction at index _xmin
	int _ymin;    // Array starts in y direction at index _ymin
	int _zmin;    // Array starts in z direction at index _zmin 

	int _nx;	// Number of working elements in x direction.
	int _ny;	// Number of working elements in y direction. 
	int _nz;	// Number of working elements in z direction. 
	int _xbegin;    // Work starts in x direction at index xbegin
	int _ybegin;    // Work starts in y direction at index ybegin
	int _zbegin;    // Work starts in z direction at index zbegin


}; 



//_____________________________________________________________________




//Constructors

template <class T> MFArray3D<T>::MFArray3D(){


_nxTotal = 0; _nyTotal = 0; _nzTotal = 0;
_nx = 0; _ny = 0; _nz = 0;
_p = 0;
_xmin = 0; _ymin = 0; _zmin = 0;
_xbegin = 0; _ybegin = 0; _zbegin = 0;


}//constructor




template <class T> MFArray3D<T>::MFArray3D(int Nz,int Ny,int Nx)
{

_nxTotal = 0; _nyTotal = 0; _nzTotal = 0;
_nx = 0; _ny = 0; _nz = 0;
_p = 0;
_xmin = 0; _ymin = 0; _zmin = 0;
_xbegin = 0; _ybegin = 0; _zbegin = 0;

	setup(0, 0, Nx - 1, Nx - 1, 0, 0, Ny - 1, Ny - 1, 0, 0,  Nz - 1, Nz - 1);

}//constructor




template <class T> MFArray3D<T>::MFArray3D(int Nzl, int Nzh, int Nyl, int Nyh, int Nxl, int Nxh)
{

_nxTotal = 0; _nyTotal = 0; _nzTotal = 0;
_nx = 0; _ny = 0; _nz = 0;
_p = 0;
_xmin = 0; _ymin = 0; _zmin = 0;
_xbegin = 0; _ybegin = 0; _zbegin = 0;

        setup(Nxl, Nxl, Nxh, Nxh, Nyl, Nyl, Nyh, Nyh, Nzl, Nzl,  Nzh, Nzh);

}//constructor


template <class T> MFArray3D<T>::MFArray3D( int Nzmin, int Nzbegin, int Nzend, int Nzmax, 
int Nymin, int Nybegin, int Nyend, int Nymax, int Nxmin, int Nxbegin, int Nxend, int Nxmax )
{

_nxTotal = 0; _nyTotal = 0; _nzTotal = 0;
_nx = 0; _ny = 0; _nz = 0;
_p = 0;
_xmin = 0; _ymin = 0; _zmin = 0;
_xbegin = 0; _ybegin = 0; _zbegin = 0;

        setup(Nxmin, Nxbegin, Nxend, Nxmax, 
			Nymin, Nybegin, Nyend, Nymax,
				Nzmin, Nzbegin, Nzend, Nzmax );


}//constructor











template <class T> MFArray3D<T>::MFArray3D(const MFArray3D<T>& A)
{
			//Copy Constructor called when MFArray3D<double> C = A;

_nxTotal = 0; _nyTotal = 0; _nzTotal = 0;
_nx = 0; _ny = 0; _nz = 0;
_p = 0;
_xmin = 0; _ymin = 0; _zmin = 0;
_xbegin = 0; _ybegin = 0; _zbegin = 0;

        if(!A.isEmpty()) 
		{

        // Allocate Memory for object being initialized using A.
        setup(A.getXmin(), A.getXbegin(), A.getXend(), A.getXmax(), 
        A.getYmin(), A.getYbegin(), A.getYend(), A.getYmax(), 
        A.getZmin(), A.getZbegin(), A.getZend(), A.getZmax()); 
        
        // Check memory allocation, exit( prev. return) if memory not allocated
        if( isEmpty() ){cerr<<" Memory could not be allocated in Copy Constructor "<<endl; exit(0); }
			//return;

	copyInnerOuter(A);

        }
	else
		{

			cerr<<" Inside Copy Constructor "<<endl;
			cerr<<" Attempt to initialize an array using another empty array  "<<endl;
			exit(0);

		
		}//if A is empty()

}//copy constructor


template <class T>
void MFArray3D<T>::copyInner(const MFArray3D<T> &A){

	//this is a private method, check inside this class 
	//if InnerIndices same before calling this function.
	
        T const * const * const * Aptr = A.address();
        T *** ThisPtr = address();

        		for(int z = A.getZbegin(); z <= A.getZend(); z++){
                	for(int y = A.getYbegin(); y <= A.getYend(); y++){
                        for(int x = A.getXbegin(); x <= A.getXend(); x++){

                                ThisPtr[z][y][x] = Aptr[z][y][x];

                        }}}


}//copyInner function


template <class T>
void MFArray3D<T>::copyInnerOuter(const MFArray3D<T> &A){

        //this is a private method, check inside this class
        //if InnerIndices same before calling this function.

        T const * const * const * Aptr = A.address();
        T *** ThisPtr = address();

                        for(int z = A.getZmin(); z <= A.getZmax(); z++){
                        for(int y = A.getYmin(); y <= A.getYmax(); y++){
                        for(int x = A.getXmin(); x <= A.getXmax(); x++){

                                ThisPtr[z][y][x] = Aptr[z][y][x];

                        }}}


}//copyInnerOuter function



template <class T> 
void MFArray3D<T>::deleteMemory()
{


	//clear memory from object  and clean it up.

   if(_p){
#ifdef DO_TIMS
{
    int n1l,n1h,n2l,n2h,n3l,n3h;
    n1l=getZmin();n1h=getZmax();
    n2l=getYmin();n2h=getYmax();
    n3l=getXmin();n3h=getXmax();
        _p[n1l][n2l]+=n3l;
        delete [] _p[n1l][n2l];
    for(int k=n1h;k>=n1l;k--){
        _p[k]+=n2l;
        delete [] _p[k];
    }
    _p+=n1l;
    delete[] _p;
}
#else
{
    int n1l,n1h,n2l,n2h,n3l,n3h;
    n1l=getZmin();n1h=getZmax();
    n2l=getYmin();n2h=getYmax();
    n3l=getXmin();n3h=getXmax();
        _p[n1l][n2l]+=n3l;
        delete [] _p[n1l][n2l];
        _p[n1l]+=n2l;
    delete [] _p[n1l];
    _p+=n1l;
    delete[] _p;
 }
#endif
        }


_nxTotal = 0; _nyTotal = 0; _nzTotal = 0;
_nx = 0; _ny = 0; _nz = 0;
_p = 0;
_xmin = 0; _ymin = 0; _zmin = 0;
_xbegin = 0; _ybegin = 0; _zbegin = 0;




}//deleteMemory











// Destructor.
template <class T> MFArray3D<T>::~MFArray3D(){

//cout<<" calling MFArray3D destructor "<<endl;

   if(_p){
#ifdef DO_TIMS
{
    int n1l,n1h,n2l,n2h,n3l,n3h;
    n1l=getZmin();n1h=getZmax();
    n2l=getYmin();n2h=getYmax();
    n3l=getXmin();n3h=getXmax();
	_p[n1l][n2l]+=n3l;
	delete [] _p[n1l][n2l];
    for(int k=n1h;k>=n1l;k--){
    	_p[k]+=n2l;
    	delete [] _p[k];
    }
    _p+=n1l;
    delete[] _p;
}
#else
{
    int n1l,n1h,n2l,n2h,n3l,n3h;
    n1l=getZmin();n1h=getZmax();
    n2l=getYmin();n2h=getYmax();
    n3l=getXmin();n3h=getXmax();
	_p[n1l][n2l]+=n3l;
	delete [] _p[n1l][n2l];
	_p[n1l]+=n2l;
    delete [] _p[n1l];
    _p+=n1l;
    delete[] _p;
 }
#endif
	}
}









//////////////////////////FUNCTION DEFINITIONS//////////////////////


// Return pointer to specified element. Null for empty Volumes.

template <class T>
inline T * MFArray3D<T>::data(int z, int y, int x ){

//no default parameters, as negative array indices possible
// and z = 0, y = 0, x = 0 is not a good default

	
	z = (z >= getZmin() ) ? z : getZmin() ;
	z = (z <= getZmax() ) ? z : getZmax();

	y = (y >= getYmin() ) ? y : getYmin();
	y = (y <= getYmax() ) ? y : getYmax() ;

	x = (x >= getXmin() ) ? x : getXmin();
	x = (x <= getXmax() ) ? x : getXmax() ;
    
	return ( _p ? &_p[z][y][x] : 0);

}

template <class T>
inline T const * MFArray3D<T>::data(int z, int y, int x) const {

	
	z = (z >= getZmin() ) ? z : getZmin() ;
	z = (z <= getZmax() ) ? z : getZmax();

	y = (y >= getYmin() ) ? y : getYmin();
	y = (y <= getYmax() ) ? y : getYmax() ;

	x = (x >= getXmin() ) ? x : getXmin();
	x = (x <= getXmax() ) ? x : getXmax() ;

    return ( _p ? &_p[z][y][x] : 0);

}

// Return 3D pointer. Null for empty Array2Ds. 

template <class T>
inline T       *       *       * MFArray3D<T>::address() 
	{ return ( _p ? _p : 0); }

template <class T>
inline T const * const * const * MFArray3D<T>::address() const 
	{ return ( _p ? _p : 0); }


// Subscripting Operators.

template <class T>
inline T       *       * MFArray3D<T>::operator [] (int i) { return _p[i]; }

template <class T>
inline T const * const * MFArray3D<T>::operator [] (int i) const { return _p[i]; }

// Access Operators.

template <class T>
inline int MFArray3D<T>::isEmpty()  const { return !(_nxTotal * _nyTotal* _nzTotal); }


template <class T>
inline int MFArray3D<T>::isInnerSizeSame(const MFArray3D<T> &A) const{


	if(getNxInner() == A.getNxInner() &&  getNyInner() == A.getNyInner() && getNzInner() ==A.getNzInner() )
	return 1; //inner size is compatible
	else
	return 0;
	
}//isInnerSizeSame

template <class T>
inline int MFArray3D<T>::isInnerOuterSizeSame(const MFArray3D<T> &A) const{


        if(
             getNxInner() == A.getNxInner() &&  getNxTotal() == A.getNxTotal()
                && getNyInner() == A.getNyInner() && getNyTotal() == A.getNyTotal()
                    &&  getNzInner() ==A.getNzInner() &&  getNzTotal() ==A.getNzTotal()
                )
        return 1; //inner outer size is compatible
        else
        return 0;

}//isInnerOuterSizeSame

template <class T>
inline int MFArray3D<T>::isInnerIndicesSame(const MFArray3D<T> &A) const{


     if(getXbegin() == A.getXbegin() && getXend() == A.getXend() 
		&& getYbegin() == A.getYbegin() && getYend() == A.getYend()
			&& getZbegin() == A.getZbegin() && getZend() == A.getZend() )

	return 1;
	
     else

	return 0;


}//isInnerIndicesSame


template <class T>
inline int MFArray3D<T>::isInnerOuterIndicesSame(const MFArray3D<T> &A) const{


     if(
	getXbegin() == A.getXbegin() && getXend() == A.getXend()   
	 && getXmin() == A.getXmin() && getXmax() == A.getXmax()   
                && getYbegin() == A.getYbegin() && getYend() == A.getYend()
                 && getYmin() == A.getYmin() && getYmax() == A.getYmax()
                        && getZbegin() == A.getZbegin() && getZend() == A.getZend() 
                         && getZmin() == A.getZmin() && getZmax() == A.getZmax() 
	)

        return 1; //all indices are compatible
        
     else

        return 0;


}//isInnerOuterIndicesSame







template <class T>
inline int MFArray3D<T>::getNelm()  const { return (_nxTotal*_nyTotal*_nzTotal); }
template <class T>
inline int MFArray3D<T>::getNelm_Inner()  const { return (_nx*_ny*_nz); }

template <class T>
inline int MFArray3D<T>::getNxTotal() const { return _nxTotal; }
template <class T>
inline int MFArray3D<T>::getNxInner() const { return _nx; }


template <class T>
inline int MFArray3D<T>::getNyTotal() const { return _nyTotal; }
template <class T>
inline int MFArray3D<T>::getNyInner() const { return _ny; }

template <class T>
inline int MFArray3D<T>::getNzTotal() const { return _nzTotal; }
template <class T>
inline int MFArray3D<T>::getNzInner() const { return _nz; }




template <class T>
inline int MFArray3D<T>::getXmin() const { return _xmin; }
template <class T>
inline int MFArray3D<T>::getYmin() const { return _ymin; }
template <class T>
inline int MFArray3D<T>::getZmin() const { return _zmin; }

template <class T>
inline int MFArray3D<T>::getXbegin() const { return _xbegin; }
template <class T>
inline int MFArray3D<T>::getYbegin() const { return _ybegin; }
template <class T>
inline int MFArray3D<T>::getZbegin() const { return _zbegin; }

template <class T>
inline int MFArray3D<T>::getXend() const { return (_xbegin + _nx - 1); }
template <class T>
inline int MFArray3D<T>::getYend() const { return (_ybegin + _ny - 1); }
template <class T>
inline int MFArray3D<T>::getZend() const { return (_zbegin + _nz - 1); }

template <class T>
inline int MFArray3D<T>::getXmax() const { return (_xmin + _nxTotal - 1); }
template <class T>
inline int MFArray3D<T>::getYmax() const { return (_ymin + _nyTotal - 1); }
template <class T>
inline int MFArray3D<T>::getZmax() const { return (_zmin + _nzTotal - 1); }


////
template <class T>
inline T MFArray3D<T>::getTrilinear(double z, double y, double x) const{
//return value of array at point(z,y,x)

if(z< getZmax()+1 && z>=getZmin() && y<getYmax()+1 && y>=getYmin() && x <getXmax()+1 && x>=getXmin())
{

	T const * const * const* InputPPPtr = address();

	int z0, y0, x0, z1, y1, x1;
	double fz, fy, fx; 
	T d000, d001, d010, d011, d100, d101, d110, d111;
	
	z0 = (int) floor(z); fz = z - z0; z1= (z0 < getZmax() )? (z0 + 1):getZbegin();
	y0 = (int) floor(y); fy = y - y0; y1 =(y0 < getYmax() )? (y0 + 1):getYbegin();
	x0 = (int) floor(x); fx = x - x0; x1 =(x0 < getXmax() )? (x0 + 1):getXbegin();
	//z0 = (int) floor(z); fz = z - z0; z1= (z0 < getZmax() )? (z0 + 1):z0;
	//y0 = (int) floor(y); fy = y - y0; y1 =(y0 < getYmax() )? (y0 + 1):y0;
	//x0 = (int) floor(x); fx = x - x0; x1 =(x0 < getXmax() )? (x0 + 1):x0;

     d000 = (double) InputPPPtr[z0][y0][x0];
     d001 = (double) InputPPPtr[z0][y0][x1];
     d010 = (double) InputPPPtr[z0][y1][x0];
     d011 = (double) InputPPPtr[z0][y1][x1];

     d100 = (double) InputPPPtr[z1][y0][x0];
     d101 = (double) InputPPPtr[z1][y0][x1];
     d110 = (double) InputPPPtr[z1][y1][x0];
     d111 = (double) InputPPPtr[z1][y1][x1];

	return (T) ( (1-fz)*( (1-fy) * ( (1-fx)*d000 + fx*d001 ) + fy * ( (1-fx)*d010 + fx*d011 )  )
                                          + fz * ( (1-fy) * ( (1-fx)*d100 + fx*d101 ) + fy * ( (1-fx)*d110 + fx*d111 )  )  );

}else{

	cout<<" Zmax = "<<getZmax()<<", Zmin = "<<getZmin()<<endl;
	cout<<" Ymax = "<<getYmax()<<", Ymin = "<<getYmin()<<endl;
	cout<<" Xmax = "<<getXmax()<<", Xmin = "<<getXmin()<<endl;

     cout<<" Point z="<<z<<", y="<<y<<", x="<<x<<" is out of array bounds "<<endl;
     cout<<show()<<endl;
     cout<<" Breaking in getTrilinear(int z, int y, int x) in MFArray3D.h"<<endl;
     exit(0);

	return 0;

}//if check inside


}//getTrilinear
/////



////
template <class T>
inline T MFArray3D<T>::getNearestNeighbour(double z, double y, double x) const{
//return value of array at point(z,y,x)

if(z< getZmax()+1 && z>=getZmin() && y<getYmax()+1 && y>=getYmin() && x <getXmax()+1 && x>=getXmin())
{

	T const * const * const* InputPPPtr = address();

	int z0, y0, x0, z1, y1, x1, zn, yn, xn;
	double fz, fy, fx; 
	T d000, d001, d010, d011, d100, d101, d110, d111;
	
	z0 = (int) floor(z); fz = z - z0; z1= (z0 < getZmax() )? (z0 + 1):getZbegin();
	y0 = (int) floor(y); fy = y - y0; y1 =(y0 < getYmax() )? (y0 + 1):getYbegin();
	x0 = (int) floor(x); fx = x - x0; x1 =(x0 < getXmax() )? (x0 + 1):getXbegin();

	if(fz < 0.5) zn = z0;
	else zn = z1;
	if(fy < 0.5) yn = y0;
	else yn = y1;
	if(fx < 0.5) xn = x0;
	else xn = x1;

	return (T) InputPPPtr[zn][yn][xn];

}else{

	cout<<" Zmax = "<<getZmax()<<", Zmin = "<<getZmin()<<endl;
	cout<<" Ymax = "<<getYmax()<<", Ymin = "<<getYmin()<<endl;
	cout<<" Xmax = "<<getXmax()<<", Xmin = "<<getXmin()<<endl;

     cout<<" Point z="<<z<<", y="<<y<<", x="<<x<<" is out of array bounds "<<endl;
     cout<<show()<<endl;
     cout<<" Breaking in getTrilinear(int z, int y, int x) in MFArray3D.h"<<endl;
     exit(0);

	return 0;

}//if check inside


}//getNearestNeighbour
/////





template <class T>
int MFArray3D<T>::show() const {

cout<<"Indices[X: "<<getXmin()<<","<<getXbegin()<<","<<getXend()<<","<<getXmax();
cout<<", Y: "<<getYmin()<<","<<getYbegin()<<","<<getYend()<<","<<getYmax();
cout<<", Z: "<<getZmin()<<","<<getZbegin()<<","<<getZend()<<","<<getZmax()<<" ]";
cout<<endl;

return 1;

}

template <class T>
int MFArray3D<T>::show(char *ch) const {

cout<<ch<<".Indices[X: "<<getXmin()<<","<<getXbegin()<<","<<getXend()<<","<<getXmax();
cout<<", Y: "<<getYmin()<<","<<getYbegin()<<","<<getYend()<<","<<getYmax();
cout<<", Z: "<<getZmin()<<","<<getZbegin()<<","<<getZend()<<","<<getZmax()<<" ]";
cout<<endl;

return 1;

}



template <class T>
int MFArray3D<T>::printXslice(int k) const {
        //SAGGITTAL
       if(k < _xmin || k > _nxTotal + _xmin -1 ) {
             cerr<<" X slice requested out of bounds :-("<<endl;
             cerr<<"             Valid range ["<<_xmin<<" , "<<_xmin+_nxTotal -1<<"]"<<endl; return 0;}

        cout<<"x="<<k<<endl;
        for ( int y = _ymin + _nyTotal - 1; y >=_ymin; y--){
        cout<<"y="<<y<<": ";
        for (int z = _zmin + _nzTotal-1; z>= _zmin; z--){
        cout<<(float) _p[z][y][k]<< " ";
        }//z
        cout<<endl;
        }//y
        cout<<endl;
return 1;
}

template <class T>
int MFArray3D<T>::printYslice(int k) const {
        //CORONAL
       if(k < _ymin || k > _nyTotal + _ymin -1 ) {
             cerr<<" Y slice requested out of bounds :-("<<endl;
             cerr<<"             Valid range ["<<_ymin<<" , "<<_ymin+_nyTotal-1<<"]"<<endl; return 0;}

        cout<<"y="<<k<<endl;
        for (int z = _zmin + _nzTotal -1; z>= _zmin; z--){
        cout<<"z="<<z<<": ";
        for ( int x = _xmin; x <=_xmin +_nxTotal -1; x++){
        cout<<(float) _p[z][k][x]<< " ";
        }//x
        cout<<endl;
        }//z
        cout<<endl;
return 1;
}


template <class T>
int MFArray3D<T>::printZslice(int k) const {
        //TRANSVERSE
	//cout.precision(2);
       if(k < _zmin || k > _nzTotal + _zmin -1 ) {
             cerr<<" Z slice requested out of bounds :-("<<endl;
             cerr<<"             Valid range ["<<_zmin<<" , "<<_zmin+_nzTotal-1<<"]"<<endl; return 0;}

        cout<<"z="<<k<<endl;
        for ( int y = getYmax(); y >=getYmin(); y--){
        cout<<"y="<<y<<": ";
        for (int x = getXmin(); x<= getXmax(); x++){
        cout<<(float) _p[k][y][x]<< " ";
        }//x
        cout<<endl;
        }//y
        cout<<endl;
return 1;
}









//-------------------------------------------------------------
//		Miscellaneous Functions
//-------------------------------------------------------------

template <class T>
void  MFArray3D<T>::getHistogram(int *histogram, int num_bins, int max_value, int min_value) const{

	cout<<" Here 1"<<endl;
	for(int i = 0; i < num_bins; i++)histogram[i] = 0; //clear the passed histogram

	cout<<" Here 2"<<endl;
	if(num_bins <1 || num_bins >256)
		{cout<<" Number of bins < 1 or > 256 "<<endl; exit(0);}

	T data_max_value = getMaxValue();
	T data_min_value = getMinValue();

	if(data_max_value > max_value || data_max_value < min_value){cout<<" Max Value in data is > 255 or < 0 "<<endl; exit(0);}
	if(data_min_value > max_value || data_min_value < min_value){cout<<" Min Value in data is > 255 or < 0 "<<endl; exit(0);}

	double bin_width = (double)(max_value - min_value + 1)/(double)num_bins;

	cout<<" bin_width = "<<bin_width<<endl;

	int binindex = 0;


        T const * const * const* dataptr = address();

        for(int z = getZbegin(); z <=getZend(); z++){
        for(int y = getYbegin(); y <= getYend(); y++){
        for(int x = getXbegin(); x <=getXend(); x++){


	binindex = (int) ( dataptr[z][y][x]/ bin_width);

	histogram[binindex]++;


        }}}


	



}//getHistogram(int num_bins)





template <class T>
T MFArray3D<T>::getMaxValue() const {

        T const * const * const* dataptr = address();
        T max = dataptr[getZbegin()][getYbegin()][getXbegin()];

	for(int z = getZbegin(); z <=getZend(); z++){
	for(int y = getYbegin(); y <= getYend(); y++){
	for(int x = getXbegin(); x <=getXend(); x++){
		if(dataptr[z][y][x] > max)max = dataptr[z][y][x];
	}}}

	return max;
}



template <class T>
T MFArray3D<T>::getMaxAbsValue() const {

        T const * const * const* dataptr = address();
        T max = dataptr[getZbegin()][getYbegin()][getXbegin()];
		//get its absolute value
	max = (max > 0)? max: -1*max;

	T absval = 0;

	
        for(int z = getZbegin(); z <=getZend(); z++){
        for(int y = getYbegin(); y <= getYend(); y++){
        for(int x = getXbegin(); x <=getXend(); x++){
	
	  absval = dataptr[z][y][x];
		//get its absolute value
	  absval = (absval > 0)? absval: -1*absval;

                if(absval > max)max = absval;
        }}}

        return max;
}






template <class T>
T MFArray3D<T>::getMaxValue(int &zpos, int &ypos, int &xpos) const {

        T const * const * const* dataptr = address();
        T max = dataptr[getZbegin()][getYbegin()][getXbegin()];

        for(int z = getZbegin(); z <=getZend(); z++){
        for(int y = getYbegin(); y <= getYend(); y++){
        for(int x = getXbegin(); x <=getXend(); x++){
                if(dataptr[z][y][x] > max){
			max = dataptr[z][y][x];
		zpos = z; ypos = y; xpos = x; }

        }}}

        return max;
}



template <class T>
T MFArray3D<T>::getMinValue() const {

 T const * const * const* dataptr = address();
        T min= dataptr[getZbegin()][getYbegin()][getXbegin()];

        for(int z = getZbegin(); z <=getZend(); z++){
        for(int y = getYbegin(); y <= getYend(); y++){
        for(int x = getXbegin(); x <=getXend(); x++){
                if(dataptr[z][y][x] < min){
                        min = dataptr[z][y][x]; }

        }}}

        return min;

}


template <class T>
T MFArray3D<T>::getMinAbsValue() const {

 T const * const * const* dataptr = address();
 T min= dataptr[getZbegin()][getYbegin()][getXbegin()]; 
 
	//get its absolute value
	min = (min > 0)? min: -1*min;

 T absval = 0;

        for(int z = getZbegin(); z <=getZend(); z++){
        for(int y = getYbegin(); y <= getYend(); y++){
        for(int x = getXbegin(); x <=getXend(); x++){

		absval = dataptr[z][y][x] ;
		//get its absolute value
		absval = (absval > 0) ? absval: -1*absval;

                if(absval  < min) min = absval; 

        }}}

        return min;

}























template <class T>
T MFArray3D<T>::getMinValue(int &zpos, int &ypos, int &xpos) const {

        T const * const * const* dataptr = address();
        T min= dataptr[getZbegin()][getYbegin()][getXbegin()];

        for(int z = getZbegin(); z <=getZend(); z++){
        for(int y = getYbegin(); y <= getYend(); y++){
        for(int x = getXbegin(); x <=getXend(); x++){
                if(dataptr[z][y][x] < min){
                        min = dataptr[z][y][x];
                zpos = z; ypos = y; xpos = x; }

        }}}

        return min;
}


template <class T>
double MFArray3D<T>::getMeanValue() const {

	T const* const * const * dataptr = address();
    	double mean = 0;

	if(!isEmpty() )
	{
	for(int z = getZbegin(); z <=getZend(); z++){
        for(int y = getYbegin(); y <= getYend(); y++){
        for(int x = getXbegin(); x <=getXend(); x++){

			mean += dataptr[z][y][x];
		}}}

		mean = mean/ _nx * _ny * _nz ;
	}

	return mean;

}


template <class T>
double MFArray3D<T>::getStdDeviation() const {

	T const* const * const * dataptr = address();
    	double mean = getMeanValue();
	double stddev =  0;
	double temp = 0;

	if(!isEmpty() )
	{
 	for(int z = getZbegin(); z <=getZend(); z++){
        for(int y = getYbegin(); y <= getYend(); y++){
        for(int x = getXbegin(); x <=getXend(); x++){

			temp = dataptr[z][y][x] - mean;
			stddev += temp * temp ; 
		}}}

		stddev = stddev/ ( _nx * _ny * _nz - 1); 
        //this is the estimate of Std. Dev. 
	//with unknown mean. If mean known, divide by N
	}

	return stddev;

}



//////////////////OPERATOR OVERLOADING///////////////////////////////////////////////////



template <class T>
void MFArray3D<T>::operator *= (T k){

        T ***ptr = address();

        if(!ptr){cerr<<" Null array being initialized in assignment"<<endl; exit(1);}

        for(int z = getZmin(); z <= getZmax(); z++){
        for(int y = getYmin(); y <= getYmax(); y++){
        for(int x = getXmin(); x <= getXmax(); x++){

        ptr[z][y][x] *= k;

        }}}


}//operator *=

template <class T>
void MFArray3D<T>::operator += (T k){

        T ***ptr = address();

        if(!ptr){cerr<<" Null array being initialized in assignment"<<endl; exit(1);}
        for(int z = getZmin(); z <= getZmax(); z++){
        for(int y = getYmin(); y <= getYmax(); y++){
        for(int x = getXmin(); x <= getXmax(); x++){

        ptr[z][y][x] += k;

        }}}


}//operator +=



template <class T>
MFArray3D<T>& MFArray3D<T>::operator = (T k){

	T ***ptr = address();

	if(!ptr){cerr<<" Null array being initialized in assignment"<<endl; exit(1);}
	for(int z = getZmin(); z <= getZmax(); z++){
	for(int y = getYmin(); y <= getYmax(); y++){
	for(int x = getXmin(); x <= getXmax(); x++){

	ptr[z][y][x] = k;
	
	}}}

	return *this;

}//operator = k


template <class T>
MFArray3D<T>& MFArray3D<T>::operator = (const MFArray3D<T> &A){

	if(!A.isEmpty()){


		if(  isEmpty()  ) 
			{

				setup(  A.getXmin(), A.getXbegin(), A.getXend(), A.getXmax(),
			        A.getYmin(), A.getYbegin(), A.getYend(), A.getYmax(),
        			A.getZmin(), A.getZbegin(), A.getZend(), A.getZmax());

			if( isEmpty() ) {cerr<< " memory Allocation failed "<<endl; exit(1);}
			copyInnerOuter(  A   );

			} 	
		else if( isInnerOuterIndicesSame( A ) )
			{

				cout<<" Inside copyInnerOuter( A ); called from Copy Operator "<<endl;

				copyInnerOuter( A );

			}
		else if( isInnerIndicesSame( A ) )
                        {

                                copyInner( A );
                        }
		else
			{
			//if dimensions are not the same and ptr of calling object is not NULL either
			// then the array exists and does not match in dimension.
			cerr<<" Inside Copy operator = "<<endl;
			cerr<<"MFArray3D<T>::= MFArray3D's have different dimensions."<<endl;
			exit(1);
			}

	}//if(!A.isEmpty())
	else
		{
			cerr<<" Inside Copy operator = "<<endl;
			cerr<<" Array to be copied from is Empty, copy not defined....bailing out "<<endl;
			exit(0);
		} //if(!A.isEmpty())



	return *this;

}//operator = Array A






template <class T>
void MFArray3D<T>::setDim(int Nz, int Ny, int Nx) {

 setup(0, 0, Nx - 1, Nx - 1, 0, 0, Ny - 1, Ny - 1, 0, 0,  Nz - 1, Nz - 1);

}


template <class T>
void MFArray3D<T>::setIndices(int Zmin, int Zbegin, int Zend, int Zmax,
                        int Ymin, int Ybegin, int Yend, int Ymax,int Xmin, int Xbegin, int Xend, int Xmax)
{


	setup(Xmin , Xbegin, Xend, Xmax,  Ymin , Ybegin, Yend, Ymax, Zmin, Zbegin, Zend, Zmax   );

}









//////////////////Allocate memory///////////////////////////////

template <class T>
void MFArray3D<T>::setup(int Xmin, int Xbegin, int Xend, int Xmax,
                        int Ymin, int Ybegin, int Yend, int Ymax, 
                        int Zmin, int Zbegin, int Zend, int Zmax){

//cout<<" Allocating memory for MFArray3D object "<<endl;

    // Delete Old Memory First.
    // ------------------------
    if ( _p ) {
#ifdef DO_TIMS
{
    int n1l,n1h,n2l,n2h,n3l,n3h;
    n1l=getZmin();n1h=getZmax();
    n2l=getYmin();n2h=getYmax();
    n3l=getXmin();n3h=getXmax();
	_p[n1l][n2l]+=n3l;
	delete [] _p[n1l][n2l];
    for(int k=n1h;k>=n1l;k--){
    	_p[k]+=n2l;
    	delete [] _p[k];
    }
    _p+=n1l;
    delete[] _p;
}
#else
    {
    int n1l,n1h,n2l,n2h,n3l,n3h;
    n1l=getZmin();n1h=getZmax();
    n2l=getYmin();n2h=getYmax();
    n3l=getXmin();n3h=getXmax();
	_p[n1l][n2l]+=n3l;
	delete [] _p[n1l][n2l];
	_p[n1l]+=n2l;
    delete [] _p[n1l];
    _p+=n1l;
    delete[] _p;
    }
#endif
        _p = 0;
        _nxTotal = 0;
        _nyTotal = 0;
        _nzTotal = 0;
        _xmin=0;
        _ymin=0;
        _zmin=0;
        _nx = 0;
        _ny = 0;
        _nz = 0;
        _xbegin=0;
        _ybegin=0;
        _zbegin=0;

    }
	
	if(Xbegin < Xmin || Xend < Xbegin || Xmax < Xend){cerr << "Dimensions not good in X "<<endl; exit(1);}
	if(Ybegin < Ymin || Yend < Ybegin || Ymax < Yend){cerr << "Dimensions not good in Y "<<endl; exit(1);}
	if(Zbegin < Zmin || Zend < Zbegin || Zmax < Zend){cerr << "Dimensions not good in Z "<<endl; exit(1);}


		T **ptr1;
		T  *ptr2;
		

		_xmin = Xmin; 
		_ymin = Ymin; 
		_zmin = Zmin;

		_nxTotal = Xmax - Xmin + 1; 
		_nyTotal = Ymax - Ymin + 1; 
		_nzTotal = Zmax - Zmin + 1;

		_xbegin = Xbegin; 
		_ybegin = Ybegin; 
		_zbegin = Zbegin;

		_nx = Xend - Xbegin + 1; 
		_ny = Yend - Ybegin + 1; 
		_nz = Zend - Zbegin + 1;
		//allocate memory here
#ifdef DO_TIMS
{
        int i,k;
        int n1l,n1h,n2l,n2h,n3l,n3h;
        n1l=Zmin;n1h=Zmax;
        n2l=Ymin;n2h=Ymax;
        n3l=Xmin;n3h=Xmax;
		_p=new T** [(n1h-n1l+1)];
        if (!_p){
                cout << "allocation failure at 0 in cube3()\n";
                exit(1);
        }
        _p-=n1l;
        for(k=n1l;k<=n1h;k++){
                  _p[k]= new T* [n2h-n2l+1];
                
                if (!_p[k]){
                        cout << "allocation failure at 1 in cube3()\n";
                        exit(1);
                }
                _p[k] -= n2l;
                        for(i=n2l;i<=n2h;i++) {
                                 if(i == n2l && k==n1l){ 
                                          _p[k][i]= new T [               (n1h-n1l+1)*
                                                                          (n2h-n2l+1)*
                                                                          (n3h-n3l+1)];
                                        if (!_p[k][i]){
                                                cout << "allocation failure at 2 in cube3()\n";
                                                exit(1);
                                        }               
                                        _p[k][i] -= n3l;
                                }
                                else {
                                        if(i == n2l) 
                                                _p[k][i]=_p[k-1][i]+(n2h-n2l+1)*(n3h-n3l+1);
                                        else 
                                        _p[k][i]=_p[k][i-1]+(n3h-n3l+1);
                                }
                } 
        }
}
#else
{
		_p   = new T** [_nzTotal];
		ptr1 = new T*  [_nzTotal*_nyTotal];
		ptr2 = new T   [_nzTotal*_nyTotal*_nxTotal];

		//Exit if memory cannot be allocated.
		// ----------------------------------------------------------
		if( !_p || !ptr1 || !ptr2 ) {
			if( _p   ) delete [] _p;
			if( ptr1 ) delete [] ptr1;
			if( ptr2 ) delete [] ptr2;
		cerr<<" Failed to allocate Memory "<<endl; exit(1);
		}



		ptr1[0] = ptr2 - _xmin; //filing the sagittal view plane of pointers to data
		for(int i=1; i< _nzTotal * _nyTotal; i++)
			ptr1[i] = ptr1[i-1] + _nxTotal;


		_p[0] = ptr1 - _ymin;
		for(int z=1; z< _nzTotal; z++)
			_p[z] = _p[z-1] + _nyTotal;

		_p  = _p - _zmin;
					//_p[_zmin][_ymin][_xmin] translates to data[0][0][0]
}
#endif

}//setup








//////////////////////////////INPUT/OUTPUT FUNCTIONS//////////////////////////////////





// Save MFArray3D<T> to file.
template <class T>
int MFArray3D<T>::save(ofstream &fout, streampos byteoffset) const{

// Check file.
        if(!fout)
        {cerr<<" File cannot be used to save, Null Stream"<<endl; return 1;}

//      Seek position from the current put pointer to #bytes in byteoffset.
        fout.seekp(byteoffset,ios::cur);
        if(!fout)
        {cerr<<"Error in seeking Position"<<endl; return 1;}

//      Write out data size in bytes first

        int size = sizeof(T);
        fout.write( (char *) &size,   sizeof(size));

//      Write Array Indices Min, Begin, End, Max first.

        int _zend, _zmax, _yend, _ymax, _xend, _xmax;
        _zend = getZend(); _zmax = getZmax();
        _yend = getYend(); _ymax = getYmax();
        _xend = getXend(); _xmax = getXmax();

        fout.write( (char *) &_zmin,   sizeof(_zmin));
        fout.write( (char *) &_zbegin, sizeof(_zbegin));
        fout.write( (char *) &_zend,   sizeof(_zend));
        fout.write( (char *) &_zmax,   sizeof(_zmax));

        fout.write( (char *) &_ymin,   sizeof(_ymin));
        fout.write( (char *) &_ybegin, sizeof(_ybegin));
        fout.write( (char *) &_yend,   sizeof(_yend));
        fout.write( (char *) &_ymax,   sizeof(_ymax));

        fout.write( (char *) &_xmin,   sizeof(_xmin));
        fout.write( (char *) &_xbegin, sizeof(_xbegin));
        fout.write( (char *) &_xend,   sizeof(_xend));
        fout.write( (char *) &_xmax,   sizeof(_xmax));

        if(!fout)
        {cerr<<" Error writing array indices to file"<<endl; return 1;}

//      Write entire contents of the Array into file

       if( _nxTotal * _nyTotal * _nzTotal)
       {
       fout.write( (char *) data(getZmin(), getYmin(), getXmin()), _nxTotal * _nyTotal *_nzTotal * sizeof(T) );
                if(!fout)
                {cerr<<" Error writing Data to file"<<endl; return 1; }
       }//if


//since called with a ostream, let the calling function do the closing
        //fout.close();

        return 1;


}//end of save 





//load MFArray3D from file
template <class T>
int MFArray3D<T>::load(ifstream &fin, streampos byteoffset){


//check stream vital signs

        if(!fin)
        {cerr<<" Invalid or Null output stream passed to load function, Can't load from this stream "<<endl;
	return 0;}

//Seek position from the current "get pointer" to #bytes in byteoffset.

        fin.seekg(byteoffset,ios::cur);
        if(!fin){cerr<<"Error in seeking Position to specified offset"<<byteoffset<<endl; return 1;}

//read size of data stored in the stream

        int size;
        fin.read( (char *) &size,    sizeof(size));

#if 1
		int ByteSwapNeeded = 0;
		if(size > 64 || size < 0){ ByteSwapNeeded = 1; classByteSwap::BYTE_SWAP_INT(&size);};
		if(ByteSwapNeeded)cout<<" Size after ByteSwapping is "<<size<<endl;
#endif

        if(size!=sizeof(T)){cerr<< "******Datatype size in stream and in Array to be loaded don't match****** "<<endl;
                           cerr<<"Size in bytes of data in stream="<<size<<", in Array="<<sizeof(T)<<endl;
                           return 0;}

//read Dimensions

	int _zmn, _zbg, _ze, _zmx, _ymn, _ybg, _ye, _ymx, _xmn, _xbg, _xe, _xmx; 

        fin.read( (char *) &_zmn,   sizeof(_zmn));
        fin.read( (char *) &_zbg, sizeof(_zbg));
        fin.read( (char *) &_ze,   sizeof(_ze));
        fin.read( (char *) &_zmx,   sizeof(_zmx));

        fin.read( (char *) &_ymn,   sizeof(_ymn));
        fin.read( (char *) &_ybg, sizeof(_ybg));
        fin.read( (char *) &_ye,   sizeof(_ye));
        fin.read( (char *) &_ymx,   sizeof(_ymx));

        fin.read( (char *) &_xmn,   sizeof(_xmn));
        fin.read( (char *) &_xbg, sizeof(_xbg));
        fin.read( (char *) &_xe,   sizeof(_xe));
        fin.read( (char *) &_xmx,   sizeof(_xmx));

#if 1
	if(ByteSwapNeeded){


		classByteSwap::BYTE_SWAP_INT(&_zmn);
		classByteSwap::BYTE_SWAP_INT(&_zbg);
		classByteSwap::BYTE_SWAP_INT(&_ze);
		classByteSwap::BYTE_SWAP_INT(&_zmx);

		classByteSwap::BYTE_SWAP_INT(&_ymn);
		classByteSwap::BYTE_SWAP_INT(&_ybg);
		classByteSwap::BYTE_SWAP_INT(&_ye);
		classByteSwap::BYTE_SWAP_INT(&_ymx);
		

		classByteSwap::BYTE_SWAP_INT(&_xmn);
		classByteSwap::BYTE_SWAP_INT(&_xbg);
		classByteSwap::BYTE_SWAP_INT(&_xe);
		classByteSwap::BYTE_SWAP_INT(&_xmx);

	}//if(ByteSwapNeeded)
#endif

//check dimensions

        int _nzTl = _zmx - _zmn + 1;
        int _nyTl = _ymx - _ymn + 1;
        int _nxTl = _xmx - _xmn + 1;

        if( (_nxTl==0) || (_nyTl==0) | (_nzTl==0) )
        {cerr<<" One of the dimensions is Zero, File Corrupted "<<endl; return 0;}

//allocate memory for the data with the indices preserved.....

     setup(_xmn , _xbg, _xe, _xmx,  _ymn , _ybg, _ye, _ymx, _zmn, _zbg, _ze, _zmx   );



//check memory allocation
        if(!_p){cerr<<" memory could not be allocated "<<endl; return 0;}


//read in data from stream and store at the given pointer

     fin.read( (char *) data(getZmin(), getYmin(), getXmin()), _nxTotal * _nyTotal * _nzTotal * sizeof(T));

    if(!fin)
        {cerr<<"Error reading data from stream "<<endl; return 0;}


#if 1
        if(ByteSwapNeeded){

		

		for(int k = getZmin(); k <= getZmax(); k++)
			for(int j = getYmin(); j <= getYmax(); j++)
				for(int i = getXmin(); i <= getXmax(); i++){

				 if(size==8)classByteSwap::BYTE_SWAP_DOUBLE( data(k, j, i) );
				 else{cerr<<"Inside MFArray3D.h:  ByteSwap not implemented to size other than 8 "<<endl; exit(0);}
				 
				}//i, j, k


        }//if(ByteSwapNeeded)

	
#endif


        //calling function closes stream.....
        //fin.close();

        return 1;


}//load(ifstream fin....)




//load MFArray3D from file
template <class T>
int MFArray3D<T>::save(const char * file, streampos byteoffset) const{

//create an output file-stream
	//ofstream fout(file, ios::out|ios::binary|ios::noreplace);
	ofstream fout(file, ios::out|ios::binary); //open a binary file by default
	if(!fout)
	{cerr<<" File could not be opened or already exists"<<endl; return 0;}

//call save to save data
    save(fout, byteoffset);

//close stream
	fout.close();

	return 1;
}//save(char * file, streampos byteoffset) const



//load MFArray3D from file
template <class T>
int MFArray3D<T>::load(const char * file, streampos byteoffset){

//create an input file-stream
	ifstream fin(file, ios::in|ios::binary);
	if(!fin)
	{cerr<<" File could not be opened "<<endl; return 0;}

//call load to load data
    load(fin, byteoffset);

//close stream
	fin.close();

	return 1;
}//load(char * file, streampos byteoffset)







//load MFArray3D from file
template <class T>
int MFArray3D<T>::saveAsText(char * file, streampos byteoffset) const{

cout<<" Inside saveAsText() "<<endl;

//create an output file-stream
        //ofstream fout(file, ios::out|ios::binary|ios::noreplace);
        ofstream fout(file, ios::out); //open a binary file by default
        if(!fout)
        {cerr<<" File could not be opened or already exists"<<endl; return 0;}

//call save to save data
	
	////////////////////////

	        for(int z = getZbegin();  z<= getZend(); z++){
                	for(int y = getYbegin(); y<= getYend(); y++){
                        	for(int x = getXbegin(); x <= getXend(); x++){

                                	//fout<<(FLT)_p[z][y][x]; 
                                	fout<<(double) _p[z][y][x]; 
					if( y!=getYend() && x!= getXend() ) fout<<",";

                                                   }//x
		
                                                }//y
				fout<<endl;

                                           }//z

	////////////////////////

//close stream
        fout.close();

        return 1;
}//save(char * file, streampos byteoffset) const


//-------------------------------------------------------------------
template <class T>
int MFArray3D<T>::saveAsVTK(char * file, streampos byteoffset) const{

        ofstream fout(file, ios::out); //open a binary file by default
        if(!fout)
        {cerr<<" File could not be opened or already exists"<<endl; return 0;}

 fout<<"# vtk DataFile Version 3.0"<<endl;
 fout<<"vtk output"<<endl;
 fout<<"ASCII"<<endl;
 fout<<"DATASET STRUCTURED_POINTS"<<endl;
 fout<<"DIMENSIONS "<<getXend()+1<<" "<<getYend()+1<<" "<<getZend()+1<<endl;
 fout<<"SPACING 1 1 1"<<endl;
 fout<<"ORIGIN 0 0 0"<<endl;
 fout<<"POINT_DATA "<<(getXend()+1)*(getYend()+1)*(getZend()+1)<<endl;
 fout<<"SCALARS Scalars float"<<endl;
 fout<<"LOOKUP_TABLE default"<<endl;

 for(int z = getZbegin();  z<= getZend(); z++){
                for(int y = getYbegin(); y<= getYend(); y++){
                 for(int x = getXbegin(); x <= getXend(); x++){
    fout<<(double) _p[z][y][x]<<endl;
                        }//x
  }//y
 }//z

        fout.close();
        return 1;
}
//-------------------------------------------------------------------


#endif
//#define __ARRAY3D_H__
