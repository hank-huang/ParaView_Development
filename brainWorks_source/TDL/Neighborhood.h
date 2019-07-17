#ifndef __NEIGHBORHOOD_H__
#define __NEIGHBORHOOD_H__
///////////////////////////////////////////////////////////////////////////
//
//  This class defines a neighborhood structure for a point in a surface
// 
///////////////////////////////////////////////////////////////////////////
#include <OS.h>
#include <ADL/Array1D.h>

class Neighborhood
{

public:
  // Default constructor
  Neighborhood() { clear(); }

  // copy constructor
  Neighborhood(Neighborhood &N) { clear(); *this = N; }

  // destructor
  virtual ~Neighborhood() {}

  // assignment operator
  const Neighborhood & operator= (const Neighborhood &N);

  int numNeighbors() const { return mNum; }

  inline int getNeighbor(int i)      const { return mNeighbor[i]; }
  inline int getNeighborDepth(int i) const { return(mNeighborDepth[i]); }

  // this should be called BEFORE any addNeighbor
  void setNeedDepth(bool tf)  { mNeedDepth = tf; }

  // only call these is previously used setNeedDepth(true)
  void setDepth(int i, int d)   { if(mNeedDepth) mNeighborDepth[i] = d;     }
  u_char  getDepth(int i) const { return((mNeedDepth)?mNeighborDepth[i]:1U);}

  // check if the vertex with index i is a neighbor of this vertex
  bool isNeighbor(int i) const { return(neighborIndex(i)>=0); }
 
  // get index of where 'i' falls in this neighbor list (<0 = not in list)
  int  neighborIndex(int i) const;

  // add a vertex to the neighborhood
  int  addNeighbor(int nbr, u_char depth=1, bool checkdup=true);
  void cleanNeighborhood(int numVerts);

  // remove a vertex from the neighborhood (true if found & removed)
  bool removeNeighbor(int);

  // clear a neighborhood
  void clear() { 
	mNum = 0; 
        mNeedDepth = false;
	mNeighbor.setDim(0); 
	mNeighborDepth.setDim(0); 
  }

  //
  // print function
  //
  void print (ostream &os) const;
private:

  int mNum;			  // number of vertices in the neighborhood
  Array1D<int>    mNeighbor;	  // list of the indices of the vertices in
				  // the neighborhood
  bool mNeedDepth;
  Array1D<u_char> mNeighborDepth; // depth of each neighbor (if reqd)
};

//
// non-member print function
//
ostream & operator<< (ostream &os, Neighborhood const &N);

#endif // __NEIGHBORHOOD_H__
