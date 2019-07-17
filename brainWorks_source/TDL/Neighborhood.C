///////////////////////////////////////////////////////////////////////////
//
// File: Neighborhood.C
//
// Author: Rob Teichman
//
// Purpose: Neighborhood class body
// 
///////////////////////////////////////////////////////////////////////////
#include <OS.h>
#include <TDL/Neighborhood.h>

//
// operator=
//
// assignment operator
//
const Neighborhood & Neighborhood::operator= (const Neighborhood &N)
{
  mNum = N.mNum;
  mNeighbor = N.mNeighbor;
  mNeedDepth = N.mNeedDepth;
  mNeighborDepth = N.mNeighborDepth;
  return(*this);
}


//
// addNeighbor
//
// Add the point with index nbr to the neighborhood
//
int Neighborhood::addNeighbor(int nbr, u_char depth, bool check)
{
  int  i;

  // if this is the first neighbor
  if(mNeighbor.isEmpty()) {
    mNeighbor.setDim(20);
    mNeighbor[0] = nbr;
    if(mNeedDepth) {
       mNeighborDepth.setDim(20);
       mNeighborDepth[0] = depth;
    }
    mNum=1;
  } 
  else {
    if(check && isNeighbor(nbr))
	return(false);

    // reallocate if necessary
    if(mNum >= mNeighbor.getNelm()) {
       Array1D<int> tempInt = mNeighbor;
       mNeighbor.setDim(mNum+20);
       for(i=0; i<mNum; i++) 
	 mNeighbor[i] = tempInt[i];
       if(mNeedDepth) {
          Array1D<u_char> tempChar = mNeighborDepth;
          mNeighborDepth.setDim(mNum+20);
          for(i=0; i<mNum; i++) 
	    mNeighborDepth[i] = tempChar[i];
       }
    }
    mNeighbor[mNum] = nbr;
    if(mNeedDepth) mNeighborDepth[mNum] = depth;
    mNum++;
  }
  return(true);
}

//
// neighborIndex
//
int Neighborhood::neighborIndex(int nbr) const
{
  for (int i=0; i<mNum; i++)
    if (mNeighbor[i] == nbr) 
	return(i);

  return(-1);
}

//
// removeNeighbor
//
// removed the specified point from this neighborhood
//
bool Neighborhood::removeNeighbor(int nbr)
{
  int i,j;

  i = neighborIndex(nbr);
  if(i<0) return(false); // not in list

  // remove neighbor from list
  Array1D<int> tempInt = mNeighbor;

  for (j=0;j<i;j++)
    mNeighbor[j] = tempInt[j];
  for (j=i;j<mNum-1;j++)
    mNeighbor[j] = tempInt[j+1];

  if(mNeedDepth) {
    Array1D<u_char> tempChar = mNeighborDepth;

    for (j=0;j<i;j++)
      mNeighborDepth[j] = tempChar[j];
    for (j=i;j<mNum-1;j++)
      mNeighborDepth[j] = tempChar[j+1];
  }
  mNum--;
  return(true);
}


void Neighborhood::cleanNeighborhood(int numVerts) {
	Array1D<int> tempNghbs;
	tempNghbs.setDim(mNum);

	Array1D<u_char> flag;
	flag.setDim(numVerts);
	flag=0;

	int unique=0;
	int i;
	for (i=0; i<mNum; ++i) {
		if (flag[mNeighbor[i]]==0) {
			tempNghbs[unique++]=mNeighbor[i];
			flag[mNeighbor[i]]=(u_char)1;
		}
	}
	mNum=unique;
	mNeighbor.setDim(unique);
	for (i=0; i<unique; ++i)
		mNeighbor[i] = tempNghbs[i];
}

//
// print out a Neighborhood list
//
void Neighborhood::print(ostream &os) const
{
  os << "{ ";
  for (int i=0; i<mNum; i++)
    os << mNeighbor[i] << " ";
  os << "}";
}

//
// operator<<
// non-member print function
//
ostream & operator<< (ostream &os, Neighborhood const &N)
{
  N.print(os);
  return os;
}
