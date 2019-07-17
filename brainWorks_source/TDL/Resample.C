///////////////////////////////////////////////////////////////////////////
//
// File: Resample.C
//
// Author: Guoling Tong
//
// Purpose: Implementation of the Resample class, which provides a uniform
//          interface for various resampling.
// 
///////////////////////////////////////////////////////////////////////////

#include <TDL/Resample.h>
#include <TDL/Convolution.h>

template <class T>
void Resample<T>::resample(typename Convolution<T>::ConvolutionType flag,T const *volIn,
			   int x1,int y1,int z1,T *volOut,int x2,int y2,int z2)
{

  Convolution<T>::convolution(flag,volIn,x1,y1,z1,volOut,x2,y2,z2);

}


#ifdef GNU_COMPILER
template class Resample<unsigned char>;
template class Resample<unsigned short>;
template class Resample<short>;
template class Resample<int>;
template class Resample<float>;
#endif

#ifdef SGI_COMPILER
#pragma instantiate Resample<unsigned char>
#pragma instantiate Resample<unsigned short>
#pragma instantiate Resample<short>
#pragma instantiate Resample<int>
#pragma instantiate Resample<float>
#endif

#ifdef RS6K
#pragma define(Resample<unsigned char>)
#pragma define(Resample<unsigned short>)
#pragma define(Resample<short>)
#pragma define(Resample<int>)
#pragma define(Resample<float>)
#endif
