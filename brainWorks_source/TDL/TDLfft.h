#ifndef _TDLfft_h_
#define _TDLfft_h_

#include <ADL/Array3D.h>

#define FFT_MAX_DIM 4096

typedef struct {
        float re;
        float im;
} bwcomplex;


class TDLfft {
public:
       enum FFTtype { Forward = -1, Reverse = 1 };

       void pad(Array3D<float> &dat);
       void unpad(Array3D<float> &dat);
       bool FFT(Array3D<float> &dat,FFTtype dir);
private:
       void fftr(int,FFTtype);
       void fftl(int,FFTtype);
       void scale(Array3D<float>&);
       bwcomplex wk[FFT_MAX_DIM+1];
};


#endif
