#include <WarpLDMM.h>
#include <TDL/ItXVolumeUtils.h>
#include "MFArray3D.h"

#ifdef GNU_COMPILER
template class MFArray3D<double>;
#else
/*#error "MFArray3D Instantiation code needed in WarpLDMM.C"*/
#endif

bool WarpLDMM::apply(Surface &S,
		const char *xf,
		const char *yf,
		const char *zf)
{
	MFArray3D<double> mx,my,mz;

	if(!mx.load(xf) || !my.load(yf) || !mz.load(zf))
		return(false);
	double nxi = mx.getNxInner();
	double nyi = my.getNyInner();
	double nzi = mz.getNzInner();

	for(int i=0;i<S.getNumVert();i++) {
		double px = S.vertices()[i].x();
		double py = S.vertices()[i].y();
		double pz = S.vertices()[i].z();

		while(px < 0.)  px += nxi; 
		while(px > nxi) px -= nxi;
		while(py < 0.)  py += nyi; 
		while(pz > nyi) py -= nyi;
		while(pz < 0.)  pz += nzi; 
		while(pz > nzi) pz -= nzi;

		double nx = mx.getTrilinear(pz,py,px);
		double ny = my.getTrilinear(pz,py,px);
		double nz = mz.getTrilinear(pz,py,px);

		S.vertices()[i].set(nx,ny,nz);
	}
	return(true);	
}
