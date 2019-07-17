#ifndef _SurfaceTracker_h_
#define _SurfaceTracker_h_

#include <OS.h>
#include <ADL/Array1D.h>
#include <TDL/IDLdefines.h>

class Surface;

class SurfaceTracker {
public:
	enum TrackType {
		Sulci, Gyri, Geodesi, SameType
		};

	 SurfaceTracker();
	~SurfaceTracker();

	void setSurface(Surface *S, TrackType t = SurfaceTracker::SameType);
	void setTrackType(TrackType t);
	ItXECode findPath(int startindx, int endindx, Array1D<int> &retindexes);
	void clean();
protected:
	Surface   *mSurface;
	TrackType  mTrackType;
	int        mNumPoints;
	int        threshold;
	void      *vnodes;
	float    **cost;
	int       *mask;
	int       *work_space1;
	int       *work_space2;
	double     maxcurv;
	double     mincurv;

	void setupCost();
private:
	};

#endif
