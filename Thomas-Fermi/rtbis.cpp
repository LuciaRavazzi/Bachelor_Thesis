#include <cmath>
#include "nr.h"
using namespace std;

DP NR::rtbis(DP func(const DP), const DP x1, const DP x2, const DP xacc, const DP N)
{
	const int JMAX=40;
	int j;
	DP dx,f,fmid,xmid,rtb;

	f=func(x1)-N;
	fmid=func(x2)-N;
	if (f*fmid >= 0.0){
	  nrerror("Root must be bracketed for bisection in rtbis");
	  cout << "f=func(x1) " << f << " fmid=func(x2) " << fmid << endl;
        }
	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	for (j=0;j<JMAX;j++) {
		fmid=func(xmid=rtb+(dx *= 0.5))-N;
		if (fmid <= 0.0) rtb=xmid;
		if (fabs(dx) < xacc || fmid == 0.0) return rtb;
	}
	nrerror("Too many bisections in rtbis");
	return 0.0;
}
