#ifndef _INTEGRALE_
#define _INTEGRALE_

#include "funzioneBase.h"
#include "nr.h"

class integrale{
	public:
	integrale(double a,double b,funzioneBase *f);
	double trapezi(Vec_DP& passo, double* n, int dim, double r);

	private:
	double _a,_b;
	funzioneBase* _f;
};
#endif
