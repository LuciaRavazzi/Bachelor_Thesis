#ifndef _INTEGRANDA_
#define _INTEGRANDA_

#include "funzioneBase.h"
#include <cmath>

//LA CLASSE È UNA POSSBILE REALIZZAZIONE DI FUNZIONE BASE.

class integranda:public funzioneBase{
	public:
	double eval(double x, double n, int a) const {return pow(x,a)*n;} //n sarà estratto da una distribuzione esponenziale.
};
#endif
