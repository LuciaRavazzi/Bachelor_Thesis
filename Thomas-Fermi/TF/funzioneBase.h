#ifndef _FUNZIONEBASE_
#define _FUNZIONEBASE_


//TALE CLASSE RAPPRESENTA UNA GENERICA FUNZIONE.

class funzioneBase {
	public:
	virtual double eval(double x, double n, int a) const=0;
};

#endif
