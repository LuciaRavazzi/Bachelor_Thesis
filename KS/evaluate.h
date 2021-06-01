
#ifndef evaluate_h
#define evaluate_h

#include "nr.h"
#include "lapack.h"
#include "funzioni.h"
#include "function.h"
#include "wave_function.h"
#include "my_nr.h"
#include "density.h"
#include "potential.h"


//Calcola la norma al quadrato di un Vec_DP
double norm2(Vec_DP v);

double Evaluate_wave_function(Density &d, Potenziale &potenziale, ofstream &out2, ofstream &out3);


#endif
