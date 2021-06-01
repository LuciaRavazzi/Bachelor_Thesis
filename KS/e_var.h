
#ifndef e_var_h
#define e_var_h

#include "funzioni.h"
#include "function.h"
#include "wave_function.h"
#include "density.h" 
#include "potential.h"
#include "V_Diretto.h"

double Evaluate_E_var1(Density &d, Potenziale &V, ofstream &out);	//Valuta l'energia variazionale

double Evaluate_E_var2(Density &d,V_Diretto &V, ofstream &out,int,int,int);

#endif
