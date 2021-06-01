#ifndef my_nr_h
#define my_nr_h

/*In questo file raccolgo le mie aggiunte alle Numerical Recipies
- interpolazione della funzione;
- accorgimenti per addattare le interfacce delle funzioni
*/


#include "nr.h"
using namespace std;

namespace NR 
{
//Spline per tabularsi i valori della primitiva
void splineP(Vec_I_DP &x, Vec_I_DP &y, Vec_I_DP &y2, Vec_O_DP &yp);

//Splint per la derivata prima e seconda (non vengono sfruttate nel mio algoritmo)
double splint1(const Vec_DP &xa, const Vec_DP &ya, const Vec_DP &y2a, const double x);
double splint2(const Vec_DP &xa, const Vec_DP &ya, const Vec_DP &y2a, const double x);
}

/*
void derivs(const double x, Vec_I_DP &y, Vec_O_DP &r);
void load(const double x, Vec_I_DP &v, Vec_O_DP &y);
void score(const double xf, Vec_I_DP &y, Vec_O_DP &f);

DP fshoot(const DP);	//funzione che "nasconde" lo shooting per inserirlo all'interno della bisezione
*/

#endif
