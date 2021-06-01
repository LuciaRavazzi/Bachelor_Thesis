#ifndef potential_h
#define potential_h

#include "nr.h"
#include "function.h"
#include "density.h"

using namespace std;

//Classe che descrive il potenziale attrattivo del nucleo Vn.

class Potenziale: public Function {
 public:
  //costruttore di potenziale nullo
  Potenziale(const Vec_DP &x);   					    
  //costruttore di potenziale del nucleo 
  Potenziale(const Density &, const Vec_DP &nucleo);        	     
  //costruttore con valori di potenziale già noti (tutto già fatto)      
  Potenziale(const Vec_DP &x, const Vec_DP &y, double yp1=0, double ypn=0); 
        
  void Evaluate(const Density &, const Vec_DP &nucleo);							
  void operator= (const Potenziale &f);
  void SetParametri();
  private:
  double p,A,alpha1,beta1,beta2,beta3,beta4;

};

Potenziale operator*(const double, Potenziale &);
Potenziale operator+(Potenziale , Potenziale );
Potenziale operator-(Potenziale , Potenziale );

//restituisce la differenza fra 2 potenziali
double Get_err(Potenziale &s, Potenziale &d); 


#endif

