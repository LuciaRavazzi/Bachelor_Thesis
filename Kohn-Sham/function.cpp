#include "function.h"
#include "my_nr.h"

#include <iostream>
#include <fstream>
#include <cmath>

Function::Function(const Vec_DP &x):X(x), Y(0.,x.size()), Y2(x.size()), Yp(x.size())
	{
	Set();
	}

Function::Function(const Vec_DP &x, const Vec_DP &y, double yp1, double ypn): X(x), Y(y), Y2(x.size()), Yp(x.size())
	{
	Set(yp1,ypn);
	}

Function::Function(const Function &f): X(f.X), Y(f.Y), Y2(f.size()), Yp(f.size())
	{
	for (unsigned int i=0; i<f.size(); ++i)
		{
		Y2[i]=f.Y2[i];	Yp[i]=f.Yp[i];
		}
	}
    
    

void Function::Set_Y(const Vec_DP &y, double yp1, double ypn)
	{
	Y=y;
	Set(yp1,ypn);
	}

void Function::Set(double yp1, double ypn)
	{
	NR::spline(X,Y,yp1,ypn,Y2);
	NR::splineP(X,Y,Y2,Yp);
	}
	
    
    
    
unsigned int Function::size() const
	{
	return X.size();
	}

Vec_DP	Function::Get_X() const
	{
	return X;
	}
 /*Vec_DP	Function::Get_Y() const
	{
	return Y;
	}
 */   
    
    
    

double Function::Value(const double x) 	const 	
	{
	double y=0; 
	NR::splint (X,Y,Y2,x,y); 
    return y;
	}
double Function::D1(const double x) 	const 	{return NR::splint1(X,Y,Y2,x);}
double Function::D2(const double x) 	const	{return NR::splint2(X,Y,Y2,x);}
double Function::Primitive(const unsigned int n) const	{return Yp[n];}





double Function::Normalize()
	{
	//La normalizzazione avviene rispetto al quadrato della funzione d'onda e considerando
	//una distribuzione di carica radiale. (Niente fattore moltiplicativo legato all'angolo solido)
	Function f(X);
	for(unsigned int i=0; i<f.size(); ++i) f[i]=Value(X[i])*Value(X[i]);
	f.Set();
	
	double val=f.Primitive(f.size()-1)-f.Primitive(0); //f.Primitive(0)=0, quindi non sevre a niente sottrarlo
	double C_n=1;  //Meglio -1 che indica errore?
	if(val!=0)
		{
		C_n=1./sqrt(val);
		for(unsigned int i=0; i<size(); ++i) Y[i]=C_n*Y[i];
		Set();
		}
	return C_n;
	}







//Stampa la mesh e i valori della funzione nei punti della mesh
void Function::Print(const char* file) const
	{
	ofstream output;	output.open(file);
	for (unsigned int i=0; i<X.size(); ++i) output<<X[i]<<"	"<<Y[i]<<endl;
	output.close();
	}


//Overloading di Print()
//Stampa usando la spline dal primo punto della mesh X[0] a max
//Il passo p è di default 0.001
void Function::Print(const char* file, double max, double p) const
	{
	ofstream output;	output.open(file);
	double x=X[0];
	do
		{
		output<<x<<" "<<Value(x)<<endl;
		x+=p;
		}while(x<max);
	output.close();
	}



double Function::operator() (const double x) const {return Value(x);}
double& Function::operator[](const unsigned int n) {return Y[n];}



void Function::operator= (Function &f)
	{
	for (unsigned int i=0; i<f.size(); ++i)
		{
		X[i]=f.X[i]; Y[i]=f.Y[i]; Y2[i]=f.Y2[i]; Yp[i]=f.Yp[i];
		}
	}


    
//Prodotto per una valore reale
Function operator*(const double val, Function &f)
	{
	Vec_DP Y(f.size());
	for (unsigned int i=0; i<f.size(); ++i) Y[i]=val*f[i];
	Function r(f.Get_X(), Y);
	return r;
	}
    
    
    
Function operator+(Function &s, Function &d)
	{
	Vec_DP Y(s.size());
	for (unsigned int i=0; i<s.size(); ++i) Y[i]=s[i]+d[i];
	Function r(s.Get_X(), Y);
	return r;
	}
    
Function operator-(Function &s, Function &d)
	{
	Vec_DP Y(s.size());
	for (unsigned int i=0; i<s.size(); ++i) Y[i]=s[i]-d[i];
	Function r(s.Get_X(), Y);
	return r;
	}
    
    
    
    
    
