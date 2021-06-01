#include "density.h"

#include <fstream>
#include <cmath>
#include <algorithm>
#include "potential.h"
#include "function.h"

using namespace std;

//|r*R(r)|^2

extern int N;
extern Dati dati;
extern double *rad_D;
extern double *rad_inv_D;
extern double *G_inv;

extern int lwork;
extern double *work;
extern int info;


extern double *matrice;
extern double *autovalori;

extern double *T;



Density::Density(const Vec_DP &x, unsigned int Z, unsigned int N): Z_(Z), N_(N), n_autofunzioni(0), Function(x)
	{
	Def_wave_function();
	Evaluate_Y();
	}
    

//funzione che in base a Z sceglie quanti orbitali creare.
//nmax= Get_periodo()+1 per memorizzare anche alcuni orbitali vuoti (posso comunque plottarli)
void Density::Def_wave_function()
{
    nmax= Get_periodo()+1;
    lmax= Get_periodo();
	
    for(unsigned int l=0; l<=lmax; ++l)
    	for(int n=l+1; n<=nmax; ++n){
            psi.push_back(Wave_function(X, n, l, 1)); //Ogni orbitale ha degenerazione 2l+1 
            psi.push_back(Wave_function(X, n, l, 0));
        }
    n_autofunzioni=psi.size();	
}	

unsigned int Density::Get_periodo() const
{
    unsigned int Z=Get_Z();
    unsigned int periodo=0;
    if (Z>118) cerr<<"Hai inserito un numero atomico troppo grande"<<endl;
    else
    {
        if (Z>86) periodo=7;
        else
        {
            if (Z>54) periodo=6;
            else
            {
                if (Z>36) periodo=5;
                else
                {
                    if (Z>18) periodo=4;
                    else
                    {
                        if (Z>10) periodo=3;
                        else
                        {
                            if (Z>2) 	periodo=2;
                            else 		periodo=1;
                        }
                    }
                }
            }
        }
    }
    return periodo;
}

void Density::Evaluate_Y_int(){

 Put_e();

 for (unsigned int j=0; j<Y.size(); ++j){
  Y[j]=exp(-X[j]);
 }

 Set();
}


//Y sono i valori della densità di probabilità totale nei punti della mesh
void Density::Evaluate_Y()
	{
	Put_e();			//Ordina le subshell e le riempie
				
    for (unsigned int j=0; j<Y.size(); ++j)
		{
		double v_el=0;	
		for (unsigned int i=0; i<n_autofunzioni; ++i)
			v_el+= ( psi[i][j]*psi[i][j] * psi[i].Get_n_pieni() ); //=|P(r)*Y(theta,phi)|^2 con P(r)=r*R(r)
            
		Y[j]=v_el;
        }
	
	Set();	//Interpola i dati
	}


//funzione che inserisce gli "elettroni" all'interno delle subshell
void Density::Put_e()
	{
	sort(psi.begin(), psi.end());		//ordine le funzioni d'onda in base all'energia
	unsigned int n_e=N_;				//variabile che conteggia il numero di elettroni ancora da inserire
	
	for (unsigned int i=0; i<n_autofunzioni;++i)
		{
		if (n_e==0) psi[i].Set_n_pieni(0);
		else
			{
			unsigned int max=2*psi[i].Get_l()+1;
			if (n_e>=max)
				{
				psi[i].Set_n_pieni(max); 
				n_e-=max;
				}
			else
				{
				psi[i].Set_n_pieni(n_e);
				n_e=0;
				}
			}
		cout << psi[i].Get_n() << " " << psi[i].Get_l() << " " << psi[i].Get_n_pieni() << endl;
		}
	}

/*
Fino a questo punto è tutto gestito dal costruttore di density
*/



void Density::Print_EigenValue(const char* file) const{
 ofstream out; 
 out.open(file);
 out<<"Ecco l'elenco degli autovalori"<<endl;
 out.close();
 for (unsigned int i=0; i<n_autofunzioni; ++i) psi[i].Print_EigenValue(file);
}
    
    
    
void Density::Print_EigenValue(ofstream &out) const {
 for (unsigned int i=0; i<n_autofunzioni; ++i)
  psi[i].Print_EigenValue(out);
}
    
    

void Density::Print_EigenFunction(unsigned int n_ciclo, double max) const
	{
	for (unsigned int i=0; i<n_autofunzioni; ++i) 
    	psi[i].Print_EigenFunction(n_ciclo, max);
	}
    
void Density::Print_EigenFunction(double max) const
	{
	for (unsigned int i=0; i<n_autofunzioni; ++i) 
    	psi[i].Print_EigenFunction(max);
	}
    
unsigned int Density::Get_Z() const
	{
	return Z_;
	}

unsigned int Density::Get_N() const
	{
	return N_;
	}
    
    
unsigned int Density::Get_n_autofunzioni() const
	{
    return n_autofunzioni;
    }


//Le due Density hanno per costruzione stessi valori di Z_, N_, n_autofunzioni
void Density::operator= (const Density &f){

	for (unsigned int i=0; i<f.size(); ++i)
		{
		X[i]=f.X[i]; Y[i]=f.Y[i]; Y2[i]=f.Y2[i]; Yp[i]=f.Yp[i];
		}
        
	for (unsigned int i=0; i<n_autofunzioni; ++i) 
    	psi[i]= f.psi[i];
	
}

Density::Density(const Vec_DP &x, const Vec_DP &y, double yp1, double ypn):Function(x,y,0,0){}

Density operator*(const double val, Density &f){
 Vec_DP Y(f.size());
 for (unsigned int i=0; i<f.size(); ++i) Y[i]=val*f[i];
 Density r(f.Get_X(), Y);
return r;
}

Density operator+(Density s, Density d){
 Vec_DP Y(s.size());
 for (unsigned int i=0; i<s.size(); ++i) Y[i]=s[i]+d[i];
 Density r(s.Get_X(), Y);
return r;
}

Density operator-(Density s, Density d){
 Vec_DP Y(s.size());
 for (unsigned int i=0; i<s.size(); ++i) Y[i]=s[i]-d[i];
 Density r(s.Get_X(), Y);
 return r;
}


unsigned int Density::Get_nmax()
	{
    return nmax;
    }
    
unsigned int Density::Get_lmax()
	{
    return lmax;
    }
