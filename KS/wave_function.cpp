#include "wave_function.h"
#include "my_nr.h"
#include "funzioni.h"

#include <cmath>
#include <fstream>
#include <stdio.h>
#include <iomanip>

using namespace std;
extern unsigned int l_;
extern Dati dati;

//r*R(r)

Wave_function::Wave_function(const Vec_DP &x, unsigned int ne, unsigned int le, int up_or_down, double Energia, unsigned int pieni): n(ne), l(le), E_(Energia), n_pieni(pieni), spin(up_or_down), Function(x), 
    Y_rad_D(x.size()) {}
	
double Wave_function::EigenValue() const
	{
	return E_;		 
	}

double Wave_function::Get_n_pieni() const
	{
	return n_pieni;
	}

unsigned int Wave_function::Get_l() const
	{
	return l;
	}
	
unsigned int Wave_function::Get_n() const
	{
	return n;
	}
unsigned int Wave_function::Get_spin() const
	{
	return spin;
	}
	
void Wave_function::Set_n_pieni(double n_e)
	{
	if (n_e <= (2*l+1) ) n_pieni=n_e;
	else 
		{
		cerr<<"Stai cercando di sovrariempire la subshell"<<endl;
		n_pieni=0;
		}
	}

void Wave_function::Print_EigenValue(const char* file) const{
 ofstream out; out.open(file, ios::app);
 out<<"n= "<<n<<" l= "<<l<<" n_contenuti= "<<n_pieni<<" spin= "<<spin-0.5<<"	E= "<< setprecision(10) <<EigenValue()<<endl;
 out.close();
}
    
	
void Wave_function::Print_EigenValue(ofstream &out) const {
 out<<"    "<< EigenValue();
}
    
    
    
    
void Wave_function::Print_EigenFunction(unsigned int n_ciclo,double max)	const
	{
	char output[50];
	sprintf (output, "risultati/autofunzioni/Ciclo_%dn_%dl_%dms_%d.dat",n_ciclo, n, l, spin);
	
	if (max==0) Print(output);
	else Print(output,max);
	}
    
    
void Wave_function::Print_EigenFunction(double max)	const
	{
	char output[50];
	sprintf (output, "risultati/autofunzioni/n_%dl_%dms_%d.dat", n, l, spin);
	
	if (max==0) Print(output);
	else Print(output,max);
	}

    
    
double Wave_function::Set_E(double E, ofstream &out3)
	{
    out3 << Get_n() << "  " << Get_l() << "  " << E_ << "  " << E <<"       ";
	double err=(E_-E)*(E_-E);
	E_=E;
	return err;
	}
	
bool operator < (const Wave_function &s, const Wave_function &d) 
	{
	return ( s.EigenValue() < d.EigenValue() );
	}
        
    
void Wave_function::Set_l(unsigned int ll){
	l=ll;
}


void Wave_function::Set_n(unsigned int nn){
	n=nn;
} 
    
    
    
void Wave_function::operator= (const Wave_function &f)
	{
	for (unsigned int i=0; i<f.size(); ++i)
		{
		X[i]=f.X[i]; Y[i]=f.Y[i]; Y2[i]=f.Y2[i]; Yp[i]=f.Yp[i];
		}
	n= f.Get_n();
    l= f.Get_l();    
    n_pieni= f.Get_n_pieni();
    E_= f.EigenValue();
    spin= f.Get_spin();
    
}


void Wave_function::Set_Y_rad_D(const Vec_DP &y){
    Y_rad_D=y;
}




        
