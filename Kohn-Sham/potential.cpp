#include "potential.h"

#include <cmath>

extern double *D;


Potenziale::Potenziale(const Vec_DP &x):Function(x){}

Potenziale::Potenziale(const Density & density, const Vec_DP &nucleo):Function(density.Get_X()){
 Evaluate(density, nucleo);
}


void Potenziale::Evaluate(const Density &density, const Vec_DP &nucleo){
 p=1.00;
 A=0.031091;
 alpha1=0.21370;
 beta1=7.5957;
 beta2=3.5876;
 beta3=1.6382;
 beta4=0.49294;    


 Vec_DP sum(0.,size());
 Vec_DP app(0.,size());
 Vec_DP app2(0.,size());
 double rs=0.;

    
//Potenziale diretto
 //cout << "Potenziale diretto" << endl;
    
    for(int i=0; i<X.size(); ++i) {
        for(int k=0; k<i; ++k) 
            sum[i]+=D[k]/X[i]*density.Y[k];    
        for(int k=i; k<X.size(); ++k)
        	sum[i]+=D[k]/X[k]*density.Y[k];
     }



//Potenziale del nucleo esteso
 for(int i=0; i<X.size(); ++i){ 
  sum[i]+=nucleo[i];
  //cout << "Potenziale diretto+nucleo " << sum[i] << endl;
 }

//Potenziale di scambio.    
 for(int i=0; i<X.size(); ++i) {
   double appo=density.Y[i]/(4.*M_PI*pow(X[i],2.));
   sum[i] += -pow((3.*appo)/M_PI,1./3.);
 }

//Potenziale di correlazione.
  double appoggio=0.;
  for(int i=0; i<X.size(); ++i){
   if(density.Y[i]!=0){
    double appo=density.Y[i]/(4.*M_PI*pow(X[i],2.));
    rs=pow(3./(4.*M_PI*appo),1./3.); //Se le densità è nulla non è definito.
    double appo1=1.+alpha1*rs;
    app[i]=2.*A*(beta1*sqrt(rs)+beta2*rs+beta3*pow(rs,3./2.)+beta4*pow(rs,p+1.));
    double appo2=1.+1./app[i];
    //cout << appo2 << endl;
    app2[i]=2.*A*((beta1/(2.*sqrt(rs)))+beta2+beta3*1.5*sqrt(rs)+beta4*(p+1.)*pow(rs,p)); 
    appoggio += appo*(-2.*A*alpha1*log(appo2)-2.*A*appo1*(1./appo2)*app2[i]*(-1./pow(app[i],2.)))*(-1./3.)*pow(3./(4.*M_PI),1./3.)*pow(appo,-4./3.)-2.*A*appo1*log(appo2); //il potenziale di correlazione nei punti della mesh.
   //cout << i << " " << pow(appo,1) << endl;
   if(pow(appo,-4./3.) > pow(10, 70)){
     sum[i] += 0.;
   } else {
     sum[i] += appoggio;
   }
   } else {
    sum[i] += 0.;
   }
  }
 

//MODIFICA	C.C.
    // sum[X.size()-1]+= abs(nucleo[0])*1000;
    //sum[X.size()-2]+= abs(nucleo[1]);
    
    Y= sum;

    Set();
}




Potenziale::Potenziale(const Vec_DP &x, const Vec_DP &y, double yp1, double ypn):Function(x,y,yp1,ypn){}



void Potenziale::operator= (const Potenziale &f){
	for (unsigned int i=0; i<f.size(); ++i) {
		X[i]=f.X[i];
		Y[i]=f.Y[i]; 
		Y2[i]=f.Y2[i]; 
		Yp[i]=f.Yp[i];
	}
}

Potenziale operator*(const double val, Potenziale &f)
	{
	Vec_DP Y(f.size());
	for (unsigned int i=0; i<f.size(); ++i) Y[i]=val*f[i];
	Potenziale r(f.Get_X(), Y);
	return r;
	}
Potenziale operator+(Potenziale s, Potenziale d)
	{
	Vec_DP Y(s.size());
	for (unsigned int i=0; i<s.size(); ++i) Y[i]=s[i]+d[i];
	Potenziale r(s.Get_X(), Y);
	return r;
	}
Potenziale operator-(Potenziale s, Potenziale d)
	{
	Vec_DP Y(s.size());
	for (unsigned int i=0; i<s.size(); ++i) Y[i]=s[i]-d[i];
	Potenziale r(s.Get_X(), Y);
	return r;
	}
	
double Get_err(Potenziale &s, Potenziale &d)
	{
	double err=0;
	for (unsigned int i=0; i<s.size(); ++i) err+=(s[i]-d[i])*(s[i]-d[i]);
	return sqrt(err/s.size());
	}
    
    
    
    
