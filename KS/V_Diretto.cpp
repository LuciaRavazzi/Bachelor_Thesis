#include "V_Diretto.h"
#include <cmath>

extern double *D;

V_Diretto::V_Diretto(const Vec_DP &x):Function(x){}

V_Diretto::V_Diretto(const Density & density):Function(density.Get_X()){
 Evaluate(density);
}

V_Diretto::V_Diretto(const Vec_DP &x, const Vec_DP &y, double yp1, double ypn):Function(x,y,yp1,ypn){}


void V_Diretto::Evaluate(const Density &density){
//Potenziale diretto nei punti della mesh radiale.
 Vec_DP sum(0.,size());

  for(int i=0; i<X.size(); ++i) {
        for(int k=0; k<i; ++k){ 
            sum[i]+=D[k]/X[i]*density.Y[k]; 
            //cout << "density.Y[i] " << density.Y[i] << "sum[i] " << sum[i] << endl; 
        }  
        for(int k=i; k<X.size(); ++k)
        	sum[i]+=D[k]/X[k]*density.Y[k];
     }

    Y=sum;

    Set();
}

void V_Diretto::operator= (const V_Diretto &f){
 for (unsigned int i=0; i<f.size(); ++i) {
  X[i]=f.X[i];
  Y[i]=f.Y[i]; 
  Y2[i]=f.Y2[i]; 
  Yp[i]=f.Yp[i];
 }
}

V_Diretto operator*(const double val, V_Diretto &f){
 Vec_DP Y(f.size());
 for (unsigned int i=0; i<f.size(); ++i) Y[i]=val*f[i];
 V_Diretto r(f.Get_X(), Y);
return r;
}

V_Diretto operator+(V_Diretto s, V_Diretto d){
 Vec_DP Y(s.size());
 for (unsigned int i=0; i<s.size(); ++i) Y[i]=s[i]+d[i];
 V_Diretto r(s.Get_X(), Y);
return r;
}

V_Diretto operator-(V_Diretto s, V_Diretto d){
 Vec_DP Y(s.size());
 for (unsigned int i=0; i<s.size(); ++i) Y[i]=s[i]-d[i];
 V_Diretto r(s.Get_X(), Y);
 return r;
}

double Get_err(V_Diretto &s, V_Diretto &d){
 double err=0;
 for (unsigned int i=0; i<s.size(); ++i) err+=(s[i]-d[i])*(s[i]-d[i]);
 return sqrt(err/s.size());
}
