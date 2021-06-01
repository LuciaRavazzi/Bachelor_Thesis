#include "my_nr.h"

void NR::splineP(Vec_I_DP &x, Vec_I_DP &y, Vec_I_DP &y2, Vec_O_DP &yp)
	{
	yp[0]=0.; //tutte le mie primitive sono nulle nell'origine
	for (unsigned int i=0; i<x.size()-1; ++i)
		{
		double h=x[i+1]-x[i];
		yp[i+1]=yp[i]+(y[i+1]+y[i])*h/2.;
		}
	}
	
double NR::splint1(const Vec_DP &xa, const Vec_DP &ya, const Vec_DP &y2a, const double x)
{
	int k;
	double h,b,a;

	int n=xa.size();
	int klo=0;
	int khi=n-1;
	while (khi-klo > 1) 
		{
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
		}
	h=xa[khi]-xa[klo];
	if (h == 0.0) nrerror("Bad xa input to routine splint");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	return  ( (ya[khi]-ya[klo])/h -(3*a*a-1)*h*y2a[klo]/6 +(3*b*b-1)*h*y2a[khi]/6 );
}

double NR::splint2(const Vec_DP &xa, const Vec_DP &ya, const Vec_DP &y2a, const double x)
{
	int k;
	double h,b,a;

	int n=xa.size();
	int klo=0;
	int khi=n-1;
	while (khi-klo > 1) 
		{
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
		}
	h=xa[khi]-xa[klo];
	if (h == 0.0) nrerror("Bad xa input to routine splint");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	return ( a*y2a[klo] +b*y2a[khi] );
}

/*
DP fshoot(const DP E)
	{
	Vec_DP v(1); 
	Vec_DP f(1); 
	v[0]=E;
	NR::shoot(v,f);
	return f[0];
	}
*/
