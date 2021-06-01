/*
//Test per LAPACK: diagonalize e inverse

double *a= new double[N*N];

double *w= new double[N];


a[0] = 1;
a[1] = 2;
a[2] = 3;
a[3] = 2;
a[4] = 5;
a[5] = -2;
a[6] = 3;
a[7] = -2;
a[8] = -2;



diagonalize(N,a,w,lwork,work,info);


cout << "INFO = " << info << endl;

cout << "Autovalori" << endl;

for(int i=0; i<N; ++i)
cout << w[i] << endl;

cout << "Autovettori" << endl;


for(int i=0; i<N; ++i){
cout << "Autovettore " << i+1 << endl;
for (int j=0; j<N; ++j)
cout << a[i*N+j] << endl;
}


double *b= new double[N*N];

//La matrice è stata sovrascritta con gli autovettori

a[0] = 1;
a[1] = 2;
a[2] = 3;
a[3] = 2;
a[4] = 5;
a[5] = -2;
a[6] = 3;
a[7] = -2;
a[8] = -2;


inverse(N, a, ipiv, b, lwork, work, info);



cout << "INFO = " << info << endl;

cout << "Matrice inversa" << endl;

for(int i=0; i<N; ++i){
cout << endl;
for (int j=0; j<N; ++j)
cout << b[i*N+j] << "  ";
}    

cout << endl << endl;


cout << "b(1,2) = " << b[ind(1,2)] << endl;
*/

