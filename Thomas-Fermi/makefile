

tesi.exe: tesi.o integrale.o funzioni.o rtbis.o
	g++ tesi.o integrale.o funzioni.o rtbis.o -o tesi.exe

tesi.o: main.cpp integrale.h funzioneBase.h integranda.h funzioni.h main.h 
	g++ -c  main.cpp -o tesi.o 

integrale.o:integrale.cpp funzioneBase.h integrale.h nr.h
	g++ -c integrale.cpp -o integrale.o

funzioni.o: funzioni.cpp funzioni.h nr.h
	g++ -c funzioni.cpp -o funzioni.o

rtbis.o: rtbis.cpp nr.h
	g++ -c rtbis.cpp -o rtbis.o













