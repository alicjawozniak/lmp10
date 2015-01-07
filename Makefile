aprox: main.o splines.o points.o aproksymator_na_bazie.o gaus/libge.a
	$(CC) -o aprox  main.o splines.o points.o aproksymator_na_bazie.o -L gaus -l ge

intrp: main.o splines.o points.o interpolator.o gaus/libge.a
	$(CC) -o intrp  main.o splines.o points.o interpolator.o -L gaus -l ge

prosta: main.o splines.o points.o prosta.o
	$(CC) -o prosta  main.o splines.o points.o prosta.o	

aproksymator_na_bazie.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c aproksymator_na_bazie.c

interpolator.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c interpolator.c

hermite: main.o splines.o points.o hermite.o gaus/libge.a
	$(CC) -o herm  main.o splines.o points.o hermite.o -L gaus -l ge 

hermite.o: makespl.h points.h gaus/piv_ge_solver.h hermite.c
	$(CC) -I gaus -c hermite.c

.PHONY: clean

clean:
	-rm *.o aprox intrp prosta herm
