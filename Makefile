all:	lsa
CCOPTS=-O2
lsa: lsa.o
	g++ $(CCOPTS) $^  -L htslib/lib -lhts -lz -o $@ 

lsa.o: lsa.cpp
	g++ $(CCOPTS) -c  -I htslib/include  lsa.cpp 

clean:
	rm -f lsa test *.o


