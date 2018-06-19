all:	lsa
CCOPTS=-O2
lsa: lsa.o
	g++  $(CCOPTS) $^  -L htslib/lib -lhts -lz -o $@ 

lsa.o: lsa.cpp MinCount.h CompareLists.h TupleOps.h Sorting.h MMIndex.h Options.h Clustering.h
	g++ $(CCOPTS) -c  -I htslib/include  lsa.cpp

clean:
	rm -f lsa test *.o


