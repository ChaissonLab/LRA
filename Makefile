all:	lsa
CCOPTS=-g -std=c++14 

htslib/lib/libhts.a:
	cd htslib && autoheader && autoconf && ./configure --prefix=$(PWD)/htslib/ && make -j 4

lsa: lsa.o
	g++ -static $(CCOPTS) $^  -L htslib/lib -lhts -lz -o $@ 

lsa.o: lsa.cpp MinCount.h CompareLists.h TupleOps.h Sorting.h MMIndex.h Options.h Clustering.h Genome.h
	g++ $(CCOPTS) -c  -I htslib/include  lsa.cpp 

clean:
	rm -f lsa test *.o


