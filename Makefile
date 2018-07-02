all:	lsa testsdp
PROF=/home/cmb-16/mjc/shared/lib/
CCOPTS=-O3 -std=c++14 
#-D _TESTING_ -lprofiler 
#  -L$(PROF) 


htslib/lib/libhts.a:
	cd htslib && autoheader && autoconf && ./configure --prefix=$(PWD)/htslib/ && make -j 4

lsa: lsa.o
	g++ -static $(CCOPTS) $^  -L htslib/lib -lhts -lz -o $@

lsa.o: lsa.cpp MinCount.h CompareLists.h TupleOps.h Sorting.h MMIndex.h Options.h Clustering.h Genome.h SparseDynamicProgramming.h Alignment.h GlobalChain.h
	g++ $(CCOPTS) -c  -I htslib/include -I seqan/include  lsa.cpp 

testsdp: TestSDP.o
	g++ $(CCOPTS) $^ -o $@

TestSDP.o: TestSDP.cpp SparseDynamicProgramming.h
	g++ $(CCOPTS) -c TestSDP.cpp

clean:
	rm -f lsa test *.o


