all:	lra 
PROF=/home/cmb-16/mjc/shared/lib/
CCOPTS=-O3 -std=c++14 
#-D _TESTING_ -lprofiler 
#  -L$(PROF) 

HEADERS=MinCount.h \
  CompareLists.h \
  TupleOps.h \
  Sorting.h \
  MMIndex.h \
  Options.h\
  Clustering.h \
  Genome.h \
  SparseDynamicProgramming.h \
  Alignment.h \
  GlobalChain.h \
  Read.h \
  MapRead.h \
  Input.h


htslib/lib/libhts.a:
	cd htslib && autoheader && autoconf && ./configure --prefix=$(PWD)/htslib/ && make -j 4

lra: lra.o
	g++ -static $(CCOPTS) $^  -L htslib/lib -lhts -lz -lpthread -o $@

lra.o: lra.cpp $(HEADERS)
	g++ $(CCOPTS) -c  -I htslib/include -I seqan/include  lra.cpp 


clean:
	rm -f lsa test *.o


