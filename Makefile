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
  Alignment.h \
  Read.h \
  MapRead.h \
  Input.h \
  SimpleGlobalChain.h \
  AffineOneGapAlign.h


htslib/lib/libhts.a:
	cd htslib && autoheader && autoconf && ./configure --disable-lzma --disable-bz2 --prefix=$(PWD)/htslib/ && make -j 4 && make install

tag: TestAffineOneGapAlign.cpp AffineOneGapAlign.h
	g++ -g TestAffineOneGapAlign.cpp -o tag 
# -D _MAT_PRINT_

lra: lra.o
	g++ -static $(CCOPTS) $^  -L htslib/lib -lhts -lz -lpthread -o $@

lra.o: lra.cpp $(HEADERS) htslib/lib/libhts.a
	g++ $(CCOPTS) -c  -I htslib/include -I seqan/include  lra.cpp 


clean:
	rm -f lra lra.o


