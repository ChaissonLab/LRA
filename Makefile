all:	lra alchemy2 tag
PROF=/home/cmb-16/mjc/shared/lib/
CCOPTS_BASE=-std=c++14 
DEBUG?=""
ifneq ($(DEBUG), "")
CCOPTS=$(CCOPTS_BASE) $(DEBUG)
else
CCOPTS=-O3 $(CCOPTS_BASE)  -DNDEBUG -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0
endif
STATIC=
ifeq ($(OPT), "1")
CCOPTS=-g $(CCOPTS_BASE) -lprofiler
STATIC=
endif


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
  AffineOneGapAlign.h \
<<<<<<< HEAD
  NaiveDP.h
=======
  MergeSplit.h
>>>>>>> 5df421ef36b52b71d17225fcdf8a9c67658749d7


htslib/lib/libhts.a:
	cd htslib && autoheader && autoconf && ./configure --disable-s3 --disable-lzma --disable-bz2 --prefix=$(PWD)/htslib/ && make -j 4 && make install

tag: TestAffineOneGapAlign.cpp AffineOneGapAlign.h
	g++ -g TestAffineOneGapAlign.cpp -o tag  -D _MAT_PRINT_

lra: lra.o
<<<<<<< HEAD
	g++ $(STATIC) $(CCOPTS) $^  -L htslib/lib -lhts -lz -lcurl -lpthread -o $@
=======
	g++ $(STATIC) $(CCOPTS) $^  -L htslib/lib -lhts -lz -lpthread -o $@ 
>>>>>>> 5df421ef36b52b71d17225fcdf8a9c67658749d7

alchemy2: Alchemy2.o
	g++ $(STATIC) $(CCOPTS) $^  -L htslib/lib -lhts -lz -lpthread -o $@

qti: QueryTime.o
	g++ $(STATIC) $(CCOPTS) $^  -L htslib/lib -lhts -lz -lpthread -o $@

lra.o: lra.cpp $(HEADERS) htslib/lib/libhts.a
	g++ $(CCOPTS) -c  -I htslib/include -I seqan/include  lra.cpp

Alchemy2.o: Alchemy2.cpp  htslib/lib/libhts.a
	g++ $(CCOPTS) -c  -I htslib/include -I seqan/include  Alchemy2.cpp

QueryTime.o: QueryTime.cpp $(HEADERS) htslib/lib/libhts.a
	g++ $(CCOPTS) -c  -I htslib/include -I seqan/include  QueryTime.cpp


clean:
	rm -f lra lra.o


