all:  lra alchemy2 tag tgc
#PROF=/home/cmb-16/mjc/shared/lib/
PROF=/home/cmb-16/mjc/jingwenr/software/lib
CCOPTS_BASE=
DEBUG?=""
OPT?=""
ifneq ($(DEBUG), "")
CCOPTS=$(CCOPTS_BASE) $(DEBUG)
else
CCOPTS=-O3 $(CCOPTS_BASE)  -DNDEBUG -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0
endif
STATIC=
ifneq ($(OPT), "")
#CCOPTS=-g $(CCOPTS_BASE) -lprofiler
#STATIC=
STATIC=-g -L $(PROF) -lprofiler
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
  Fragment.h\
  BasicEndpoint.h\
  PrioritySearchTree.h\
  AffineOneGapAlign.h \
  GlobalChain.h \
  NaiveDP.h \
  SparseDP.h

CPP=g++ -std=c++14 

htslib/lib/libhts.a:
	cd htslib && autoheader && autoconf && ./configure --disable-s3 --disable-lzma --disable-bz2 --prefix=$(PWD)/htslib/ && make -j 4 && make install

tag: TestAffineOneGapAlign.cpp AffineOneGapAlign.h
	$(CPP) -g TestAffineOneGapAlign.cpp -o tag 
# -D _MAT_PRINT_

tgc: TestGlobalChain.cpp GlobalChain.h Fragment.h BasicEndpoint.h PrioritySearchTree.h
	$(CPP) -g TestGlobalChain.cpp -o tgc

lra: lra.o
	$(CPP) $(STATIC) $(CCOPTS) $^  -L htslib/lib -lhts -lz -Wl,-rpath=$(PWD)/htslib/lib -lpthread -o $@

alchemy2: Alchemy2.o
	$(CPP) $(STATIC) $(CCOPTS) $^  -L htslib/lib -lhts -lz -lpthread -o $@

qti: QueryTime.o
	$(CPP) $(STATIC) $(CCOPTS) $^  -L htslib/lib -lhts -lz -lpthread -o $@

lra.o: lra.cpp $(HEADERS) htslib/lib/libhts.a
	$(CPP) $(CCOPTS) -c  -I htslib/include -I seqan/include  lra.cpp

Alchemy2.o: Alchemy2.cpp  htslib/lib/libhts.a
	$(CPP) $(CCOPTS) -c  -I htslib/include -I seqan/include  Alchemy2.cpp

QueryTime.o: QueryTime.cpp $(HEADERS) htslib/lib/libhts.a
	$(CPP) $(CCOPTS) -c  -I htslib/include -I seqan/include  QueryTime.cpp


clean:
	rm -f lra lra.o


