all:  lra 
PROF=/home/cmb-16/mjc/shared/lib/
#PROF=/home/cmb-16/mjc/jingwenr/software/lib
CCOPTS_BASE=
DEBUG?=""
OPT?=""
ifneq ($(DEBUG), "")
CCOPTS=$(CCOPTS_BASE) $(DEBUG)
else
CCOPTS=-O2 $(CCOPTS_BASE)  -DNDEBUG 
endif
STATIC=
ifneq ($(OPT), "")
#CCOPTS=-g $(CCOPTS_BASE) -lprofiler
#STATIC=
STATIC=-g -L $(PROF) -lprofiler
CCOPTS=-g
endif

#-D _TESTING_ -lprofiler 
#  -L$(PROF) 

HEADERS=MinCount.h \
  SeqUtils.h \
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
  Read.h \
  SplitClusters.h \
  SparseDP.h \
  Timing.h \

CXX=g++
# -std=c++14 

tag: TestAffineOneGapAlign.cpp AffineOneGapAlign.h
	$(CXX) -g TestAffineOneGapAlign.cpp -o tag 
# -D _MAT_PRINT_

tgc: TestGlobalChain.cpp GlobalChain.h Fragment.h BasicEndpoint.h PrioritySearchTree.h
	$(CXX) -g TestGlobalChain.cpp -o tgc

lra: lra.o
	$(CXX) $(STATIC) $(CCOPTS) $^ -I $(CONDA_PREFIX)/include -L $(CONDA_PREFIX)/lib  -lhts -lz -lpthread -ldeflate -lbz2 -lhts -o $@ 

alchemy2: Alchemy2.o
	$(CXX) $(STATIC) $(CCOPTS) $^  -L htslib/lib  -lhts -lz -lpthread -o $@  -Wl,-rpath,$(PWD)/htslib/lib

qti: QueryTime.o
	$(CXX) $(STATIC) $(CCOPTS) $^  -L htslib/lib -lhts -lz -lpthread -o $@

lra.o: lra.cpp $(HEADERS) 
	$(CXX) $(CCOPTS) -c  -I $(CONDA_PREFIX)/include/  lra.cpp 
#  $(CXX) $(CCOPTS) -c  -I htslib/include  lra.cpp 

Alchemy2.o: Alchemy2.cpp Genome.h htslib/lib/libhts.a
	$(CXX) $(CCOPTS) -c  -I htslib/include Alchemy2.cpp

QueryTime.o: QueryTime.cpp $(HEADERS) htslib/lib/libhts.a
	$(CXX) $(CCOPTS) -c  -I htslib/include QueryTime.cpp

writeblock: WriteBlock.cpp
	$(CXX) $(CCOPTS) WriteBlock.cpp -o writeblock

IndexInformativeKmers.o: IndexInformativeKmers.cpp 
	$(CXX) $(CCOPTS) $^ -c -I htslib -I bwa

bwa/bwa.o:
	cd bwa && make

iik: IndexInformativeKmers.o bwa/bwa.o bwa/kstring.o bwa/utils.o bwa/kthread.o bwa/kstring.o bwa/ksw.o bwa/bwt.o bwa/bntseq.o bwa/bwamem.o bwa/bwamem_pair.o bwa/bwamem_extra.o bwa/malloc_wrap.o	bwa/QSufSort.o bwa/bwt_gen.o bwa/rope.o bwa/rle.o bwa/is.o bwa/bwtindex.o
	$(CXX) $(CCOPTS) $^ -o iik -L htslib/lib/libhts.a -lz -lpthread -I htslib

clean:
	rm -f lra lra.o iik IndexInformativeKmers.o
