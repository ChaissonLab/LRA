all:  lra 
PROF=/home/cmb-16/mjc/shared/lib/
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
CCOPTS=-O -g $(CCOPTS_BASE) -lprofiler
#STATIC=
STATIC=-g -L $(PROF) -lprofiler
#CCOPTS=-g -fsanitize=undefined
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
  Read.h \
  SparseDP.h \
  Timing.h \
  IndelRefine.h \
  Mapping_ultility.h \
  SparseDP_Forward.h \
  DivideSubByCol1.h \
  DivideSubByCol2.h \
  DivideSubByRow1.h \
  DivideSubByRow2.h \
  AlignmentBlock.h \
  ChainRefine.h \
  Chain.h \
  ClusterRefine.h \
  Fragment_Info.h \
  LinearExtend.h \
  LocalRefineAlignment.h \
  LogLookUpTable.h \
  Map_highacc.h \
  Map_lowacc.h \
  SubRountine.h \
  Types.h \
  SubProblem.h  \
  SplitClusters.h 

CXX=g++ -std=c++11 

# tag: TestAffineOneGapAlign.cpp AffineOneGapAlign.h
# 	$(CXX) -g TestAffineOneGapAlign.cpp -o tag 
# # -D _MAT_PRINT_

# tgc: TestGlobalChain.cpp GlobalChain.h Fragment.h BasicEndpoint.h PrioritySearchTree.h
# 	$(CXX) -g TestGlobalChain.cpp -o tgc

# tir: TestIndelRefine.cpp IndelRefine.h
# 	$(CXX) -g TestIndelRefine.cpp  -I $(CONDA_PREFIX)/include -L $(CONDA_PREFIX)/lib  -lhts -o tir -lbz2 -lz

lra: lra.o
	$(CXX) $(STATIC) $(CCOPTS) $^ -I $(CONDA_PREFIX)/include -L $(CONDA_PREFIX)/lib  -lhts -lz -lpthread -ldeflate -lbz2  -o $@ 

# alchemy2: Alchemy2.o
# 	$(CXX) $(STATIC) $(CCOPTS) $^  -L $(CONDA_PREFIX)/lib  -lhts -lz -lpthread -o $@  -Wl,-rpath,$(PWD)/htslib/lib

# qti: QueryTime.o
# 	$(CXX) $(STATIC) $(CCOPTS) $^  -L $(CONDA_PREFIX)/lib -lhts -lz -lpthread -o $@

lra.o: lra.cpp $(HEADERS) 
	$(CXX) $(CCOPTS) -c  -I $(CONDA_PREFIX)/include  lra.cpp 
#  $(CXX) $(CCOPTS) -c  -I htslib/include  lra.cpp 

# Alchemy2.o: Alchemy2.cpp Genome.h htslib/lib/libhts.a
Alchemy2.o: Alchemy2.cpp Genome.h
	$(CXX) $(CCOPTS) -c  -I $(CONDA_PREFIX)/include Alchemy2.cpp

# QueryTime.o: QueryTime.cpp $(HEADERS) htslib/lib/libhts.a
QueryTime.o: QueryTime.cpp $(HEADERS) htslib/lib/libhts.a
	$(CXX) $(CCOPTS) -c  -I $(CONDA_PREFIX)/include QueryTime.cpp

writeblock: WriteBlock.cpp
	$(CXX) $(CCOPTS) WriteBlock.cpp -o writeblock

IndexInformativeKmers.o: IndexInformativeKmers.cpp 
	$(CXX) $(CCOPTS) $^ -c -I $(CONDA_PREFIX)/include/htslib -I bwa

# bwa/bwa.o:
# 	cd bwa && make

# iik: IndexInformativeKmers.o bwa/bwa.o bwa/kstring.o bwa/utils.o bwa/kthread.o bwa/kstring.o bwa/ksw.o bwa/bwt.o bwa/bntseq.o bwa/bwamem.o bwa/bwamem_pair.o bwa/bwamem_extra.o bwa/malloc_wrap.o	bwa/QSufSort.o bwa/bwt_gen.o bwa/rope.o bwa/rle.o bwa/is.o bwa/bwtindex.o
# 	$(CXX) $(CCOPTS) $^ -o iik -L $(CONDA_PREFIX)/lib/libhts.a -lz -lpthread -I $(CONDA_PREFIX)/include

clean:
	rm -f lra lra.o 
#   rm -f lra lra.o iik IndexInformativeKmers.o alchemy2 Alchemy2.o
