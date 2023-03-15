PROG=   lra
PROG_EXTRA= alchemy2 qti
LIBS=   -lz -lpthread -lhts
PROF=/home/cmb-16/mjc/shared/lib/
DEBUG?=""
OPT?=""
STATIC=
CXX=g++ -std=c++14 
CFLAGS=-g
asan?=""
tsan?=""

ifneq ($(DEBUG), "")
CFLAGS=-g
else
CFLAGS=-O2  -DNDEBUG 
endif

ifneq ($(asan), "")
  CFLAGS+=-fsanitize=address
  LIBS+=-fsanitize=address
endif

ifneq ($(tsan), "")
  CFLAGS+=-fsanitize=thread
  LIBS+=-fsanitize=thread
endif

ifneq ($(OPT), "")
#STATIC=-L $(PROF) -lprofiler
endif

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
  SplitClusters.h \
  RefineBreakpoint.h

all:$(PROG) 

# tag: TestAffineOneGapAlign.cpp AffineOneGapAlign.h
# 	$(CXX) -g TestAffineOneGapAlign.cpp -o tag 
# # -D _MAT_PRINT_

# tgc: TestGlobalChain.cpp GlobalChain.h Fragment.h BasicEndpoint.h PrioritySearchTree.h
# 	$(CXX) -g TestGlobalChain.cpp -o tgc

# tir: TestIndelRefine.cpp IndelRefine.h
# 	$(CXX) -g TestIndelRefine.cpp  -I $(CONDA_PREFIX)/include -L $(CONDA_PREFIX)/lib  -lhts -o tir -lbz2 -lz

lra: lra.o
	$(CXX) $(CFLAGS) $(STATIC) $^ -I  -L/usr/lib64  -L $(CONDA_PREFIX)/lib $(LIBS) -o $@ -Wl,-rpath-link=$(CONDA_PREFIX)/lib 

alchemy2: Alchemy2.o
	$(CXX) $(CFLAGS) $(STATIC)  $^  -L $(CONDA_PREFIX)/lib  $(LIBS) -o $@  -Wl,-rpath,$(CONDA_PREFIX)/lib

qti: QueryTime.o
	$(CXX) $(CFLAGS) $(STATIC) $^  -L $(CONDA_PREFIX)/lib $(LIBS) -o $@

lra.o: lra.cpp $(HEADERS) 
	$(CXX) $(CFLAGS) -c  -I $(CONDA_PREFIX)/include  lra.cpp 

Alchemy2.o: Alchemy2.cpp Genome.h
	$(CXX) $(CFLAGS) -c  -I $(CONDA_PREFIX)/include Alchemy2.cpp

QueryTime.o: QueryTime.cpp $(HEADERS) $(CONDA_PREFIX)/lib/libhts.a
	$(CXX) $(CFLAGS) -c  -I $(CONDA_PREFIX)/include QueryTime.cpp
clean:
	rm -f $(PROG) $(PROG_EXTRA) *.o 
