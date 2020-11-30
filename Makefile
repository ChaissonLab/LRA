CC=g++
CFLAGS = -O2 -DNDEBUG 
INCLUDES = -I $(CONDA_PREFIX)/include  
LFLAGS = -L $(CONDA_PREFIX)/lib 
LIBS = -lhts -lz -lpthread -ldeflate -lbz2 

SRCS = lra.cpp
OBJS = $(SRCS:.c=.o)
MAIN = lra


.PHONY: depend clean

all:    $(MAIN)

$(MAIN): $(OBJS) 
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it

