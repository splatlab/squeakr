TARGETS= squeakr-count squeakr-query squeakr-inner-prod

ifdef D
	DEBUG=-g
	OPT=
else
	DEBUG=
	OPT=-Ofast
endif

ifdef NH
	ARCH=
else
	ARCH=-msse4.2 -D__SSE4_2_
endif

ifdef P
	PROFILE=-pg -no-pie # for bug in gprof.
endif

CXX = g++ -std=c++11
CC = g++ -std=c++11
LD= g++ -std=c++11

CXXFLAGS = -Wall $(DEBUG) $(PROFILE) $(OPT) $(ARCH) -m64 -I. -Wno-unused-result -Wno-strict-aliasing -Wno-unused-function -Wno-sign-compare

LDFLAGS = $(DEBUG) $(PROFILE) $(OPT) -lpthread -lssl -lcrypto -lboost_system -lboost_thread -lm -lbz2 -lz

#
# declaration of dependencies
#

all: $(TARGETS)

# dependencies between programs and .o files

squeakr-count:                  main.o 								 hashutil.o threadsafe-gqf/gqf.o
squeakr-query: 					 kmer_query.o 					 hashutil.o threadsafe-gqf/gqf.o
squeakr-inner-prod: 			 kmer_inner_prod.o 			 hashutil.o threadsafe-gqf/gqf.o

# dependencies between .o files and .h files

main.o: 								 									threadsafe-gqf/gqf.h hashutil.h chunk.h kmer.h reader.h
kmer_query.o: 					 									threadsafe-gqf/gqf.h hashutil.h chunk.h kmer.h
kmer_inner_prod.o: 			 									threadsafe-gqf/gqf.h hashutil.h
hashutil.o: 																									 hashutil.h

# dependencies between .o files and .cc (or .c) files

%.o: %.cc
threadsafe-gqf/gqf.o: threadsafe-gqf/gqf.c threadsafe-gqf/gqf.h

#
# generic build rules
#

$(TARGETS):
	$(LD) $^ $(LDFLAGS) -o $@

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDE) $< -c -o $@

%.o: %.c
	$(CC) $(CXXFLAGS) $(INCLUDE) $< -c -o $@

clean:
	rm -f *.o threadsafe-gqf/gqf.o $(TARGETS)

