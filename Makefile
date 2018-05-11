TARGETS= squeakr

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
CC = gcc -std=gnu11
LD= g++ -std=c++11

LOC_INCLUDE=include
LOC_SRC=src
OBJDIR=obj

CXXFLAGS += -Wall $(DEBUG) $(PROFILE) $(OPT) $(ARCH) -m64 -I. -I$(LOC_INCLUDE) \
-Wno-unused-result -Wno-strict-aliasing -Wno-unused-function -Wno-sign-compare

CFLAGS += -Wall $(DEBUG) $(PROFILE) $(OPT) $(ARCH) -m64 -I. -I$(LOC_INCLUDE)\
-Wno-unused-result -Wno-strict-aliasing -Wno-unused-function -Wno-sign-compare \
-Wno-implicit-function-declaration

LDFLAGS += $(DEBUG) $(PROFILE) $(OPT) -lpthread -lboost_system \
-lboost_thread -lm -lz -lrt

#
# declaration of dependencies
#

all: $(TARGETS)

# dependencies between programs and .o files
squeakr:					$(OBJDIR)/count.o $(OBJDIR)/query.o $(OBJDIR)/innerprod.o $(OBJDIR)/list.o $(OBJDIR)/hashutil.o $(OBJDIR)/kmer.o $(OBJDIR)/util.o

# dependencies between .o files and .h files

$(OBJDIR)/count.o: 			$(LOCAL_INCLUDE)/cqf.h $(LOCAL_INCLUDE)/hashutil.h $(LOCAL_INCLUDE)/chunk.h $(LOCAL_INCLUDE)/kmer.h $(LOCAL_INCLUDE)/reader.h
$(OBJDIR)/query.o: 			$(LOCAL_INCLUDE)/cqf.h $(LOCAL_INCLUDE)/hashutil.h $(LOCAL_INCLUDE)/chunk.h $(LOCAL_INCLUDE)/kmer.h
$(OBJDIR)/innerprod.o: $(LOCAL_INCLUDE)/cqf.h $(LOCAL_INCLUDE)/hashutil.h
$(OBJDIR)/list.o: 		 $(LOCAL_INCLUDE)/cqf.h $(LOCAL_INCLUDE)/hashutil.h
$(OBJDIR)/hashutil.o: 	$(LOCAL_INCLUDE)/hashutil.h
$(OBJDIR)/kmer.o: 			$(LOC_SRC)/kmer.cc $(LOC_INCLUDE)/kmer.h
$(OBJDIR)/util.o: 			$(LOC_SRC)/utill.cc $(LOC_INCLUDE)/util.h

# dependencies between .o files and .cc (or .c) files
$(OBJDIR)/gqf.o: $(LOC_SRC)/cqf/gqf.c $(LOC_INCLUDE)/cqf/gqf.h

#
# generic build rules
#

$(TARGETS):
	$(LD) $^ $(LDFLAGS) -o $@

$(OBJDIR)/%.o: $(LOC_SRC)/%.cc | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

$(OBJDIR)/%.o: $(LOC_SRC)/%.c | $(OBJDIR)
	$(CC) $(CFLAGS) $(INCLUDE) -c -o $@ $<

$(OBJDIR)/%.o: $(LOC_SRC)/cqf/%.c | $(OBJDIR)
	$(CC) $(CFLAGS) $(INCLUDE) -c -o $@ $<

$(OBJDIR):
	@mkdir -p $(OBJDIR)

clean:
	rm -f $(OBJDIR)/*.o core $(TARGETS)

