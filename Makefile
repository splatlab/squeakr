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

CXXFLAGS += -Wall $(DEBUG) $(PROFILE) $(OPT) $(ARCH) -m64 -I. -I$(LOC_INCLUDE)

CFLAGS += -Wall $(DEBUG) $(PROFILE) $(OPT) $(ARCH) -m64 -I. -I$(LOC_INCLUDE)

LDFLAGS += $(DEBUG) $(PROFILE) $(OPT) -lpthread -lboost_system \
-lboost_thread -lm -lbz2 -lz -lrt

#
# declaration of dependencies
#

all: $(TARGETS)

# dependencies between programs and .o files
squeakr:					$(OBJDIR)/kmer.o $(OBJDIR)/hashutil.o $(OBJDIR)/util.o \
									$(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
									$(OBJDIR)/SqueakrFS.o \
									$(OBJDIR)/count.o $(OBJDIR)/query.o $(OBJDIR)/innerprod.o \
									$(OBJDIR)/list.o $(OBJDIR)/squeakr.o

# dependencies between .o files and .h files

$(OBJDIR)/squeakr.o:		$(LOC_SRC)/squeakr.cc
$(OBJDIR)/count.o: 			$(LOC_INCLUDE)/gqf_cpp.h $(LOC_INCLUDE)/chunk.h \
												$(LOC_INCLUDE)/kmer.h \
												$(LOC_INCLUDE)/reader.h $(LOC_INCLUDE)/util.h \
												$(LOC_INCLUDE)/SqueakrFS.h
$(OBJDIR)/query.o: 			$(LOC_INCLUDE)/gqf_cpp.h $(LOC_INCLUDE)/kmer.h \
												$(LOC_INCLUDE)/util.h
$(OBJDIR)/innerprod.o: 	$(LOC_INCLUDE)/gqf_cpp.h
$(OBJDIR)/list.o: 		 	$(LOC_INCLUDE)/gqf_cpp.h $(LOC_INCLUDE)/kmer.h \
												$(LOC_INCLUDE)/util.h
$(OBJDIR)/kmer.o: 			$(LOC_SRC)/kmer.cc $(LOC_INCLUDE)/kmer.h
$(OBJDIR)/util.o: 			$(LOC_SRC)/util.cc $(LOC_INCLUDE)/util.h

# dependencies between .o files and .cc (or .c) files
$(OBJDIR)/gqf.o: 				$(LOC_SRC)/gqf/gqf.c $(LOC_INCLUDE)/gqf/gqf.h
$(OBJDIR)/gqf_file.o: 	$(LOC_SRC)/gqf/gqf_file.c $(LOC_INCLUDE)/gqf/gqf_file.h
$(OBJDIR)/hashutil.o: 	$(LOC_INCLUDE)/gqf/hashutil.h

#
# generic build rules
#

$(TARGETS):
	$(LD) $^ $(LDFLAGS) -o $@

$(OBJDIR)/%.o: $(LOC_SRC)/%.cc | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

$(OBJDIR)/%.o: $(LOC_SRC)/%.c | $(OBJDIR)
	$(CXX) $(CFLAGS) $(INCLUDE) -c -o $@ $<

$(OBJDIR)/%.o: $(LOC_SRC)/gqf/%.c | $(OBJDIR)
	$(CXX) $(CFLAGS) $(INCLUDE) -c -o $@ $<

$(OBJDIR):
	@mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR) core $(TARGETS)

