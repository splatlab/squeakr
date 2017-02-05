CXX = g++ -std=c++11

#CXXFLAGS = -Wall -g -I. -pthread -Wno-unused-result -Wno-strict-aliasing
#CXXFLAGS = -Wall -Ofast -m64 -I. -Wno-unused-result -Wno-strict-aliasing -DLOG_WAIT_TIME -DLOG_CLUSTER_LENGTH
CXXFLAGS = -Wall -Ofast -m64 -I. -Wno-unused-result -Wno-strict-aliasing

LDFLAGS = -lpthread -lssl -lcrypto -lboost_system -lboost_thread
LIBS = libs/libbz2.a libs/libz.a

TARGET_MAIN	= main
MAIN_SRC = main.cc hashutil.cc threadsafe-gqf/gqf.c

TARGET_QUERY	= query
QUERY_SRC = kmer_query.cc hashutil.cc threadsafe-gqf/gqf.c

TARGET_INNER_PROD = inner-prod
INNER_PROD_SRC = kmer_inner_prod.cc hashutil.cc threadsafe-gqf/gqf.c

$(TARGET_MAIN): $(MAIN_SRC)
	$(CXX) $(CXXFLAGS) $(MAIN_SRC) $(INCLUDE) $(LDFLAGS) $(LIBS) -o $@

$(TARGET_QUERY): $(QUERY_SRC)
	$(CXX) $(CXXFLAGS) $(QUERY_SRC) $(INCLUDE) $(LDFLAGS) $(LIBS) -o $@

$(TARGET_INNER_PROD): $(INNER_PROD_SRC)
	$(CXX) $(CXXFLAGS) $(INNER_PROD_SRC) $(INCLUDE) $(LDFLAGS) $(LIBS) -o $@

clean:
	rm -f $(TARGET_MAIN) $(TARGET_QUERY) $(TARGET_INNER_PROD) *.o core
