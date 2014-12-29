CXX = g++
CXXFLAGS = -Iseqan/core/include -Iseqan/extras/include -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros
CXXFLAGS += -O3 -DNDEBUG -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1
LDFLAGS = -lz -lbz2 -lrt

OBJS = main.o

bam2graph: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)
