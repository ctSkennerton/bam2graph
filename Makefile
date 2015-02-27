CXX = g++
CXXFLAGS = -Iseqan/include -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros
CXXFLAGS += -O3 -DNDEBUG -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1
LDFLAGS = -lz -lbz2 -lrt

PACKAGE_VERSION = 0.1
PACKAGE_DATE = "2014-12-30"

OBJS = main.o


bam2graph: version.h $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

ifneq "$(wildcard .git)" ""
PACKAGE_VERSION := $(shell git log -n 1 --pretty="%h" )
PACKAGE_DATE := $(shell git log -n 1 --pretty="%ai" )

# Force version.h to be remade if $(PACKAGE_VERSION) has changed.
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))
endif

version.h:
	echo '#define PACKAGE_VERSION "$(PACKAGE_VERSION)"' > $@
	echo '#define PACKAGE_DATE "$(PACKAGE_DATE)"' >> $@

force:

.PHONY: force
