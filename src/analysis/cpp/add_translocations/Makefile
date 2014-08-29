BOOST_ROOT ?= /home/craczy/boost_1_49_0
#SEQAN_ROOT ?= /home/oschulz-trieglaff/code/seqan-release/seqan-1.3.1/
SEQAN_ROOT ?= /home/oschulz-trieglaff/code/seqan-trunk/core/include/

#LDFLAGS := -L${BOOST_ROOT}/lib -lboost_program_options -lboost_filesystem 
LDFLAGS := -L${BOOST_ROOT}/lib

CPPFLAGS := -I${BOOST_ROOT}/include -I$(SEQAN_ROOT) -I.
CPPFLAGS += -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
CPPFLAGS += -O3 -g
#CPPFLAGS += -g
CPPFLAGS += -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros

CPP = /usr/bin/g++

TARGETS    = AddTranslocations
ADD_TR_OBJ = AddTranslocations.o VcfVariant.o

#ALL_OBJ = $(SAMPLE_INV_OBJ) 

all: ${TARGETS}

%.o: %.cpp %.h
	$(CPP) -c $(CPPFLAGS) $<

AddTranslocations: $(ADD_TR_OBJ)
	$(CPP) $(CPPFLAGS) $(LDFLAGS) -o $@ $^

clean: 
	/bin/rm -f *.o *~ ${TARGETS}

