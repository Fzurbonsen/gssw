CC:=gcc
CFLAGS+= -Wall -O3 -g
CPPFLAGS+= -std=c++11
CWD=$(shell pwd)
OBJ_DIR:=$(CWD)/obj
BIN_DIR:=$(CWD)/bin
SRC_DIR:=$(CWD)/src
LIB_DIR:=$(CWD)/lib
INC_DIR:=$(CWD)/include

OBJ=$(OBJ_DIR)/gssw.o
OBJ+=$(OBJ_DIR)/vg_gwfa_pipeline.o
EXE=gssw_example
EXEADJ=gssw_example_adj
EXETEST=gssw_test
EXE_PIPELINE_TEST=vg_gwfa_pipeline_test

LIB_FLAGS= -lz -lm -lstdc++ -ledlib -lcsswl -lgwfa
LDFLAGS = -L$(LIB_DIR)

INCFLAGS= -I$(INC_DIR)

.PHONY:all clean cleanlocal test

LIB=$(LIB_DIR)/libedlib.a
LIB+=$(LIB_DIR)/libcsswl.a
LIB+=$(LIB_DIR)/libgwfa.a

INCLUDE=$(SRC_DIR)/vg_gwfa_pipeline.hpp
INCLUDE+=$(SRC_DIR)/vg_gwfa_pipeline_wrapper.h
INCLUDE+=$(SRC_DIR)/gwfa/*.h
INCLUDE+=$(SRC_DIR)/edlib/edlib/include/edlib.h

all:$(BIN_DIR)/$(EXE) $(BIN_DIR)/$(EXEADJ) $(BIN_DIR)/$(EXETEST) $(BIN_DIR)/$(EXE_PIPELINE_TEST) $(LIB_DIR)/libgssw.a

$(BIN_DIR)/$(EXE):$(OBJ) $(SRC_DIR)/example.c
	# Make dest directory
	@mkdir -p $(@D)
	$(CC) $(LDFLAGS) $(CPPFLAGS) $(CFLAGS) $(SRC_DIR)/example.c -o $@ $(OBJ) $(LIB_FLAGS)

$(BIN_DIR)/$(EXEADJ):$(OBJ) $(SRC_DIR)/example_adj.c
	@mkdir -p $(@D)
	$(CC) $(LDFLAGS) $(CPPFLAGS) $(CFLAGS) $(SRC_DIR)/example_adj.c -o $@ $(OBJ) $(LIB_FLAGS)

$(BIN_DIR)/$(EXETEST):$(OBJ) $(SRC_DIR)/gssw_test.c
	@mkdir -p $(@D)
	$(CC) $(LDFLAGS) $(CPPFLAGS) $(CFLAGS) $(SRC_DIR)/gssw_test.c -o $@ $(OBJ) $(LIB_FLAGS)

$(BIN_DIR)/$(EXE_PIPELINE_TEST): $(OBJ) $(SRC_DIR)/vg_gwfa_pipeline_test.c
	@mkdir -p $(@D)
	$(CC) $(LDFLAGS) $(CPPFLAGS) $(CFLAGS) $(SRC_DIR)/vg_gwfa_pipeline_test.c -o $@ $(OBJ) $(LIB_FLAGS)

# $(OBJ):$(INC_DIR) $(LIB) $(SRC_DIR)/gssw.c $(SRC_DIR)/vg_gwfa_pipeline.cpp
# 	@mkdir -p $(@D)
# 	$(CC) $(LDFLAGS) $(CPPFLAGS) $(CFLAGS) -c -o $@ $(SRC_DIR)/gssw.c $(INCFLAGS)

$(OBJ_DIR)/gssw.o: $(LIB) $(SRC_DIR)/gssw.c $(INC_DIR)
	@mkdir -p $(@D)
	$(CC) $(LDFLAGS) $(CPPFLAGS) $(CFLAGS) -c $(SRC_DIR)/gssw.c -o $@ $(INCFLAGS)

$(OBJ_DIR)/vg_gwfa_pipeline.o: $(LIB) $(SRC_DIR)/vg_gwfa_pipeline.cpp $(INC_DIR)
	@mkdir -p $(@D)
	$(CC) $(LDFLAGS) $(CPPFLAGS) $(CFLAGS) -c $(SRC_DIR)/vg_gwfa_pipeline.cpp -o $@ $(INCFLAGS)

$(LIB_DIR)/libedlib.a:
	@mkdir -p $(LIB_DIR)/
	+ cd $(SRC_DIR)/edlib && cd build && cmake -D CMAKE_BUILD_TYPE=Release .. && $(MAKE) && cp lib/libedlib.a $(LIB_DIR)

$(LIB_DIR)/libcsswl.a:
	@mkdir -p $(LIB_DIR)/
	+ cd $(SRC_DIR)/csswl/src && $(MAKE) && ar rcs libcsswl.a *.o && cp libcsswl.a $(LIB_DIR)

$(LIB_DIR)/libgwfa.a:
	@mkdir -p $(LIB_DIR)/
	+ cd $(SRC_DIR)/gwfa && $(MAKE) && ar rcs libgwfa.a gfa-base.o gfa-io.o gfa-sub.o gwf-ed.o kalloc.o && cp libgwfa.a $(LIB_DIR)

$(INC_DIR): $(INCLUDE)
	@mkdir -p $(INC_DIR)/
	@mkdir -p $(INC_DIR)/gwfa
	@mkdir -p $(INC_DIR)/edlib
	@mkdir -p $(INC_DIR)/csswl
	+ cp -r $(SRC_DIR)/vg_gwfa_pipeline_wrapper.h $(INC_DIR)
	+ cp -r $(SRC_DIR)/vg_gwfa_pipeline.hpp $(INC_DIR)
	+ cp -r $(SRC_DIR)/gssw.h $(INC_DIR)
	+ cp -r $(SRC_DIR)/gwfa/*.h $(INC_DIR)/gwfa
	+ cp -r $(SRC_DIR)/edlib/edlib/include/edlib.h $(INC_DIR)/edlib
	+ cp -r $(SRC_DIR)/csswl/src/*.h $(INC_DIR)/csswl

$(LIB_DIR)/libgssw.a:$(OBJ)
	@mkdir -p $(@D)
	ar rvs $@ $^
	
test:$(BIN_DIR)/$(EXETEST)
	$(BIN_DIR)/$(EXETEST)

cleanlocal:
	$(RM) -r lib/
	$(RM) -r bin/
	$(RM) -r obj/
	$(RM) -r include/
	cd $(SRC_DIR)/gwfa && $(MAKE) clean && rm -f libgwfa.a
	cd $(SRC_DIR)/edlib && $(MAKE) clean
	cd $(SRC_DIR)/csswl/src $(MAKE) clean && rm -f libcsswl.a

clean:cleanlocal




