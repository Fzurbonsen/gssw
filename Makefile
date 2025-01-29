CC:=gcc
CFLAGS+=-Wall -O3 -g
OBJ_DIR:=obj
BIN_DIR:=bin
SRC_DIR:=src
LIB_DIR:=lib
INC_DIR:=include

CWD:=$(shell pwd)

INCLUDE_FLAGS :=-I$(CWD)/$(INC_DIR)
LD_LIB_DIR_FLAGS :=-L$(CWD)/$(LIB_DIR)

OBJ=gssw.o
EXE=gssw_example
EXEADJ=gssw_example_adj
EXETEST=gssw_test

EDLIB_DIR:=src/edlib
S_GWFA_DIR:=src/s_gwfa

LIB_DEPS=
LIB_DEPS+=$(LIB_DIR)/libedlib.a
LIB_DEPS+=$(LIB_DIR)/libs_gwfa.a

.PHONY:all clean cleanlocal test

all:$(BIN_DIR)/$(EXE) $(BIN_DIR)/$(EXEADJ) $(BIN_DIR)/$(EXETEST) $(LIB_DIR)/libgssw.a

$(BIN_DIR)/$(EXE):$(OBJ_DIR)/$(OBJ) $(SRC_DIR)/example.c $(LIB_DEPS)
	# Make dest directory
	@mkdir -p $(@D)
	$(CC) $(INCLUDE_FLAGS) $(LD_LIB_DIR_FLAGS) $(LDFLAGS) $(CPPFLAGS) $(CFLAGS) $(SRC_DIR)/example.c -o $@ $< -ls_gwfa -ledlib -lm -lz

$(BIN_DIR)/$(EXEADJ):$(OBJ_DIR)/$(OBJ) $(SRC_DIR)/example_adj.c $(LIB_DEPS)
	@mkdir -p $(@D)
	$(CC) $(INCLUDE_FLAGS) $(LD_LIB_DIR_FLAGS) $(LDFLAGS) $(CPPFLAGS) $(CFLAGS) $(SRC_DIR)/example_adj.c -o $@ $< -ls_gwfa -ledlib -lm -lz

$(BIN_DIR)/$(EXETEST):$(OBJ_DIR)/$(OBJ) $(SRC_DIR)/gssw_test.c $(LIB_DEPS)
	@mkdir -p $(@D)
	$(CC) $(INCLUDE_FLAGS) $(LD_LIB_DIR_FLAGS) $(LDFLAGS) $(CPPFLAGS) $(CFLAGS) $(SRC_DIR)/gssw_test.c -o $@ $< -ls_gwfa -ledlib -lm -lz

$(OBJ_DIR)/$(OBJ):$(SRC_DIR)/gssw.h $(SRC_DIR)/gssw.c $(LIB_DEPS)
	@mkdir -p $(@D)
	$(CC) $(INCLUDE_FLAGS) $(LD_LIB_DIR_FLAGS) $(LDFLAGS) $(CPPFLAGS) $(CFLAGS) -c -o $@ $(SRC_DIR)/gssw.c

$(LIB_DIR)/libgssw.a:$(OBJ_DIR)/$(OBJ) $(LIB_DEPS)
	@mkdir -p $(@D)
	ar rvs $@ $<

$(LIB_DIR)/libedlib.a: $(EDLIB_DIR)/edlib/src/edlib.cpp $(EDLIB_DIR)/edlib/include/edlib.h
	@mkdir -p $(LIB_DIR)
	@mkdir -p $(INC_DIR)
	+cd $(EDLIB_DIR) && cd build && cmake -D CMAKE_BUILD_TYPE=Release .. && $(MAKE) $(FILTER) && cp lib/libedlib.a $(CWD)/$(LIB_DIR)/ && cp ../edlib/include/edlib.h $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libs_gwfa.a: $(S_GWFA_DIR)/src/s_gwfa.c $(S_GWFA_DIR)/src/s_gwfa.h
	@mkdir -p $(LIB_DIR)
	@mkdir -p $(INC_DIR)
	+cd $(S_GWFA_DIR) && $(MAKE) clean && $(MAKE) $(FILTER) && cp lib/libs_gwfa.a $(CWD)/$(LIB_DIR)/ && cp src/s_gwfa.h $(CWD)/$(INC_DIR)/

	
test:$(BIN_DIR)/$(EXETEST)
	$(BIN_DIR)/$(EXETEST)

cleanlocal:
	$(RM) -r lib/
	$(RM) -r bin/
	$(RM) -r obj/
	cd $(DEP_DIR) && cd s_gwfa && $(MAKE) clean
	cd $(DEP_DIR) && cd edlib && $(MAKE) clean && rm -rf build && mkdir build

clean:cleanlocal




