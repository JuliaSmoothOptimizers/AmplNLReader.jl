ASL_LIB_DIR = $(shell brew --prefix asl)/lib
ASL_INC_DIR = $(shell brew --prefix asl)/include
JULIA_LIB_DIR = $(shell brew --prefix julia)/lib/julia
JULIA_INC_DIR = $(shell brew --prefix julia)/include

LIBDIRS = -L$(ASL_LIB_DIR) -L$(JULIA_LIB_DIR)
INCDIRS = -I$(ASL_INC_DIR) -I$(JULIA_INC_DIR)

SRC = jampl.c
LIB = jampl.dylib

all: $(LIB)

$(LIB): $(SRC)
	clang -std=c99 -dynamiclib -fPIC $< $(INCDIRS) $(LIBDIRS) -lasl -ljulia -o $@

clean:
	rm -f *.o

mrclean: clean
	rm -f $(LIB)
