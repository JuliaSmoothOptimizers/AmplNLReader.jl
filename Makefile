ASL_LIB_DIR = $(shell brew --prefix asl)/lib
ASL_INC_DIR = $(shell brew --prefix asl)/include
JULIA_LIB_DIR = $(shell brew --prefix julia)/lib/julia
JULIA_INC_DIR = $(shell brew --prefix julia)/include

LIBDIRS = -L$(ASL_LIB_DIR) -L$(JULIA_LIB_DIR)
INCDIRS = -I. -I$(ASL_INC_DIR) -I$(JULIA_INC_DIR)

HDR = jampl.h
SRC = jampl.c
LIB = jampl.dylib

all: $(LIB)

$(LIB): $(HDR) $(SRC)
	clang -std=c99 -dynamiclib -fPIC $(SRC) $(INCDIRS) $(LIBDIRS) -lasl -ljulia -o $@

clean:
	rm -f ${SRC:.c=.o}

mrclean: clean
	rm -f $(LIB)
