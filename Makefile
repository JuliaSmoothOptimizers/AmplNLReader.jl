ASL_LIB_DIR = $(shell brew --prefix asl)/lib
ASL_INC_DIR = $(shell brew --prefix asl)/include

LIBDIRS = -L$(ASL_LIB_DIR)
INCDIRS = -I. -I$(ASL_INC_DIR)

DEBUG = #-DDEBUG_AMPL_JL

HDR = jampl.h
SRC = jampl.c
LIB = jampl.dylib

all: $(LIB)

$(LIB): $(HDR) $(SRC)
	clang -g -std=c99 -dynamiclib $(DEBUG) -fPIC $(SRC) $(INCDIRS) $(LIBDIRS) -lasl -o $@

clean:
	rm -f ${SRC:.c=.o}

mrclean: clean
	rm -f $(LIB)
