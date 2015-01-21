ASL_LIB_DIR = $(shell brew --prefix asl)/lib
ASL_INC_DIR = $(shell brew --prefix asl)/include

LIBDIRS = -L$(ASL_LIB_DIR)
INCDIRS = -I. -I$(ASL_INC_DIR)

ifeq ($(strip $(shell uname)),Darwin)
    SO = dylib
    SHARED = -dynamiclib
else
    SO = so
    SHARED = -shared -Wl,-soname,$@
endif

HDR = jampl.h
SRC = jampl.c
LIB = libjampl.$(SO)

all: $(LIB)

$(LIB): $(HDR) $(SRC)
	cc -g -std=c99 $(SHARED) -fPIC $(SRC) $(INCDIRS) $(LIBDIRS) -lasl -lm -ldl -o $@

clean:
	rm -f ${SRC:.c=.o}

mrclean: clean
	rm -f $(LIB)
