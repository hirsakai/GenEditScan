UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	CC := gcc
	CFLAGS := -O3 -Wall
else
	CC := /opt/homebrew/opt/llvm/bin/clang
	CFLAGS := -O3 -Wno-deprecated-non-prototype
endif

CPROB := ./cprob
CMATH := ./cmath

INCS := $(CPROB)/mconf.h

COBJS := $(CPROB)/chdtr.o $(CPROB)/mtherr.o $(CPROB)/igam.o $(CPROB)/igami.o \
$(CPROB)/gamma.o $(CPROB)/const.o $(CPROB)/ndtri.o $(CPROB)/polevl.o \
$(CMATH)/isnan.o

MAIN := libprob.a

%.o: %.c
	$(CC) -c $(CFLAGS) -o $@ $<

# all
all: $(MAIN)

# cprob
$(MAIN): $(COBJS) $(INCS)
	ar rv $(MAIN) $(COBJS)
	ranlib $(MAIN)

# clean
clean:
	-rm -f $(COBJS) $(MAIN) *~

.PHONY: all clean
