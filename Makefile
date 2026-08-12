CC      ?= cc
CFLAGS  ?= -O3 -Wall -Wextra
CFLAGS  += $(shell pkg-config --cflags htslib)
LDLIBS   = $(shell pkg-config --libs htslib) -lm -lpthread

.PHONY: all clean

all: hapbsa

hapbsa: hapbsa.c
	$(CC) $(CFLAGS) -o $@ hapbsa.c $(LDLIBS)

clean:
	rm -f hapbsa
