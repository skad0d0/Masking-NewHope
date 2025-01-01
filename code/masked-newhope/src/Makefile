CC = /usr/bin/gcc
CFLAGS = -w -Wall -Wextra -g -O3 -fomit-frame-pointer -march=native
NISTFLAGS = -O3 -fomit-frame-pointer -march=native -fPIC -no-pie

HEADERS = params.h poly.h randombytes.h ntt.h reduce.h verify.h cpapke.h gadgets.h utils.h masked_kem.h ccakem.h cpucycles.h
SOURCES = poly.c reduce.c fips202.c verify.c cpapke.c ntt.c precomp.c masked_poly.c random.c gadgets.c newhope_gadgets.c utils.c masked_kem.c polycomp.c ccakem.c cpucycles.c

ORDER=3
RNG=1
MACRO = -D NEWHOPE_MASKING_ORDER=$(ORDER) -D RNG_MODE=$(RNG)

all : test_main
			
test_main: $(HEADERS) $(SOURCES) main.c randombytes.c
	$(CC) $(CFLAGS) $(SOURCES) -DNEWHOPE_N=512 $(MACRO) randombytes.c main.c -o main_test

.PHONY: clean

clean:
	-rm main_test
