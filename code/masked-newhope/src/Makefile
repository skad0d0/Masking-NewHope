CC = /usr/bin/gcc
CFLAGS = -w -Wall -Wextra -g -O3 -fomit-frame-pointer -march=native
NISTFLAGS = -O3 -fomit-frame-pointer -march=native -fPIC -no-pie

HEADERS = params.h poly.h randombytes.h ntt.h reduce.h verify.h cpapke.h gadgets.h utils.h masked_kem.h ccakem.h
SOURCES = poly.c reduce.c fips202.c verify.c cpapke.c ntt.c precomp.c masked_poly.c random.c gadgets.c newhope_gadgets.c utils.c masked_kem.c polycomp.c ccakem.c
# SOURCES_GADGETS = gadget.c newhope_gadget.c
# HEADERS_GADGETS = gadget.h
CCAKEM  = ccakem.c
# CPAKEM  = cpakem.c

# all:	test_ccakem512 \
# 			speed_ccakem512 \
# 				test_main \

all : test_main \
			

test_ccakem512: $(HEADERS) $(SOURCES) $(CCAKEM) test_newhope.c randombytes.c
	ln -sf ccakem.h api.h
	$(CC) $(CFLAGS) $(SOURCES) -DNEWHOPE_N=512 $(CCAKEM) randombytes.c test_newhope.c -o $@

speed_ccakem512: $(HEADERS) $(SOURCES) $(CCAKEM) speed.c randombytes.c cpucycles.h cpucycles.c
	ln -sf ccakem.h api.h
	$(CC) $(CFLAGS) $(SOURCES) -DNEWHOPE_N=512 $(CCAKEM) randombytes.c cpucycles.c speed.c -o $@

test_main: $(HEADERS) $(SOURCES) main.c randombytes.c
	$(CC) $(CFLAGS) $(SOURCES) -DNEWHOPE_N=512  randombytes.c main.c -o main_test

.PHONY: clean

clean:
	-rm PQCgenKAT_cpakem512
	-rm PQCgenKAT_cpakem1024
	-rm PQCgenKAT_ccakem512
	-rm PQCgenKAT_ccakem1024
	-rm test_cpakem512
	-rm test_cpakem1024
	-rm test_ccakem512
	-rm test_ccakem1024
	-rm speed_cpakem512
	-rm speed_cpakem1024
	-rm speed_ccakem512
	-rm speed_ccakem1024
	# -rm api.h
	-rm PQCkemKAT_*.req
	-rm PQCkemKAT_*.rsp
