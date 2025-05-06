TAR=test_kem_hqc192_vlen1024
OBJ=code.o code_custom.o fft_custom.o fips202.o fips202_custom.o gf.o gf_custom.o hqc.o hqc_custom.o kem.o kem_custom.o parsing_custom.o reed_muller_custom.o reed_solomon_custom.o shake_ds.o shake_prng.o vector.o vector_custom.o
# OBJ+=test_gf.o
# OBJ+=test_vector.o
# OBJ+=test_reed_solomon.o
# OBJ+=test_reed_muller.o
OBJ+=test_kem.o

GCC:=/your-xuantie-gnu-path/bin/riscv64-unknown-elf-gcc
OBJDUMP:=/your-xuantie-gnu-path/bin/riscv64-unknown-elf-objdump
LINKSCRIPT:=riscv64-virt.ld

TARGET_ISA=riscv
GEM5_PATH=/your-gem5-path

CFLAGS=-O3 -march=rv64imafdcv -fno-tree-vectorize -Wall -I$(GEM5_PATH)/include
LDFLAGS=-L$(GEM5_PATH)/util/m5/build/$(TARGET_ISA)/out -lm5 -static

.DEFAULT_GOAL := default

.PHONY: default
default: $(addsuffix .riscv,$(TAR))

.PHONY: dumps
dumps: $(addsuffix .dump,$(TAR))

%.o:%.c
	$(GCC) $(CFLAGS) -c $^ -o $@

%.riscv:$(OBJ)
	$(GCC) $(LDFLAGS) $^ -o $@

%.dump: %.riscv
	$(OBJDUMP) -D $< > $@

.PHONY:
clean:
	rm -f *.riscv *.o *.dump
