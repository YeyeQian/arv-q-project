TAR=test_kem
OBJ=cbd.o cbd_custom.o fips202.o fips202_custom.o kem.o kem_custom.o ntt.o ntt_custom.o pack_unpack.o pack_unpack_custom.o poly_mul.o poly_mul_custom.o poly.o poly_custom.o reduce.o SABER_indcpa.o SABER_indcpa_custom.o verify.o verify_custom.o
OBJ+=test_kem.o

GCC:=/home/data/userhome/zhaoyifan/MyToolChains/xuantie-riscv-gnu-c908e/bin/riscv64-unknown-elf-gcc
OBJDUMP:=/home/data/userhome/zhaoyifan/MyToolChains/xuantie-riscv-gnu-c908e/bin/riscv64-unknown-elf-objdump
LINKSCRIPT:=riscv64-virt.ld

TARGET_ISA=riscv
GEM5_PATH=/your-gem5-path

CFLAGS=-O3 -march=rv64imafdcv -fno-tree-vectorize -Wall -I$(GEM5_PATH)/include
LDFLAGS=-L$(GEM5_PATH)/util/m5/build/$(TARGET_ISA)/out -lm5 -static -lm

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