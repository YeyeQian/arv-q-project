TAR=test_sign
OBJ=fips202.o ntt.o ntt_custom.o packing.o packing_custom.o poly.o polyvec.o polyvec_custom.o reduce.o rounding.o sign.o sign_custom.o symmetric-shake.o symmetric-shake-custom.o poly_custom.o fips202_custom.o
OBJ+=sign_custom_asm.o
OBJ+=poly_custom_rorg_asm.o
OBJ+=sign_custom_optimal_c908.o
# OBJ+=test_symmetric_shake.o
OBJ+=test_sign.o
# OBJ+=test_ntt.o

GCC:=/home/data/userhome/zhaoyifan/MyToolChains/xuantie-riscv-gnu-c908e/bin/riscv64-unknown-elf-gcc
OBJDUMP:=/home/data/userhome/zhaoyifan/MyToolChains/xuantie-riscv-gnu-c908e/bin/riscv64-unknown-elf-objdump
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