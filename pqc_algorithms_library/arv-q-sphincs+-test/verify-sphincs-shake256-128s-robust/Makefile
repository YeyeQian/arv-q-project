TAR=test_spx_custom

HASH=shake256
THASH=robust
OBJ=address.o wots.o utils.o fors.o sign.o hash_$(HASH).o thash_$(HASH)_$(THASH).o fips202.o fips202_custom.o fors_custom.o hash_$(HASH)_custom.o sign_custom.o thash_$(HASH)_$(THASH)_custom.o wots_custom.o utils_custom.o
#OBJ+=test_fors.o
OBJ+=test_spx_custom.o

GCC:=/home/data/userhome/zhaoyifan/MyToolChains/xuantie-riscv-gnu-c908e/bin/riscv64-unknown-elf-gcc
OBJDUMP:=/home/data/userhome/zhaoyifan/MyToolChains/xuantie-riscv-gnu-c908e/bin/riscv64-unknown-elf-objdump
# GCC:=/home/data/userhome/zhaoyifan/MyToolChains/my-riscv-vector-gnu/bin/riscv64-unknown-elf-gcc
# OBJDUMP:=/home/data/userhome/zhaoyifan/MyToolChains/my-riscv-vector-gnu/bin/riscv64-unknown-elf-objdump
LINKSCRIPT:=riscv64-virt.ld

TARGET_ISA=riscv
GEM5_PATH=/your-gem5-path

CFLAGS=-O3 -march=rv64imafdcv -fno-tree-vectorize -Wall -I$(GEM5_PATH)/include
LDFLAGS=-L$(GEM5_PATH)/util/m5/build/$(TARGET_ISA)/out -lm5 -static

SRC_DIRS:=test
CUR_DIRS=$(shell pwd)
SRC_DIRS:=$(patsubst %,$(CUR_DIRS)/%,$(SRC_DIRS))
VPATH:=$(SRC_DIRS)


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
