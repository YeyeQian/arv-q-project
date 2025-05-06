TAR=test_nist
OBJ=codec.o common.o fft.o fpr.o keygen.o nist.o rng.o shake.o sign.o vrfy.o codec_custom.o common_custom.o fft_custom.o rng_custom.o shake_custom.o sign_custom_no_rvvf.o sign_custom_with_rvvf.o utility_custom.o vrfy_custom.o nist_custom.o
#OBJ+=test_fft.o
#OBJ+=test_vrfy.o
#OBJ+=test_common.o
#OBJ+=test_sign.o
OBJ+=test_nist.o

#GCC:=/home/data/userhome/zhaoyifan/MyToolChains/xuantie-riscv-gnu-c908e/bin/riscv64-unknown-elf-gcc
#OBJDUMP:=/home/data/userhome/zhaoyifan/MyToolChains/xuantie-riscv-gnu-c908e/bin/riscv64-unknown-elf-objdump
GCC:=/your-xuantie-gnu-path/bin/riscv64-unknown-elf-gcc
OBJDUMP:=/your-xuantie-gnu-path/bin/riscv64-unknown-elf-objdump
LINKSCRIPT:=riscv64-virt.ld

TARGET_ISA=riscv

CFLAGS=-O3 -march=rv64imafdcv -fno-tree-vectorize -Wall
LDFLAGS=-static
LDLIBS=-lm

.DEFAULT_GOAL := default

.PHONY: default
default: $(addsuffix .riscv,$(TAR))

.PHONY: dumps
dumps: $(addsuffix .dump,$(TAR))

%.o:%.c
	$(GCC) $(CFLAGS) -c $^ -o $@

%.riscv:$(OBJ)
	$(GCC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

%.dump: %.riscv
	$(OBJDUMP) -D $< > $@

.PHONY:
clean:
	rm -f *.riscv *.o *.dump
