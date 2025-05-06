TAR=test_decrypt
OBJ=benes.o benes_custom.o bm.o bm_custom.o controlbits.o crypto_hash.o crypto_hash_custom.o decrypt.o decrypt_custom.o encrypt.o encrypt_custom.o gf.o gf_custom.o operations.o operations_custom.o pk_gen.o pk_gen_custom.o root.o root_custom.o sk_gen.o sk_gen_custom.o synd.o synd_custom.o transpose.o transpose_custom.o util.o pk_sk_cons.o
OBJ+=test_decrypt.o

GCC:=/your-xuantie-gnu-path/bin/riscv64-unknown-elf-g++
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
