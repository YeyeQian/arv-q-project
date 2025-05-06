<<<<<<< HEAD
TAR=test_kem
OBJ=cbd.o fips202.o indcpa.o kem.o ntt.o poly.o polyvec.o reduce.o sha256.o sha512.o symmetric-shake.o verify.o fips202_custom.o poly_custom.o symmetric-shake-custom.o cbd_custom.o polyvec_custom.o indcpa_custom.o kem_custom.o ntt_custom.o poly_reorg_custom.o poly_custom_redundant.o poly_custom_shuffling.o
=======
TAR=test_kem_kyber512_vlen1024
OBJ=cbd.o fips202.o indcpa.o kem.o ntt.o poly.o polyvec.o reduce.o sha256.o sha512.o symmetric-shake.o verify.o fips202_custom.o poly_custom.o symmetric-shake-custom.o cbd_custom.o polyvec_custom.o indcpa_custom.o kem_custom.o ntt_custom.o poly_reorg_custom.o
>>>>>>> ebdbc1a76475082974a90356871ef1a4d5f927e8
# OBJ+=test_ntt.o
# OBJ+=test_poly_compress.o
# OBJ+=test_poly_mul.o
# OBJ+=test_fips202.o
# OBJ+=test_rejsample.o
# OBJ+=test_indcpa_custom.o
# OBJ+=test_poly.o
# OBJ+=test_symmetric-shake.o
# OBJ+=test_polyvec.o
# OBJ+=test_poly_reorg.o
# OBJ+=test_indcpa_custom_asm.o
OBJ+=test_kem.o
<<<<<<< HEAD
# OBJ+=test_poly_against_sca.o
=======
# OBJ+=poly_custom_redundant.o
# OBJ+=ntt_custom_sca.o
>>>>>>> ebdbc1a76475082974a90356871ef1a4d5f927e8

GCC:=/home/data/userhome/zhaoyifan/MyToolChains/xuantie-riscv-gnu-c908e/bin/riscv64-unknown-elf-gcc
OBJDUMP:=/home/data/userhome/zhaoyifan/MyToolChains/xuantie-riscv-gnu-c908e/bin/riscv64-unknown-elf-objdump
# GCC:=/home/data/userhome/zhaoyifan/MyToolChains/my-riscv-vector-gnu/bin/riscv64-unknown-elf-gcc
# OBJDUMP:=/home/data/userhome/zhaoyifan/MyToolChains/my-riscv-vector-gnu/bin/riscv64-unknown-elf-objdump
LINKSCRIPT:=riscv64-virt.ld

TARGET_ISA=riscv
GEM5_PATH=/your-gem5-path

CFLAGS=-O3 -march=rv64imafdcv -fno-tree-vectorize -Wall -I$(GEM5_PATH)/include
# CFLAGS=-O3 -march=rv64imafdcv -Wall -I$(GEM5_PATH)/include
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