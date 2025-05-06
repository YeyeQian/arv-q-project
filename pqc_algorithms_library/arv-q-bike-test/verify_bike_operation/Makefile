TAR=PQCgenKAT_kem_bike5
OBJ=conversions.o conversions_custom.o decode.o decode_custom.o extGCD_polyinv.o extGCD_polyinv64.o hash_wrapper.o hash_wrapper_custom.o kem.o poly_op_helper.o poly_op_util.o poly_op_util_custom.o sampling.o sampling_custom.o shake_prng.o shake_prng_custom.o threshold.o gf2x_port.o gf2x_custom.o
# OBJ+=temp_test_main.o
# OBJ+=test_gf2x.o
# OBJ+=test_conversion.o
# OBJ+=test_poly_op_util.o
OBJ+=PQCgenKAT_kem.o

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

%.o:%.cpp
	$(GCC) $(CFLAGS) -c $^ -o $@

%.riscv:$(OBJ)
	$(GCC) $(LDFLAGS) $^ -o $@

%.dump: %.riscv
	$(OBJDUMP) -D $< > $@

.PHONY:
clean:
	rm -f *.riscv *.o *.dump
