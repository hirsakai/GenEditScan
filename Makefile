# Makefile for GenEditScan
# Copyright (C) 2018 National Agriculture and Food Research Organization (NARO)
#
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	CC := g++
else
	CC := /opt/homebrew/opt/llvm/bin/clang++
endif
#
LDFLAGS := -L./cprob

CFLAGS := -std=c++17 -O3 -Wall -fopenmp

COBJS := bitwise_operation.o complementary.o fastq_extension.o fastq_match.o gtest.o kmer_extension.o kmer_match.o \
		main.o statistics_file.o vector_sequence.o

LIBS := -lz -lprob

MAIN := geneditscan

%.o: %.cpp
	$(CC) -c $(CFLAGS) -o $@ $<

# all
all: $(MAIN)

# kmer_analysis
$(MAIN): $(COBJS) ./cprob/libprob.a
	$(CC) $(CFLAGS) -o $(MAIN) $(COBJS) $(LIBS) $(LDFLAGS)

./cprob/libprob.a:
	cd $(@D); $(MAKE)

# clean
clean:
	-rm -f $(COBJS) $(MAIN) *~
	cd ./cprob; $(MAKE) clean

.PHONY: all clean

# dependencies (g++ -MM source.cpp)
bitwise_operation.o: bitwise_operation.cpp bitwise_operation.h options.h
complementary.o: complementary.cpp complementary.h
fastq_extension.o: fastq_extension.cpp fastq_extension.h \
 bitwise_operation.h options.h
fastq_match.o: fastq_match.cpp fastq_match.h bitwise_operation.h \
 options.h
gtest.o: gtest.cpp gtest.h options.h
kmer_extension.o: kmer_extension.cpp kmer_extension.h bitwise_operation.h \
 options.h statistics_file.h gtest.h outside_data.h fastq_extension.h \
 complementary.h
kmer_match.o: kmer_match.cpp kmer_match.h bitwise_operation.h options.h \
 statistics_file.h gtest.h outside_data.h fastq_match.h vector_sequence.h
main.o: main.cpp bitwise_operation.h options.h statistics_file.h gtest.h \
 outside_data.h kmer_match.h fastq_match.h kmer_extension.h \
 fastq_extension.h
statistics_file.o: statistics_file.cpp statistics_file.h gtest.h \
 options.h outside_data.h complementary.h
vector_sequence.o: vector_sequence.cpp vector_sequence.h \
 bitwise_operation.h options.h complementary.h
