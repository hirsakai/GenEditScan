# GenEditScan: Gene Editing and Scanning tools
This repository provides programs for the k-mer analysis to detect unintended external DNA in a host genome.

GenEditScan is based upon kmer (https://github.com/taitoh1970/kmer).

## Installation

### Linux
Move to the GenEditScan src directory and type `make`.

### Mac
First of all, install Homebrew. See the following page for installation.

https://brew.sh/

Next, install llvm.

`brew install llvm`

In the GenEditScan src directory, just type `make`. If you don't have `make` in your computer, you also have to install it.

## How to use the programs
`./geneditscan kmer -v vector.fasta -m mutant_read1.fastq.gz,mutant_read2.fastq.gz -w wildtype_read1.fastq.gz,wildtype_read2.fastq.gz -k kmer -o out_prefix`

The input files and options to be specified are as follows:

`vector.fasta`: Vector sequence (FASTA format)  
`read1.fastq.gz`: R1 (forward) read  
`read2.fastq.gz`: R2 (reverse) read  
`kmer`: # of k (This value must be 8 or more and 20 is recommended)  
`out_prefix`: Names used as a prefix of output files  

The output files are as flollows:

`out_prefix.statistics.txt`: K-mer sequence detection results file
`out_prefix.outside.txt`: Analysis results file around the detected k-mer sequences
`out_prefix.mutant.merFreq.txt`: Mutant's mer frequency file
`out_prefix.wildtype.merFreq.txt`: Wild type's mer frequency file

## All options
Usage : ./geneditscan kmer [options]

[required]
-v | --vector   : Vector file
-m | --mutant   : Mutant files (connect with comma)
-w | --wild     : Wild type files (connect with comma)

[optional]
-k | --kmer     : K-mer (20)
-f | --fdr      : Threshold by FDR (0.01)
-b | --bases    : Number of bases on each side (5)
-o | --out      : Output prefix (out_prefix)
-t | --threads  : Number of threads (all threads)
-l | --length   : Maximum read length (512)
-r | --read     : Number of lines of Fastq file to be read in memory (10000000)
-i | --interval : Log output interval (1000000)
-h | --help     : Print this menu

## Dependencies
Netlib Cephes library (cprob and cmath)
https://netlib.org/cephes/

## Citation
If you're using GenEditScan in your work, please cite:

Itoh T, Onuki R, Tsuda M, Oshima M, Endo M, Sakai H, Tanaka T, Ohsawa R, Tabei Y. Foreign DNA detection by high-throughput sequencing to regulate genome-edited agricultural products. Sci Rep. 2020 Mar 18;10(1):4914. doi: 10.1038/s41598-020-61949-5. PMID: 32188926; PMCID: PMC7080720.
https://www.nature.com/articles/s41598-020-61949-5
