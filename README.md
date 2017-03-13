# Subsample genotype data

In genome wide association (GWA) studies a recurring problem is differences in statistical power due to sample sizes and rare sequence variants. One way to account for this is to subsample the data in order to "lock" the minor allele frequency (MAF) in the data set, *i.e.* all SNPs will have the same MAF.

This is a tool that will help with this subsampling.

## Requirements

- CMake
- htslib
- argp

On Mac, argp can can be installed through Homebrew with the argp-standalone formula.

## Build

```bash
git clone --recursive https://github.com/maehler/vcfsubsample
cd vcfsubsample
mkdir build
cd build
cmake ..
make
make install
```

## Usage

```
> vcfsubsample --help
Usage: vcfsubsample [OPTION...] VCF
vcfsubsample -- subsample a VCF file in order to fix the minor allele frequency
across SNPs

     --exact-samples        Subsample to exactly --min-samples
     --keep-missing         Keep SNPs for which subsampling is not possible
     --maf=FLOAT            Target MAF to aim for
     --margin=FLOAT         Allow target MAF within this margin
     --max-mgf=FLOAT        Maximum genotype frequency to allow for each SNP
     --min-samples=N        Minimum number of samples to allow
     --samplenames          Output names of suggested samples to use for
                            subsampling instead of the number of genotypes of
                            each class
 -?, --help                 Give this help list
     --usage                Give a short usage message
 -V, --version              Print program version

Report bugs to <niklas.mahler@gmail.com>.
```

Genotype data can be read from stdin by using `-` for the VCF filename.
