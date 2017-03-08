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
```
