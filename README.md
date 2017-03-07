# Subsample genotype data

In genome wide association (GWA) studies a recurring problem is differences in statistical power due to sample sizes and rare sequence variants. One way to account for this is to subsample the data in order to "lock" the minor allele frequency (MAF) in the data set, *i.e.* all SNPs will have the same MAF.

This is a tool that will help with this subsampling.

## Build

The build requires [htslib](https://github.com/samtools/htslib) and CMake. htslib is included as a submodule, so the easiest is to clone the repository recursively.

```bash
git clone --recursive ... && cd vcfsubsample
mkdir build && cd build
cmake ..
make
```
