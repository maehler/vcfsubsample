#ifndef SUBSAMPLE_H
#define SUBSAMPLE_H

#include "genotype.h"

#define SUBSAMPLE_OK 0
#define SUBSAMPLE_E_MGF 1
#define SUBSAMPLE_E_SAMPLES 2
#define SUBSAMPLE_E_UNKNOWN 3

int subsample_genotype(struct o_genotype *ogt, double target_maf, double margin,
  double max_mgf, int min_samples);

#endif
