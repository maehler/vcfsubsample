#ifndef SUBSAMPLE_H
#define SUBSAMPLE_H

#include "genotype.h"

int subsample_genotype(struct o_genotype *ogt, double target_maf, double margin,
  double max_mgf, int min_samples);

#endif
