#include <math.h>

#include "genotype.h"

int count_samples(struct genotype *genotype) {
  return genotype->hom_ref + genotype->hom_alt + genotype->het;
}

double gt_maf(struct genotype *genotype) {
  double n_samples = count_samples(genotype);
  double ref_freq = (2 * genotype->hom_ref + genotype->het) / (2 * n_samples);
  return fmin(ref_freq, 1 - ref_freq);
}

double gt_mgf(struct genotype *genotype) {
  double n_samples = count_samples(genotype);
  return fmax(fmax(genotype->hom_ref, genotype->het), genotype->hom_alt) / n_samples;
}
