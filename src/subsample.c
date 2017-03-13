#include <math.h>
#include <stdio.h>

#include "subsample.h"
#include "genotype.h"

int subsample_genotype(struct o_genotype *ogt, double target_maf, double margin,
    double max_mgf, int min_samples, unsigned int exact_samples) {
  double maf;
  double mgf;

  maf = ogt_maf(ogt);
  mgf = ogt_mgf(ogt);

  if (mgf > max_mgf) {
    return SUBSAMPLE_E_MGF;
  }

  if (ogt_count_samples(ogt) < min_samples) {
    return SUBSAMPLE_E_SAMPLES;
  }

  if (maf > target_maf - margin && maf < target_maf + margin) {
    if (!exact_samples || (exact_samples && ogt_count_samples(ogt) == min_samples)) {
      return SUBSAMPLE_OK;
    }
  }

  if (maf < target_maf && ogt->hom_major > 0) {
    ogt->hom_major--;
    return subsample_genotype(ogt, target_maf, margin, max_mgf, min_samples, exact_samples);
  }
  if (maf > target_maf && ogt->hom_minor > 0) {
    ogt->hom_minor--;
    return subsample_genotype(ogt, target_maf, margin, max_mgf, min_samples, exact_samples);
  }
  if (maf > target_maf && ogt->het > 0) {
    ogt->het--;
    return subsample_genotype(ogt, target_maf, margin, max_mgf, min_samples, exact_samples);
  }
  if (maf == target_maf && exact_samples) {
    ogt->hom_major--;
    return subsample_genotype(ogt, target_maf, margin, max_mgf, min_samples, exact_samples);
  }
  return SUBSAMPLE_E_UNKNOWN;
}
