#include <math.h>

#include "genotype.h"

int gt_count_samples(struct genotype *gt) {
  return gt->hom_ref + gt->hom_alt + gt->het;
}

int ogt_count_samples(struct o_genotype *ogt) {
  return ogt->hom_major + ogt->hom_minor + ogt->het;
}

int gt_to_ogt(struct genotype *gt, struct o_genotype *ogt) {
  ogt->het = gt->het;
  if (2 * gt->hom_ref + gt->het > 2 * gt->hom_alt + gt->het) {
    ogt->hom_major = gt->hom_ref;
    ogt->hom_minor = gt->hom_alt;
    return GENOTYPE_SAME;
  } else {
    ogt->hom_major = gt->hom_alt;
    ogt->hom_minor = gt->hom_ref;
    return GENOTYPE_SWITCHED;
  }
}

double gt_maf(struct genotype *gt) {
  double n_samples = gt_count_samples(gt);
  double ref_freq = (2 * gt->hom_ref + gt->het) / (2 * n_samples);
  return fmin(ref_freq, 1 - ref_freq);
}

double ogt_maf(struct o_genotype *ogt) {
  double n_samples = ogt_count_samples(ogt);
  return (2 * ogt->hom_minor + ogt->het) / (2 * n_samples);
}

double gt_mgf(struct genotype *gt) {
  double n_samples = gt_count_samples(gt);
  return fmax(fmax(gt->hom_ref, gt->het), gt->hom_alt) / n_samples;
}

double ogt_mgf(struct o_genotype *ogt) {
  double n_samples = ogt_count_samples(ogt);
  return fmax(fmax(ogt->hom_major, ogt->het), ogt->hom_minor) / n_samples;
}
