#include <math.h>
#include <stdio.h>

#include "subsample.h"
#include "genotype.h"

int subsample_genotype(struct o_genotype *ogt, double target_maf, double margin, double max_mgf) {
  double maf;
  double mgf;

  maf = ogt_maf(ogt);
  mgf = ogt_mgf(ogt);

  if (mgf > max_mgf) {
    return 1;
  }

  if (maf > target_maf - margin & maf < target_maf + margin) {
    return 0;
  }

  if (maf < target_maf & ogt->hom_major > 0) {
    ogt->hom_major--;
    return subsample_genotype(ogt, target_maf, margin, max_mgf);
  }
  if (maf > target_maf & ogt->hom_minor > 0) {
    ogt->hom_minor--;
    return subsample_genotype(ogt, target_maf, margin, max_mgf);
  }
  if (maf > target_maf & ogt->het > 0) {
    ogt->het--;
    return subsample_genotype(ogt, target_maf, margin, max_mgf);
  }
  return 1;
}
