#ifndef GENOTYPE_H
#define GENOTYPE_H

struct genotype {
  unsigned int hom_ref, hom_alt, het;
};

double gt_maf(struct genotype*);
double gt_mgf(struct genotype*);

#endif
