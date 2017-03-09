#ifndef GENOTYPE_H
#define GENOTYPE_H

// A genotype struct that is based on reference and alternative alleles
struct genotype {
  unsigned int hom_ref, hom_alt, het;
};

// A genotype struct that is based on major and minor alleles
struct o_genotype {
  unsigned int hom_major, hom_minor, het;
};

int gt_to_ogt(struct genotype *gt, struct o_genotype *ogt);
double gt_maf(struct genotype *gt);
double ogt_maf(struct o_genotype *ogt);
double gt_mgf(struct genotype *gt);
double ogt_mgf(struct o_genotype *ogt);

#endif
