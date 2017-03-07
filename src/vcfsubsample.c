#include <stdio.h>
#include <stdlib.h>
#include <argp.h>

#include <htslib/vcf.h>

int main(int argc, char *argv[]) {
  fprintf(stdout, "Subsample that VCF!\n");

  if (argc != 2) {
    fprintf(stderr, "usage: %s <vcf_file>\n", argv[0]);
    return 1;
  }

  htsFile *inf = bcf_open(argv[1], "r");
  if (inf == NULL) {
    return EXIT_FAILURE;
  }

  bcf_hdr_t *hdr = bcf_hdr_read(inf);
  fprintf(stdout, "file %s contains %i samples\n", argv[1], bcf_hdr_nsamples(hdr));

  bcf_hdr_destroy(hdr);
  bcf_close(inf);

  return 0;
}
