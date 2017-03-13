#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include <htslib/vcf.h>
#include <stdbool.h>

#include "vcfsubsample.h"
#include "subsample.h"
#include "genotype.h"

const char *argp_program_version = vcfsubsample_VERSION_STRING;
const char *argp_program_bug_address = vcfsubsample_BUG_ADDRESS;

// Documentation
static char doc[] =
  "vcfsubsample -- subsample a VCF file in order to fix the minor allele frequency across SNPs";
static char args_doc[] = "VCF";

// Keys for arguments without short options
#define OPT_MAF 1
#define OPT_MARGIN 2
#define OPT_MAX_MGF 3
#define OPT_MIN_SAMPLES 4
#define OPT_SAMPLE_NAMES 5
#define OPT_KEEP_MISSING 6

// Options
static struct argp_option options[] = {
  {"maf",         OPT_MAF,          "FLOAT", 0, "Target MAF to aim for"},
  {"margin",      OPT_MARGIN,       "FLOAT", 0, "Allow target MAF within this margin"},
  {"max-mgf",     OPT_MAX_MGF,      "FLOAT", 0, "Maximum genotype frequency to allow for each SNP"},
  {"min-samples", OPT_MIN_SAMPLES,  "N",     0, "Minimum number of samples to allow"},
  {"samplenames", OPT_SAMPLE_NAMES,  0,      0, "Output names of suggested samples to use for "
                                               "subsampling instead of the number of genotypes of each class"},
  {"keep-missing", OPT_KEEP_MISSING, 0,      0, "Keep SNPs for which subsampling is not possible"},
  {0}
};

struct arguments {
  char *args[1];
  double maf, margin, max_mgf;
  unsigned int min_samples, samplenames, keep_missing;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
  char *rest;
  struct arguments *arguments = state->input;

  switch(key) {
    case OPT_MAF:
      arguments->maf = strtod(arg, &rest);
      if (*rest != 0 | arguments->maf <= 0 | arguments->maf > 0.5) {
        fprintf(stderr, "Error: margin must be a number between 0 (exclusive) and 0.5\n");
        return ARGP_ERR_UNKNOWN;
      }
      break;
    case OPT_MARGIN:
      arguments->margin = strtod(arg, &rest);
      if (*rest != 0 | arguments->margin < 0 | arguments->margin > 0.5) {
        fprintf(stderr, "Error: margin must be a number beteween 0 and 0.5\n");
        return ARGP_ERR_UNKNOWN;
      }
      break;
    case OPT_MAX_MGF:
      arguments->max_mgf = strtod(arg, &rest);
      if (*rest != 0 | arguments->max_mgf <= 0 | arguments->max_mgf > 1) {
        fprintf(stderr, "Error: max-mgf must be a number between 0 (exclusive) and 1\n");
        return ARGP_ERR_UNKNOWN;
      }
      break;
    case OPT_MIN_SAMPLES:
      arguments->min_samples = strtol(arg, &rest, 10);
      if (*rest != 0 | arguments->min_samples < 1) {
        fprintf(stderr, "Error: min-samples must be a positive integer\n");
        return ARGP_ERR_UNKNOWN;
      }
      break;
    case OPT_SAMPLE_NAMES:
      arguments->samplenames = true;
      break;
    case OPT_KEEP_MISSING:
      arguments->keep_missing = true;
      break;
    case ARGP_KEY_ARG:
      if (state->arg_num >= 1) {
        argp_usage(state);
      }

      arguments->args[state->arg_num] = arg;
      break;
    case ARGP_KEY_END:
      if (state->arg_num < 1) {
        argp_usage(state);
      }
      break;
    default:
      return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

static struct argp argp = {options, parse_opt, args_doc, doc};

int main(int argc, char *argv[]) {
  struct arguments arguments;

  int nsnps = 0;

  int ngt;
  int ngt_arr;
  int *gt_arr = NULL;
  int *orig_gt_arr = NULL;
  int n_subset_hom_minor;
  int n_subset_het;
  int n_subset_hom_major;
  int n_subset_samples;

  int sstatus;
  int skipped_mgf = 0;
  int skipped_samples = 0;
  int nseq;
  int gt1, gt2;
  int gt_switch;

  const char **seqnames = NULL;

  arguments.maf = DEFAULT_MAF;
  arguments.margin = DEFAULT_MARGIN;
  arguments.max_mgf = DEFAULT_MAX_MGF;
  arguments.min_samples = DEFAULT_MIN_SAMPLES;
  arguments.samplenames = DEFAULT_SAMPLENAMES;
  arguments.keep_missing = DEFAULT_KEEP_MISSING;

  argp_parse(&argp, argc, argv, 0, 0, &arguments);

  fprintf(stderr, "Arguments:\n"
                  "\tInput file:  %s\n"
                  "\tTarget MAF:  %.3f\n"
                  "\tMAF margin:  %.3f\n"
                  "\tMax MGF:     %.3f\n"
                  "\tMin samples: %i\n",
          strcmp(arguments.args[0], "-") == 0 ? "stdin" : arguments.args[0],
          arguments.maf, arguments.margin,
          arguments.max_mgf, arguments.min_samples);

  htsFile *vcf = bcf_open(arguments.args[0], "r");
  if (vcf == NULL) {
    return EXIT_FAILURE;
  }

  bcf_hdr_t *hdr = bcf_hdr_read(vcf);
  fprintf(stderr, "file %s contains %i samples\n", arguments.args[0], bcf_hdr_nsamples(hdr));
  orig_gt_arr = malloc(bcf_hdr_nsamples(hdr) * sizeof(int));

  seqnames = bcf_hdr_seqnames(hdr, &nseq);

  bcf1_t *rec = bcf_init();

  // Output format
  printf("chrom\tpos\toriginal_maf\tsubsampled_maf\tnsamples\t");
  if (arguments.samplenames) {
    printf("samples\n");
  } else {
    printf("n_hom_major\tn_het\tn_hom_minor\n");
  }

  while (bcf_read(vcf, hdr, rec) == 0) {
    if (!bcf_is_snp(rec) | (rec->n_allele != 2)) {
      // Only support biallelic SNPs
      continue;
    }
    nsnps++;

    struct genotype gt;
    struct o_genotype ogt;

    gt.hom_ref = 0;
    gt.hom_alt = 0;
    gt.het = 0;

    ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);

    for (int i = 0; i < ngt; i += 2) {
      gt1 = bcf_gt_allele(gt_arr[i]);
      gt2 = bcf_gt_allele(gt_arr[i + 1]);
      if (gt1 == gt2 && gt1 == 0) {
        gt.hom_ref++;
        orig_gt_arr[i / 2] = GENOTYPE_HOM_REF;
      } else if (gt1 == gt2) {
        gt.hom_alt++;
        orig_gt_arr[i / 2] = GENOTYPE_HOM_ALT;
      } else {
        gt.het++;
        orig_gt_arr[i / 2] = GENOTYPE_HET;
      }
    }

    gt_switch = gt_to_ogt(&gt, &ogt);

    sstatus = subsample_genotype(&ogt, arguments.maf, arguments.margin,
        arguments.max_mgf, arguments.min_samples);

    if (sstatus != SUBSAMPLE_OK) {
      switch (sstatus) {
        case SUBSAMPLE_E_MGF:
          if (arguments.keep_missing) {
            printf("%s\t%i\t%f\tNA\t0", seqnames[rec->rid], rec->pos + 1, gt_maf(&gt));
            if (arguments.samplenames) {
              printf("\tNA\n");
            } else {
              printf("\tNA\tNA\tNA\n");
            }
          }
          skipped_mgf++;
          break;
        case SUBSAMPLE_E_SAMPLES:
          if (arguments.keep_missing) {
            printf("%s\t%i\t%f\tNA\t0", seqnames[rec->rid], rec->pos + 1, gt_maf(&gt));
            if (arguments.samplenames) {
              printf("\tNA\n");
            } else {
              printf("\tNA\tNA\tNA\n");
            }
          }
          skipped_samples++;
          break;
        case SUBSAMPLE_E_UNKNOWN:
          fprintf(stderr, "Error: an unknown error occured when subsampling");
          return EXIT_FAILURE;
      }
      continue;
    }

    // bcf1_t::pos is zero-based
    printf("%s\t%i\t%f\t%f\t%i",
      seqnames[rec->rid], rec->pos + 1,
      gt_maf(&gt), ogt_maf(&ogt),
      ogt_count_samples(&ogt));

    if (arguments.samplenames) {
      n_subset_hom_major = 0;
      n_subset_hom_minor = 0;
      n_subset_het = 0;
      n_subset_samples = 0;

      printf("\t");

      for (int i = 0; i < bcf_hdr_nsamples(hdr); i++) {
        if (n_subset_hom_major < ogt.hom_major &&
            ((gt_switch == GENOTYPE_SAME && orig_gt_arr[i] == GENOTYPE_HOM_REF) ||
              (gt_switch == GENOTYPE_SWITCHED && orig_gt_arr[i] == GENOTYPE_HOM_ALT))) {
          if (n_subset_samples > 0) {
            printf(",");
          }
          printf("%s", hdr->samples[i]);
          n_subset_hom_major++;
          n_subset_samples++;
        }

        if (n_subset_hom_minor < ogt.hom_minor &&
            ((gt_switch == GENOTYPE_SAME && orig_gt_arr[i] == GENOTYPE_HOM_ALT) ||
              (gt_switch == GENOTYPE_SWITCHED && orig_gt_arr[i] == GENOTYPE_HOM_REF))) {
          if (n_subset_samples > 0) {
            printf(",");
          }
          printf("%s", hdr->samples[i]);
          n_subset_hom_minor++;
          n_subset_samples++;
        }

        if (n_subset_het < ogt.het) {
          if (n_subset_samples > 0) {
            printf(",");
          }
          printf("%s", hdr->samples[i]);
          n_subset_het++;
          n_subset_samples++;
        }
      }
    } else {
      printf("\t%i\t%i\t%i", ogt.hom_major, ogt.het, ogt.hom_minor);
    }
    printf("\n");
  }

  fprintf(stderr, "%i biallelic SNPs in file\n", nsnps);
  fprintf(stderr, "Downsampling impossible for %i SNPs:\n",
    skipped_mgf + skipped_samples);
  fprintf(stderr, "\t%i due to too high MGF\n", skipped_mgf);
  fprintf(stderr, "\t%i due to too few samples remaining\n", skipped_samples);

  bcf_hdr_destroy(hdr);
  bcf_close(vcf);
  bcf_destroy(rec);

  free(seqnames);
  free(orig_gt_arr);

  return EXIT_SUCCESS;
}
