#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include <htslib/vcf.h>

#include "vcfsubsample.h"

const char *argp_program_version = vcfsubsample_VERSION_STRING;

// Documentation
static char doc[] =
  "vcfsubsample -- subsample a VCF file in order to fix the minor allele frequency across SNPs";
static char args_doc[] = "VCF";

// Keys for arguments without short options
#define OPT_MAF 1
#define OPT_MARGIN 2
#define OPT_MAX_MGF 3
#define OPT_MIN_SAMPLES 4

// Options
static struct argp_option options[] = {
  {"maf",         OPT_MAF,         "FLOAT", 0, "MAF to aim for"},
  {"margin",      OPT_MARGIN,      "FLOAT", 0, "MAF +/- margin is ok"},
  {"max-mgf",     OPT_MAX_MGF,     "FLOAT", 0, "Maximum genotype frequency to allow for each SNP"},
  {"min-samples", OPT_MIN_SAMPLES, "N",     0, "Minimum number of samples to allow"},
  {0}
};

struct arguments {
  char *args[1];
  double maf, margin, max_mgf;
  unsigned int min_samples;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
  char *rest;
  struct arguments *arguments = state->input;

  switch(key) {
    case OPT_MAF:
      arguments->maf = strtod(arg, &rest);
      break;
    case OPT_MARGIN:
      arguments->margin = strtod(arg, &rest);
      if (*rest != 0 | arguments->margin < 0) {
        fprintf(stderr, "Error: margin must be a number beteween 0 and 1\n");
        return ARGP_ERR_UNKNOWN;
      }
      break;
    case OPT_MAX_MGF:
      arguments->max_mgf = strtod(arg, &rest);
      if (*rest != 0 | arguments->max_mgf < 0 | arguments->max_mgf > 1) {
        fprintf(stderr, "Error: max-mgf must be a number between 0 and 1\n");
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

  arguments.maf = DEFAULT_MAF;
  arguments.margin = DEFAULT_MARGIN;
  arguments.max_mgf = DEFAULT_MAX_MGF;
  arguments.min_samples = DEFAULT_MIN_SAMPLES;

  argp_parse(&argp, argc, argv, 0, 0, &arguments);

  htsFile *inf = bcf_open(arguments.args[0], "r");
  if (inf == NULL) {
    return EXIT_FAILURE;
  }

  bcf_hdr_t *hdr = bcf_hdr_read(inf);
  fprintf(stdout, "file %s contains %i samples\n", arguments.args[0], bcf_hdr_nsamples(hdr));

  bcf_hdr_destroy(hdr);
  bcf_close(inf);

  return 0;
}
