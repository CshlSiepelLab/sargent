#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <tskit.h>
#include <phast/misc.h>
#include "sargent.help"

#define CHECK(ret) do {                                                 \
    if ((ret) < 0) {                                                    \
        fprintf(stderr, "tskit error at %s:%d: %s\n",                  \
                __FILE__, __LINE__, tsk_strerror(ret));                 \
        exit(EXIT_FAILURE);                                             \
    }                                                                   \
} while (0)

static void print_summary(const tsk_treeseq_t *ts) {
  int ret;
  tsk_size_t n = tsk_treeseq_get_num_samples(ts);
  const tsk_id_t *samples = tsk_treeseq_get_samples(ts);
  double L = tsk_treeseq_get_sequence_length(ts);

  printf("samples        : %lld\n", (long long) n);
  printf("nodes          : %lld\n", (long long) tsk_treeseq_get_num_nodes(ts));
  printf("edges          : %lld\n", (long long) tsk_treeseq_get_num_edges(ts));
  printf("trees          : %lld\n", (long long) tsk_treeseq_get_num_trees(ts));
  printf("sites          : %lld\n", (long long) tsk_treeseq_get_num_sites(ts));
  printf("mutations      : %lld\n", (long long) tsk_treeseq_get_num_mutations(ts));
  printf("sequence length: %.6g\n", L);

  if (tsk_treeseq_get_num_sites(ts) > 0) {
    double pi;
    double windows[2] = {0.0, L};
    tsk_size_t sizes[1] = {n};
    ret = tsk_treeseq_diversity(
        ts, 1, sizes, samples, 1, windows, TSK_STAT_SITE, &pi);
    CHECK(ret);
    printf("diversity (pi) : %.6g\n", pi);
  }
}

int main(int argc, char *argv[]) {
  signed char c;
  int opt_idx;
  int silent = 0;

  struct option long_opts[] = {
    {"silent",  0, 0, 'z'},
    {"version", 0, 0, 'v'},
    {"help",    0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "zvh", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'z':
      silent = 1;
      break;
    case 'v':
      printf("sargent version %s\n", SARGENT_VERSION);
      exit(0);
    case 'h':
      printf("%s", HELP);
      exit(0);
    case '?':
      fprintf(stderr, "Bad argument.  Try 'sargent -h'.\n");
      exit(EXIT_FAILURE);
    }
  }

  if (optind != argc - 1) {
    fprintf(stderr, "ERROR: input file required.  Try 'sargent -h'.\n");
    exit(EXIT_FAILURE);
  }

  const char *infile = argv[optind];
  if (!silent)
    fprintf(stderr, "Reading %s...\n", infile);

  /* temporary: exercise phast_fopen to verify phast linking */
  FILE *f = phast_fopen(infile, "r");
  fclose(f);

  tsk_treeseq_t ts;
  int ret = tsk_treeseq_load(&ts, infile, 0);
  CHECK(ret);

  print_summary(&ts);

  tsk_treeseq_free(&ts);

  if (!silent)
    fprintf(stderr, "Done.\n");

  return 0;
}
