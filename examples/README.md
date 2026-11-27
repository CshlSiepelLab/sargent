# Example input files

A small simulated dataset that shows the input formats `sargent` expects:
200 training individuals and 50 held-out test individuals, on a 500 kb
region.

| File | Contents |
|------|----------|
| `toy.phen` | Phenotypes for the training individuals, one value per line, in the same order as the individuals in `toy.trees` and `toy.fam`. |
| `toy.trees` | ARG for the training individuals, in [tskit](https://tskit.dev/) `.trees` format.  Diploid individuals are taken from the tree sequence's individual table, in order. |
| `toy.bed`, `toy.bim`, `toy.fam` | Typed genotypes for the training individuals, in PLINK binary format.  Only a subset of the variants in the ARG are typed. |
| `toy.full.trees` | ARG for the training and test individuals together, with the training individuals first (needed for `--predict`). |
| `toy.test.bed`, `toy.test.bim`, `toy.test.fam` | Typed genotypes for the test individuals (needed for `--predict`). |
| `toy.test.phen` | Phenotypes for the test individuals, used to assess prediction accuracy (needed for `--predict`). |

Fit the model:

```bash
sargent --stage1-gibbs --joint --uniform-arg -o toy_fit \
    examples/toy.phen examples/toy.trees examples/toy
```

Fit the model and predict phenotypes for the test individuals:

```bash
sargent --stage1-gibbs --joint --uniform-arg --predict -o toy_pred \
    examples/toy.phen examples/toy.trees examples/toy
```

With `--predict`, the test-set files are found from the genotype prefix
(`examples/toy`) by default.  Each command takes a few seconds.  The first run
also writes a cache file (`examples/toy.gall`) next to the `.trees` file,
which later runs reuse.

These files were generated with the repository's simulator:

```bash
python scripts/simulate.py --n-indiv 200 --n-test 50 --length 5e5 \
    --bed-coverage 0.25 --h2-snp 0.3 --h2-arg 0.3 --seed 1 --out toy
```
