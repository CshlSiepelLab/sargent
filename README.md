# SARGENT: Sparse ARG-based ENcoding of Traits

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24606224.svg)](https://doi.org/10.5281/zenodo.24606224)

SARGENT fits Bayesian sparse linear mixed models for complex traits that
combine effects of typed SNPs with a random effect derived from an ancestral
recombination graph (ARG), capturing heritability not tagged by the typed SNPs.

## Installing with Homebrew (macOS and Linux)

```bash
brew install cshlsiepellab/tools/sargent
```

The example files are installed in `$(brew --prefix)/share/sargent/examples`.

## Building from source

### Dependencies

#### Required

- C compiler (clang or gcc)
- CMake >= 3.16

#### Optional

- OpenMP — used automatically if available; on macOS install via Homebrew:
  ```bash
  brew install libomp
  ```


### Building

```bash
git clone --recurse-submodules https://github.com/CshlSiepelLab/sargent.git
cd sargent
mkdir build && cd build
cmake ..
make
```

The [tskit](https://tskit.dev/) C library (C_1.3.1) is obtained as a git submodule from the tskit repository and compiled into `sargent`; no separate tskit installation is needed.  (Note that GitHub's automatically generated source archives do not include submodules, so clone the repository as shown.) OpenMP is detected and used automatically if available.

The `sargent` binary will be at `build/src/progs/sargent`.

### Installation

To install the binary to a custom prefix (e.g. `~/sargent`):

```bash
cmake .. -DCMAKE_INSTALL_PREFIX=~/sargent
make && make install
```

Add to your shell profile to put `sargent` on your PATH:

```bash
export PATH="$HOME/sargent/bin:$PATH"
```

## Usage

```bash
sargent [OPTIONS] <phenotypes> <input.trees> [<geno_prefix>]
```

Run `sargent --help` for a full description of arguments, options, and
output files.  The `examples/` directory contains a small simulated dataset
that illustrates the input formats; see `examples/README.md`.

## License

SARGENT is distributed under the BSD 3-Clause License (see `LICENSE`).  It
includes the [tskit](https://tskit.dev/) C library (MIT License) as a
submodule; see `AUTHORS`.

## Contact

Please report bugs or questions via
[GitHub issues](https://github.com/CshlSiepelLab/sargent/issues).

## Citation

If you use SARGENT, please cite the software as

> Kuru N, Siepel A. SARGENT: Sparse ARG-based ENcoding of Traits.
> Zenodo. https://doi.org/10.5281/zenodo.24606224

(see also `CITATION.cff`).
