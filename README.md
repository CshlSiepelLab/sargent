# SARGENT: Sparse ARG-based ENcoding of Traits

## Dependencies

### Required

- C compiler (clang or gcc)
- CMake >= 3.16
- [Meson](https://mesonbuild.com/) and [Ninja](https://ninja-build.org/) (to build the bundled tskit library)

### Optional

- OpenMP — used automatically if available; on macOS install via Homebrew:
  ```bash
  brew install libomp
  ```

On macOS with conda (for cmake/meson/ninja):
```bash
conda install -c conda-forge cmake meson ninja
```

## Building

```bash
git clone --recurse-submodules https://github.com/asiepel/sargent.git
cd sargent
mkdir build && cd build
cmake ..
make
```

The build system automatically compiles [tskit](https://tskit.dev/) C_1.3.1 from the bundled submodule — no separate tskit installation is needed. OpenMP is detected and used automatically if available.

The `sargent` binary will be at `build/src/progs/sargent`.

## Installation

To install the binary to a custom prefix (e.g. `~/sargent`):

```bash
cmake .. -DCMAKE_INSTALL_PREFIX=~/sargent
make && make install
```

Add to your shell profile to put `sargent` on your PATH:

```bash
export PATH="$HOME/sargent/bin:$PATH"
```
