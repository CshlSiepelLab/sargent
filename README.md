# sargent

## Dependencies

- C compiler (clang or gcc)
- CMake >= 3.16
- [Meson](https://mesonbuild.com/) and [Ninja](https://ninja-build.org/) (used to build the bundled tskit library)

On macOS with conda:
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

The build system automatically compiles [tskit](https://tskit.dev/) C_1.3.1 from the bundled submodule — no separate tskit installation is needed.

The `sargent` binary will be at `build/src/progs/sargent`.

## Optional: PHAST support

If you have [PHAST](http://compgen.cshl.edu/phast/) installed, pass its prefix to cmake:

```bash
cmake .. -DPHAST_ROOT=/path/to/phast
```
