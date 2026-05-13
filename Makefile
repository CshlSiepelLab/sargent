BUILD_DIR  := build
PREFIX     := $(HOME)/sargent
TSKIT_ROOT ?= /usr/local/Caskroom/miniforge/base/envs/sargent
PHAST_ROOT ?=

CMAKE_OPTS := -DTSKIT_ROOT=$(TSKIT_ROOT) \
              -DCMAKE_C_COMPILER=/opt/homebrew/opt/llvm/bin/clang \
              -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
ifneq ($(PHAST_ROOT),)
  CMAKE_OPTS += -DPHAST_ROOT=$(PHAST_ROOT)
endif

.PHONY: all install clean

all:
	cmake -B $(BUILD_DIR) $(CMAKE_OPTS)
	cmake --build $(BUILD_DIR)
	ln -sf $(abspath $(BUILD_DIR))/compile_commands.json compile_commands.json

install: all
	cmake --install $(BUILD_DIR) --prefix $(PREFIX)

clean:
	rm -rf $(BUILD_DIR)
