#!/usr/bin/env python3
"""
analyze_profile.py — resolve and rank hot functions from a samply Firefox-profiler
JSON(.gz) profile.

Usage:
    python3 analyze_profile.py profile.json.gz [profile2.json.gz ...]
    python3 analyze_profile.py --top 20 profile.json.gz

Symbol resolution:
  samply --save-only stores function names as hex offset strings (e.g. "0x3abc").
  These offsets are RELATIVE TO EACH LIBRARY'S LOAD BASE, which samply normalises
  to 0x0 (field libs[i].start == 0x0 for all libraries).

  For shared libraries (dylibs): nm addresses already start near 0 (PIC code),
  so the profile offset == nm offset.  Direct lookup.

  For the main Mach-O executable (sargent, on macOS ARM64): the binary is linked
  with a preferred load address of 0x100000000.  nm therefore shows symbols at
  0x100000000 + offset.  So we add 0x100000000 before the nm lookup.

  These rules are auto-detected from nm: if the lowest text symbol in the nm
  output is ≥ 0x100000000 we classify it as a main-exe binary and add that base;
  otherwise we treat it as a dylib and use the nm address directly.
"""

import argparse
import collections
import gzip
import json
import subprocess
import sys


# ---------------------------------------------------------------------------
# Symbol table helpers
# ---------------------------------------------------------------------------

def load_nm_symbols(binary_path: str):
    """Return (sorted_symbols, base_offset) where sorted_symbols is a list of
    (abs_addr, name) tuples sorted by address, and base_offset is the amount
    to ADD to profile addresses before binary-searching sorted_symbols.
    """
    result = subprocess.run(
        ["nm", "-n", binary_path],
        capture_output=True, text=True
    )
    syms = []
    for line in result.stdout.splitlines():
        parts = line.split()
        if len(parts) == 3 and parts[1] in ("T", "t", "D", "d"):
            try:
                addr = int(parts[0], 16)
                name = parts[2].lstrip("_")
                syms.append((addr, name))
            except ValueError:
                pass
    syms.sort()
    if not syms:
        return [], 0

    # Auto-detect whether this is a main executable (linked at 0x100000000+)
    # or a position-independent dylib (linked near 0x0).
    min_addr = syms[0][0]
    base_offset = 0x100000000 if min_addr >= 0x100000000 else 0
    return syms, base_offset


def lookup_symbol(syms, abs_addr: int) -> str:
    """Binary-search sorted (addr, name) list for the symbol containing abs_addr."""
    lo, hi = 0, len(syms) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        if syms[mid][0] <= abs_addr:
            lo = mid + 1
        else:
            hi = mid - 1
    if hi >= 0:
        return syms[hi][1]
    return "?(unknown)"


# ---------------------------------------------------------------------------
# Profile analysis
# ---------------------------------------------------------------------------

def analyze(profile_path: str, top_n: int = 20, verbose: bool = False):
    print(f"\n{'='*72}")
    print(f"Profile: {profile_path}")
    print(f"{'='*72}")

    open_fn = gzip.open if profile_path.endswith(".gz") else open
    with open_fn(profile_path) as f:
        data = json.load(f)

    # Build per-library symbol tables
    lib_paths = [lib.get("path", "") for lib in data["libs"]]
    lib_syms   = {}   # lib_index -> (sorted_symbols, base_offset)
    for i, path in enumerate(lib_paths):
        if path and path not in lib_syms.values():
            syms, base = load_nm_symbols(path)
            lib_syms[i] = (syms, base)
            if verbose:
                print(f"  lib[{i}] {path.split('/')[-1]}: {len(syms)} symbols, "
                      f"base={hex(base)}")

    # Unwrap the (single) main thread
    threads = data.get("threads", [])
    if not threads:
        print("  (no threads in profile)")
        return

    all_wall_s = 0.0
    for thread_idx, t in enumerate(threads):
        if not t.get("isMainThread", False) and len(threads) > 1:
            continue  # analyse only the main thread for single-threaded programs

        strings           = t["stringArray"]
        func_name_indices = t["funcTable"]["name"]
        func_resources    = t["funcTable"]["resource"]
        stack_prefix      = t["stackTable"]["prefix"]
        stack_frame       = t["stackTable"]["frame"]
        frame_func        = t["frameTable"]["func"]
        sample_stacks     = t["samples"]["stack"]
        resource_lib      = t["resourceTable"]["lib"]

        total = sum(1 for s in sample_stacks if s is not None)
        if total == 0:
            print("  (0 samples)")
            continue

        # Assume 1000 Hz sampling rate (samply default)
        wall_s = total / 1000.0
        all_wall_s += wall_s
        print(f"\nThread '{t.get('name','?')}': {total} samples, "
              f"wall ≈ {wall_s:.1f} s")

        leaf_counts  = collections.Counter()
        by_lib_counts = collections.Counter()

        for s in sample_stacks:
            if s is None:
                continue
            frame_idx = stack_frame[s]
            func_idx  = frame_func[frame_idx]
            addr_str  = strings[func_name_indices[func_idx]]

            res_idx   = func_resources[func_idx]
            lib_idx   = resource_lib[res_idx] if (
                res_idx is not None and res_idx >= 0 and
                res_idx < len(resource_lib)) else None
            lib_path  = lib_paths[lib_idx] if lib_idx is not None else ""
            lib_short = lib_path.split("/")[-1] if lib_path else "?"

            if not addr_str.startswith("0x"):
                leaf_counts[f"{lib_short}::{addr_str}"] += 1
                by_lib_counts[lib_short] += 1
                continue

            addr = int(addr_str, 16)
            if lib_idx is not None and lib_idx in lib_syms:
                syms, base = lib_syms[lib_idx]
                name = lookup_symbol(syms, addr + base)
                leaf_counts[f"{lib_short}::{name}"] += 1
            else:
                leaf_counts[f"{lib_short}::0x{addr:x}"] += 1
            by_lib_counts[lib_short] += 1

        # --- Top-N by self-time ---
        print(f"\n  Top {top_n} functions by self-time:")
        print(f"  {'%wall':>6}  {'samples':>8}  function")
        print(f"  {'------':>6}  {'-------':>8}  --------")
        for name, cnt in leaf_counts.most_common(top_n):
            pct = 100.0 * cnt / total
            print(f"  {pct:6.1f}%  {cnt:8d}  {name}")

        # --- By library ---
        print(f"\n  Time by library:")
        for lib, cnt in by_lib_counts.most_common():
            pct = 100.0 * cnt / total
            print(f"  {pct:6.1f}%  {cnt:8d}  {lib}")

    return all_wall_s


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Resolve and rank hot functions from samply Firefox-profiler JSON(s)."
    )
    parser.add_argument("profiles", nargs="+", metavar="profile.json[.gz]")
    parser.add_argument("--top", type=int, default=20,
                        help="Number of top functions to show (default: 20)")
    parser.add_argument("--verbose", action="store_true",
                        help="Print library symbol counts")
    args = parser.parse_args()

    for path in args.profiles:
        analyze(path, top_n=args.top, verbose=args.verbose)


if __name__ == "__main__":
    main()
