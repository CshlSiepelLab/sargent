import msprime

ts = msprime.sim_ancestry(
    samples=50, sequence_length=5e6,
    recombination_rate=1e-8, population_size=1e4, random_seed=42,
)
ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=42)
ts.dump("example.trees")
print(f"wrote example.trees: {ts.num_samples} samples, "
      f"{ts.num_trees} trees, {ts.num_sites} sites")
