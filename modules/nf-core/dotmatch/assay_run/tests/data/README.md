# Workflow Test Fixtures

These tiny CRISPR guide-counting fixtures are shared by local workflow examples
and wrapper tests. They are intentionally small enough for Planemo, nf-test, and
CI smoke checks.

`sample_a.fastq` covers the core DotMatch outcomes:

- `ambiguous`: `ACGT` exactly matches `guide_a`, but `guide_b` is inside the
  configured one-edit radius.
- `ambiguous`: `ACGG` is one Hamming edit from both `guide_a` and `guide_b`.
- `unmatched`: `CCCC` is outside the configured edit threshold.
- `invalid`: `AC` is shorter than the configured guide window.

`sample_b.fastq` adds a second sample with one unique exact `guide_c`
assignment and one unmatched read so MAGeCK-style multi-sample output is
exercised.
