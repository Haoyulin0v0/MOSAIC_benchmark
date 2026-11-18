# Validation scope

The automated scientific tests independently check:

- All168 combinations of four lengths,14 identities and three models: alignment
  path replay, M/X/I/D consumption, C hotspot budget and same-seed determinism.
- Different RNG label tuples remain distinct; fixed floating-point expression
  order is preserved for the150bp/85% rounding boundary.
- Distinct, non-overlapping source loci across all four lengths; perfect-model
  deduplication; exact short/native target bytes and sample target counts.
- A tampered generated input is rejected by its frozen manifest.
- Reverse-strand indel paths agree across PAF, MAF, SAM and BLAST aligned strings.
- Region90 clips both axes, requires the correct strand, excludes gaps, and applies
  the paired-length gate per whole native HSP. Small HSPs cannot be pooled to pass.
- Bounding-box-only PAF and inconsistent path endpoints are rejected.
- Successful and timed-out native subprocesses retain receipts and partial output,
  charge the shared budget, and refuse unaffordable work before launching it.

The six-tool synthetic integration job then exercises the actual single workflow,
indexes, perfect gates, native adapters, long-background not-supported state for
original Bowtie, paired summaries, and separate timing. GitHub Actions runs these
checks and uploads its logs/summary/toolchain as reviewable evidence.

This is software and interface validation. It does not certify the scientific
adequacy of1000 loci, natural evolutionary realism of U/C/I, family-independent
test sources, suitability of a100bp/90% threshold, or full-size compute feasibility.
Those require the declared production inputs, threshold/resource qualification,
and interpretation of real measurements. Demo recovery must not enter paper
benchmark results. Local validation details are recorded in the PR description;
CI statuses apply only to the commit tested.
