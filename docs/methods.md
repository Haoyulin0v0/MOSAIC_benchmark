# Reproduction and measurement protocol

## Experimental units and source sampling

The production source quota is1000 per reference and length, not1000 mutation
replicates per source. Length cohorts use distinct starts and, by default, mutually
non-overlapping intervals even across lengths. Four cohorts occupy1,900,000bp per
reference before exclusions. Every model, identity and context derived from a locus
shares its source identifier. Different positions do not prove family independence.
An optional `exclude_bed` per reference prevents overlap with registered development
or validation regions. BED coordinates must be0-based half-open; no tool-derived
recovery is consulted for sampling or exclusion.

For a requested length L, a record contributes max(0,record_length−L+1) candidate
starts. An integer is sampled uniformly over the union of these starts, then
converted to record and offset by cumulative weights. Reject non-ACGT windows,
duplicate starts, overlaps with previously accepted intervals, and declared BED
exclusions. Case-insensitive ACGT permits masked repeat sequence; the source mask
count is retained. The record order and length order are the declared FASTA/config
orders. Exhaustion of the100,000-candidate bound is an input qualification failure,
not an aligner miss. No outcome-dependent reseeding or replacement is allowed.

Each reference/length has500 forward and500 reverse loci, shuffled by its named
orientation stream. Reverse-complementing a query does not add another locus.
The five complete assemblies retain all reference records; human must contain only
NC_000022.11. Compressed input SHA256 values are pinned in `config/references.json`.
Decompression preserves all bytes. Human extraction copies only its complete
matching FASTA record, without rewrapping or renaming.

## Random streams and exact replay

Master seed2026111201. Python3.11 is the reference interpreter; its exact version is
recorded in each manifest. Stream seed payload is UTF-8:

`json.dumps([master_seed, stage, *labels], ensure_ascii=False, separators=(',', ':'))`

Hash with SHA256; use the first8 bytes as a big-endian integer for `random.Random`.
Store the complete digest. The implementation's ordered stages/labels are:

- `source`, `final`, reference_id, length, locus_ordinal (1…1000).
- `orientation`, reference_id, length.
- `mutation`, source_id, source_length, model, nominal_identity,1.
- `background_selection`, reference_id, length.
- `background`, reference_id, length, four-source_group_ordinal (0-based).

The final1 denotes one mutation version. It is not a1000-element replicate axis.
Mutation RNG states are independent of donor sampling. Background blocks are reused
across model/identity variants within a group; the target bytes are identical to
their short counterparts. Query groups are selected once and never by tool results.
`statistical_seed=2026111202` is reserved for future explicitly defined resampling;
the default summaries do not fabricate confidence intervals.

## U, C and I generation

Let n be source length and p=nominal_identity/100. Python `round` uses ties-to-even,
and the expression order is contractual. For U/C, X=round(n−p*n), not n*(1−p).
At n150 and identity85 the former gives22 substitutions. Select X distinct source
positions; each selected base is replaced uniformly by one of the other three.
U draws positions uniformly without replacement.

C always uses h50/a80/w10: H=round(n*50/100), W=max(1,round(n*10/100)),
X_hot=round(80*X/100). Draw contiguous candidate windows until their union contains
exactly H positions. Each next window has width min(W,H−current_union_size); its
start is uniform among fitting starts, and windows may overlap. Stop at the exact
union size or raise at10,000 attempts. Uniformly select X_hot hotspot positions and
X−X_hot cold positions without replacement. Save every window, actual hotspot
position, substitution position, and integer count. h describes union coverage;
w is the candidate-window scale, not a guarantee of disjoint equal-width hotspots.
At100% identity no edit is made; the three models alias one P100 record.

I preserves the previous controlled edit rule: total indel bases
G=round(n*min(.04,.25*(1−p))), insertion bases I=floor(G/2), deletion bases D=G−I,
and substitutions X=round(n−D−p*(n+I)). Delete non-overlapping runs and place inserted
runs at distinct anchors0…n. Run lengths follow a geometric stopping process with
stop probability0.5, capped at20 and truncated to the remaining budget. Substitutions
are uniform among undeleted positions; inserted bases are uniform A/C/G/T. The
10,000-attempt placement limits are explicit. These are synthetic edit processes,
not fitted evolutionary histories.

Truth stores forward-reference aligned strings `qaln/taln`, M/X/I/D counts, actual
query length, and realized identity M/(M+X+I+D). Nominal identity is discretized by
integer budgets and need not be exact. Query-strand changes are applied after
constructing the generating path; all reference coordinates remain forward,
0-based, half-open. Insertions consume query only, deletions reference only.

## Two query conditions

All isolated target sequences are individual FASTA records, possibly batched in
one invocation. They are never implicitly concatenated for the short task.

From each reference's150/1000bp cohort, a fixed random permutation selects200 loci,
then partitions them into50 consecutive groups of four. The same source groups
are used at60/80/94% under U/C/I. Draw five native Arabidopsis blocks, each2000bp,
and construct B0–q1–B1–q2–B2–q3–B3–q4–B4. Raw background is uppercased in queries;
original donor/reference bytes and masks are preserved and fingerprinted. Save
donor record/coordinates, group seed, and all embedded target positions.

Background length B is10,000bp, not total query length. S=B+sum(actual_target_length):
10,600bp for four150bp substitution targets and14,000bp for four1000bp targets.
Indels adjust S using actual query lengths. There are5400 samples and21,600 target
evaluations per tool. The short counterpart already exists in the primary bank.

Donor is not added to the search reference. Different species do not certify the
absence of homology. There is no new randomized-background condition. Consequently
the short/native difference measures the complete embedding task's effect; it
cannot isolate background composition, ordering, query length, or internal filtering.

## Native policies and tool boundaries

All policies are fixed across source identity and outcome. The exact argv are
generated in `native.POLICIES` and retained per invocation:

- minimap2: asm20, k11/w5, base-level PAF, A1/B1, N0, O4,8/E2,1, s20/m20/n2/p0.8.
- LAST: YASS reference index with lowercase masking; lastal j3/fMAF/m100.
- LASTZ: default policy, pairwise MAF, complete multi-record reference.
- BLASTN: task blastn, word7, reward1/penalty−1, gap2/1, dust/soft_masking enabled,
  both strands, evalue1e−6, max_target_seqs500. Aligned strings are exported; bounding
  boxes alone cannot recover paired positions.
- FastGA: named pinned compatibility build; k40 GIX, i0.7/f10/c85/s100/l100. This is
  a declared short-target policy, distinct from historical s1000 genome comparisons.
  Query GDB/index/readback, mapping and ALN-to-PAF conversion are separately timed.
- Original Bowtie1.3.1: FASTA input, SAM, v3, all best-stratum alignments. It is
  ungapped and accepts at most1024bp per read. Composite records exceeding that
  limit are NOT_SUPPORTED; no splitting or Bowtie2 fallback. Valid supported
  searches returning no alignments count as zero recovery.

minimap2/LAST/BLASTN/FastGA/Bowtie use the configured thread count; LASTZ is single
threaded. Native calls serialize with `native_slot=1` and a filesystem lock. Perfect
source controls must recover at every declared length before mutated jobs for that
tool/reference are admitted. Failure is a technical gate failure, not an opportunity
to remove difficult identity strata. Runtime version/hash inventories and newly
measured costs must not be confused with historical policy-selection results.

Authoritative tool references: [Bowtie manual](https://github.com/BenLangmead/bowtie/blob/v1.3.1/MANUAL.markdown),
[FastGA pinned source](https://github.com/thegenemyers/FASTGA/tree/178afcbba2f79abc1686e3bd2b8a916969467e83),
[minimap2 manual](https://lh3.github.io/minimap2/minimap2.html),
[LAST](https://gitlab.com/mcfrith/last), [LASTZ](https://lastz.github.io/lastz/),
[BLAST command-line manual](https://www.ncbi.nlm.nih.gov/books/NBK279690/).

## Outcome definitions and native format checks

Normalize every native alignment to forward reference coordinates and forward
whole-query coordinates, with an explicit relative strand and base-level operation
path. PAF requires `cg`; MAF must be pairwise; SAM honors reverse flags and clipping;
BLAST reconstructs operations from qseq/sseq, including reverse subject hits.
Check names, lengths, bounds, and path consumption. Do not infer paired coverage
from an alignment bounding box.

For each HSP independently count M/=/X positions (mismatches remain paired); gaps
are excluded. Discard HSPs with fewer than100 paired positions. For a target, require
the correct record and relative strand. Intersect every remaining paired segment
with BOTH the target's query interval and its reference-source interval. Take
unions separately on the two axes, then divide by full target-query length and
full source length. Region90 is1 iff both fractions are≥0.90. Multiple qualifying
HSPs can contribute; multiple sub-threshold HSPs cannot be combined to pass100bp.
Native output ordering and overlapping HSPs do not change the unions.

Strict generating-pair recall is a separate diagnostic: a query coordinate must
be paired to its exact generated reference coordinate. Its denominator is the
number of generated M/X pairs. It does not assume an indel representation is the
only biologically acceptable alignment. No natural precision, indel-event F1, or
equivalent-repeat-copy metric is claimed without external truth.

FAILED/TIMEOUT supported positive searches contribute0 and remain in denominators.
NOT_SUPPORTED contributes missing recovery and an explicit count. Missing inputs,
index failure, conversion failure and failed perfect controls prevent a complete
report rather than being silently replaced by zero. P100 controls are separate.

## Costs, provenance and reporting

GNU `/usr/bin/time -v` measures each native command's wall time, user+system CPU,
and maximum resident size in KiB. The Python wrapper also records elapsed wall
time, affinity, CPU model, load, timeout, RLIMIT_AS, output-file cap, argv, return
code, stdout/raw, stderr and content hashes. A separate single-threaded worker
sets resource limits before exec; no scheduler-thread `preexec_fn` is used.
Timeouts terminate the entire invocation process group, retaining partial output.
Missing GNU statistics remain null, never0 or the configured timeout.

Mapping, reference indexing, FastGA query preparation/readback, and conversion
are distinct. CPU counts and RAM scopes are not interchangeable; peak RSS is
GNU time's invocation statistic, not a sum of simultaneous processes. A hardware
snapshot describes its collection time, not historical cache/frequency conditions.

The shared budget admits each timeout reservation, charges measured wrapper elapsed
time after completion, and retains unsettled reservations on supervisor failure.
Reservations and worst-case bounds are never labeled measured runtime. Partial
attempts are immutable numbered directories. Do not delete Snakemake metadata to
force apparent success or overwrite an interrupted native attempt.

Outputs give per-target data, per-reference/length/identity/model/context strata,
equal-cell macro recovery excluding P100, and paired background gain/loss. Pairing
uses identical fragment IDs, not unrelated cohort averages. Source derivatives
and four-source composite groups share dependence. Default summaries are
descriptive; no independent-query binomial intervals or natural-family independence
claims are manufactured from the number of generated queries.
