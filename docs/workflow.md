# Workflow operation

All commands are run from the repository root. `mosaic plan` is side-effect-free;
`prepare` produces project-owned inputs; `run --dry-run` constructs the whole DAG;
`run` submits an independent systemd user service; `run --foreground` is suitable
inside CI or a scheduler allocation. There is one `workflow/Snakefile`.

## Stages

1. Explicitly fetch checksum-pinned public materials, or configure already registered
   complete FASTA files. Human is one complete chromosome; mouse/zebrafish are full.
2. Generate distinct loci, U/C/I variants, exact path truth, canonical perfect
   controls, and native-background layouts. A complete manifest is written last.
3. Inventory executable/helper hashes, version text, full policy and hardware.
4. Build indexes per reference/tool, measuring each invocation. LASTZ has no separate
   reusable index. FastGA performs GDB/GIX construction and full sequence readback.
5. Map/evaluate P100 controls and enforce a gate for each tool/reference.
6. Run mutated short and intact composite query batches; mark original Bowtie's
   >1024bp composites as unsupported without launching invalid searches.
7. Normalize native base-level paths and evaluate both coverage axes.
8. Produce TSV/JSON summaries; optional plotting is offline and does not rerun tools.

The configured short batch64 and long batch1 are defaults requiring resource
qualification on the actual workstation. The inventory alone is not a perfect
control or a six-reference resource certification. The production72h cap is a hard
deployment limit, not an estimated runtime or authorization to exceed an existing
external budget. Failure to admit the next reservation stops the workflow.

## Files and continuation

`mosaic logs --config ...` prints absolute locations. The default production log is
`work/benchmark/execution.log` under the checkout; service/submission/exit JSONs
sit next to it. The exit record includes process return code; only a complete
summary plus corresponding input/control audits establishes matrix completion.
It does not automatically establish scientific acceptance.

`mosaic resume --config ...` runs the same DAG with `--rerun-incomplete`. Successful
outputs are reused. Every actual rerun gets another `attempt_NNNN`; no failed raw
file is overwritten. Configuration, package source and input fingerprints must
still match. Scientific changes require another explicit work_dir and protocol.
Never rename an old fixed-source replicate dataset as the new distinct-locus set.

If a supervisor is interrupted during a command, `budget.json` retains its
reservation. Do not clear it blindly. First confirm the recorded service/invocation
has stopped, inspect the attempt's command/time/raw/receipt, and reconcile the
ledger using actual evidence. If elapsed time cannot be reconstructed, keep it
unknown and account a clearly labeled conservative bound in a reviewed new budget;
do not invent a GNU runtime or resume an unaccounted invocation automatically.

Background submission uses a stable run-directory identity and checks the prior
unit state. A supervisor lock also rejects concurrent writers. Foreground mode
retains the same logs and exit records. A machine without a working user systemd
bus must use explicit foreground/scheduler operation, not a session-dependent
asynchronous command masquerading as an independent service.

## Environment and release

Use `environment.yaml` or `pip install '.[workflow]'` plus six installed toolchains.
FastGA's build script checks out its pinned commit and applies the named patch;
do not reapply it to a modified tree. Preserve the build directory/patch hash and
runtime `toolchain.json`. Native tools are external dependencies, not vendored
genome/binary assets in this Git repository. An exact rerun additionally needs the
same solver resolution/builds; save a conda explicit export or container digest.

Public repository outputs intentionally exclude local absolute reference paths,
credentials, multi-GB genomes, historical raw data, and workstation-specific caches.
Production accessions/checksums and relative output paths are sufficient to fetch
the declared public material and reconstruct the controlled test inputs.
