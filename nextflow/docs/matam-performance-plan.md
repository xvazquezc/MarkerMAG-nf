# MATAM performance and graph-parallelisation plan

## Current status (2026-11-11)

Most of the MATAM side of this plan is now implemented in the fork
(`xvazquezc/matam`). The sections below are the original plan, kept for
reference; where they say MATAM 1.6.0 or a serial graph builder, they describe
the starting point.

- **ovgraphbuild is multi-threaded.** `matam_assembly.py --cpu` is passed as
  `ovgraphbuild --threads` (at most 24 workers), with byte-identical output.
  Speed-up is ~9x at 24 threads.
- **Exact pair sharding exists.** `ovgraphbuild --pair-shard-index/-count`
  gives pair-count-balanced, contiguous ranges, and merged shards rebuild the
  serial edge order.
- **There is a stop point after alignment filtering.**
  `matam_assembly.py --stop_after alignments_filtering` keeps what the next
  steps read, and `--resume_from` continues the run. MarkerMAG-nf already
  uses `--stop_after scaffolding`.
- **Production-size benchmark (synthetic).** Simulated reads were mapped
  against a SILVA 144 subset, then the graph was built with 24 threads:
  - 169,624 reads (14.4 G candidate pairs): 10 min;
  - 180,320 reads in a skewed community (16.3 G pairs): 17 min;
  - serial estimates 1.5 h and 2.7 h.

  These data are 40-65x cheaper per pair than the production run (~100 h
  serial for 170k reads). Full numbers are in the fork's `PLAN.md`.
- **Decision for now: no scatter/merge orchestration in MarkerMAG-nf.**
  - With the ~9x thread speed-up, a 100 h serial graph takes ~11 h in one
    24-CPU task, within a 48 h walltime. Give MATAM_ASSEMBLE 24 CPUs
    (`--matam_threads 24`) for large samples.
  - Splitting the graph into separate jobs is only worth it for graphs still
    over the walltime at 24 threads.
- **Candidate-pair pruning is implemented in MATAM (2026-11-12).**
  - `ovgraphbuild` now compares each read only with reads on a shared
    reference: 3.5-5.6% of the pairs on the synthetic sets.
  - Output and statistics are byte-identical. `--all-pairs` keeps the old
    loop.
  - The gain is smaller than the pair count suggests. At 24 threads it is
    1.28-1.40× wall time (1.35-1.59× CPU), and 1.4-2.3× on single-thread
    shards. The skipped pairs were only 25-40% of the CPU time; most of the
    rest is the de novo re-alignment of compatible pairs.
  - Details are in the fork's `PLAN.md`.
- **Shards and thread chunks are balanced by expected work (2026-11-20,
  MATAM fork `4f54d4e`).**
  - Each read is weighted by its comparisons, not its pair count. Output is
    unchanged.
  - At 24 threads: 1.10-1.32× wall time.
  - For 16 single-thread shards, the slowest took 1.45× the mean CPU time
    (was 2.24-3.04×). This matters most for scatter/merge: a job finishes
    with its slowest shard.
  - What is left follows the compatible pairs, whose de novo re-alignment is
    about 75% of the CPU time.
- **A one-off 47 min stall is fixed.** SeqAn's daily update check (`wget`)
  is no longer compiled into `ovgraphbuild`. It kept the graph step's output
  open, and on HPC nodes without internet access it could stall a task the
  same way.
- **Map once, subsample afterwards is implemented (2026-11-20, MATAM fork
  `e4260a8`, here `d70d5f2`).** Each sample's 16S reads are mapped once. Every depth takes its read
  mapping from those alignments (`--read_alignments`), and
  `--matam_map_once false` restores the old way. The assemblies differ less
  than they do between two MATAM seeds; see the section below.
- **The de novo re-alignment is skipped for perfect overlaps (2026-11-20,
  MATAM fork `52ed46c`).** When a vectorized score-only pass proves that SeqAn's
  alignment can only be the ungapped all-match overlap (97-99% of the
  re-aligned pairs), its statistics are filled in directly; everything else
  still goes through SeqAn. Output is byte-identical (checked pair by pair
  on both synthetic sets). At 24 threads the graph now takes 36 s instead
  of 304 s (mixed) and 69 s instead of 480 s (skewed): 6.9-8.4× on top of
  the changes above. Details are in the fork's `PLAN.md`.
- **Production benchmark done (2026-11-23): no scatter/merge needed.** The
  filtered SAMs of the `ferrari/dnasip` run on Gadi (MATAM 1.6.0, whose
  50-100% subsamples hit the 24 h walltime) were rebuilt with the current
  `ovgraphbuild`, outputs byte-identical to 1.6.0.
  - The slow graphs had 0.68-1.97 M reads, not 170k, at a normal cost per
    pair. Heavy 100%: 22.2 h serial on Gadi, now 16 min at 16 threads on a
    workstation. Bulk 50%: 37.3 h, now 23 min. Bulk 100% (timed out): 1 h
    38 min.
  - Estimated on Gadi with 8 CPUs: up to 4-5 h for the graph of Bulk 100%.
  - **What now limits MATAM_ASSEMBLE:** `componentsearch` uses about 50
    bytes per edge (14.6 GB for Heavy 100%'s 294 M edges). Bulk 100%
    (1.417 G edges) needs about 70 GB, more than the 33-37 GB requested.
    SortMeRNA mapping (7.5 h for Bulk 100% at 8 CPUs) is the other large
    step; map once runs it once per sample.
  - Details are in the fork's `PLAN.md`.
- **`componentsearch` rewritten (2026-11-24, MATAM fork `9fc3aa4`).** Same
  output as upstream for a given seed (byte-identical on Heavy 100% and
  Bulk 50%), at 8 bytes per edge instead of about 50. Heavy 100%: 2 min,
  2.3 GB (was 9 min, 14.8 GB). Bulk 100%: 5.5 min, 11.1 GB (was about
  70 GB). MATAM_ASSEMBLE's memory request no longer needs to cover
  `componentsearch` beyond the graph-building step's own needs.
- **`componentsearch` edge filter fixed (2026-11-24, fork `c0a73f3`).**
  Upstream's edge filter was given `-N` instead of `-E` and, through a
  counting bug, removed almost nothing (none of the weight-1 edges on the
  production graphs). It now removes meta edges of at most `-E` read links,
  from both sides.
- **MATAM default `-E 1` again (2026-11-25).** `-E 1` was the authors'
  default (v0.9.3), applied by the Java ComponentSearch until v1.0.0. On a
  synthetic community and on the Heavy production pairs (SILVA 144, 5 seeds),
  `-E 1` gave fewer chimeric scaffolds in every seed (synthetic 7.2 → 5.4,
  Heavy 45.6 → 42.6) with unchanged recovery and coverage. MarkerMAG-nf
  uses MATAM's defaults, so its assemblies change slightly from 1.6.0 runs.
  MATAM's `--min_overlap_edge 0` reproduces them, but the pipeline does not
  expose that option (`markermag_matam_assemble_pct.py` builds the MATAM
  command). Details are in the fork's `PLAN.md` (item 6).
- **On/off options read as booleans (2026-11-25).** Nextflow 26.10 passes
  command-line values as strings, and Groovy treats `"false"` as true:
  `--skip_matam false` skipped MATAM, `--reconstruct_only false` ran in
  reconstruct-only mode, and `--skip_cp_num false` skipped copy numbers.
  Same for a quoted `"false"` in a params file. `lib/Flags.groovy` now reads
  `help`, `reconstruct_only`, `matam_map_once`, `skip_matam` and
  `skip_cp_num` (`true`/`false` in any case, a bare flag is `true`), and
  `validate_params` stops on any other value. `matam_map_once` had its own
  check; it uses the helper now. Checked with the test profile before and
  after: `--reconstruct_only false` without reads (ran UCLUST before, now
  the reads error), `--skip_matam false` without 16S or database (16S
  error before, now the `--matam_db` error), `--skip_matam maybe` (now
  rejected). `nextflow lint`: no new warnings.
- **MATAM memory request 34 → 24 GB (2026-11-25).** `matam_mem_mb`
  30000 → 20480, so MATAM jobs request 24 / 28 / 32 GB over the three
  attempts (`matam_mem_mb` + 4 GB × attempt), within the 32 GB that 8 CPUs
  carry on Gadi's normal queue. The first attempt was 33.3 GB, just above
  it. Measured peaks: scaffolding against the complete SILVA 144 index
  15.5 GB (the largest; SILVA 138.2 is smaller), SortMeRNA read mapping
  ~8 GB, ovgraphbuild 11.4 GB and componentsearch 11.1 GB on Bulk 100%.
  16 GB would leave under 0.5 GB for scaffolding with SILVA 144, and a
  `std::bad_alloc` exits 1, which `errorStrategy` does not retry.

## Status and scope

This document records the implementation plan for addressing long-running
`ASSEMBLE_16S:MATAM_ASSEMBLE` tasks in MarkerMAG-nf.

The immediate priority is exact parallelisation of MATAM's compatibility-graph
construction. Reusing a single mapping across all subsampling percentages is a
promising later improvement, but is deliberately deferred.

The following will **not** be changed as part of this work:

- database clustering;
- the default MATAM subsampling percentages;
- MATAM's coverage threshold or its default value.

## Findings

MATAM 1.6.0 passes `--cpu` to SortMeRNA, GNU `sort`, and component assembly, so
those stages can use multiple CPUs. Its `ovgraphbuild` compatibility-graph
stage does not use the supplied CPU count. The implementation contains a serial
nested loop over all mapped-read pairs and exposes no thread or pair-range
option.

For `n` mapped reads, the number of candidate pairs is:

```text
n * (n - 1) / 2
```

For 170,175 reads this is 14,479,680,225 pairs. Consequently, graph runtime
grows approximately quadratically with the number of reads. Allocating more
CPUs to the current MATAM command will help alignment but will leave those CPUs
idle during graph construction.

Relevant upstream material:

- [MATAM 1.6.0 pipeline and resume points](https://github.com/bonsai-team/matam/blob/v1.6.0/scripts/matam_assembly.py)
- [Serial compatibility-graph loop](https://github.com/bonsai-team/matam/blob/v1.6.0/ovgraphbuild/src/compatibilityGraphBuilding.cpp)
- [Upstream report that `--cpu` does not parallelise graph construction](https://github.com/bonsai-team/matam/issues/118)

## Ownership boundary

The scalability limitation belongs to MATAM, specifically `ovgraphbuild`, rather
than to MarkerMAG or Nextflow. Exact parallelisation requires changing MATAM's
graph builder or supplying a compatible replacement.

MarkerMAG-nf still has an integration role:

1. split a MATAM run at the filtered-alignment/graph boundary;
2. scatter exact pair ranges as independent HPC jobs;
3. merge and validate graph shards;
4. resume MATAM at graph compaction;
5. expose configuration and stage-specific resources;
6. preserve provenance for the patched MATAM binary.

The preferred ownership model is a minimal MATAM 1.6.0 fork containing the
generic graph-sharding capability, with MarkerMAG-nf containing only workflow
orchestration and wrappers. The MATAM change should be suitable for an upstream
pull request.

## Immediate design: exact graph sharding

### Required MATAM changes

Create a small fork based on the pinned MATAM 1.6.0 tag and extend
`ovgraphbuild` with options such as:

```text
--pair-shard-index INDEX
--pair-shard-count COUNT
```

Each invocation must:

- read the complete filtered SAM file;
- construct the same ordered read/vertex table as an unsharded run;
- process a mutually exclusive subset of outer-loop read indices;
- compare every assigned `(i, j)` pair where `i < j`;
- write an edge CSV shard using the original global vertex identifiers;
- write or checksum the common node CSV;
- report the assigned and completed candidate-pair counts.

Shards must be balanced by **pair count**, not by assigning equal numbers of
outer-loop reads. The loop is triangular: early `i` values have many more
partners than late values. Boundaries can be computed from the cumulative pair
count:

```text
pairs_before(i) = i * (2*n - i - 1) / 2
```

Contiguous ranges of `i` should be used. Merging shards in range order then
reproduces the stock edge ordering and avoids introducing graph-algorithm
differences caused only by record order.

The current graph object is not used as a downstream in-memory result;
`componentsearch` consumes the emitted node and edge CSV files. This makes
file-level sharding practical.

### MATAM stage boundary

Graph sharding also requires a clean stopping point after alignment filtering.
The preferred MATAM-fork change is a generic option such as:

```text
--stop_after alignments_filtering
```

The pre-graph run must retain:

- the subsampled FASTQ;
- SortMeRNA's aligned FASTQ and SAM;
- the score-filtered SAM used by `ovgraphbuild`;
- MATAM's output/work directory and parameter manifest.

MATAM already supports `--resume_from graph_compaction`, which should be used
after the merged graph has been installed under the filenames expected by the
original pipeline.

If adding `--stop_after` upstream is undesirable, a MarkerMAG-nf wrapper can
execute the same mapping and filtering commands directly. That is less
desirable because it duplicates MATAM's internal filename and parameter logic.

### MarkerMAG-nf processes

Replace the monolithic per-percentage `MATAM_ASSEMBLE` execution with:

```text
MATAM_PRE_GRAPH
    -> MATAM_GRAPH_SHARD
    -> MATAM_MERGE_GRAPH
    -> MATAM_FINISH
```

#### `MATAM_PRE_GRAPH`

- Retains the current sample/percentage subsampling behaviour.
- Runs MATAM through alignment filtering.
- Emits the filtered SAM, subsampled reads, and resumable MATAM state.
- Uses a mapping-oriented resource label with multiple CPUs.
- Becomes independently cached by Nextflow.

#### `MATAM_GRAPH_SHARD`

- Scatters over `(sample, percentage, graph_shard)`.
- Runs the patched `ovgraphbuild` for one exact pair range.
- Requests one CPU initially because each shard remains serial.
- Reads the complete filtered SAM through a staged symlink.
- Emits an edge shard, node checksum, range metadata, and pair counts.

Using many single-CPU jobs avoids reserving idle CPUs and provides checkpointing
at shard granularity. A later internally threaded implementation can be added if
benchmarks show that fewer multi-CPU jobs are more efficient on the HPC system.

#### `MATAM_MERGE_GRAPH`

- Confirms that all expected shard indices are present exactly once.
- Confirms identical node-file checksums across shards.
- Confirms non-overlapping, complete pair ranges.
- Confirms that completed pair counts sum to `n*(n-1)/2`.
- Retains a single node CSV.
- Concatenates edge records in pair-range order under one header.
- Writes a manifest recording binary version, parameters, read count, pair
  count, shard count, checksums, and elapsed times.

#### `MATAM_FINISH`

- Restores or links the pre-graph artifacts and merged node/edge CSV under
  MATAM's expected work-directory names.
- Invokes `matam_assembly.py --resume_from graph_compaction` with the same
  parameters used before the split.
- Produces the same per-percentage scaffold FASTA expected by the existing
  `COMBINE_MATAM` process.

### Configuration

Add explicit parameters without changing existing biological defaults:

```groovy
params.matam_graph_shards      = 1
params.matam_graph_target_time = '12.h'
params.matam_graph_max_shards  = 128
```

The first implementation should require an explicit shard count while retaining
`1` as a compatibility default. After calibration, an `auto` mode can estimate
the shard count from the filtered unique-read count and measured pairs per
second:

```text
estimated_hours = pair_count / calibrated_pairs_per_second / 3600
shards = ceil(estimated_hours / target_hours)
```

The result should be clamped to `matam_graph_max_shards` and recorded in the
manifest. Automatic selection must use performance observations from Gadi, not
an assumed linear speed-up.

Introduce separate resource labels:

```text
process_matam_pre_graph
process_matam_graph
process_matam_merge_graph
process_matam_finish
```

This permits alignment to request multiple CPUs while graph shards request one.
Walltimes and retry behaviour can then be set per stage. A graph timeout will
rerun only the failed shard rather than repeating alignment and all completed
graph work.

### Packaging and provenance

Keep MATAM pinned to 1.6.0. Package the minimal fork or replacement
`ovgraphbuild` binary at a pinned commit and record that commit in:

- `nextflow/environment.yml`;
- task `versions.yml` or an equivalent version output;
- the merged graph manifest;
- the README and changelog.

Do not rely on an untracked binary installed manually on Gadi. The exact binary
must be reproducible through the Conda environment or another versioned build
mechanism.

Because MATAM is AGPL-licensed, confirm that the fork and distributed source
remain compliant with its licence.

## Expected runtime effects

For a graph currently taking approximately 100 hours:

| Shards | Ideal compute time per shard | Practical expectation |
|---:|---:|---|
| 4 | 25 h | Longer due to imbalance and repeated SAM/reference loading |
| 8 | 12.5 h | Likely suitable as an initial benchmark |
| 16 | 6.25 h | Better walltime, with more filesystem and scheduler load |
| 32 | 3.1 h | Useful only if shard overhead remains modest |

These values are idealised. Total CPU work is not reduced; it is distributed.
Every shard rereads the full SAM and reference, so excessive shard counts can
stress shared storage and reduce efficiency.

Graph sharding alone may still be impractical for percentage/sample combinations
whose serial projections are hundreds of hours. The first benchmarks must
measure throughput and establish a reasonable maximum graph size and shard
count, without changing the agreed biological defaults.

## Validation requirements

### Unit and algorithm tests

- Verify the pair-boundary function for small `n` and many shard counts.
- Ensure every pair occurs exactly once and no pair is omitted.
- Cover empty shards and `n < shard_count`.
- Confirm global vertex identifiers are identical across shards.
- Confirm malformed or incomplete shard sets fail before graph merge.

### Exact-equivalence tests

For the MATAM example data and a medium representative dataset:

1. run stock `ovgraphbuild` once;
2. run the sharded implementation with 1, 2, 4, and 8 shards;
3. require byte-identical node CSVs;
4. require byte-identical merged edge CSVs, or identical records after a
   documented canonical ordering;
5. compare graph-compaction outputs;
6. compare final scaffold sequences after normalising identifiers.

Any sequence difference must be treated as a correctness failure until its
cause is understood.

### Failure and resume tests

- Terminate one graph shard and verify only that shard is rerun with `-resume`.
- Remove one edge shard and verify merge fails clearly.
- Corrupt a node checksum and verify merge fails.
- Verify retries do not rerun `MATAM_PRE_GRAPH`.
- Verify a completed merged graph allows `MATAM_FINISH` to resume independently.

### HPC benchmarks

Record, per sample and percentage:

- filtered unique reads;
- theoretical pair count;
- shard boundaries and pair counts;
- SAM and reference sizes;
- graph walltime and CPU utilisation;
- peak memory;
- read bandwidth;
- queue wait time;
- merge time and output size.

Benchmark at least 1, 4, 8, 16, and 32 shards before selecting a recommended
Gadi value.

## Map once, subsample afterwards (implemented 2026-11-20, `d70d5f2`)

**What is done.**

- The MATAM fork has a new option: `matam_assembly.py --read_alignments SAM`.
  - It takes the read mapping from alignments made earlier, instead of running
    SortMeRNA.
  - It selects the pairs with at least one alignment, plus their alignments, in
    the input order. The rest of the run is unchanged.
  - The SAM's `@PG` line must show the same database, `--best` and `--evalue`.
  - Resuming works as after a mapping: there is a completion record.
  - The clustered database index is not needed.
- Every run now also publishes `marker_alignments.sam.gz`, beside
  `marker_reads_R1/R2`.
- In MarkerMAG-nf:
  - `MATAM_FILTER` also emits `<id>_16S.sam.gz`.
  - Samples that enter at the extracted-read checkpoint are mapped once by the
    new `MATAM_MAP_16S`. It writes only the alignments, so it never writes over
    a checkpoint file named like a published one.
  - `MATAM_ASSEMBLE` passes the alignments as `--read_alignments`.
  - Subsamples are drawn as before, so each percentage gets the same read
    pairs.
  - `--matam_map_once false` brings back a mapping per subsample.
- The standalone `MarkerMAG matam_16s` does the same; `-map_each_subsample`
  restores the old behaviour.

**It is not exactly the same as mapping each subsample.**

- SortMeRNA's minimum alignment score depends on the total length of all the
  reads it maps together. With `-e 1e-5`, the minimum was:
  - 65 for 1% of a 80,840-pair 16S set;
  - 68 for 10%;
  - 72 for the whole set;
  - 84 at about 1,000× that size, emulated with `-e 1e-8`, a scale that
    represents raw metagenome reads.
- `--num_alignments` keeps the first alignments that reach that score. A read
  near the threshold, or with more candidate references than `--best`, can
  therefore get other alignments, including good ones.
- After MATAM's alignment filtering, the alignments differed for:
  - 1.2–1.4% of a subsample's reads, when the mapping was of the 16S set;
  - 2.7%, when it was of a raw-read-sized set.
- SortMeRNA itself is deterministic: the same input gives identical output at
  1 or 24 threads.
- Map-once makes a read's alignments the same at every depth. With one
  mapping per subsample, they depend on which other reads were drawn.

**Effect on the assemblies.** This was measured on the mixed SILVA-subset
set (80,840 pairs of 39 known 16S sequences). Each depth was assembled to
scaffolding (`--seed 42`) four ways:

- **per-subsample:** one mapping per subsample, as before;
- **seed 7:** the same as per-subsample, with only MATAM's seed changed;
- **16S map-once:** taken from a mapping of the whole 16S set (as for
  checkpoint samples);
- **raw map-once:** taken from a mapping of the raw-read size, emulated as
  above (as for raw-read samples).

Scaffolds identical to the per-subsample assembly:

| Depth | Per-subsample scaffolds | Seed 7 | 16S map-once | Raw map-once |
|---:|---:|---:|---:|---:|
| 1% | 27 | 8 of 27 | 14 of 24 | 14 of 23 |
| 5% | 59 | 24 of 64 | 37 of 57 | 26 of 57 |
| 10% | 62 | 16 of 62 | 29 of 65 | 32 of 61 |
| 25% | 80 | 22 of 76 | 54 of 80 | 46 of 77 |
| 50% | 83 | 24 of 86 | 74 of 84 | 43 of 85 |

True sequences recovered (a scaffold at ≥ 99% identity over ≥ 95% of its
length):

| Depth | Per-subsample | Seed 7 | 16S map-once | Raw map-once |
|---:|---:|---:|---:|---:|
| 1% | 17 | 15 | 15 | 14 |
| 5% | 28 | 27 | 27 | 28 |
| 10% | 27 | 26 | 27 | 27 |
| 25% | 30 | 31 | 30 | 30 |
| 50% | 28 | 30 | 28 | 29 |

Changing only MATAM's seed changes the assemblies more than map-once does,
and the recovery of the true sequences stays within the seed-to-seed range.
The differences are therefore within MATAM's own variability.

**Cost.** The selection replaces SortMeRNA at each depth. It took 1.7–3.2 s,
against 1.4–17.8 s of mapping for 1–50% on this small database at 24
threads. With full SILVA, each depth would also have to load the clustered
index. Each depth's job also no longer needs that index.

**Original proposal (kept for reference).** In the implementation, the
selected SAM keeps the full mapping's record order, which is the input order
the mate tagging needs, so no re-sorting was required. The full mapping is
not split.

The current workflow maps each per-percentage read subset independently. With
the default percentages, this maps an aggregate 266% of the extracted-read set.
A future design can map the complete extracted-read set once per sample and then
derive percentage-specific FASTQ and SAM files from that mapping.

Proposed future stages:

```text
MATAM_MAP_FULL
    -> optional MATAM_MAP_SHARD scatter
    -> MATAM_MERGE_MAPPING
    -> MATAM_SUBSAMPLE_MAPPING per percentage
    -> alignment filtering and graph stages
```

For raw-read samples, the aligned FASTQ and SAM already created by
`MATAM_FILTER` should be retained and reused. For samples entering through the
`16s_reads_r1`/`16s_reads_r2` checkpoint, one full mapping should be introduced.

Percentage selection must keep mates together and retain all SAM alignments for
every selected read. The percentage-specific SAM must then be canonically
sorted before MATAM filtering. Before adoption, compare independently mapped
subsets with subsets derived from a full mapping to demonstrate equivalent
filtered alignments and final assemblies.

If full mapping itself exceeds the desired walltime, it can be split by read
group, mapped independently with SortMeRNA, and merged. This is safe only after
an explicit equivalence test and should avoid excessive database reloads and
shared-filesystem contention.

## Implementation order

1. Create and pin a minimal MATAM 1.6.0 fork.
2. Add a stop point after alignment filtering.
3. Add exact pair-range support to `ovgraphbuild`.
4. Add pair-partition unit tests in the MATAM fork.
5. Split MarkerMAG-nf into pre-graph, graph-shard, merge, and finish processes.
6. Add manifests, checksums, stage-specific resources, and failure handling.
7. Demonstrate exact equivalence on small and medium datasets.
8. Benchmark shard counts on Gadi.
9. Document the recommended HPC configuration.
10. Submit the generic MATAM changes upstream.
11. ~~Revisit the deferred single-mapping design~~: done (2026-11-20), see
    "Map once, subsample afterwards".

