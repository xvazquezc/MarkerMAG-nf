# MATAM performance and graph-parallelisation plan

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

## Deferred plan: map once, subsample afterwards

This improvement is intentionally not part of the immediate graph work, but is
recorded for later implementation.

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
11. Revisit the deferred single-mapping design after graph sharding is stable.

