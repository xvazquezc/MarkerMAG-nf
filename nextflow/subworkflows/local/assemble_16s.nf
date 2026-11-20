// ASSEMBLE_16S subworkflow
//
// Accepts either raw paired reads or already extracted paired 16S reads. Raw
// reads first pass through MATAM --filter_only. Both routes then scatter over all
// configured subsample percentages — each pct runs as an independent HPC job
// via MATAM_ASSEMBLE — before collecting, combining, clustering, and polishing.
//
// With params.matam_map_once (default) each sample's 16S reads are mapped
// once: raw-read samples reuse the alignments of MATAM_FILTER, and extracted
// reads are mapped by MATAM_MAP_16S. Every percentage takes its read mapping
// from those alignments instead of mapping its subsample again.
//
// Inputs:
//   ch_reads            : [meta, r1, r2] for MATAM extraction
//   ch_extracted_pairs  : [meta, reads_16s_r1, reads_16s_r2]
//   ch_matam_db         : value channel — absolute MATAM DB prefix
//
// Emits:
//   seqs_16s    : [meta, polished_fasta]

include { MATAM_FILTER   } from '../../modules/local/matam_filter/main'
include { MATAM_MAP_16S  } from '../../modules/local/matam_map_16s/main'
include { MATAM_ASSEMBLE } from '../../modules/local/matam_assemble/main'
include { COMBINE_MATAM  } from '../../modules/local/combine_matam/main'
include { UCLUST_16S     } from '../../modules/local/uclust_16s/main'
include { POLISH_16S     } from '../../modules/local/polish_16s/main'

workflow ASSEMBLE_16S {

    take:
    ch_reads            // [meta, r1, r2]
    ch_extracted_pairs  // [meta, reads_16s_r1, reads_16s_r2]
    ch_matam_db         // value: absolute MATAM DB prefix

    main:
    def pct_list = params.matam_pcts
        .tokenize(',')
        .collect { it.trim() }   // keep as strings — Nextflow filename glob must match
                                  // what markermag_matam_assemble_pct.py produces

    //
    // Extract 16S reads once for raw-read samples, then merge the result with
    // samples that entered at the extracted-reads checkpoint.
    //
    MATAM_FILTER ( ch_reads, ch_matam_db )
    ch_filtered_pairs = MATAM_FILTER.out.reads_16s
    ch_checkpoint_pairs = ch_extracted_pairs

    // From the command line, `--matam_map_once false` arrives as the string
    // "false", which Groovy treats as true.
    def map_once = params.matam_map_once.toString().toLowerCase()
    if ( !(map_once in ['true', 'false']) ) {
        error "ERROR: --matam_map_once must be true or false, not '${params.matam_map_once}'"
    }
    if ( map_once == 'true' ) {
        // Subsamples keep being drawn from the checkpoint reads themselves
        // (as without map-once); MATAM_MAP_16S only provides their alignments.
        MATAM_MAP_16S ( ch_checkpoint_pairs, ch_matam_db )
        ch_reads_16s = ch_filtered_pairs
            .join(MATAM_FILTER.out.alignments, failOnDuplicate: true, failOnMismatch: true)
            .mix(ch_checkpoint_pairs
                .join(MATAM_MAP_16S.out.alignments, failOnDuplicate: true, failOnMismatch: true))
    } else {
        ch_reads_16s = ch_filtered_pairs
            .mix(ch_checkpoint_pairs)
            .map { meta, forward, reverse -> tuple(meta, forward, reverse, []) }
    }

    //
    // Scatter: one MATAM_ASSEMBLE job per (sample × pct) combination.
    // Each job receives its own process_matam resource allocation,
    // running as an independent HPC submission.
    //
    ch_pct_inputs = ch_reads_16s
        .flatMap { meta, forward, reverse, alignments ->
            pct_list.collect { pct -> tuple(meta, pct, forward, reverse, alignments) }
        }

    MATAM_ASSEMBLE ( ch_pct_inputs, ch_matam_db )

    //
    // Collect all pct assemblies per sample, then concatenate
    //
    MATAM_ASSEMBLE.out.assembly
        .groupTuple()
        | COMBINE_MATAM

    //
    // Cluster combined assemblies then polish with Barrnap
    //
    UCLUST_16S ( COMBINE_MATAM.out.combined, params.cluster_iden )

    POLISH_16S ( UCLUST_16S.out.clustered )

    emit:
    seqs_16s = POLISH_16S.out.polished   // [meta, polished_fasta]
}
