// ASSEMBLE_16S subworkflow
//
// Accepts either raw paired reads or already extracted 16S reads. Raw reads
// first pass through MATAM --filter_only. Both routes then scatter over all
// configured subsample percentages — each pct runs as an independent HPC job
// via MATAM_ASSEMBLE — before collecting, combining, clustering, and polishing.
//
// Inputs:
//   ch_reads            : [meta, r1, r2] for MATAM extraction
//   ch_extracted_reads  : [meta, reads_16s] that bypass extraction
//   ch_db_dir           : value channel — path to parent dir of MATAM DB files
//   ch_db_name          : value channel — basename of the MATAM DB prefix
//
// Emits:
//   seqs_16s    : [meta, polished_fasta]

include { MATAM_FILTER   } from '../../modules/local/matam_filter/main'
include { MATAM_ASSEMBLE } from '../../modules/local/matam_assemble/main'
include { COMBINE_MATAM  } from '../../modules/local/combine_matam/main'
include { UCLUST_16S     } from '../../modules/local/uclust_16s/main'
include { POLISH_16S     } from '../../modules/local/polish_16s/main'

workflow ASSEMBLE_16S {

    take:
    ch_reads            // [meta, r1, r2]
    ch_extracted_reads  // [meta, reads_16s]
    ch_db_dir           // value: path to DB parent directory
    ch_db_name          // value: DB basename

    main:
    def pct_list = params.matam_pcts
        .tokenize(',')
        .collect { it.trim() }   // keep as strings — Nextflow filename glob must match
                                  // what markermag_matam_assemble_pct.py produces

    //
    // Extract 16S reads once for raw-read samples, then merge the result with
    // samples that entered at the extracted-reads checkpoint.
    //
    MATAM_FILTER ( ch_reads, ch_db_dir, ch_db_name )
    ch_reads_16s = MATAM_FILTER.out.reads_16s.mix( ch_extracted_reads )

    //
    // Scatter: one MATAM_ASSEMBLE job per (sample × pct) combination.
    // Each job receives its own process_matam resource allocation,
    // running as an independent HPC submission.
    //
    ch_pct_inputs = ch_reads_16s
        .flatMap { meta, reads ->
            pct_list.collect { pct -> tuple(meta, pct, reads) }
        }

    MATAM_ASSEMBLE ( ch_pct_inputs, ch_db_dir, ch_db_name )

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
