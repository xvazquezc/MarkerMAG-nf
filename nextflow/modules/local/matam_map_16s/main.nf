// MATAM_MAP_16S — maps already extracted 16S read pairs once per sample
// (matam_assembly.py --filter_only) and emits only MATAM's alignments of them.
//
// Samples entering at the 16s_reads_r1/r2 checkpoint keep their own reads for
// subsampling; every subsample assembly takes its read mapping from these
// alignments (params.matam_map_once). The selected reads are not copied, so
// nothing is written under the names a previous run published them with.

process MATAM_MAP_16S {
    tag   "${meta.id}"
    label 'process_matam'

    input:
    tuple val(meta), path(r1), path(r2)
    val   matam_db              // absolute MATAM DB prefix on shared storage

    output:
    tuple val(meta), path("${meta.id}_16S.sam*"), emit: alignments

    script:
    """
    markermag_matam_filter.py \\
        --r1       ${r1}                           \\
        --r2       ${r2}                           \\
        --prefix   ${meta.id}                      \\
        --matam_db '${matam_db}'                    \\
        --cpu      ${task.cpus}                    \\
        --mem_mb   ${params.matam_mem_mb}           \\
        --matam    '${params.matam_executable}'     \\
        --alignments_only
    """
}
