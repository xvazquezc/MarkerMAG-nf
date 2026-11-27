// MATAM_FILTER — runs matam_assembly.py --filter_only once per sample.
// Keeps R1/R2 separate, extracts 16S reads, and emits a paired-file tuple,
// plus MATAM's alignments of those reads: every subsample assembly takes its
// read mapping from them instead of mapping again (params.matam_map_once).
//
// It maps every raw read, so it has its own label (process_matam_filter:
// params.matam_filter_threads CPUs, 48 h); the memory is that of the other
// MATAM jobs (matam_mem_mb plus 4 GB per attempt).

process MATAM_FILTER {
    tag   "${meta.id}"
    label 'process_matam_filter'

    input:
    tuple val(meta), path(r1), path(r2)
    val   matam_db              // absolute MATAM DB prefix on shared storage

    output:
    tuple val(meta), path("${meta.id}_16S_R1*"), path("${meta.id}_16S_R2*"),
          emit: reads_16s
    tuple val(meta), path("${meta.id}_16S.sam*"),
          emit: alignments

    script:
    """
    markermag_matam_filter.py \\
        --r1       ${r1}                           \\
        --r2       ${r2}                           \\
        --prefix   ${meta.id}                      \\
        --matam_db '${matam_db}'                    \\
        --cpu      ${task.cpus}                    \\
        --mem_mb   ${params.matam_mem_mb}           \\
        --matam    '${params.matam_executable}'
    """
}
