// MATAM_FILTER — runs matam_assembly.py --filter_only once per sample.
// Merges R1/R2, extracts 16S reads, emits them for scatter over pct values.
//
// Uses the same process_matam label as MATAM_ASSEMBLE so the HPC job
// requests exactly the memory that MATAM is told it may use (matam_mem_mb).

process MATAM_FILTER {
    tag   "${meta.id}"
    label 'process_matam'

    input:
    tuple val(meta), path(r1), path(r2)
    val   matam_db              // absolute MATAM DB prefix on shared storage

    output:
    tuple val(meta), path("${meta.id}_16S_reads*"), emit: reads_16s

    script:
    """
    markermag_matam_filter.py \\
        --r1       ${r1}                           \\
        --r2       ${r2}                           \\
        --prefix   ${meta.id}                      \\
        --matam_db '${matam_db}'                    \\
        --cpu      ${task.cpus}                    \\
        --mem_mb   ${params.matam_mem_mb}
    """
}
