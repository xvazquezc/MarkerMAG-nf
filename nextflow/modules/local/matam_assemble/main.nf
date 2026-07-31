// MATAM_ASSEMBLE — runs one MATAM assembly for a single subsample percentage.
//
// Each percentage from params.matam_pcts is dispatched as an independent
// Nextflow process (separate HPC job), directly solving the CPU-efficiency
// problem of the original sequential loop in matam_16s.py.
//
// The wrapper always emits one per-depth FASTA. It is empty when the selected
// subset has no reads or MATAM succeeds without producing scaffolds.

process MATAM_ASSEMBLE {
    tag   "${meta.id}:pct${pct}"
    label 'process_matam'

    input:
    tuple val(meta), val(pct), path(reads_16s)
    val   matam_db              // absolute MATAM DB prefix on shared storage

    output:
    tuple val(meta), path("${meta.id}_subsample_*_scaffolds.fasta"),
          emit: assembly

    script:
    """
    markermag_matam_assemble_pct.py \\
        --reads_16s ${reads_16s}                    \\
        --pct       ${pct}                           \\
        --prefix    ${meta.id}                       \\
        --matam_db  '${matam_db}'                      \\
        --cpu       ${task.cpus}                     \\
        --mem_mb    ${params.matam_mem_mb}
    """
}
