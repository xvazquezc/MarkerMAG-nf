// COMBINE_MATAM — concatenates all per-pct MATAM assembly fastas into one file.
// Downstream UCLUST_16S then clusters the combined sequences.

process COMBINE_MATAM {
    tag   "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(assemblies)   // list of per-pct scaffold fastas

    output:
    tuple val(meta), path("${meta.id}_assembled_16S_combined.fasta"), emit: combined

    script:
    """
    cat ${assemblies} > ${meta.id}_assembled_16S_combined.fasta
    """
}
