// UCLUST_16S — clusters 16S sequences with USEARCH cluster_fast.
// Used after COMBINE_MATAM (MATAM path) and optionally on user-provided 16S.
// Emits the representative (longest per cluster) fasta plus the .uc table.

process UCLUST_16S {
    tag   "${meta.id}"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/uclust_16s" }, mode: 'copy'

    input:
    tuple val(meta), path(seqs_in)
    val   iden                          // clustering identity (0–1)

    output:
    tuple val(meta), path("${meta.id}_16S_clustered.fasta"),       emit: clustered
    tuple val(meta), path("${meta.id}_16S_clustered.uc"),          emit: uc
    tuple val(meta), path("${meta.id}_16S_clustered.uc.reorganised.txt"), emit: membership

    script:
    """
    MarkerMAG uclust_16s \\
        -in  ${seqs_in}                       \\
        -i   ${iden}                           \\
        -out ${meta.id}_16S_clustered.fasta
    """
}
