// POLISH_16S — runs Barrnap on assembled/clustered 16S sequences and removes
// any sequences that Barrnap does not identify as 16S rRNA.

process POLISH_16S {
    tag   "${meta.id}"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/polish_16s" }, mode: 'copy'

    input:
    tuple val(meta), path(seqs_in)

    output:
    tuple val(meta), path("${meta.id}_16S_polished.fasta"), emit: polished

    script:
    """
    MarkerMAG polish_16s \\
        -in  ${seqs_in}                      \\
        -out ${meta.id}_16S_polished.fasta
    """
}
