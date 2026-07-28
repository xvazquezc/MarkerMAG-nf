process RENAME_READS {
    tag   "${meta.id}"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/renamed_reads" }, mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_R1.*"), path("${meta.id}_R2.*"), emit: reads

    script:
    def fq_flag = (r1.name =~ /\.f(ast)?q(\.gz)?$/) ? '-fq' : ''
    """
    MarkerMAG rename_reads \\
        -r1 ${r1}        \\
        -r2 ${r2}        \\
        -p  ${meta.id}   \\
        -t  ${task.cpus} \\
        ${fq_flag}
    """
}
