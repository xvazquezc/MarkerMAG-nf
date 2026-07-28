// LINK_16S — core MarkerMAG linking step.
//
// Maps paired-end reads to both 16S sequences and MAGs, identifies
// read pairs that bridge them, validates linkages via mini-assembly
// (Round 2), and optionally estimates 16S copy number.
//
// -no_polish and -no_cluster are always passed because POLISH_16S and
// UCLUST_16S have already been run as separate upstream processes.

process LINK_16S {
    tag   "${meta.id}"
    label 'process_high'

    publishDir { "${params.outdir}/${meta.id}/link_16s" }, mode: 'copy',
               saveAs: { filename ->
                   filename.contains('_MarkerMAG_wd') ? filename : null
               }

    input:
    tuple val(meta), path(r1), path(r2), path(mag_dir), val(mag_ext), path(marker_16s)

    output:
    tuple val(meta), path("${meta.id}_MarkerMAG_wd/${meta.id}_linkages_by_genome.txt"),
          emit: linkages_gnm
    tuple val(meta), path("${meta.id}_MarkerMAG_wd/${meta.id}_linkages_by_contig.txt"),
          optional: true, emit: linkages_ctg
    tuple val(meta), path("${meta.id}_MarkerMAG_wd/${meta.id}_copy_num_by_16S.txt"),
          optional: true, emit: cp_num_16s
    tuple val(meta), path("${meta.id}_MarkerMAG_wd/${meta.id}_copy_num_by_MAG.txt"),
          optional: true, emit: cp_num_mag
    tuple val(meta), path("${meta.id}_MarkerMAG_wd/*.log"),
          emit: log

    script:
    def skip_cn_flag = params.skip_cp_num ? '-skip_cn' : ''
    """
    MarkerMAG link                   \\
        -p       ${meta.id}          \\
        -r1      ${r1}               \\
        -r2      ${r2}               \\
        -marker  ${marker_16s}       \\
        -mag     ${mag_dir}          \\
        -x       ${mag_ext}          \\
        -t       ${task.cpus}        \\
        -min_link    ${params.min_link}    \\
        -min_16s_len ${params.min_16s_len} \\
        -max_16s_div ${params.max_16s_div} \\
        -mismatch    ${params.mismatch}    \\
        -aln_len     ${params.aln_len}     \\
        -aln_pct     ${params.aln_pct}     \\
        -no_polish                   \\
        -no_cluster                  \\
        ${skip_cn_flag}
    """
}
