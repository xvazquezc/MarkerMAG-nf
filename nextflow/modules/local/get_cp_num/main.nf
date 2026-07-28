// GET_CP_NUM — standalone copy-number estimation.
//
// Can be run independently of LINK_16S (e.g. with different --subsample_pct)
// or when -skip_cn was passed to LINK_16S.

process GET_CP_NUM {
    tag   "${meta.id}"
    label 'process_high'

    publishDir { "${params.outdir}/${meta.id}/get_cp_num" }, mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2), path(mag_dir), val(mag_ext),
          path(marker_16s), path(linkages_gnm)

    output:
    tuple val(meta), path("${meta.id}_get_16S_cp_num_wd/${meta.id}_copy_num_by_16S.txt"),
          optional: true, emit: cp_num_16s
    tuple val(meta), path("${meta.id}_get_16S_cp_num_wd/${meta.id}_copy_num_by_MAG.txt"),
          optional: true, emit: cp_num_mag

    script:
    """
    MarkerMAG get_cp_num             \\
        -p           ${meta.id}      \\
        -r1          ${r1}           \\
        -r2          ${r2}           \\
        -marker      ${marker_16s}   \\
        -mag         ${mag_dir}      \\
        -x           ${mag_ext}      \\
        -linkages    ${linkages_gnm} \\
        -t           ${task.cpus}    \\
        -subsample_pct ${params.subsample_pct}
    """
}
