// MARKERMAG_CORE subworkflow
//
// Runs the main linkage analysis (LINK_16S) and optionally the standalone
// copy-number estimation (GET_CP_NUM).
//
// Inputs:
//   ch_samples : [meta, r1, r2, mag_dir, mag_ext, marker_16s]
//
// Emits:
//   linkages_gnm : [meta, linkages_by_genome.txt]
//   cp_num_16s   : [meta, copy_num_by_16S.txt]   (may be empty if skipped)
//   cp_num_mag   : [meta, copy_num_by_MAG.txt]   (may be empty if skipped)

include { LINK_16S    } from '../../modules/local/link_16s/main'
include { GET_CP_NUM  } from '../../modules/local/get_cp_num/main'

workflow MARKERMAG_CORE {

    take:
    ch_samples  // [meta, r1, r2, mag_dir, mag_ext, marker_16s]

    main:
    //
    // Core linking — polishing and clustering already done upstream,
    // so -no_polish and -no_cluster are always passed from within LINK_16S.
    //
    LINK_16S ( ch_samples )

    //
    // Optional standalone copy-number estimation (useful for re-running
    // with different --subsample_pct without re-running the full link step).
    //
    ch_cp_num_16s = Channel.empty()
    ch_cp_num_mag = Channel.empty()

    if ( !params.skip_cp_num ) {
        ch_get_cp_inputs = ch_samples
            .join( LINK_16S.out.linkages_gnm, by: 0 )
            .map { meta, r1, r2, mag_dir, mag_ext, marker, linkages ->
                tuple(meta, r1, r2, mag_dir, mag_ext, marker, linkages)
            }

        GET_CP_NUM ( ch_get_cp_inputs )

        ch_cp_num_16s = GET_CP_NUM.out.cp_num_16s
        ch_cp_num_mag = GET_CP_NUM.out.cp_num_mag
    }

    emit:
    linkages_gnm = LINK_16S.out.linkages_gnm   // [meta, linkages_by_genome.txt]
    cp_num_16s   = ch_cp_num_16s
    cp_num_mag   = ch_cp_num_mag
}
