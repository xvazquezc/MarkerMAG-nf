#!/usr/bin/env nextflow
// DSL2 is the default in Nextflow >=22.03; the explicit enable line is not needed.

// ── Subworkflows ───────────────────────────────────────────────────────────────
include { ASSEMBLE_16S   } from './subworkflows/local/assemble_16s'
include { MARKERMAG_CORE } from './subworkflows/local/markermag_core'

// ── Standalone modules (used in direct paths) ──────────────────────────────────
include { RENAME_READS } from './modules/local/rename_reads/main'
include { UCLUST_16S   } from './modules/local/uclust_16s/main'
include { POLISH_16S   } from './modules/local/polish_16s/main'

// ──────────────────────────────────────────────────────────────────────────────
// Help
// ──────────────────────────────────────────────────────────────────────────────
def help_message() {
    log.info """
    ╔══════════════════════════════════════════════════════════════════════════╗
    ║                           M a r k e r M A G                            ║
    ║        Link MAGs to 16S rRNA marker genes via paired-end reads          ║
    ╚══════════════════════════════════════════════════════════════════════════╝

    Usage:
        nextflow run main.nf [options]

    ── Input / output ────────────────────────────────────────────────────────────
      --input          CSV samplesheet  [required]
                       Columns: sample, r1, r2, mag_dir, mag_ext,
                                16s_reads, 16s_reads_r1, 16s_reads_r2,
                                16s_fasta
                       16s_fasta: use completed 16S sequences.
                       16s_reads:  single-end extracted candidate reads.
                       16s_reads_r1/r2: paired extracted candidate reads.
                       Leave all extracted-read fields empty to run extraction.
      --outdir         Output directory  [default: ${params.outdir}]
      --reconstruct_only
                       Stop after 16S clustering and polishing; skip MAG linking.
                       MAGs are optional, and r1/r2 are only required when
                       neither extracted 16S reads nor 16s_fasta is supplied.

    ── MATAM 16S assembly ────────────────────────────────────────────────────────
      --matam_db       Indexed MATAM DB directory or full prefix       [required
                       for extracted-read and raw-read assembly routes]
                       Build: index_default_ssu_rrna_db.py -d \$DBDIR
      --matam_pcts     Subsample percentages, comma-separated           [default: ${params.matam_pcts}]
                       Each value spawns an independent HPC job.
      --matam_threads  CPUs per MATAM assembly job                      [default: ${params.matam_threads}]
                       (tool default: 1; MATAM scales well to ~8–16)
      --matam_mem_mb   Memory per MATAM assembly job (MB)               [default: ${params.matam_mem_mb}]
                       (tool default: 10000; also sets --max_memory in MATAM)
      --skip_matam     Disable MATAM; all samples must provide 16s_fasta [default: ${params.skip_matam}]
      --matam_executable
                       MATAM fork executable                            [default:
                       ${params.matam_executable}]

    ── 16S QC ────────────────────────────────────────────────────────────────────
      --cluster_iden   UCLUST identity threshold for 16S clustering     [default: ${params.cluster_iden}]
                       (uclust_16s.py default: 1.0; matam_16s.py default: 0.999)
      --min_16s_len    Minimum 16S length (bp) accepted for linking     [default: ${params.min_16s_len}]
                       (tool default: 1200)

    ── Linking (MarkerMAG link) ──────────────────────────────────────────────────
      --min_link       Min paired-end linkages to report a MAG–16S hit  [default: ${params.min_link}]
                       (tool default: 9)
      --max_16s_div    Max genetic divergence (%) between linked 16S     [default: ${params.max_16s_div}]
                       (tool default: 1)
      --mismatch       Max mismatch (%) in read-to-16S alignments        [default: ${params.mismatch}]
                       (tool default: 2)
      --aln_len        Min read alignment length (bp)                    [default: ${params.aln_len}]
                       (tool default: 45)
      --aln_pct        Min read alignment percentage (%)                 [default: ${params.aln_pct}]
                       (tool default: 35)

    ── Copy number estimation ────────────────────────────────────────────────────
      --skip_cp_num    Skip copy number estimation inside LINK_16S       [default: ${params.skip_cp_num}]
      --subsample_pct  % of reads used for MAG coverage estimation       [default: ${params.subsample_pct}]
                       (tool default: 25)

    ── Preprocessing ─────────────────────────────────────────────────────────────
      --rename_reads   Rename reads to MarkerMAG format before analysis  [default: ${params.rename_reads}]

    ── Resource ceilings ─────────────────────────────────────────────────────────
      --max_cpus       Maximum CPUs per process                          [default: ${params.max_cpus}]
      --max_memory     Maximum memory per process                        [default: ${params.max_memory}]
      --max_time       Maximum wall time per process                     [default: ${params.max_time}]

    ── Profiles ──────────────────────────────────────────────────────────────────
      -profile standard      Local execution (default)
      -profile hpc           HPC scheduler (PBS/SLURM — edit conf/hpc.config)
      -profile docker        Docker container
      -profile singularity   Singularity container
      -profile conda         Conda environment
      --conda_cache_dir      Shared Conda environment cache
                             [default: ${params.conda_cache_dir}]

    Example:
        nextflow run main.nf \\
            -profile singularity,hpc \\
            --input assets/samplesheet.csv \\
            --matam_db /path/to/SILVA_138_SSURef_NR95 \\
            --matam_threads 12 \\
            --outdir results
    """.stripIndent()
}

// ──────────────────────────────────────────────────────────────────────────────
// Parameter validation
// ──────────────────────────────────────────────────────────────────────────────
def validate_params() {
    if ( !params.input ) {
        error "ERROR: --input samplesheet is required.\n" +
              "Example: --input assets/samplesheet.csv"
    }
}

// Resolve either a MATAM database directory or an explicit database prefix.
def resolve_matam_db_prefix(db_arg) {
    def db_path = file(db_arg).toAbsolutePath()

    if ( db_path.toFile().isDirectory() ) {
        def clustered_fastas = db_path.toFile().listFiles()
            .findAll { it.isFile() && it.name.endsWith(".clustered.fasta") }
            .sort { it.name }

        if ( clustered_fastas.size() == 0 ) {
            error "ERROR: No *.clustered.fasta MATAM database found in: ${db_path}"
        }
        if ( clustered_fastas.size() > 1 ) {
            def candidates = clustered_fastas.collect { it.name }.join(", ")
            error "ERROR: Multiple MATAM databases found in ${db_path}: ${candidates}. " +
                  "Provide the full prefix for the intended database."
        }

        return clustered_fastas[0].absolutePath
            .replaceFirst(/\.clustered\.fasta$/, "")
    }

    def prefix = db_path.toString()
    def clustered_fasta = file("${prefix}.clustered.fasta")
    if ( !clustered_fasta.toFile().isFile() ) {
        error "ERROR: MATAM database file not found: ${clustered_fasta}"
    }
    return prefix
}

// ──────────────────────────────────────────────────────────────────────────────
// Helper: parse samplesheet row → channel tuple
// ──────────────────────────────────────────────────────────────────────────────
def parse_samplesheet_row(row) {
    def meta = [ id: row.sample ]

    def r1        = row.r1      ? file(row.r1,      checkIfExists: true) : []
    def r2        = row.r2      ? file(row.r2,      checkIfExists: true) : []
    def mag_dir   = row.mag_dir ? file(row.mag_dir, checkIfExists: true) : []
    def mag_ext   = row.mag_ext ?: 'fa'
    def reads_16s_single = row['16s_reads']
        ? file(row['16s_reads'], checkIfExists: true)
        : []
    def reads_16s_r1 = row['16s_reads_r1']
        ? file(row['16s_reads_r1'], checkIfExists: true)
        : []
    def reads_16s_r2 = row['16s_reads_r2']
        ? file(row['16s_reads_r2'], checkIfExists: true)
        : []
    if ( (reads_16s_r1 == []) != (reads_16s_r2 == []) ) {
        error "ERROR: Sample '${meta.id}' must provide both 16s_reads_r1 and 16s_reads_r2, or neither."
    }
    if ( reads_16s_single != [] && reads_16s_r1 != [] ) {
        error "ERROR: Sample '${meta.id}' cannot combine 16s_reads with 16s_reads_r1/16s_reads_r2."
    }
    def reads_16s = reads_16s_single != []
        ? [reads_16s_single]
        : (reads_16s_r1 != [] ? [reads_16s_r1, reads_16s_r2] : [])
    def seqs_16s  = row['16s_fasta'] ? file(row['16s_fasta'], checkIfExists: true) : []

    if ( (r1 == []) != (r2 == []) ) {
        error "ERROR: Sample '${meta.id}' must provide both r1 and r2, or neither."
    }

    if ( !params.reconstruct_only ) {
        if ( r1 == [] || r2 == [] ) {
            error "ERROR: Sample '${meta.id}' must provide r1 and r2 for MAG–16S linking."
        }
        if ( mag_dir == [] ) {
            error "ERROR: Sample '${meta.id}' must provide mag_dir for MAG–16S linking."
        }
    }

    if ( seqs_16s == [] ) {
        if ( params.skip_matam ) {
            error "ERROR: Sample '${meta.id}' must provide 16s_fasta when --skip_matam is set."
        }
        if ( !params.matam_db ) {
            error "ERROR: Sample '${meta.id}' requires MATAM assembly, but --matam_db was not provided."
        }
    }

    if ( seqs_16s == [] && reads_16s == [] && (r1 == [] || r2 == []) ) {
        error "ERROR: Sample '${meta.id}' must provide r1 and r2 when MATAM extraction is required."
    }

    return tuple(meta, r1, r2, mag_dir, mag_ext, reads_16s, seqs_16s)
}

// ──────────────────────────────────────────────────────────────────────────────
// Main workflow
// ──────────────────────────────────────────────────────────────────────────────
workflow {

    if ( params.help ) {
        help_message()
        exit 0
    }

    validate_params()

    // ── Samplesheet → channel ────────────────────────────────────────────────
    ch_raw = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row -> parse_samplesheet_row(row) }

    // ── Optional read renaming ────────────────────────────────────────────────
    if ( params.rename_reads ) {
        ch_raw.branch {
                meta, r1, r2, mag_dir, mag_ext, reads_16s, s16 ->
            with_paired_reads:     r1 != [] && r2 != []
            without_paired_reads:  true
        }.set { ch_rename_branched }

        RENAME_READS (
            ch_rename_branched.with_paired_reads.map {
                    meta, r1, r2, mag_dir, mag_ext, reads_16s, s16 ->
                tuple(meta, r1, r2)
            }
        )
        ch_reads_renamed_with_pairs = RENAME_READS.out.reads
            .join(
                ch_rename_branched.with_paired_reads.map {
                        meta, r1, r2, mag_dir, mag_ext, reads_16s, s16 ->
                    tuple(meta, mag_dir, mag_ext, reads_16s, s16)
                },
                by: 0
            )
            .map { meta, r1, r2, mag_dir, mag_ext, reads_16s, s16 ->
                tuple(meta, r1, r2, mag_dir, mag_ext, reads_16s, s16)
            }

        ch_reads_renamed = ch_reads_renamed_with_pairs.mix(
            ch_rename_branched.without_paired_reads
        )
    } else {
        ch_reads_renamed = ch_raw
    }

    // ── Branch: completed 16S, extracted 16S reads, or raw reads ─────────────
    ch_reads_renamed.branch {
            meta, r1, r2, mag_dir, mag_ext, reads_16s, s16 ->
        has_16s:              s16 != []
        has_extracted_reads:  reads_16s != []
        needs_matam_filter:   true
    }.set { ch_branched }

    // ── Path A: user provided 16S ─────────────────────────────────────────────
    // Cluster and polish the provided sequences before linking.
    UCLUST_16S (
        ch_branched.has_16s.map { meta, r1, r2, mag_dir, mag_ext, reads_16s, s16 ->
            tuple(meta, s16)
        },
        params.cluster_iden
    )
    POLISH_16S ( UCLUST_16S.out.clustered )

    ch_provided_16s_ready = POLISH_16S.out.polished

    // ── Path B: run MATAM 16S assembly ──────────────────────────────────────
    if ( !params.skip_matam && params.matam_db ) {
        // Discover the index prefix when a database directory is supplied,
        // then pass its absolute value without staging or renaming it.
        ch_matam_db = Channel.value(resolve_matam_db_prefix(params.matam_db))
        ASSEMBLE_16S (
            ch_branched.needs_matam_filter.map {
                    meta, r1, r2, mag_dir, mag_ext, reads_16s, s16 ->
                tuple(meta, r1, r2)
            },
            ch_branched.has_extracted_reads.map {
                    meta, r1, r2, mag_dir, mag_ext, reads_16s, s16 ->
                tuple(meta, reads_16s)
            }.filter { meta, reads -> reads.size() == 1 }
             .map { meta, reads -> tuple(meta, reads[0]) },
            ch_branched.has_extracted_reads.map {
                    meta, r1, r2, mag_dir, mag_ext, reads_16s, s16 ->
                tuple(meta, reads_16s)
            }.filter { meta, reads -> reads.size() == 2 }
             .map { meta, reads -> tuple(meta, reads[0], reads[1]) },
            ch_matam_db
        )

        ch_matam_16s_ready = ASSEMBLE_16S.out.seqs_16s
    } else {
        ch_matam_16s_ready = Channel.empty()
    }

    // ── Merge both paths into a single 16S channel ──────────────────────────
    ch_16s_ready = ch_provided_16s_ready.mix( ch_matam_16s_ready )

    if ( !params.reconstruct_only ) {
        // ── Reconstruct full sample tuples for MARKERMAG_CORE ───────────────
        ch_samples_for_link = ch_reads_renamed
            .map { meta, r1, r2, mag_dir, mag_ext, reads_16s, s16 ->
                tuple(meta, r1, r2, mag_dir, mag_ext)
            }
            .join( ch_16s_ready, by: 0 )
            .map { meta, r1, r2, mag_dir, mag_ext, marker ->
                tuple(meta, r1, r2, mag_dir, mag_ext, marker)
            }

        // ── Core analysis ────────────────────────────────────────────────────
        MARKERMAG_CORE ( ch_samples_for_link )

        // ── Summary reporting ────────────────────────────────────────────────
        MARKERMAG_CORE.out.linkages_gnm
            .map { meta, f ->
                "${meta.id}\t${params.outdir}/${meta.id}/link_16s/${meta.id}_MarkerMAG_wd/${f.name}"
            }
            .collectFile(name: 'linkages_summary.tsv', newLine: true,
                         storeDir: "${params.outdir}")
    }
}
