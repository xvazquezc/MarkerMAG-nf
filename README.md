
## MarkerMAG-nf: a Nextflow implementation of MarkerMAG

This repository extends MarkerMAG with a reproducible Nextflow DSL2 workflow
while retaining the original MarkerMAG command-line implementation.

[![pypi licence](https://img.shields.io/pypi/l/MarkerMAG.svg)](https://opensource.org/licenses/gpl-3.0.html)
[![pypi version](https://img.shields.io/pypi/v/MarkerMAG.svg)](https://pypi.python.org/pypi/MarkerMAG) 
[![Bioconda](https://img.shields.io/conda/vn/bioconda/markermag.svg?color=green)](https://anaconda.org/bioconda/markermag)
[![DOI](https://img.shields.io/static/v1.svg?label=DOI&message=10.1093/bioinformatics/btac398&color=orange)](https://doi.org/10.1093/bioinformatics/btac398)


Publication
---

+ **Weizhi Song**, Shan Zhang, Torsten Thomas*, MarkerMAG: linking metagenome-assembled genomes (MAGs) with 16S rRNA marker genes using paired-end short reads, Bioinformatics, 2022, btac398, [https://doi.org/10.1093/bioinformatics/btac398](https://doi.org/10.1093/bioinformatics/btac398)
+ Contact: Dr. Weizhi Song (songwz03@gmail.com), Prof. Torsten Thomas (t.thomas@unsw.edu.au)
+ Center for Marine Science & Innovation, University of New South Wales, Sydney, Australia


Updates
---

+ 2022-05-08 - MarkerMAG is now available on Bioconda, please refers to "**How to install**" for details.
+ 2022-03-12 - A [demo dataset](https://doi.org/10.5281/zenodo.6466784) (together with command) has now been provided! You can use it to check if MarkerMAG is installed successfully on your system.


MarkerMAG modules
---

1. Main module

    + `link`: linking MAGs with 16S rRNA marker genes
    
1. Supplementary modules

    + `matam_16s`: assemble 16S rRNA genes with Matam ([manual](doc/README_matam_16s.md))
    + `uclust_16s`: cluster 16S rRNA gene sequences with VSEARCH
    + `polish_16s`: trim non-16S sequence ends with Barrnap
    + `barrnap_16s`: identify 16S rRNA genes from genomes/MAGs with Barrnap ([manual](doc/README_barrnap_16s.md))


How to install
---

+ MarkerMAG is implemented in [python3](https://www.python.org), 
  It has been tested on Linux and MacOS, but NOT on Windows.

Download the repository of MarkerMAG-nf and MATAM (from the [xvazquezc's fork](https://github.com/xvazquezc/matam)):

```bash
git clone https://github.com/xvazquezc/MarkerMAG-nf
git clone --recursive https://github.com/xvazquezc/matam
```

Create and activate the `markermag-nf` environment. This is an explicit user
step; Nextflow will not create or modify the environment:

```bash
cd MarkerMAG-nf
conda env create --file nextflow/environment.yml
conda activate markermag-nf
```

Build the MATAM fork with that environment and expose its commands from the
same Conda prefix:

```bash
cd ../matam
python build.py
for executable in "$PWD"/bin/*; do
    ln -sfn "$executable" "$CONDA_PREFIX/bin/$(basename "$executable")"
done
```

Verify the installation before launching the workflow:

```bash
command -v MarkerMAG
command -v matam_assembly.py
```



How to run
---

+ MarkerMAG’s input consists of 
   1. A set of user-provided MAGs
   2. A set of 16S rRNA gene sequences (either user-provided or generated with the `matam_16s` module) 
   3. Input reads need to be **quality-filtered** and in fasta format (no quality score).
   
+ MarkerMAG is designed for paired short-read data such as Illumina R1/R2
  files. Read identifiers are preserved. Mate identity is obtained from paired
  file position and SAM first/second-mate flags rather than name suffixes.

+ Although you can use your preferred tool to reconstruct 16S rRNA gene sequences from the metagenomic dataset, 
   MarkerMAG does have a supplementary module (`matam_16s`) to reconstruct 16S rRNA genes. 
   Please refer to the manual [here](doc/README_matam_16s.md) if you want to give it a go.

+ Link 16S rRNA gene sequences with MAGs ([demo dataset](https://doi.org/10.5281/zenodo.6466784)): 

      MarkerMAG link -p Demo -r1 demo_R1.fasta -r2 demo_R2.fasta -marker demo_16S.fasta -mag demo_MAGs -x fa -t 12


Nextflow implementation
---

An experimental [Nextflow DSL2](https://www.nextflow.io/) implementation is
available in the [`nextflow`](nextflow) directory. It coordinates read
preprocessing, 16S assembly or preparation, MAG–16S linking, and optional copy
number estimation as reproducible processes.

The workflow has three per-sample entry points:

1. If `16s_fasta` is supplied, the sequences are clustered with `uclust_16s`,
   polished with Barrnap, and passed to `MarkerMAG link`.
2. If the `16s_reads_r1`/`16s_reads_r2` pair is supplied, MATAM extraction is
   skipped and those reads enter assembly.
3. If all optional 16S fields are empty, MATAM extracts candidate reads from
   `r1` and `r2`.

For the two assembly routes, each value in `--matam_pcts` is assembled as an
independent Nextflow task, making the expensive assemblies suitable for
parallel HPC submission. The assemblies are then combined, clustered,
polished, and linked.

MATAM paired-read subsampling is deterministic. Use
`--matam_subsample_seed` to choose the random seed (default: `42`); use the
same seed to reproduce a previous set of subsamples.

### Requirements

+ Nextflow 23.10 or newer.
+ The user-created and activated `markermag-nf` Conda environment described
  under “How to install”, containing all workflow dependencies.
+ VSEARCH available on `PATH` from the activated environment.
+ MATAM and an indexed MATAM reference database when using the assembly route.

The workflow never creates, updates, or activates Conda environments. Activate
`markermag-nf` before starting Nextflow and ensure scheduler jobs inherit its
`PATH`.

### Samplesheet

Pass a CSV file with `--input`. Relative paths are resolved from the directory
where Nextflow is launched.

| Column | Description |
|:---|:---|
| `sample` | Unique sample identifier |
| `r1` | Forward reads; optional in `--reconstruct_only` mode when extracted 16S reads or `16s_fasta` are supplied |
| `r2` | Reverse reads; required under the same conditions as `r1` |
| `mag_dir` | Directory containing the sample's MAG files; optional in `--reconstruct_only` mode |
| `mag_ext` | MAG filename extension without the leading dot, such as `fa` |
| `16s_reads_r1` | Optional forward candidate 16S reads; requires `16s_reads_r2` |
| `16s_reads_r2` | Optional reverse candidate 16S reads; requires `16s_reads_r1` |
| `16s_fasta` | Optional prepared 16S FASTA; leave empty to run MATAM |

See [`nextflow/assets/samplesheet.csv`](nextflow/assets/samplesheet.csv) for a
template containing examples of all three checkpoints. If optional reads and
sequences are both populated, `16s_fasta` takes precedence.

### Run with supplied 16S sequences

    cd nextflow
    nextflow run main.nf \
        --input assets/samplesheet.csv \
        --skip_matam \
        --outdir results

All samples must provide `16s_fasta` when `--skip_matam` is used.

### Reconstruct 16S sequences only

Use `--reconstruct_only` to stop after 16S clustering and Barrnap polishing,
without running MAG–16S linking or copy-number estimation:

    cd nextflow
    nextflow run main.nf \
        --reconstruct_only \
        --input samplesheet.csv \
        --matam_db /path/to/indexed/SILVA_138_SSURef_NR95 \
        --outdir results

Input requirements depend on the selected checkpoint:

+ With `16s_fasta`, `r1`, `r2`, and `mag_dir` may all be empty.
+ With `16s_reads_r1`/`16s_reads_r2`, `r1`, `r2`,
  and `mag_dir` may all be empty, but `--matam_db` remains required.
+ With no optional 16S input, `r1` and `r2` are required for MATAM extraction;
  `mag_dir` remains optional.

The final sequences are written to
`<outdir>/<sample>/polish_16s/<sample>_16S_polished.fasta`.

### Resume from extracted 16S reads

Put paired candidate reads in `16s_reads_r1` and `16s_reads_r2`. Interleaved
and single-end reads are not supported. Leave `16s_fasta` empty. The original
`r1` and `r2` files are still required for the later MarkerMAG linking stage.

    cd nextflow
    nextflow run main.nf \
        --input assets/samplesheet.csv \
        --matam_db /path/to/indexed/SILVA_138_SSURef_NR95 \
        --matam_pcts 1,5,10,25,50,75,100 \
        --outdir results

This checkpoint bypasses `MATAM_FILTER`; the supplied reads enter
`MATAM_ASSEMBLE` directly. A MATAM database is still required for assembly.

Paired files are subsampled together by record position before each MATAM task.

### Assemble 16S sequences with MATAM

    cd nextflow
    nextflow run main.nf \
        --input assets/samplesheet.csv \
        --matam_db /path/to/indexed/SILVA_138_SSURef_NR95 \
        --matam_pcts 1,5,10,25,50,75,100 \
        --matam_threads 8 \
        --matam_mem_mb 30000 \
        --outdir results

`--matam_db` may be a directory containing exactly one `*.clustered.fasta`
MATAM index, or the full database prefix when a directory contains multiple
indexes. The resolved prefix is passed to MATAM as an absolute path and must
be visible from every compute node.

The database must already have its SortMeRNA indexes (`<prefix>.clustered.smr/idx`
and `<prefix>.complete.smr/idx`), as built by `matam_db_preprocessing.py` without
`--skip_index`. MarkerMAG-nf checks for them before starting any task, and MATAM
uses only these indexes: it never builds its own in each task (hours for SILVA). To index an existing database once, run
`index_ref_db.py -i <prefix> --max_memory 30000`.

Use `--matam_executable /path/to/matam/bin/matam_assembly.py` to select
the fork explicitly.

### Execution profiles

| Profile | Purpose |
|:---|:---|
| `standard` | Local execution using tools already on `PATH` |
| `test` | Small local resources and the supplied synthetic test |
| `hpc` | PBS Pro by default; edit [`conf/hpc.config`](nextflow/conf/hpc.config) for the local scheduler and queue |
| `docker` | Run with the local `markermag:latest` image |
| `singularity` | Run from `docker://markermag:latest` |

Use `-profile hpc` for scheduler execution. Docker and Singularity remain
optional, explicitly built container runtimes. Build the Docker image from the
repository root:

    docker build -f nextflow/Dockerfile -t markermag:latest .

### Nextflow outputs

Per-sample outputs are written beneath `--outdir/<sample>/`:

+ `uclust_16s/`: clustered 16S FASTA, VSEARCH `.uc` table, and membership table.
+ `polish_16s/`: Barrnap-polished 16S FASTA.
+ `link_16s/`: MarkerMAG linkage tables, logs, plots, and working outputs.
+ `get_cp_num/`: standalone copy-number outputs when enabled.

The top-level `linkages_summary.tsv` maps each sample to its persistent
genome-level linkage table.

### Quick test

The test generator creates a small deterministic dataset from the bundled SILVA
reference:

    cd nextflow
    python test/generate_test_data.py
    nextflow run main.nf \
        -profile test \
        --input test/samplesheet.csv \
        --outdir test/results

The test follows the supplied-16S route and does not require a MATAM database.
Nextflow runtime state is written to ignored paths and can be removed with
`nextflow clean -f` after a run.


Output files
---

1. Summary of identified linkages at genome level:

    | Marker | MAG | Linkage | Round |
    |:---:|:---:|:---:|:---:|
    | matam_16S_7   | MAG_6 | 181| Rd1 |
    | matam_16S_12  | MAG_9 | 102| Rd1 |
    | matam_16S_6   | MAG_59| 55 | Rd2 |

2. Summary of identified linkages at contig level (with figure):

    |Marker___MAG (linkages)	|Contig	        |Round_1	|Round_2	|
    |:---:|:---:|:---:|:---:|
    |matam_16S_7___MAG_6(181)	            |Contig_1799	|176	    |0          |
    |matam_16S_7___MAG_6(181)	            |Contig_1044	|5	        |0          |
    |matam_16S_12___MAG_9(102)	            |Contig_840	    |102	    |0          |
    |matam_16S_6___MAG_59(39)	            |Contig_171	    |0	        |55         |

   ![linkages](doc/images/linkages_plot_2.png)

3. Copy number of linked 16S rRNA genes.


4. Visualization of individual linkage.
  
   MarkerMAG supports the visualization of identified linkages (needs [Tablet](https://ics.hutton.ac.uk/tablet/)). 
   Output files for visualization can be found in the [Prefix]_linkage_visualization_rd1/2 folders.
   You can visualize how the linking reads are aligned to MAG contig and 16S rRNA gene by double-clicking the corresponding ".tablet" file. 
   Fifty Ns are added between the linked MAG contig and 16S rRNA gene.
 
   ![linkages](doc/images/linking_reads.png)
 
   *If you saw error message from Tablet that says input files format can not be understood, 
   please refer to [here](https://github.com/cropgeeks/tablet/issues/15) for a potential solution.
