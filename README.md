# dragenflow

A Nextflow DSL2 pipeline for running DRAGEN alignment and variant calling workflows. Built on the nf-core framework, it supports multiple analysis modes (germline, somatic, RNA-seq, methylation) with configurable compute profiles for on-premise DRAGEN hardware and AWS.

## Table of Contents

- [Pipeline Overview](#pipeline-overview)
- [Analysis Modes](#analysis-modes)
- [Samplesheet Formats](#samplesheet-formats)
- [Run Command](#run-command)
- [Compute Profiles](#compute-profiles)
- [Testing with Example Data](#testing-with-example-data)
- [Key Parameters](#key-parameters)
- [Outputs](#outputs)
- [Reference Files](#reference-files)

---

## Pipeline Overview

All analysis modes share these common steps:

1. **Parse samplesheet** — validates and routes samples to alignment or demux entry points
2. **Gather alignment samples** — collects FASTQs from reads, fastq_lists, demux paths, or converts CRAMs/BAMs back to FASTQ
3. **Create fastq_list** — generates per-sample DRAGEN fastq_list CSV
4. **DRAGEN alignment** — runs `dragen` with mode-specific arguments
5. **Variant annotation** — VEP annotation of SNVs, SVs, and CNVs → TSV tables
6. **MultiQC** — aggregates QC metrics

Mode-specific downstream steps are described in [Analysis Modes](#analysis-modes) below.

---

## Analysis Modes

Each mode is activated by passing the corresponding profile. The `workflow` parameter is set automatically by the profile.

### `alignonly`

Alignment only — no variant calling.

```bash
nextflow run dhslab/dragenflow -r main \
  -profile dhslab,ris2,dragen4,alignonly \
  --input samplesheet.csv \
  --outdir results/
```

### `germline`

Germline SNV, SV, and CNV calling. Outputs VEP-annotated VCFs and TSV tables for each variant type.

```bash
nextflow run dhslab/dragenflow -r main \
  -profile dhslab,ris2,dragen4,germline \
  --input samplesheet.csv \
  --outdir results/
```

DRAGEN flags enabled: `--enable-variant-caller`, `--enable-sv`, `--enable-cnv`

### `somaticheme`

Somatic tumor/normal calling for **hematologic malignancies** (liquid tumor mode). Requires a samplesheet with matched tumor and normal samples. Includes SNV, SV, and CNV calling with systematic noise filtering.

```bash
nextflow run dhslab/dragenflow -r main \
  -profile dhslab,ris2,dragen4,somaticheme \
  --input somatic_samplesheet.csv \
  --outdir results/
```

Default noise files (hg38):

- SNV: `IDPF_WGS_hg38_v.2.0.0_systematic_noise.snv.bed.gz`
- SV: `WGS_FF_Heme_hg38_v3.1.0_systematic_noise.sv.bedpe.gz`

### `somaticsolid`

Somatic tumor/normal calling for **solid tumors**. Similar to `somaticheme` but without liquid tumor UMI options; includes TMB estimation.

```bash
nextflow run dhslab/dragenflow -r main \
  -profile dhslab,ris2,dragen4,somaticsolid \
  --input somatic_samplesheet.csv \
  --outdir results/
```

### `tumoronlyheme`

Tumor-only variant calling for **hematologic malignancies** — no matched normal required. Enables ploidy estimation, DUX4 fusion caller, and CNV calling against a population VCF.

```bash
nextflow run dhslab/dragenflow -r main \
  -profile dhslab,ris2,dragen4,tumoronlyheme \
  --input samplesheet.csv \
  --outdir results/
```

Default population VCF: `1000G_phase1.snps.high_confidence.hg38.vcf.gz`

### `rnaseq`

RNA-seq alignment with quantification and gene fusion detection. Downstream steps annotate gene/transcript expression tables.

```bash
nextflow run dhslab/dragenflow -r main \
  -profile dhslab,ris2,dragen4,rnaseq \
  --input samplesheet.csv \
  --outdir results/
```

DRAGEN flags enabled: `--enable-rna`, `--enable-rna-quantification`, `--enable-rna-gene-fusion`, `--rrna-filter-enable`

Default annotation GTF: `Homo_sapiens.GRCh38.105.chr.sorted.gtf.gz`

### `bsseq`

Bisulfite sequencing — whole-genome methylation calling (5mC). Downstream steps generate CpG methylation BED and BigWig files.

```bash
nextflow run dhslab/dragenflow -r main \
  -profile dhslab,ris2,dragen4,bsseq \
  --input samplesheet.csv \
  --outdir results/
```

Uses a methylation-specific DRAGEN reference directory (`dragen_hg38_5mCv4.3.6`).

### `download_nirvana`

Not a samplesheet-driven analysis mode — downloads the Nirvana annotation data bundle used by `--use_nirvana`/`--nirvana_path` instead of running alignment or variant calling. Set with `--workflow download_nirvana` directly (there is no dedicated profile) and no `--input` is required:

```bash
nextflow run dhslab/dragenflow -r main \
  -profile dhslab,ris2,dragen4 \
  --workflow download_nirvana \
  --outdir /path/to/nirvana_annotation_data
```

The downloaded data lands in `<outdir>/nirvana_annotation_data` and can be pointed to via `--nirvana_path` in subsequent runs. `--nirvana_assembly` (default `GRCh38`) selects which assembly's annotation bundle is fetched.

### `idtumi`

UMI-aware adapter processing using IDT UDI-UMI (10x19) format. Disables duplicate marking (UMI deduplication is used instead). Combine with another analysis profile for variant calling.

```bash
nextflow run dhslab/dragenflow -r main \
  -profile dhslab,ris2,dragen4,idtumi,tumoronlyheme \
  --input samplesheet.csv \
  --outdir results/
```

---

## Samplesheet Formats

The pipeline accepts several CSV formats. The first column is always `id` (or a pair-level ID for wide somatic format).

### Standard samplesheet

For most modes (germline, alignonly, rnaseq, bsseq, tumoronly):

```csv
id,read1,read2,fastq_list,cram,bam
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,,,
sample2,,,/path/to/sample2_fastq_list.csv,,
sample3,,,,/path/to/sample3.cram,
```

Columns are optional — provide whichever input type is available. Multiple input types for the same sample are merged (e.g., existing CRAM + new FASTQs will both be realigned together), as in `assets/stub/cram_fastq_wgs_mastersheet.csv`, where `10929-BM-lib1` is realigned from its existing CRAM plus an additional lane of FASTQs:

```csv
id,cram,read1,read2
10940-BM-lib1,assets/stub/crams/10940-BM-lib1_tumor.cram,,
10929-BM-lib1,assets/stub/crams/10929-BM-lib1_input.cram,,
10929-BM-lib1,,assets/stub/demux_fastq/10929-BM-lib1_S1_L007_R1_001.fastq.gz,assets/stub/demux_fastq/10929-BM-lib1_S1_L007_R2_001.fastq.gz
```

A BAM can be substituted for a CRAM the same way — see `assets/stub/bam_wgs_test_mastersheet.csv`:

```csv
id,bam
10929-BM-lib1,assets/stub/crams/10929-BM-lib1_input.bam
```

### DRAGEN fastq_list

A pre-built DRAGEN-format fastq_list can be referenced from the `fastq_list` column, as in `assets/stub/fastq_list_wgs_mastersheet.csv`:

```csv
id,fastq_list
10929-BM-lib1,assets/stub/demux_fastq/fastq_list.csv
```

where the referenced file is a standard DRAGEN fastq_list (one row per lane):

```csv
RGID,RGSM,RGLB,Lane,Read1File,Read2File
TAATGTGTCT.TATGCCTTAC.5,10929-BM-lib1,UnknownLibrary,7,assets/stub/demux_fastq/10929-BM-lib1_S1_L007_R1_001.fastq.gz,assets/stub/demux_fastq/10929-BM-lib1_S1_L007_R2_001.fastq.gz
TAATGTGTCT.TATGCCTTAC.6,10929-BM-lib1,UnknownLibrary,8,assets/stub/demux_fastq/10929-BM-lib1_S1_L008_R1_001.fastq.gz,assets/stub/demux_fastq/10929-BM-lib1_S1_L008_R2_001.fastq.gz
```

### Somatic — long format (recommended)

Tumor and normal samples are listed as separate rows sharing the same pair-level `id`, distinguished by `sample_id` and `sample_type`. See `assets/stub/somatic_reads_mastersheet.csv`:

```csv
id,sample_id,sample_type,read1,read2
KDM7A-KO-CART,CD34-CART-DNA,normal,assets/stub/demux_fastq/CD34-CART-DNA_S6_L003_R1_001.fastq.gz,assets/stub/demux_fastq/CD34-CART-DNA_S6_L003_R2_001.fastq.gz
KDM7A-KO-CART,KDM7A-KO-CART-DNA,tumor,assets/stub/demux_fastq/KDM7A-KO-CART-DNA_S3_L005_R1_001.fastq.gz,assets/stub/demux_fastq/KDM7A-KO-CART-DNA_S3_L005_R2_001.fastq.gz
```

### Somatic — wide format

Tumor and normal can also be specified in a single row using prefixed columns. See `assets/stub/somatic_reads_wide_mastersheet.csv`:

```csv
id,tumor_id,normal_id,tumor_read1,tumor_read2,normal_read1,normal_read2
KDM7A-KO-CART,KDM7A-KO-CART-tumor,KDM7A-KO-CART-normal,assets/stub/demux_fastq/KDM7A-KO-CART-DNA_S3_L005_R1_001.fastq.gz,assets/stub/demux_fastq/KDM7A-KO-CART-DNA_S3_L005_R2_001.fastq.gz,assets/stub/demux_fastq/CD34-CART-DNA_S6_L003_R1_001.fastq.gz,assets/stub/demux_fastq/CD34-CART-DNA_S6_L003_R2_001.fastq.gz
```

Or referencing a samplemap. See `assets/stub/somatic_samplemap_wide_mastersheet.csv`:

```csv
id,tumor_id,normal_id,samplemap
KDM7A-KO-CART,KDM7A-KO-CART-DNA,CD34-CART-DNA,assets/stub/demux_fastq/Samplemap2.somatic.csv
```

### MGI Samplemap format

MGI sequencer output (`Samplemap2.csv`) is also accepted directly via the `samplemap` column. The `Library Name` field in the samplemap must match the sample identifier — `sample_id` when both `sample_id` and `sample_type` are present (long-format somatic, see `assets/stub/somatic_samplemap_mastersheet.csv`), otherwise `id` (see `assets/stub/samplemap_idtumi_mastersheet.csv`):

```csv
id,samplemap
GM24385_unedited,assets/stub/demux_fastq/Samplemap2.umi.csv
```

```csv
id,sample_id,sample_type,samplemap
KDM7A-KO-CART,CD34-CART-DNA,normal,assets/stub/demux_fastq/Samplemap2.somatic.csv
KDM7A-KO-CART,KDM7A-KO-CART-DNA,tumor,assets/stub/demux_fastq/Samplemap2.somatic.csv
```

---

## Run Command

```bash
nextflow run dhslab/dragenflow -r main \
  -profile dhslab,ris2,dragen4,<analysis-profile> \
  --input /path/to/samplesheet.csv \
  --outdir /path/to/results/
```

### Profiles (stack as comma-separated list)

| Profile           | Purpose                              |
| ----------------- | ------------------------------------ |
| `dhslab`          | dhslab group settings (queue, paths) |
| `ris`             | WashU compute1, LSF scheduler        |
| `ris2`            | WashU compute2, SLURM scheduler      |
| `dragen4`         | On-premise DRAGEN 4.3.6 hardware     |
| `dhslabdragenaws` | AWS Batch with DRAGEN 4.4.6          |
| `alignonly`       | Alignment-only workflow              |
| `germline`        | Germline variant calling             |
| `somaticheme`     | Somatic heme tumor/normal            |
| `somaticsolid`    | Somatic solid tumor/normal           |
| `tumoronlyheme`   | Tumor-only heme                      |
| `rnaseq`          | RNA-seq                              |
| `bsseq`           | Bisulfite/methylation                |
| `idtumi`          | IDT UMI adapter mode                 |

---

## Compute Profiles

### WashU RIS compute1 (LSF)

```bash
-profile dhslab,ris,dragen4,<analysis-profile>
```

Requires `--user_group`, `--queue`, and `--job_group_name` (set in `dhslab` profile or on command line).

### WashU RIS compute2 (SLURM)

```bash
-profile dhslab,ris2,dragen4,<analysis-profile>
```

Uses `condo-dspencer` partition by default. Override with `--slurm_partition`.

### AWS Batch

```bash
nextflow secrets set AWS_ACCESS_KEY <key>
nextflow secrets set AWS_SECRET_KEY <secret>
nextflow secrets set AWS_DRAGEN_USER <dragen_user>
nextflow secrets set AWS_DRAGEN_PASSWORD <dragen_password>

nextflow run dhslab/dragenflow -r main \
  -profile dhslab,dhslabdragenaws,<analysis-profile> \
  --input samplesheet.csv \
  --outdir s3://your-bucket/results/ \
  -bucket-dir s3://dhslab-dragen-data/work/
```

Uses DRAGEN 4.4.6 on AWS with Wave/Fusion enabled.

---

## Testing with Example Data

`assets/stub/` contains small real FASTQ/CRAM/BAM files plus pre-built samplesheets ("mastersheets") that exercise every input pathway described above, so the pipeline can be exercised end-to-end without production data.

| File                                     | Input type                              | Pairs with                                 |
| ---------------------------------------- | --------------------------------------- | ------------------------------------------ |
| `reads_wgs_mastersheet.csv`              | Paired FASTQs                           | `alignonly`, `germline`                    |
| `fastq_list_wgs_mastersheet.csv`         | DRAGEN fastq_list CSV                   | `alignonly`, `germline`                    |
| `bam_wgs_test_mastersheet.csv`           | BAM for realignment                     | `alignonly`, `germline`                    |
| `cram_fastq_wgs_mastersheet.csv`         | CRAM + FASTQ merged per sample          | `alignonly`, `germline`                    |
| `reads_rnaseq_mastersheet.csv`           | Paired FASTQs                           | `rnaseq`                                   |
| `reads_idtumi_mastersheet.csv`           | Paired UMI FASTQs                       | `idtumi` (+ `tumoronlyheme`/`somaticheme`) |
| `samplemap_idtumi_mastersheet.csv`       | MGI samplemap                           | `idtumi` (+ `tumoronlyheme`/`somaticheme`) |
| `somatic_reads_mastersheet.csv`          | Tumor/normal FASTQs, long format        | `somaticheme`, `somaticsolid`              |
| `somatic_reads_wide_mastersheet.csv`     | Tumor/normal FASTQs, wide format        | `somaticheme`, `somaticsolid`              |
| `somatic_samplemap_mastersheet.csv`      | Tumor/normal MGI samplemap, long format | `somaticheme`, `somaticsolid`              |
| `somatic_samplemap_wide_mastersheet.csv` | Tumor/normal MGI samplemap, wide format | `somaticheme`, `somaticsolid`              |

> `cram_wgs_test_mastersheet.csv` currently points to a CRAM filename that isn't present in `assets/stub/crams/` — use `cram_fastq_wgs_mastersheet.csv` or `bam_wgs_test_mastersheet.csv` for a working CRAM/BAM realignment example instead.

`assets/stub/dragen_path/` holds a complete set of pre-computed DRAGEN outputs for `10929-BM-lib1` (VCFs, metrics, CRAM), useful as a reference for expected output file naming/contents.

### `-profile stub` / `-stub-run`

For quick validation of pipeline logic (channel wiring, samplesheet parsing, module I/O) without DRAGEN hardware or a license, add the `stub` profile and Nextflow's `-stub-run` flag. Every process returns canned output instead of actually running:

```bash
nextflow run main.nf -profile stub,alignonly -stub-run \
  --input assets/stub/cram_fastq_wgs_mastersheet.csv \
  --outdir results_test
```

> Nextflow's `-stub-run`/`-stub` flag has a known quirk where it also populates the pipeline's positional-argument list, producing a harmless `WARN: nf-core pipelines do not accept positional arguments. The positional argument \`true\` has been detected.` on every stub run. This is a Nextflow CLI parsing artifact, not a real usage mistake, and can be ignored.

### Example runs against real DRAGEN hardware

```bash
# Germline WGS from paired FASTQs
nextflow run main.nf -profile dhslab,ris2,dragen4,germline \
  --input assets/stub/reads_wgs_mastersheet.csv \
  --outdir results_test

# Somatic heme tumor/normal via MGI samplemap (wide format)
nextflow run main.nf -profile dhslab,ris2,dragen4,somaticheme \
  --input assets/stub/somatic_samplemap_wide_mastersheet.csv \
  --outdir results_somatic

# RNA-seq
nextflow run main.nf -profile dhslab,ris2,dragen4,rnaseq \
  --input assets/stub/reads_rnaseq_mastersheet.csv \
  --outdir results_rnatest

# IDT UMI adapter mode, tumor-only heme calling, restricted to a target BED
nextflow run main.nf -profile dhslab,ris2,dragen4,idtumi,tumoronlyheme \
  --input assets/stub/reads_idtumi_mastersheet.csv \
  --target_bedfile /path/to/targets.bed \
  --outdir results_umi_test
```

Add `--run_dragen false` to any run to execute everything up to (but not including) the DRAGEN call — useful for checking that inputs were parsed and staged correctly before committing compute time.

---

## Key Parameters

| Parameter                     | Default                | Description                                                |
| ----------------------------- | ---------------------- | ---------------------------------------------------------- |
| `--input`                     | —                      | Path to samplesheet CSV (required)                         |
| `--outdir`                    | —                      | Output directory (required)                                |
| `--refdir`                    | hg38 dragen_hg38v4.3.6 | DRAGEN hash-table reference directory                      |
| `--fasta`                     | hg38_mgi_patch.fa      | FASTA reference (for CRAM decoding, VEP)                   |
| `--target_bedfile`            | null                   | Target BED for exome/panel; enables exome SV/CNV calling   |
| `--extra_dragen_args`         | null                   | Additional raw arguments passed to `dragen`                |
| `--mark_duplicates`           | true                   | Enable duplicate marking                                   |
| `--alignment_file_format`     | CRAM                   | Output format: CRAM or BAM                                 |
| `--umi`                       | null                   | UMI library type (e.g., `random-simplex`)                  |
| `--readfamilysize`            | 3                      | Minimum UMI family size                                    |
| `--liquid_tumor`              | false                  | Enable liquid tumor UMI mode                               |
| `--solid_tumor`               | false                  | Enable solid tumor UMI mode                                |
| `--dux4caller`                | false                  | Enable DRAGEN DUX4 fusion caller                           |
| `--use_nirvana`               | true                   | Use Nirvana for variant annotation (instead of dbSNP only) |
| `--vepcache`                  | VEP113_cache           | VEP cache directory                                        |
| `--hotspot_vcf`               | null                   | Somatic hotspot VCF                                        |
| `--hotspot_bed`               | null                   | BED file to generate hotspot VCF on the fly                |
| `--snv_noisefile`             | null                   | Systematic noise BED for SNV filtering                     |
| `--sv_noisefile`              | null                   | Systematic noise BEDPE for SV filtering                    |
| `--cnv_population_vcf`        | null                   | Population BAF VCF for CNV calling                         |
| `--dragen_cnv_filter_length`  | null                   | Minimum CNV length filter (bp)                             |
| `--dragen_cnv_merge_distance` | null                   | CNV merge distance (bp)                                    |
| `--downsample_rna`            | null                   | Downsample RNA reads (reads count)                         |
| `--run_dragen`                | true                   | Set to false to skip DRAGEN (for testing)                  |

---

## Outputs

Results are written per-sample to `<outdir>/<sample_id>/`. DRAGEN output files are copied directly from the DRAGEN output directory.

| Mode               | Key outputs                                                                                                |
| ------------------ | ---------------------------------------------------------------------------------------------------------- |
| All                | `<id>.cram` / `<id>.bam`, alignment QC metrics, `pipeline_info/`                                           |
| germline / somatic | `<id>.hard-filtered.vcf.gz` (SNVs), `<id>.sv.vcf.gz`, `<id>.cnv.vcf.gz`, VEP-annotated VCFs and TSV tables |
| rnaseq             | `<id>.quant.genes.sf`, `<id>.quant.sf`, `<id>.fusion_candidates.final`, annotated expression tables        |
| bsseq              | `<id>.methylation_call_file.gz`, methylation BED, BigWig                                                   |
| All                | `multiqc_report.html`                                                                                      |

---

## Reference Files

Default paths (hg38, dhslab group storage):

| Resource               | Path                                                                                                                 |
| ---------------------- | -------------------------------------------------------------------------------------------------------------------- |
| DRAGEN reference (DNA) | `/storage2/fs1/dspencer/Active/shared/refdata/hg38/dragen_hg38v4.3.6`                                                |
| DRAGEN reference (5mC) | `/storage2/fs1/dspencer/Active/shared/refdata/hg38/dragen_hg38_5mCv4.3.6`                                            |
| FASTA                  | `/storage2/fs1/dspencer/Active/shared/refdata/hg38/sequence/hg38_mgi_patch.fa`                                       |
| dbSNP                  | `/storage2/fs1/dspencer/Active/shared/refdata/hg38/dragenfiles/dbsnp.vcf.gz`                                         |
| VEP cache              | `/storage2/fs1/dspencer/Active/shared/refdata/hg38/VEP113_cache`                                                     |
| Nirvana annotation     | `/storage2/fs1/dspencer/Active/shared/refdata/hg38/dragenfiles/nirvana_annotation_data_323`                          |
| SNV noise (heme WGS)   | `/storage2/fs1/dspencer/Active/shared/refdata/hg38/dragenfiles/IDPF_WGS_hg38_v.2.0.0_systematic_noise.snv.bed.gz`    |
| SV noise (heme WGS)    | `/storage2/fs1/dspencer/Active/shared/refdata/hg38/dragenfiles/WGS_FF_Heme_hg38_v3.1.0_systematic_noise.sv.bedpe.gz` |
| Population BAF VCF     | `/storage2/fs1/dspencer/Active/shared/refdata/hg38/dragenfiles/1000G_phase1.snps.high_confidence.hg38.vcf.gz`        |
| Ensembl GTF (RNA)      | `/storage2/fs1/dspencer/Active/shared/refdata/hg38/ensemble/Homo_sapiens.GRCh38.105.chr.sorted.gtf.gz`               |
| Cytobands              | `assets/data/hg38.cytoBandIdeo.bed.gz`                                                                               |
