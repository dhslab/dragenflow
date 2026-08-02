# dragenflow

A Nextflow DSL2 pipeline for running DRAGEN alignment and variant calling workflows. Built on the nf-core framework, it supports multiple analysis modes (germline, somatic, RNA-seq, methylation) with configurable compute profiles for on-premise DRAGEN hardware and AWS.

## Table of Contents

- [Pipeline Overview](#pipeline-overview)
- [Analysis Modes](#analysis-modes)
- [Samplesheet Formats](#samplesheet-formats)
- [Run Command](#run-command)
- [Compute Profiles](#compute-profiles)
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

Columns are optional — provide whichever input type is available. Multiple input types for the same sample are merged (e.g., existing CRAM + new FASTQs will both be realigned together).

### DRAGEN fastq_list

A pre-built DRAGEN-format fastq_list can be referenced from the `fastq_list` column:

```csv
RGID,RGSM,RGLB,Lane,Read1File,Read2File
TAATGTGTCT.TATGCCTTAC.5,sample1,UnknownLibrary,5,/path/R1.fastq.gz,/path/R2.fastq.gz
```

### Somatic — long format (recommended)

Tumor and normal samples are listed as separate rows with `individual_id` and `sample_type`:

```csv
id,individual_id,sample_type,read1,read2
PatientA-tumor,PatientA,tumor,/path/tumor_R1.fastq.gz,/path/tumor_R2.fastq.gz
PatientA-normal,PatientA,normal,/path/normal_R1.fastq.gz,/path/normal_R2.fastq.gz
```

### Somatic — wide format

Tumor and normal can also be specified in a single row using prefixed columns:

```csv
id,tumor_id,normal_id,tumor_read1,tumor_read2,normal_read1,normal_read2
PatientA,PatientA-tumor,PatientA-normal,/path/tumor_R1.fastq.gz,/path/tumor_R2.fastq.gz,/path/normal_R1.fastq.gz,/path/normal_R2.fastq.gz
```

Or referencing a samplemap:

```csv
id,tumor_id,normal_id,samplemap
PatientA,PatientA-tumor,PatientA-normal,/path/to/Samplemap2.csv
```

### MGI Samplemap format

MGI sequencer output (Samplemap2.csv) is also accepted directly via the `samplemap` column. The `Library Name` field in the samplemap must match the sample `id`.

---

## Run Command

```bash
nextflow run dhslab/dragenflow -r main \
  -profile dhslab,ris2,dragen4,<analysis-profile> \
  --input /path/to/samplesheet.csv \
  --outdir /path/to/results/
```

### Profiles (stack as comma-separated list)

| Profile | Purpose |
|---|---|
| `dhslab` | dhslab group settings (queue, paths) |
| `ris` | WashU compute1, LSF scheduler |
| `ris2` | WashU compute2, SLURM scheduler |
| `dragen4` | On-premise DRAGEN 4.3.6 hardware |
| `dhslabdragenaws` | AWS Batch with DRAGEN 4.4.6 |
| `alignonly` | Alignment-only workflow |
| `germline` | Germline variant calling |
| `somaticheme` | Somatic heme tumor/normal |
| `somaticsolid` | Somatic solid tumor/normal |
| `tumoronlyheme` | Tumor-only heme |
| `rnaseq` | RNA-seq |
| `bsseq` | Bisulfite/methylation |
| `idtumi` | IDT UMI adapter mode |

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

## Key Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | — | Path to samplesheet CSV (required) |
| `--outdir` | — | Output directory (required) |
| `--refdir` | hg38 dragen_hg38v4.3.6 | DRAGEN hash-table reference directory |
| `--fasta` | hg38_mgi_patch.fa | FASTA reference (for CRAM decoding, VEP) |
| `--target_bedfile` | null | Target BED for exome/panel; enables exome SV/CNV calling |
| `--extra_dragen_args` | null | Additional raw arguments passed to `dragen` |
| `--mark_duplicates` | true | Enable duplicate marking |
| `--alignment_file_format` | CRAM | Output format: CRAM or BAM |
| `--umi` | null | UMI library type (e.g., `random-simplex`) |
| `--readfamilysize` | 3 | Minimum UMI family size |
| `--liquid_tumor` | false | Enable liquid tumor UMI mode |
| `--solid_tumor` | false | Enable solid tumor UMI mode |
| `--dux4caller` | false | Enable DRAGEN DUX4 fusion caller |
| `--use_nirvana` | true | Use Nirvana for variant annotation (instead of dbSNP only) |
| `--vepcache` | VEP113_cache | VEP cache directory |
| `--hotspot_vcf` | null | Somatic hotspot VCF |
| `--hotspot_bed` | null | BED file to generate hotspot VCF on the fly |
| `--snv_noisefile` | null | Systematic noise BED for SNV filtering |
| `--sv_noisefile` | null | Systematic noise BEDPE for SV filtering |
| `--cnv_population_vcf` | null | Population BAF VCF for CNV calling |
| `--dragen_cnv_filter_length` | null | Minimum CNV length filter (bp) |
| `--dragen_cnv_merge_distance` | null | CNV merge distance (bp) |
| `--downsample_rna` | null | Downsample RNA reads (reads count) |
| `--run_dragen` | true | Set to false to skip DRAGEN (for testing) |

---

## Outputs

Results are written per-sample to `<outdir>/<sample_id>/`. DRAGEN output files are copied directly from the DRAGEN output directory.

| Mode | Key outputs |
|---|---|
| All | `<id>.cram` / `<id>.bam`, alignment QC metrics, `pipeline_info/` |
| germline / somatic | `<id>.hard-filtered.vcf.gz` (SNVs), `<id>.sv.vcf.gz`, `<id>.cnv.vcf.gz`, VEP-annotated VCFs and TSV tables |
| rnaseq | `<id>.quant.genes.sf`, `<id>.quant.sf`, `<id>.fusion_candidates.final`, annotated expression tables |
| bsseq | `<id>.methylation_call_file.gz`, methylation BED, BigWig |
| All | `multiqc_report.html` |

---

## Reference Files

Default paths (hg38, dhslab group storage):

| Resource | Path |
|---|---|
| DRAGEN reference (DNA) | `/storage2/fs1/dspencer/Active/shared/refdata/hg38/dragen_hg38v4.3.6` |
| DRAGEN reference (5mC) | `/storage2/fs1/dspencer/Active/shared/refdata/hg38/dragen_hg38_5mCv4.3.6` |
| FASTA | `/storage2/fs1/dspencer/Active/shared/refdata/hg38/sequence/hg38_mgi_patch.fa` |
| dbSNP | `/storage2/fs1/dspencer/Active/shared/refdata/hg38/dragenfiles/dbsnp.vcf.gz` |
| VEP cache | `/storage2/fs1/dspencer/Active/shared/refdata/hg38/VEP113_cache` |
| Nirvana annotation | `/storage2/fs1/dspencer/Active/shared/refdata/hg38/dragenfiles/nirvana_annotation_data_323` |
| SNV noise (heme WGS) | `/storage2/fs1/dspencer/Active/shared/refdata/hg38/dragenfiles/IDPF_WGS_hg38_v.2.0.0_systematic_noise.snv.bed.gz` |
| SV noise (heme WGS) | `/storage2/fs1/dspencer/Active/shared/refdata/hg38/dragenfiles/WGS_FF_Heme_hg38_v3.1.0_systematic_noise.sv.bedpe.gz` |
| Population BAF VCF | `/storage2/fs1/dspencer/Active/shared/refdata/hg38/dragenfiles/1000G_phase1.snps.high_confidence.hg38.vcf.gz` |
| Ensembl GTF (RNA) | `/storage2/fs1/dspencer/Active/shared/refdata/hg38/ensemble/Homo_sapiens.GRCh38.105.chr.sorted.gtf.gz` |
| Cytobands | `assets/data/hg38.cytoBandIdeo.bed.gz` |
