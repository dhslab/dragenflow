# ![nf-core/dragenflow](docs/images/nf-core-dragenflow_logo_light.png#gh-light-mode-only) ![nf-core/dragenflow](docs/images/nf-core-dragenflow_logo_dark.png#gh-dark-mode-only)

## Introduction

**nf-core/dragenflow** is a bioinformatics pipeline that runs a variety of dragen commands and workflows for downstream analysis

## Pipeline Summary:
all workflows:
- samplesheet check
- make fasqlists (from reads/crams/bams)
- concatenate fastqlists
- run dragen

rna downstream analysis:
- get sizes file and strandedness
- annotate rnaseq
- bedtools genomecov
- ucsc bedclip, bedgraph to bigwig

tumor downstream analysis
- annotate small variants

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

## Usage

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->

### Samplesheet Format

>**If running with mgi samplesheet:**
>
>Pass flag --mgi true, and use Samplemap2.csv with the following columns:
>```csv
>FASTQ Path - Read 1,FASTQ Path - Read 2,Flowcell ID,Index Sequence,Flowcell Lane,ESP ID,Pool Name,Species,Illumina Sample Type,Library Type,Library Name,Date Complete,Total Reads,Total Bases,Avg >Q Score Read 1,Avg Q Score Read 2,% >Q30 Read 1,% >Q30 Read 2,PhiX Error Rate Read 1,PhiX Error Rate Read 2,% Pass Filter Clusters Read 1,% Pass Filter Clusters Read 2
>```
>
>**If running with custom samplesheet:**
>
>First column should be id, remaining columns are data type, or a combination of data types (read1,read2/bam/cram)
>
>Examples:
>```csv
>id,read1,read2
>sample1,sample1_R1.fastq.gz,sample1_R1.fastq.gz
>sample2,sample2_R1.fastq.gz,sample2_R1.fastq.gz
>```
>
>```csv
>id,bam,cram
>sample1,,sample1.cram
>sample2,sample2.bam,
>```

### Run Command

>```bash
>nextflow run dhslab/dragenflow -r dev \
   >-profile ris,<dragen2/dragen4/dragenaws> \
   >--input /path/to/samplesheet \
   >--outdir <OUTDIR> \ 
   >--workflow <rna/5mc/align/somatic/tumor/idtumis>
>```

### Optional Parameters
> --dragen_args \<dragen arguments> : provides additional arguments in dragen command
>
> --mgi true : pass if mgi samplesheet is used
>
> -bucket-dir \<path to s3 bucket dir> : pass if aws is used

### AWS
> run the following commands:
>
> nextflow secrets set AWS_ACCESS_KEY \<aws access key>
>
> nextflow secrets set AWS_SECRET_KEY \<aws secret key>
>
> export DRAGEN_USERNAME \<dragen username>
>
> export DRAGEN_PASSWORD \<dragen password>
>
> to check if secrets are set/exist in NXF_HOME, run:
>
> nextflow secrets list
