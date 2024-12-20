# ![nf-core/dragenflow](docs/images/nf-core-dragenflow_logo_light.png#gh-light-mode-only) ![nf-core/dragenflow](docs/images/nf-core-dragenflow_logo_dark.png#gh-dark-mode-only)

## Introduction

**nf-core/dragenflow** is a bioinformatics pipeline that runs a variety of dragen commands and workflows for downstream analysis

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

:::warning
Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).
:::

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/dragenflow/usage) and the [parameter documentation](https://nf-co.re/dragenflow/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/dragenflow/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/dragenflow/output).

## Credits

nf-core/dragenflow was originally written by Nidhi.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#dragenflow` channel](https://nfcore.slack.com/channels/dragenflow) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/dragenflow for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
