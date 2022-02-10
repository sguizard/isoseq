# nf-core/isoseq: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* [CCS](#ccs) - Generate CCS sequences
* [LIMA](#lima) - Remove primer sequences from CCS
* [ISOSEQ REFINE](#isoseq-refine) - Select sequences with polyA tails and remove it
* [BAMTOOLS CONVERT](#bamtools-convert) - Convert bam file into fasta file
* [TAMA POLYA CLEAN UP](#tama-polya-clean-up) - Remove remaining polyA tails
* [ULTRA or MINIMAP2](#ultra-minimap2) - Map selected reads on genome
* [SAMTOOLS SORT](#samtools-sort) - Sort alignment and convert sam file into bam file
* [TAMA FILE LIST](#tama-file-list) - Prepare list file for TAMA collapse
* [TAMA COLLAPSE](#tama-collapse) - Clean gene models
* [TAMA MERGE](#tama-merge) - Merge all annotations into one for each sample with TAMA merge
* [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### CCS

<details markdown="1">
<summary>Output files</summary>

* `01_PBCCS/`
    * `70dpf_Liver.chunk9.bam`: The CCS sequences
    * `70dpf_Liver.chunk9.bam.pbi`: The Pacbio index of CCS files
    * `70dpf_Liver.chunk9.metrics.json.gz`: Statistics for each zmws
    * `70dpf_Liver.chunk9.report.json`: General statistics about generated CCS sequences in json format
    * `70dpf_Liver.chunk9.report.txt`: General statistics about generated CCS sequences in txt format

</details>

[CCS](https://github.com/PacificBiosciences/ccs) generate a Circular Consensus Sequence from subreads. It reports the number of selected and discarded zmws and the reason why.

### LIMA

<details markdown="1">
<summary>Output files</summary>

* `02_LIMA/`
    * `70dpf_Liver.chunk9_flnc.json`: Metadata about generated xml file
    * `70dpf_Liver.chunk9_flnc.lima.clips`: Clipped sequences
    * `70dpf_Liver.chunk9_flnc.lima.counts`: Statistics about detected primers pairs
    * `70dpf_Liver.chunk9_flnc.lima.guess`: Statistics about detected primers pairs
    * `70dpf_Liver.chunk9_flnc.lima.report`: Detailed statistics on primers pairs for each sequence
    * `70dpf_Liver.chunk9_flnc.lima.summary`: General statistics about selected and rejected sequences
    * `70dpf_Liver.chunk9_flnc.primer_5p--primer_3p.bam`: Selected sequences
    * `70dpf_Liver.chunk9_flnc.primer_5p--primer_3p.bam.pbi`: Pacbio index of selected sequences
    * `70dpf_Liver.chunk9_flnc.primer_5p--primer_3p.consensusreadset.xml`: Selected sequences metadata

</details>

[LIMA](https://github.com/pacificbiosciences/barcoding/) clean generated CCS. It selects sequences containing valid pairs of primers and removed it.

### ISOSEQ REFINE

<details markdown="1">
<summary>Output files</summary>

* `03_ISOSEQ3_REFINE/`
    * `70dpf_Liver.chunk9.bam`: Sequences sequences
    * `70dpf_Liver.chunk9.bam.pbi`: Pacbio index of selected sequences
    * `70dpf_Liver.chunk9.consensusreadset.xml`: Metadata
    * `70dpf_Liver.chunk9.filter_summary.json`: Number of Full Length, Full Length Non Chimeric, Full Length Non Chimeric PolyA
    * `70dpf_Liver.chunk9.report.csv`: Primers and insert length of each read

</details>

[ISOSEQ REFINE](https://github.com/PacificBiosciences/IsoSeq) select reads with polyA tails and discard polyA tail.

### BAMTOOLS CONVERT

<details markdown="1">
<summary>Output files</summary>

* `04_BAMTOOLS_CONVERT/`
    * `70dpf_Liver.chunk9.fasta`: The reads in fasta format.

</details>

[BAMTOOLS CONVERT](https://github.com/pezmaster31/bamtools) convert reads in BAM format into fasta format.

### TAMA POLYA CLEAN UP

<details markdown="1">
<summary>Output files</summary>

* `/`
    * `70dpf_Liver.chunk9_tama.fa`: The polyA tail free reads.
    * `70dpf_Liver.chunk9_polya_flnc_report.txt`: Length of removed tails.
    * `70dpf_Liver.chunk9_tama_tails.fa`: Sequence of removed tails.

</details>

[GSTAMA_POLYACLEANUP](https://github.com/GenomeRIK/tama) TAMA cleanup remove polyA tails from the selected reads.

### ULTRA or MINIMAP2

<details markdown="1">
<summary>Output files</summary>

* `06_ULTRA/` or `06_MINIMAP2/`
    * `70dpf_Liver.chunk9.sam`: The aligned reads.

</details>

[`MINIMAP2`](https://github.com/lh3/minimap2) or [`uLTRA`](https://github.com/ksahlin/ultra) aligns reads ont the genome.

### SAMTOOLS SORT

<details markdown="1">
<summary>Output files</summary>

* `08_SAMTOOLS_SORT/`
    * `70dpf_Liver.chunk9.bam`: The sorted aligned reads.

</details>

[SAMTOOLS SORT](http://www.htslib.org/doc/samtools-sort.html) sort the aligned reads and convert the sam file in bam file.

### TAMA COLLAPSE

<details markdown="1">
<summary>Output files</summary>

* `09_GSTAMA_COLLAPSE/`
    * `70dpf_Lung.chunk9_collapsed.bed`: This is a bed12 format file containing the final collapsed version of your transcriptome
    * `70dpf_Lung.chunk9_local_density_error.txt`: This file contains the log of filtering for local density error around the splice junctions
    * `70dpf_Lung.chunk9_polya.txt`: This file contains the reads with potential poly A truncation
    * `70dpf_Lung.chunk9_read.txt`: This file contains information for all mapped reads from the input SAM/BAM file.
    * `70dpf_Lung.chunk9_strand_check.txt`: This file shows instances where the sam flag strand information contrasted the GMAP strand information.
    * `70dpf_Lung.chunk9_trans_read.bed`: This file uses bed12 format to show the transcript model for each read based on the mapping prior to collapsing.This file uses bed12 format to show the transcript model for each read based on the mapping prior to collapsing.
    * `70dpf_Lung.chunk9_trans_report.txt`: This file contains collapsing information for each transcript
    * `70dpf_Lung.chunk9_varcov.txt`: This file contains the coverage information for each variant detected.
    * `70dpf_Lung.chunk9_variants.txt`: This file contains the variants called

</details>

[TAMA COLLAPSE](https://github.com/GenomeRIK/tama/wiki/Tama-Collapse) TAMA Collapse is a tool that allows you to collapse redundant transcript models in your Iso-Seq data.

### TAMA FILE LIST

<details markdown="1">
<summary>Output files</summary>

* `10_GSTAMA_FILELIST/`
    * `70dpf_Lung.tsv`: A tsv listing bed files to merge with TAMA merge

</details>

TAMA FILELIST is a home script for generating input file list for TAMA merge.

### TAMA MERGE

<details markdown="1">
<summary>Output files</summary>

* `11_GSTAMA_MERGE/`
    * `70dpf_Lung.bed`: This is the main merged annotation file.
    * `70dpf_Lung_gene_report.txt`: This contains a report of the genes from the merged file.
    * `70dpf_Lung_merge.txt`: This contains a bed12 format file which shows the coordinates of each input transcript matched to the merged transcript ID.
    * `70dpf_Lung_trans_report.txt`: This contains the source information for each merged transcript.

</details>

[TAMA MERGE](https://github.com/GenomeRIK/tama/wiki/Tama-Merge) TAMA Merge is a tool that allows you to merge multiple transcriptomes while maintaining source information.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

* `multiqc/`
    * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
    * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
    * `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

* `pipeline_info/`
    * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
    * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
    * Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
