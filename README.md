# Paired-End Exome Alignment and Variant Calling
Simple exome workflow for the Xiao lab.

There are 2 separate ways to perform this analysis:
1. Cluster (Drexel's Proteus)
2. Snakemake pipeline

## Cluster Workflow
* For this workflow, see: Cluster_Workflow.md

## Snakemake Pipeline
* For this workflow, the Snakefile and config.yaml are provided.
* Edit config.yaml for fastq files as input and reference fasta files.
* A subset of fastq files mapping to chr2 are provided, as is a chr2:200000-700000 fasta with BWA and Samtools indices.
* This pipeline assumes paired-end, Illumina reads and BWA index has been run on your reference of choice.
