#!/bin/bash
t=10
genome_dir="/mnt/e/tanmay_lele_19thmay"
trimmed_dir="trimmed_reads/ancestral"
alignment_dir="bismark_alignment/ancestral"
fastqc_dir="trimmed_fastqc_reports/ancestral"

mkdir -p "$trimmed_dir"
mkdir -p "$alignment_dir"
mkdir -p "$fastqc_dir"

for fq in raw_fastq/ancestral/*/*.fastq.gz; do
    sample=$(basename "$fq" .fastq.gz)
    sample_dir=$(dirname "$fq")
    echo "Processing $sample..."
    trim_galore --rrbs --gzip --cores 5 -o "$trimmed_dir" "$fq" #directional by default; reads had frequent CGG motifs and library is Msp1 digested library 
    trimmed_file="${trimmed_dir}/${sample}_trimmed.fq.gz"
    echo "Trimming done $sample..."
    fastqc -t $t -o "$fastqc_dir" "$trimmed_file"
    echo "QC done done $sample..."
    bismark --genome "$genome_dir" --bowtie2 -p $t -o "$alignment_dir" "$trimmed_file" #does end-to-end alignment by default, the library is directional as non-directional mode gives similar results
    echo "Alignment done $sample..."
    #do bismark methylation extractor step separately
done

