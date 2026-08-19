#!/bin/bash
bam_dir="/mnt/e/tanmay_lele_19thmay/RRBS/fastq/bismark_alignment/bam_for_coverage_test"
meth_dir="/mnt/e/tanmay_lele_19thmay/RRBS/fastq/bismark_meth_extractor/11_09_25"

for bam_file in "$bam_dir"/*.bam; do
    sample_name=$(basename "$bam_file" .bam)
    echo "Processing $sample_name..."
    bismark_methylation_extractor -s --merge_non_CpG --comprehensive --bedGraph --parallel 3 --gzip -o "$meth_dir" "$bam_file"
    echo "Extraction completed for $sample_name"
    #rm "$meth_dir/Non_CpG_context_${sample_name}.txt" #Remove the non-cpg context file
    echo "Deleted Non_CpG_context_${sample_name}.txt"
done

