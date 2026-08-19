#/bin/bash
t=10
#find . -type f -name "*.fastq.gz" -exec fastqc -t $t -o fastqc_reports {} +
multiqc fastqc_reports -o multiqc_summary