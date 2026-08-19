#/bin/bash
#script is for preprocessing and aligning the atacseq reads
raw_fastq_dir="/mnt/e/tanmay_lele_19thmay/ATAC/Ancestral/fastq"
bowtie2_index="/mnt/d/sejal/GRCh38_noalt_as" #bowtie2 index human genome
trim_dir="/mnt/e/tanmay_lele_19thmay/ATAC/Ancestral/trimmed"
fastqc_dir="/mnt/e/tanmay_lele_19thmay/ATAC/Ancestral/fastqc"
bam_dir="/mnt/e/tanmay_lele_19thmay/ATAC/Ancestral/bam"
peak_dir="/mnt/e/tanmay_lele_19thmay/ATAC/Ancestral/peaks"
bw_dir="/mnt/e/tanmay_lele_19thmay/ATAC/Ancestral/bw"

num_t=10
for r1 in "$trim_dir"/*_R1_trim.fastq.gz; do
    baseprefix=$(basename "$r1" | sed 's/_R1_trim.fastq.gz//')
    r2="$raw_fastq_dir/${baseprefix}_R2_trim.fastq.gz" #change this to trim directory if needs to be run again
    echo "Processing $r1 and $r2"
    #fastp -i "$r1" -I "$r2" -o "$trim_dir/${baseprefix}_R1_trim.fastq.gz" -O "$trim_dir/${baseprefix}_R2_trim.fastq.gz" -w $num_t
    echo "Trimming done for $r1 and $r2"
    #rm "$raw_fastq_dir/${baseprefix}_R1_001.fastq.gz" "$raw_fastq_dir/${baseprefix}_R2_001.fastq.gz" #remove the raw fastq files
    #fastqc --threads $num_t -o "$fastqc_dir" "$trim_dir/${baseprefix}_R1_trim.fastq.gz" "$trim_dir/${baseprefix}_R2_trim.fastq.gz"
    bowtie2 --local -k 1 --no-mixed --no-discordant -x "${bowtie2_index}/GRCh38_noalt_as" -1 "$trim_dir/${baseprefix}_R1_trim.fastq.gz" -2 "$trim_dir/${baseprefix}_R2_trim.fastq.gz" --threads $num_t | samtools view -@ $num_t -bS - > "$bam_dir/${baseprefix}.bam"
    echo "Alignment done"
    samtools sort "$bam_dir/${baseprefix}.bam" -@ $num_t -o "$bam_dir/${baseprefix}_sorted.bam" 
    echo "Sorted the aligned file"
    samtools index "$bam_dir/${baseprefix}_sorted.bam" "$bam_dir/${baseprefix}_sorted.bam.bai" -@ $num_t
    rm "$bam_dir/${baseprefix}.bam" #remove the aligned bam file
    echo "Indexing done"
    samtools view -h "$bam_dir/${baseprefix}_sorted.bam" | grep -v chrM | samtools sort -@ $num_t -O bam -o "$bam_dir/${baseprefix}.rmChrM.bam" -T . 
    rm "$bam_dir/${baseprefix}_sorted.bam" "$bam_dir/${baseprefix}_sorted.bam.bai" #remove the sorted bam file
    samtools index -@ $num_t "$bam_dir/${baseprefix}.rmChrM.bam" "$bam_dir/${baseprefix}.rmChrM.bam.bai"
    picard AddOrReplaceReadGroups INPUT="$bam_dir/${baseprefix}.rmChrM.bam" OUTPUT="$bam_dir/${baseprefix}.rg.bam" RGID="${baseprefix}" RGLB="lib1" RGPL="ILLUMINA" RGPU="unit1" RGSM="${baseprefix}" VALIDATION_STRINGENCY=LENIENT
    picard MarkDuplicates QUIET=true INPUT="$bam_dir/${baseprefix}.rg.bam" OUTPUT="$bam_dir/${baseprefix}.marked.bam" METRICS_FILE="$bam_dir/${baseprefix}.dup.metrics" REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=.
    echo "Marked PCR and optical duplicates"
    rm "$bam_dir/${baseprefix}.rmChrM.bam" "$bam_dir/${baseprefix}.rmChrM.bam.bai" "$bam_dir/${baseprefix}.rg.bam"
    samtools view -@ $num_t -h -b -f 2 -F 1548 -q 30 "$bam_dir/${baseprefix}.marked.bam" | samtools sort -@ $num_t -o "$bam_dir/${baseprefix}.filtered.bam"
    samtools index -@ $num_t "$bam_dir/${baseprefix}.filtered.bam" "$bam_dir/${baseprefix}.filtered.bam.bai"
    echo "Removed duplicates, filtered, sorted and indexed"
    rm "$bam_dir/${baseprefix}.marked.bam" "$bam_dir/${baseprefix}.marked.bai"
    bedtools intersect -nonamecheck -v -abam "$bam_dir/${baseprefix}.filtered.bam" -b "$bam_dir/hg38-blacklist.v2.bed" > "$bam_dir/${baseprefix}.tmp.bam" #Remove reads within the blacklist regions
    echo "Removed blacklisted regions and created a tmp bam file"
    rm "$bam_dir/${baseprefix}.filtered.bam" "$bam_dir/${baseprefix}.filtered.bam.bai"
    samtools sort -@ $num_t "$bam_dir/${baseprefix}.tmp.bam" -o "$bam_dir/${baseprefix}.blacklist-filtered.bam"  #Sort and index the bam file
    samtools index -@ $num_t "$bam_dir/${baseprefix}.blacklist-filtered.bam" "$bam_dir/${baseprefix}.blacklist-filtered.bam.bai"
    rm "$bam_dir/${baseprefix}.tmp.bam"
    echo "Final sorted blacklist-filtered bam file ready!"
    bamCoverage --bam "$bam_dir/${baseprefix}.blacklist-filtered.bam" -o "$bw_dir/${baseprefix}.blacklist-filtered.bw" --numberOfProcessors $num_t --binSize 50 --normalizeUsing CPM --smoothLength 200 --effectiveGenomeSize 2913022398
    echo "Big-wig coverage files ready!"
    macs3 callpeak -t "$bam_dir/${baseprefix}.blacklist-filtered.bam" -f BAMPE -n "${baseprefix}" -q 0.01 --outdir "$peak_dir"
done


