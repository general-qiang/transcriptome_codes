#!/bin/bash
set -e  # Exit on error

# Change to the main directory
cd /Users/leitchlab/Documents

output_dir="trimmed_reads"
r1_output_dir="$output_dir/R1_trimmed"
r2_output_dir="$output_dir/R2_trimmed"

# Create output directories
# mkdir -p $r1_output_dir
# mkdir -p $r2_output_dir

# list of sample names
samples=("Blood-cell-pellet_S93_L006" "Brain_ant_S94_L006" "Brain_post_S95_L006" "Branchial-Pouch_S104_L006" "Egg_S97_L006" "Gut-ant_S100_L006" "Gut-post_S99_L006" "Heart_S101_L006" "Notochord_S103_L006" "Olfactory-bulb_S102_L006" "Sensory-tentacles_S96_L006" "Slime-Gland_S98_L006")

# loop through each sample
for sample in "${samples[@]}"; do
    echo "Processing sample: $sample"
    
    R1="forwardreads/${sample}_R1_001.fastq"  # Forward read file
    R2="reversereads/${sample}_R2_001.fastq"  # Reverse read file
    
    r1_paired="$output_dir/${sample}_R1_trimmed.fastq"
    r2_paired="$output_dir/${sample}_R2_trimmed.fastq"

    
#run trimmomatic
    java -jar /Users/leitchlab/Documents/Trimmomatic-0.40/dist/JAR/trimmomatic-0.40-rc1.jar PE \
    -threads 4 \
    $R1 $R2 \
    $r1_paired \
    $r2_paired \
    ILLUMINACLIP:/Users/leitchlab/Documents/Trimmomatic-0.40/adapters/TruSeq3-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 MINLEN:36

    
    echo "finished processing sample: $sample"
done

echo "All samples processed."
