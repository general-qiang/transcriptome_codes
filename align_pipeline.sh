#!/bin/bash

# Change to the directory containing the reference genome
cd /Users/leitchlab/Documents/genome

# Unzip the reference genome if it's compressed
gunzip Pata_MYAQb_rnm.fa.gz

# Build the HISAT2 index for the reference genome
hisat2-build Pata_MYAQb_rnm.fa Pata_MYAQb_rnm_index

# Change back to the main directory
cd /Users/leitchlab/Documents

# List of sample names (without _R1/_R2 or file extensions)
samples=("Blood-cell-pellet_S93_L006" "Brain_ant_S94_L006" "Brain_post_S95_L006" "Branchial-Pouch_S104_L006" "Egg_S97_L006" "Gut-ant_S100_L006" "Gut-post_S99_L006" "Heart_S101_L006" "Notochord_S103_L006" "Olfactory-bulb_S102_L006" "Sensory-tentacles_S96_L006" "Slime-Gland_S98_L006")

# Loop through each sample
for sample in "${samples[@]}"; do
    echo "Processing sample: $sample"

    # Define input and output file names
    R1="forwardreads/${sample}_R1_001.fastq"  # Forward read file
    R2="reversereads/${sample}_R2_001.fastq"  # Reverse read file
    SAM="${sample}_aligned.sam"               # Output SAM file
    BAM="${sample}_aligned.bam"               # Output BAM file
    SORTED_BAM="${sample}_sorted_aligned.bam" # Sorted BAM file

    # Run HISAT2 alignment
    hisat2 -p 8 -x genome/Pata_MYAQb_rnm_index -1 $R1 -2 $R2 -S $SAM

    # Convert SAM to BAM
    samtools view -Sb $SAM > $BAM

    # Sort the BAM file
    samtools sort -o $SORTED_BAM $BAM

    # Index the sorted BAM file
    samtools index $SORTED_BAM

    # Clean up intermediate files (optional)
    rm $SAM $BAM

    echo "Finished processing sample: $sample"
done

echo "All samples processed."
