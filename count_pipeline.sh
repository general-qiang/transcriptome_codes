#!/bin/bash

# Change to the main directory
cd /Users/leitchlab/Documents

ANNOTATION="/Users/leitchlab/Documents/genome/Pata_MYAQb_rnm_ah2p.gtf"

# List of sample names (without _R1/_R2 or file extensions)
samples=("Blood-cell-pellet_S93_L006" "Brain_ant_S94_L006" "Brain_post_S95_L006" "Branchial-Pouch_S104_L006" "Egg_S97_L006" "Gut-ant_S100_L006" "Gut-post_S99_L006" "Heart_S101_L006" "Notochord_S103_L006" "Olfactory-bulb_S102_L006" "Sensory-tentacles_S96_L006" "Slime-Gland_S98_L006")

# Loop through each sample
for sample in "${samples[@]}"; do
    echo "Processing sample: $sample"

    # Define input and output file names
    SORTED_BAM="aligned_output/${sample}_sorted_aligned.bam" # Sorted BAM file
    COUNT_FILE="aligned_output/${sample}_counts.txt"         # Output count file
    
    # Run htseq
   featureCounts -T 4 -s 0 -a $ANNOTATION -t exon -g gene_id -o $COUNT_FILE -p $SORTED_BAM


    echo "Finished processing sample: $sample"
done



echo "All samples processed."
