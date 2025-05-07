#!/bin/bash

# List of sample names (without _R1/_R2 or file extensions)
samples=("Blood-cell-pellet_S93_L006" "Brain_ant_S94_L006" "Brain_post_S95_L006" "Branchial-Pouch_S104_L006" "Egg_S97_L006" "Gut-ant_S100_L006" "Gut-post_S99_L006" "Heart_S101_L006" "Notochord_S103_L006" "Olfactory-bulb_S102_L006" "Sensory-tentacles_S96_L006" "Slime-Gland_S98_L006")

cd /Users/leitchlab/Documents/count_output

# Loop through each sample
for sample in "${samples[@]}"; do
    echo "Processing sample: $sample"

    COUNT="${sample}_counts.txt.summary"
    
    cat $COUNT
    
done

