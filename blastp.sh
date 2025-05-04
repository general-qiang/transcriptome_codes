#!/bin/bash -x

# blast against the entire database
# makeblastdb -in uniprot_sprot.fasta -dbtype prot -out uniprot_sprot
# makeblastdb -in uniprot_trembl.fasta -dbtype prot -out uniprot_trembl

# blastp -query Pata_ah2p.pep.fa -db uniprot_sprot -outfmt 6 -evalue 1e-5 -out blast_results.txt

# get a file with only gene ids and search result
# awk '{print $1 "\t" $2}' blast_results.txt | sort | uniq > gene2uniprot.txt


# awk 'NR==FNR {ids[$1]; next} ($1 in ids)' filtered_id.txt id_map_p.txt > filtered_map.txt





# blast against a selected database consisting of only ion-channels.

makeblastdb -in uniprotkb_ion_channel_AND_model_organis_2025_04_08.fasta -dbtype prot -out uniprot_ion_channel_db

blastp -query Pata_ah2p.pep.fa -db uniprot_ion_channel_db -out pat_ah2p_vs_uniprot_ion_channel.txt -evalue 1e-5 -outfmt 6
