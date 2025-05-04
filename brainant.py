import csv
import numpy as np
from collections import defaultdict
import os
import re

# Clean and process UniProt IDs
def clean_uniprot_id(uniprot_id):
    """Remove unwanted characters from UniProt IDs and extract the primary ID."""
    if not uniprot_id or len(uniprot_id) < 3:  # Skip empty or too-short IDs
        return None
    # Remove unwanted characters and extract the primary ID
    uniprot_id = uniprot_id.replace('c("', '').replace('"', '').replace('(', '').strip()
    if '|' in uniprot_id:  # Handle IDs like "tr|A0A8M9QNH1|A0A8M9QNH1_DANRE"
        parts = uniprot_id.split('|')
        if len(parts) > 1:
            cleaned = parts[1]  # Extract the primary ID (e.g., A0A8M9QNH1)
        else:
            cleaned = parts[0]
    else:
        cleaned = uniprot_id
    return cleaned if cleaned and len(cleaned) > 2 else None

# Parse UniProt text format from local file
def parse_uniprot_entry(entry_text):
    """Parse a single UniProt entry from text into annotations."""
    annotations = {
        'entry_name': '',
        'accession': '',
        'protein_name': '',
        'gene_name': '',
        'function': [],
        'subcellular_location': [],
        'go_terms': []
    }

    # Split the entry into lines and process each line
    for line in entry_text.split('\n'):
        line = line.strip()  # Remove leading/trailing whitespace
        if line.startswith('ID  '):
            annotations['entry_name'] = line[5:].split()[0]
        elif line.startswith('AC  '):
            # Extract the first accession number
            annotations['accession'] = line[5:].split(';')[0].strip()
        elif line.startswith('DE   RecName: Full='):
            # Extract the full recommended name
            annotations['protein_name'] = line.split('Full=')[1].split(';')[0].strip()
        elif line.startswith('GN   Name='):
            # Extract the gene name
            annotations['gene_name'] = line.split('Name=')[1].split(';')[0].strip()
        elif line.startswith('CC   -!- FUNCTION:'):
            # Extract function descriptions
            function_text = line.split('FUNCTION:')[1].strip()
            annotations['function'].append(function_text)
        elif line.startswith('CC   -!- SUBCELLULAR LOCATION:'):
            # Extract subcellular location descriptions
            location_text = line.split('SUBCELLULAR LOCATION:')[1].strip()
            annotations['subcellular_location'].append(location_text)
        elif line.startswith('DR   GO;'):
            # Extract Gene Ontology (GO) terms
            parts = [p.strip() for p in line[5:].split(';')]
            if len(parts) >= 3:
                go_id = parts[1]
                go_term = parts[2].split(':')[1] if ':' in parts[2] else parts[2]
                annotations['go_terms'].append(f"{go_id}: {go_term}")

    return annotations

# Create a local UniProt database index
def create_uniprot_index(uniprot_file_path):
    """Create an index of accession numbers to file positions for quick lookup."""
    index = {}
    with open(uniprot_file_path, 'r') as f:
        current_pos = 0  # Track the position manually
        current_id = None
        while True:
            line = f.readline()
            if not line:  # End of file
                break
            if line.startswith('ID  '):
                current_id = line[5:].split()[0]
            elif line.startswith('AC  '):
                accessions = [acc.strip() for acc in line[5:].split(';') if acc.strip()]
                for acc in accessions:
                    index[acc] = current_pos
            elif line.startswith('//'):  # End of entry marker
                current_id = None
            current_pos = f.tell()  # Update position after reading the line
    return index

# Load UniProt database and create index
UNIPROT_DB_PATH = 'uniprotkb_ion_channel_AND_model_organis_2025_04_14.txt'
uniprot_index = create_uniprot_index(UNIPROT_DB_PATH)

# Fetch annotations from local UniProt database
def get_uniprot_annotations(uniprot_id):
    """Fetch functional annotations from local UniProt database for a given ID."""
    # Clean the ID first
    uniprot_id = clean_uniprot_id(uniprot_id)
    if not uniprot_id:
        print(f"Skipping invalid UniProt ID: {uniprot_id}")
        return None
    
    try:
        # Handle TrEMBL entries (tr|...) by extracting primary ID
        if uniprot_id.startswith('tr|'):
            primary_id = uniprot_id.split('|')[1]
            uniprot_id = primary_id
        
        # Look up in our index
        if (uniprot_id not in uniprot_index):
            print(f"UniProt ID not found in local database: {uniprot_id}")
            return None
            
        with open(UNIPROT_DB_PATH, 'r') as f:
            f.seek(uniprot_index[uniprot_id])
            entry_lines = []
            for line in f:
                if line.startswith('//'):  # End of entry
                    break
                entry_lines.append(line)
            entry_text = ''.join(entry_lines)
            return parse_uniprot_entry(entry_text)
            
    except Exception as e:
        print(f"Error processing {uniprot_id}: {str(e)}")
        return None

# Process BLAST results
blast_results = defaultdict(dict)

with open('pat_ah2p_vs_uniprot_ion_channel.txt', 'r', encoding='utf-8') as file:
    for line_num, line in enumerate(file, 1):
        row = line.strip().split()
        if len(row) < 12:
            print(f"Skipping incomplete BLAST row {line_num}")
            continue
        subject = clean_uniprot_id(row[1])  # Clean the subject ID
        if not subject:
            print(f"Skipping invalid subject ID in row {line_num}: {row[1]}")
            continue
        blast_results[subject]['custom_id'] = row[0]
        blast_results[subject]['perc_identity'] = float(row[2])
        blast_results[subject]['evalue'] = float(row[10])
        blast_results[subject]['bit_score'] = float(row[11])

print(f"First 10 BLAST results: {list(blast_results.items())[:10]}")  # Debugging

# Calculate BLAST index
def blast_index(perc_identity, evalue, bit_score, 
                w_identity=0.5, w_evalue=0.3, w_bit=0.2):
    try:
        e_transform = -np.log10(evalue + 1e-100)
        norm_identity = perc_identity / 100
        norm_evalue = e_transform / 0.00001  
        norm_bit = bit_score / 3439  
        index = (w_identity * norm_identity + 
                 w_evalue * norm_evalue + 
                 w_bit * norm_bit)
        print(f"Calculated BLAST index: {index} (perc_identity={perc_identity}, evalue={evalue}, bit_score={bit_score})")  # Debugging
        return index
    except Exception as e:
        print(f"Error calculating BLAST index: {e}")
        return None

for subject in blast_results:
    data = blast_results[subject]
    data['index'] = blast_index(data['perc_identity'], 
                               data['evalue'], 
                               data['bit_score'])

# Process DE data
DE_brainant = defaultdict(dict)

with open('DE_ic_sens_brainant.csv', 'r') as csvfile:
    csvreader = csv.reader(csvfile) # Skip header row
    try:
        header = next(csvreader)
    except StopIteration:
        print("Error: Empty DE file")
        exit(1)
    
    for row_num, row in enumerate(csvreader, 1):
        if len(row) < 8:
            print(f"Skipping incomplete DE row {row_num}")
            continue
        
        gene_id = row[1].strip('"\' ')
        logFC_str = row[2].strip('"\' ')
        pval_str = row[5].strip('"\' ')
        subject_ids = row[7].strip('"\' ')
        
        if not subject_ids:
            print(f"Row {row_num}: No subject_ids found for gene {gene_id}")  # Debugging
            continue
        
        try:
            DE_brainant[gene_id]['logFC'] = float(logFC_str)
            DE_brainant[gene_id]['pval'] = float(pval_str)
            cleaned_ids = []
            for s in subject_ids.split(','):
                cleaned_id = clean_uniprot_id(s.strip())
                cleaned_ids.append(cleaned_id)
            DE_brainant[gene_id]['subject_ids'] = [sid for sid in cleaned_ids if sid]
        except ValueError as e:
            print(f"Skipping gene '{gene_id}' (row {row_num}) due to invalid logFC/pval: {e}")
            print(f"Problematic values: logFC='{logFC_str}', pval='{pval_str}'")
            continue

# Add BLAST index and fetch UniProt annotations
for gene_id in DE_brainant:
    DE_brainant[gene_id]['blast_results'] = []
    DE_brainant[gene_id]['uniprot_annotations'] = []
    
    for subject_id in DE_brainant[gene_id]['subject_ids']:
        # Debugging: Check subject ID and BLAST keys
        print(f"Checking subject ID: {subject_id}")  # Debugging
        print(f"Available BLAST keys: {list(blast_results.keys())[:10]}")  # Debugging
        
        # Add BLAST index
        if subject_id in blast_results:
            blast_index_value = blast_results[subject_id]['index']
            print(f"Appending BLAST index for subject ID {subject_id}: {blast_index_value}")  # Debugging
            DE_brainant[gene_id]['blast_results'].append(blast_index_value)
        else:
            print(f"No BLAST result for subject ID: {subject_id}")  # Debugging
            DE_brainant[gene_id]['blast_results'].append(None)
        
        # Add UniProt annotations from local database
        annotations = get_uniprot_annotations(subject_id)
        DE_brainant[gene_id]['uniprot_annotations'].append(annotations)

# Debugging UniProt annotations
for gene_id in DE_brainant:
    for i, subject_id in enumerate(DE_brainant[gene_id]['subject_ids']):
        annotations = DE_brainant[gene_id]['uniprot_annotations'][i]
        if not annotations:
            print(f"No annotations for subject ID: {subject_id}")  # Debugging

# Prepare output data
output_data = []
for gene_id in DE_brainant:
    for i, subject_id in enumerate(DE_brainant[gene_id]['subject_ids']):
        annotations = DE_brainant[gene_id]['uniprot_annotations'][i]
        blast_index_value = DE_brainant[gene_id]['blast_results'][i]
        # print(f"Gene ID: {gene_id}, Subject ID: {subject_id}, BLAST Index: {blast_index_value}")  # Debugging
        if annotations:
            output_data.append({
                'gene_id': gene_id,
                'subject_id': subject_id,
                'logFC': DE_brainant[gene_id]['logFC'],
                'pval': DE_brainant[gene_id]['pval'],
                'blast_index': blast_index_value,
                'protein_name': annotations.get('protein_name', ''),
                'gene_name': annotations.get('gene_name', ''),
                'function': ' | '.join(annotations.get('function', [])),
                'subcellular_location': ' | '.join(annotations.get('subcellular_location', [])),
                'go_terms': ' | '.join(annotations.get('go_terms', []))
            })

# Save to CSV
output_fields = ['gene_id', 'subject_id', 'logFC', 'pval', 'blast_index', 
                'protein_name', 'gene_name', 'function', 
                'subcellular_location', 'go_terms']

with open('sens_brainant_complete_chart.csv', 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=output_fields)
    writer.writeheader()
    writer.writerows(output_data)

print("Processing complete. Results saved to sens_brainant_complete_chart.csv.")