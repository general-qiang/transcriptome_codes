# Transcriptome Codes

The workflow is designed as such: 
The raw reads from illumina RNA sequencing are first checked with fastQC and cleaned using trimmomatic to ge rid of adapters and poor quality reads. 
The trimmed reads are aligned to a genome (usually in fa.gz format) using HISAT2. The HISAT2 output would include a reading of alignment rate that reflects the overall quality of your sequence. 
The output BAM files are then processed with featurecounts to get counts of alignment. The counting output (txt or csv) is then used for differential expression analysis and gene ontology. 


The provided code is for using edgeR for DE analysis and sorting the output. 
