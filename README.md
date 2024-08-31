# CerevisiAnnotate

CerevisiAnnotate is a Python script that refines Transdecoder pan-transcriptome annotations using a rule-based approach. It processes FASTA files containing Transdecoder outputs, performs UniProt queries, and conducts BLAST searches against a Saccharomyces cerevisiae S288C database to provide comprehensive annotations. It was originally designed for annotation of a Saccharomyces pan-transcriptome, but can with minor edits be applied to other systems as well.

CerevisiAnnotate was originally designed as part of a project to streamline Transdecoder annotation, assigning one gene name to one gene-level annotation (as defined by Transdecoder clustering). CerevisiAnnotate outputs metadata in a tabulated format to provide full traceability and transparency on how gene annotation is performed.

## Features

- Parses Transdecoder FASTA outputs
- Selects representative proteins for each gene
- Queries UniProt REST API for detailed protein information
- Performs BLAST searches against S288C database for additional annotations
- Generates a TSV file with consolidated annotation information

## Requirements

- Python 3.x
- Biopython
- requests
- simplejson
- BLAST+ (blastp command-line tool)

## Usage

```
python CerevisiAnnotate.py <input_fasta> <blast_db_file> <output_tsv>
```

- `<input_fasta>`: Path to the input FASTA file containing Transdecoder outputs
- `<blast_db_file>`: Path to the BLAST database file for S288C
- `<output_tsv>`: Path for the output TSV file containing annotations

## Output

The script generates a TSV file with the following columns:

- TrinityGene
- ORFtype
- Transcript_Length
- Annotation_Score
- UniRef
- Used_Annotation
- Original_Accession
- Original_Organism
- Original_Evidence_Level
- Original_Description
- Original_Shortname
- Original_OLN_KEGG
- S288C_Accession
- S288C_Organism
- S288C_Evidence_Level
- S288C_Description
- S288C_Shortname
- S288C_OLN
- BLAST_Query
- BLAST_Reference
- BLAST_Evalue
- BLAST_Length
- BLAST_Identity
- Number_Alternative_Proteins
- Number_Votes
- Number_Alternative_Protein_Votes
- Gene_Name_Orig
- Gene_Name_S288C

## Note

This script was written by Patrick Rynkiewicz and is licensed under GNU GPLv3.

## Last Updated

21 August 2024
