# Refines Transdecoder pan-transcriptome annotations using a rule-based approach

# Written 24 June 2022
# Updated 04 March 2022
# Written by Patrick Rynkiewicz
# License: GNU GPLv3

from statistics import mean
from Bio import SeqIO
import requests
import sys
import re
import simplejson
import time
import subprocess
import os.path

# Stores metadata associated with transdecoder sequence records


class Transdecoder_Record:
    def __init__(self, id, gene, isoform,
                 protein, ORFtype, tx_length, score,
                 uniref, seq):
        self.id = id                                        # Entire ID
        # Just Gene name (TRINITY_DN000_g0)
        self.gene = gene
        self.isoform = isoform                              # Just isoform (i0)
        self.protein = protein                              # Just protein (p0)
        # complete, internal, 5_prime_partial etc.
        self.ORFtype = ORFtype
        self.tx_length = tx_length                          # Transcript length
        # Transdecoder Annotation Score
        self.score = score
        # UniRef seed (UniRef100_P0000)
        self.uniref = uniref
        self.seq = seq

# Stores metadata associated with uniprot records


class UniProt_Record:
    def __init__(self, accession, organism, evidence_level, description, shortname, oln, gene_name):
        self.accession = accession                          # UniProt Accession
        self.organism = organism                            # Organism Name
        self.evidence_level = evidence_level                # Evidence Level per UniProt
        self.description = description                      # Description per UniProt
        self.shortname = shortname                          # Short Name for Protein
        self.oln = oln                                      # Ordered Locus Name
        self.gene_name = gene_name

# Stores BLAST alignment information


class BLAST_Alignment:
    def __init__(self, query, ref, evalue, length, identity):
        self.query = query          # Qname
        self.ref = ref              # Ref
        self.evalue = evalue        # Evalue if BLASTed
        self.length = length        # Qlen if BLASTed
        self.identity = identity    # Identity % if BLASTed

# Stores annotation record information, including rule-based annotation
# voting metadata.


class Annotation_Record:
    def __init__(self, TD_Record, UP_Record,
                 used_annotation,
                 num_alt, num_votes, num_alt_votes,
                 orig_protein_name, s288c_protein_name):
        self.TD_Record = TD_Record                          # TransDecoder Record
        self.UP_Record = UP_Record                          # UniProt Record
        self.used_annotation = used_annotation              # Original or S288C_BLAST
        # Number of alternative annotations
        self.num_alt = num_alt
        # Tally of annotations for representative protein
        self.num_votes = num_votes
        # Tally of annotations for alternative protein
        self.num_alt_votes = num_alt_votes
        self.orig_protein_name = orig_protein_name          # Original Protein Name
        # S288C Protein Name if BLASTed
        self.s288c_protein_name = s288c_protein_name


"""
DEBUG: Prints transdecoder record information in tabulated format

Inputs:
record
    - Transdecoder_Record object

Outputs:
    - Prints table
"""


def formatted_print_transdecoder_record(record):
    print(f"""{record.gene}\t{record.isoform}\t{record.protein}\
    \t{record.ORFtype}\t{record.tx_length}\t{record.score}\t{record.uniref}\t{record.seq}""")


"""
Parses FASTA record and extracts Trandecoder_Record information

Inputs:
fasta_record
    - FASTA record from SeqIO.parse() call

Outputs:
- Populated Transdecoder_Record for the protein record.
"""


def parse_fasta_record(fasta_record):
    desc = fasta_record.description
    whitespace_split = desc.split(" ")
    try:
        # Parse record from object description and extract metadata to
        # Transdecoder record.
        this_id = whitespace_split[0]
        this_gene = re.findall("(TRINITY\w+g[0-9]+)", this_id)[0]
        this_isoform = re.findall("(i[0-9]+)", this_id)[0]
        this_protein = re.findall("(p[0-9]+)", this_id)[0]
        this_ORFtype = whitespace_split[4][5:]
        this_length = whitespace_split[5][4:]
        this_score = whitespace_split[6].split(",")[1][6:]
        this_uniref = re.findall("(UniRef\w+)", whitespace_split[6])[0]
        this_seq = fasta_record.seq
        return(Transdecoder_Record(id=this_id,
                                   gene=this_gene,
                                   isoform=this_isoform,
                                   protein=this_protein,
                                   ORFtype=this_ORFtype,
                                   tx_length=this_length,
                                   score=this_score,
                                   uniref=this_uniref,
                                   seq=this_seq))
    except Exception:
        return("N/A")


"""
Parses FASTA file and assembles dictionary of Transdecoder_Records

Inputs:
fasta_file
    - String path to FASTA file

Outputs:
gene_dictionary
    - Dictionary wherein
        {
            TRINITY Gene Name : [Transdecoder_Record objects]
        }
"""


def pull_fasta_info(fasta_file):
    print("Parsing FASTA info...")
    gene_dictionary = {}
    num_records = 0
    num_parsed = 0
    # Parse FASTA record info as Transdecoder_Record objects.
    for record in SeqIO.parse(fasta_file, "fasta"):
        num_records = num_records + 1
        transdecoder_record = parse_fasta_record(record)
        # If there is a UniRef record associated with the annotation
        if not transdecoder_record == "N/A":
            num_parsed = num_parsed + 1
            # New gene, add dictionary record
            if transdecoder_record.gene not in gene_dictionary:
                gene_dictionary[transdecoder_record.gene] = [
                    transdecoder_record]
            # Gene exists in dictionary, add to record list values
            elif transdecoder_record.gene in gene_dictionary:
                gene_dictionary[transdecoder_record.gene].append(
                    transdecoder_record)
    print(f"\t{num_records} Total Records in FASTA File")
    print(f"\t{num_parsed} Parsed Records with UniProt Annotation")
    return(gene_dictionary)


"""
Collapses alternative gene annotations from a list into one representative

Inputs:
protein list 
    - list of Transdecoder_Record objects associated with a gene

Outputs:
output metadata
    - list containing
        - number of alternative protein annotations for the gene
        - number of votes (counts for the annotation)
        - number of votes for alternative annotations
protein_record
    - Transdecoder_Record object representing the best annotation for the gene

"""


def choose_representative_protein(protein_list):
    annotations = {}
    # Assemble annotations dictionary of {uniref : [scores]}
    for protein in protein_list:
        if protein.uniref not in annotations:
            annotations[protein.uniref] = [float(protein.score)]
        elif protein.uniref in annotations:
            score_list = annotations[protein.uniref]
            score_list.append(float(protein.score))
            annotations[protein.uniref] = score_list
    # One annotation present more than once. Return best.
    if len(annotations) == 1:
        for item in annotations:
            uniprot_id = item
            best_score = max(annotations[item])
            num_alt = 0
            num_votes = len(annotations)
            num_alt_votes = 0
            output_metadata = [num_alt, num_votes, num_alt_votes]
            for protein_record in protein_list:
                if protein_record.uniref == uniprot_id \
                        and protein_record.score == str(format(best_score, ".2f")):
                    return(output_metadata, protein_record)
    # More than one competing annotation
    elif len(annotations) > 1:
        # First and second best annotations found
        first = {}
        first_ref = ""
        second = {}
        second_ref = ""
        for uniref in annotations:
            # Retrieve annotation score list for this potential annotation
            score_list = annotations[uniref]
            # 1) No indicated most common annotation yet, automatically promote
            if len(first) == 0:
                first = {uniref: score_list}
                first_ref = uniref
            # 2) Length of score list is longer than best annotation, replace
            elif len(score_list) > len(first[first_ref]):
                # 2.1) No second, move first to second.
                if len(second) == 0:
                    second.clear()
                    second = {first_ref: first[first_ref]}
                    second_ref = first_ref
                # 2.2) First is longer than second, swap
                elif len(first[first_ref]) > len(second[second_ref]):
                    second.clear()
                    second = {first_ref: first[first_ref]}
                    second_ref = first_ref
                # 2.3) Current first is better than current second but has same length, swap
                elif len(first[first_ref]) == len(second[second_ref]) and mean(first[first_ref]) > mean(second[second_ref]):
                    second.clear()
                    second = {first_ref: first[first_ref]}
                    second_ref = first_ref
                first.clear()
                first = {uniref: score_list}
                first_ref = uniref
            # 3) Same length as first but better score, make first
            elif len(score_list) == len(first[first_ref]) and mean(score_list) > mean(first[first_ref]):
                # 3.1) No second, move first to second.
                if len(second) == 0:
                    second.clear()
                    second = {first_ref: first[first_ref]}
                    second_ref = first_ref
                # 3.2) First is longer than second, swap
                elif len(first[first_ref]) > len(second[second_ref]):
                    second.clear()
                    second = {first_ref: first[first_ref]}
                    second_ref = first_ref
                # 3.3) Current first is better than current second but has same length, swap
                elif len(first[first_ref]) == len(second[second_ref]) and mean(first[first_ref]) > mean(second[second_ref]):
                    second.clear()
                    second = {first_ref: first[first_ref]}
                    second_ref = first_ref
                first.clear()
                first = {uniref: score_list}
                first_ref = uniref
            # 4) Length of score list is shorter than best annotation, but no second exists, so make it second
            elif len(score_list) <= len(first[first_ref]) and len(second) == 0:
                second = {uniref: score_list}
                second_ref = uniref
            # 5) Length of score list is shorter than best annotation, but longer than second annotation. replace second.
            elif len(score_list) < len(first[first_ref]) and len(score_list) > len(second[second_ref]):
                second.clear()
                second = {uniref: score_list}
                second_ref = uniref
            # 6) Length of score list is the same as best annotation, and average score is better. replace best annotation
            elif len(score_list) == len(first[first_ref]) and mean(score_list) > mean(first[first_ref]):
                first.clear()
                first = {uniref: score_list}
                first_ref = uniref
            # 7) Length of score list is shorter than best annotation, but the same as second annotation.
            #    If average score is better, replace second annotation
            elif len(score_list) == len(second[second_ref]) and mean(score_list) > mean(second[second_ref]):
                second.clear()
                second = {uniref: score_list}
                second_ref = uniref
            # 8) Length of score list is longer than second best annotation, replace
            elif len(score_list) > len(second[second_ref]):
                second.clear()
                second = {uniref: score_list}
                second_ref = uniref
        overall_best = {}
        overall_best_ref = ""
        num_alt_votes = 0
        # If second ref has a better average annotation score but is one record shorter, make it the best.
        # Otherwise make first the best
        if mean(second[second_ref]) > mean(first[first_ref]) and (len(first[first_ref]) - 1) == len(second[second_ref]):
            overall_best = second
            overall_best_ref = second_ref
            num_alt_votes = len(first[first_ref])
        else:
            overall_best = first
            overall_best_ref = first_ref
            num_alt_votes = len(second[second_ref])
        # Assemble metadata
        num_alt = len(annotations) - 1
        num_votes = len(overall_best[overall_best_ref])
        best_score = max(overall_best[overall_best_ref])
        output_metadata = [num_alt, num_votes, num_alt_votes]
        # Output record
        for protein_record in protein_list:
            if protein_record.uniref == overall_best_ref \
                    and protein_record.score == str(format(best_score, ".2f")):
                return(output_metadata, protein_record)


"""
Collapses gene dictionary items to contain one annotation per gene

Inputs:
gene_dict
    - Dictionary wherein
        {
            TRINITY Gene Name : [Transdecoder_Record objects]
        }

Outputs:
collapsed_dict
    - Consolidated dictionary wherein one gene name is associated with one 
      Transdecoder annotation record and additional metadata can be stored
"""


def collapse_annotations(gene_dict):
    print("Choosing best annotations...")
    unannotated_misformatted_count = 0
    singletons = 0
    multiple = 0
    collapsed_dict = {}
    for gene in gene_dict:
        protein_list = gene_dict[gene]
        # Only one annotated protein, check UniProt seed protein for OLN
        if len(protein_list) == 1:
            protein_record = protein_list[0]
            # Metadata holds num_alt, num_votes, and num_alt_votes
            metadata = [0, 1, 0]
            singletons = singletons + 1
            collapsed_dict[protein_record] = metadata
        # Multiple potential annotations, rule-based approach
        elif len(protein_list) > 1:
            metadata, protein_record = choose_representative_protein(
                protein_list)
            multiple = multiple + 1
            collapsed_dict[protein_record] = metadata
        # No annotation, misformatted record
        else:
            unannotated_misformatted_count = unannotated_misformatted_count + 1
    # Print summary statistics and return dictionary
    print(f"\t{len(collapsed_dict)} / {len(gene_dict)} genes successfully annotated with one UniProt record")
    print(f"\t{singletons} Unambiguous Annotations")
    print(f"\t{multiple} Collapsed Annotations")
    print(f"\t{unannotated_misformatted_count} Misformatted and Unannotated Proteins")
    return(collapsed_dict)


"""
Maps chosen best protein to KEGG ID by obtaining an exact or proximate 
orderedLocusName (OLN) for the gene.

Inputs: 
protein_td_record
    - Transdecoder_Record object of chosen best representative protein for the
      Trinity gene

Outputs:
Dictionary with populated annotation information and metadata
"""


def map_uniprot(protein_td_record):
    uniref_seed_protein = protein_td_record.uniref.split("_")[1]
    # Try to annotate original accession
    out_dict = {}
    up_record_orig = uniprot_query(uniref_seed_protein)
    if len(up_record_orig.oln) == 0:
        alignment_obj = blast_s288c(protein_td_record)
        # If the alignment failed, just return the original object
        if alignment_obj == "":
            out_dict["Original_UniProt_Record"] = up_record_orig
        # If the alignment succeeded, return alignment object and both UP Objects
        else:
            up_record_s288c = uniprot_query(alignment_obj.ref)
            out_dict["BLAST_Alignment"] = alignment_obj
            out_dict["Original_UniProt_Record"] = up_record_orig
            out_dict["S288C_UniProt_Record"] = up_record_s288c
            return(out_dict)
    # If the original record had an associated OLN, return it.
    else:
        out_dict["Original_UniProt_Record"] = up_record_orig
        return(out_dict)
    # Failing case, return empty string
    return("")


"""
Queries UniProt REST API for gene OLN given seed protein ID

Inputs:
protein_id
    - UniProt ID derived from UniRef record. Seed protein for the
      UniRef100 cluster
type
    - "original" or "s288c_blast" based on query

Outputs:
UniProt_Record with as many completed fields as could be retrieved from API
"""


def uniprot_query(protein_id):
    header = "https://rest.uniprot.org/uniprotkb/"
    url = f"{header}{protein_id}"
    response_text = requests.get(url).text
    print(response_text)
    response_json = simplejson.loads(response_text)
    # UniProt Record Data Fields
    accession = ""
    organism = ""
    evidence_level = ""
    description = ""
    shortname = ""
    oln = ""
    gene_name = ""
    # Try to get fields
    try:
        accession = response_json["primaryAccession"]
    except KeyError as err:
        this_err = err
    try:
        organism = response_json["organism"]["scientificName"]
    except KeyError as err:
        this_err = err
    try:
        evidence_level = response_json["proteinExistence"]
    except KeyError as err:
        this_err = err
    try:
        description = response_json["proteinDescription"]["recommendedName"]["fullName"]["value"]
    except KeyError as err:
        this_err = err
    try:
        shortname = response_json["proteinDescription"]["recommendedName"]["shortNames"][0]["value"]
    except KeyError as err:
        this_err = err
    try:
        oln = response_json["genes"][0]["orderedLocusNames"][0]["value"]
    except KeyError as err:
        this_err = err
    try:
        gene_name = response_json["genes"][0]["geneName"]["value"]
    except KeyError as err:
        this_err = err
    return(UniProt_Record(accession=accession,
                          organism=organism,
                          evidence_level=evidence_level,
                          description=description,
                          shortname=shortname,
                          oln=oln,
                          gene_name=gene_name))


"""
Performs local BLASTP alignment to save S288C annotation information for
each annotated gene where applicable, since S288C tends to have the most robust
gene annotations

Inputs:
td_record
    - Transdecoder_Record for processing with annotation information

Outputs:
BLAST_Alignment record with alignment information and best match to S288C gene
"""


def blast_s288c(td_record):
    query_file = "annotation_query.fasta"
    if os.path.exists(query_file):
        os.remove(query_file)
    # Temporarily store query as a file for input to BLASTP
    with open(query_file, "w") as query_fh:
        query_fh.write(f">{td_record.id}\n")
        query_fh.write(str(td_record.seq))
    # Assemble blastp alignment command
    command_str = f"""
    blastp -query {query_file} -db {db_file} -max_target_seqs 1 -outfmt 6 -evalue 1e-5
    """
    #blast_outfmt6 = subprocess.check_output(command_str, shell=True)
    #blast_str = blast_outfmt6.decode('utf-8')

    with open(os.devnull, 'w') as devnull:
        blast_output = subprocess.check_output(command_str, shell=True, stderr=devnull, text=True)
        blast_str = blast_output.strip()

    # Parse BLAST alignment info to BLAST_Alignment record and output.
    if len(blast_str) > 0:
        blast_groups = blast_str.split("\t")
        ref = blast_groups[1]
        #uniprot_id = ref.split("|")[1]
        uniprot_id = ref
        print(ref)
        query = blast_groups[0]
        evalue = blast_groups[10]
        length = blast_groups[3]
        identity = blast_groups[2]
        return(BLAST_Alignment(query=query,
                               ref=uniprot_id,
                               evalue=evalue,
                               length=length,
                               identity=identity))
    else:
        return("")


"""
Compiles and outputs annotation information including all fields to TSV

Inputs:
collapsed_dict
    - dictionary containing annotation information for best annotation per gene
      as well as extended annotation information from BLAST/UniProt search

Outputs:
TSV file with all annotation information and chosen single best annotation per
gene

"""


def extend_annotations(collapsed_dict):
    print("Extending Annotations Using UniProtKB REST API...")
    original_annotations = 0
    s288c_BLAST_annotations = 0
    with open(outfile, "w") as outfh:
        # Output header to file
        outfh.write("""\
TrinityGene\tORFtype\tTranscript_Length\tAnnotation_Score\t\
UniRef\tUsed_Annotation\tOriginal_Accession\tOriginal_Organism\t\
Original_Evidence_Level\tOriginal_Description\tOriginal_Shortname\t\
Original_OLN_KEGG\tS288C_Accession\tS288C_Organism\tS288C_Evidence_Level\t\
S288C_Description\tS288C_Shortname\tS288C_OLN\tBLAST_Query\tBLAST_Reference\t\
BLAST_Evalue\tBLAST_Length\tBLAST_Identity\tNumber_Alternative_Proteins\t\
Number_Votes\tNumber_Alternative_Protein_Votes\tGene_Name_Orig\tGene_Name_S288C\n\
""")
        for annotation in collapsed_dict:
            td_record = annotation
            metadata = collapsed_dict[annotation]
            annotation_dict = map_uniprot(td_record)
            """
            ["BLAST_Alignment"]
            ["Original_UniProt_Record"]
            ["S288C_UniProt_Record"]
            """
            # If uniprot annotation didn't fail completely
            if not annotation_dict == "":
                # Output annotation with BLAST and both UP Records
                if len(annotation_dict) == 3:
                    s288c_BLAST_annotations = s288c_BLAST_annotations + 1
                    aln = annotation_dict["BLAST_Alignment"]
                    up_orig = annotation_dict["Original_UniProt_Record"]
                    up_S288C = annotation_dict["S288C_UniProt_Record"]
                    outfh.write(f"""\
{td_record.gene}\t{td_record.ORFtype}\t{td_record.tx_length}\t{td_record.score}\t\
{td_record.uniref}\tS288C\t{up_orig.accession}\t{up_orig.organism}\t{up_orig.evidence_level}\t\
{up_orig.description}\t{up_orig.shortname}\t{up_orig.oln}\t{up_S288C.accession}\t{up_S288C.organism}\t\
{up_S288C.evidence_level}\t{up_S288C.description}\t{up_S288C.shortname}\t{up_S288C.oln}\t\
{aln.query}\t{aln.ref}\t{aln.evalue}\t{aln.length}\t{aln.identity}\t{metadata[0]}\t{metadata[1]}\t\
{metadata[2]}\t{up_orig.gene_name}\t{up_S288C.gene_name}\n\
""")
                elif len(annotation_dict) == 1:
                    original_annotations = original_annotations + 1
                    up_orig = annotation_dict["Original_UniProt_Record"]
                    outfh.write(f"""\
{td_record.gene}\t{td_record.ORFtype}\t{td_record.tx_length}\t{td_record.score}\t\
{td_record.uniref}\tOriginal\t{up_orig.accession}\t{up_orig.organism}\t{up_orig.evidence_level}\t\
{up_orig.description}\t{up_orig.shortname}\t{up_orig.oln}\t\t\t\t\t\t\t\t\t\t\t\t\
{metadata[0]}\t{metadata[1]}\t{metadata[2]}\t{up_orig.gene_name}\t\n\
""")
    print(f"\t{original_annotations} Original Annotations Output")
    print(f"\t{s288c_BLAST_annotations} BLAST S288C Annotations Output")


def main():
    # Read in command-line arguments and define global variables
    global db_file
    global outfile
    fasta_file = sys.argv[1]
    db_file = sys.argv[2]
    outfile = sys.argv[3]
    # Extract FASTA record information from file
    gene_dict = pull_fasta_info(fasta_file)
    # Collapse annotations to best representative annotation per gene
    collapsed_annotations = collapse_annotations(gene_dict)
    # Compile and output TSV with best annotations and metadata
    extend_annotations(collapsed_annotations)
    # Cleanup -- remove temporary query file used for BLAST
    if os.path.exists("annotation_query.fasta"):
        os.remove("annotation_query.fasta")


if __name__ == "__main__":
    main()
