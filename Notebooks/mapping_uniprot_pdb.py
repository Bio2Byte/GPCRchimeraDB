#using CIF
import requests
import json
import os
import numpy as np
from Bio import pairwise2, SeqIO, PDB
import json

Three_to_One_AA = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
    'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def pad_line(line):
    """Helper function to pad line to 80 characters in case it is shorter"""
    size_of_line = len(line)
    if size_of_line < 80:
        padding = 80 - size_of_line + 1
        line = line.strip('\n') + ' ' * padding + '\n'
    return line[:81]  # 80 + newline character

def extract_residues_pos_pdb(pdb_filename, chain_id):
    """Extract residue names and positions for a given chain in a PDB file."""
    
    # Initialize PDB parser
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_filename)
    
    residues = {}
    
    # Iterate through the structure to find the specified chain
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    if PDB.is_aa(residue, standard=True):  # Ignore non-amino acids
                        res_name = Three_to_One_AA[residue.get_resname().upper()]
                        res_id = residue.get_id()[1]  # Residue position
                        residues[res_id] = res_name
    return residues

    
def get_first_pos_pdb(pdb_file,chain_interest):
    #SIFT mapping identifies 1st residue from PDB in uniprot. But that 1st residue is not necessarly the first residue that is resolved in the PDB
    #Here check if the first residue that is resolved corresponds to a residue at the same position in sequence
    #if yes, no need to renumber

    fhandle = open(pdb_file, 'r')
    _pad_line = pad_line
    for line in fhandle:
        line = _pad_line(line)
        if line.startswith("ATOM") and line[21]==chain_interest:
            pos = int(line[24:27])
            return pos

    
def extract_pdb_sequence(pdb_filename, chain_id):
    
    """Extract the amino acid sequence of a given chain from a PDB file."""
    pdb_parser = PDB.PDBParser(QUIET=True)
    structure = pdb_parser.get_structure("pdb_structure", pdb_filename)

    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                sequence = []
                for residue in chain:
                    if PDB.is_aa(residue, standard=True):  # Ignore non-amino acid residues
                        sequence.append(Three_to_One_AA[residue.get_resname()])
                return "".join(sequence)

def count_aligned_residues(seq1, seq2):
    assert len(seq1) == len(seq2), "Sequences must be of equal length"
    well_aligned_count = sum(1 for a, b in zip(seq1, seq2) if a != '-' and b != '-' and a == b)
    return well_aligned_count
      
def align_uniprot_pdb(uniprot_seq, pdb_seq):
    """Align the UniProt and PDB sequences and find mismatches."""

    #depending on the sequences to align, glocal or local alignment can be better
    
    alignments_global = pairwise2.align.globalms(
        uniprot_seq, pdb_seq, match=2, mismatch=-1,
        open=-2, extend=-1,
        one_alignment_only=True
    )
    aligned_uniprot, aligned_pdb = alignments_global[0].seqA, alignments_global[0].seqB
    aligned_count_global = count_aligned_residues(aligned_uniprot, aligned_pdb)

    alignments_local = pairwise2.align.localms(
        uniprot_seq, pdb_seq, match=2, mismatch=-1,
        open=-2, extend=-1,
        one_alignment_only=True
    )
    aligned_uniprot, aligned_pdb = alignments_local[0].seqA, alignments_local[0].seqB
    aligned_count_local= count_aligned_residues(aligned_uniprot, aligned_pdb)

    if aligned_count_global>aligned_count_local:
        aligned_uniprot, aligned_pdb = alignments_global[0].seqA, alignments_global[0].seqB
    else:
        aligned_uniprot, aligned_pdb = alignments_local[0].seqA, alignments_local[0].seqB

    return aligned_uniprot, aligned_pdb

def find_mismatches_alignment_start(aligned_uniprot, aligned_pdb):
    found = False
    start_uniprot = 0
    start_pdb = 0
    for i in range(len(aligned_uniprot)):
        if not found:
            if aligned_uniprot[i] != "-":
                start_uniprot +=1
            if aligned_pdb[i] != "-":
                start_pdb +=1
            if aligned_uniprot[i] == aligned_pdb[i] and aligned_uniprot[i] != "-":
                found = True
                for j in range(1,6):
                    if aligned_uniprot[i+j] == aligned_pdb[i+j] and aligned_uniprot[i+j] != "-":
                        continue
                    else:
                        found = False
                        break
                if found:
                    idx_match_alignment = i
                    return start_uniprot, start_pdb, idx_match_alignment
                    
    return None, None, None  # No match found


def find_mismatches_alignment_stop(aligned_uniprot, aligned_pdb):
    found = False
    stop_uniprot = 0
    stop_pdb = 0
    for i in range(len(aligned_uniprot)-1,0,-1):
        if not found:
            if aligned_uniprot[i] != "-":
                stop_uniprot +=1
            if aligned_pdb[i] != "-":
                stop_pdb +=1
            if aligned_uniprot[i] == aligned_pdb[i] and aligned_uniprot[i] != "-":
                found = True
                for j in range(5,0,-1):
                    if aligned_uniprot[i-j] == aligned_pdb[i-j] and aligned_uniprot[i-j] != "-":
                        continue
                    else:
                        found = False
                        break
                if found:
                    idx_match_alignment = i
                    return stop_uniprot, stop_pdb, idx_match_alignment
                    
    return None, None, None  # No match found
        
def SIFT_mapping(pdb_id,uniprot_id):
    #find equivalence start uniprot and pdb 
    # https://www.ebi.ac.uk/pdbe/api/doc/sifts.html
    response = requests.get(f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}")
    data = response.json()
    uniprot_pdb_start = data[pdb_id.lower()]['UniProt'][uniprot_id]['mappings'][0]['unp_start'] #residue represented in PDB
    pdb_start = data[pdb_id.lower()]['UniProt'][uniprot_id]['mappings'][0]['start']['author_residue_number'] #corresponding residue in PDB
    pdb_idx = data[pdb_id.lower()]['UniProt'][uniprot_id]['mappings'][0]['start']['residue_number']
    return uniprot_pdb_start, pdb_start, pdb_idx

def map_MSA_seq(sequence_aligned):
    previous = 0
    translate = {}
    sequence_nogaps = sequence_aligned.replace("-","")
    for res in range(len(sequence_nogaps)):
        idx_msa = previous + sequence_aligned[previous:].index(sequence_nogaps[res])
        translate[idx_msa+1]=res+1
        previous = idx_msa + 1
    return translate

def match_uniprot_pdb_numbering(pdb_file,pdb_id,chain,uniprot_id,sequence,json_folder):
    #align sequence uniprot and PDB
    #if match found 5 residues in a row, renumber
    #else check if SIFT can help to find match
    modify = False
    add_uniprot =0
    pdb_seq = extract_pdb_sequence(pdb_file, chain)
    #sequence in pdb with residue number
    numbering_pdb = extract_residues_pos_pdb(pdb_file, chain)
    #start of uniprot seq in pdb and position 1st residue of PDB in uniprot seq
    aligned_uniprot, aligned_pdb=align_uniprot_pdb(sequence, pdb_seq)
    start_uniprot, start_pdb, idx_start_match_alignment = find_mismatches_alignment_start(aligned_uniprot, aligned_pdb)
    stop_uniprot, stop_pdb, idx_stop_match_alignment=find_mismatches_alignment_stop(aligned_uniprot, aligned_pdb)
    if start_uniprot == None: #when alignment is too messy (too many inserts into sequence) SIFT can sometimes help
        start_uniprot_sift, start_uniprot_in_pdb, idx_start_uniprot_in_pdb = SIFT_mapping(pdb_id,uniprot_id) #check if SIFT can help with renumbering
        try: #SIFT gives start PDB according to Uniprot nut also for unmodelled parts so residue not necessarly really in PDB
            pdb_seq_matching = pdb_seq[idx_start_uniprot_in_pdb:]
            aligned_uniprot, aligned_pdb=align_uniprot_pdb(sequence, pdb_seq_matching)
            start_uniprot, start_pdb, idx_start_match_alignment = find_mismatches_alignment_start(aligned_uniprot, aligned_pdb)
            if start_uniprot == None:
                raise
            stop_uniprot, stop_pdb, idx_stop_match_alignment=find_mismatches_alignment_stop(aligned_uniprot, aligned_pdb)
        except:
            #try if 1st residue of PDB is start of matching section
            difference = np.abs(start_uniprot_sift-idx_start_uniprot_in_pdb)
            idx_start_uniprot_in_pdb = 0
            start_uniprot_in_pdb=list(numbering_pdb)[0]
            aligned_uniprot, aligned_pdb=align_uniprot_pdb(sequence[start_uniprot_in_pdb-difference:], pdb_seq)
            start_uniprot, start_pdb, idx_start_match_alignment = find_mismatches_alignment_start(aligned_uniprot, aligned_pdb)
            start_uniprot = start_uniprot_in_pdb
            add_uniprot=start_uniprot_in_pdb-1
            stop_uniprot, stop_pdb, idx_stop_match_alignment=find_mismatches_alignment_stop(aligned_uniprot, aligned_pdb)
    else:
        #position in PDB of 1st residue uniprot
        start_uniprot_in_pdb = list(numbering_pdb.keys())[start_pdb-1]
        idx_start_uniprot_in_pdb=list(numbering_pdb).index(start_uniprot_in_pdb)

    #start_uniprot contains the position of the 1st residue that is elucidated in PDB
    #start_uniprot_in_pdb is the position of start_uniprot in the PDB
    uniprot_pdb_dict = {}
    idx_alignment = idx_start_match_alignment
    translate_MSA_seq = map_MSA_seq(aligned_uniprot)
    for uniprot_res_aligned, pdb_res_aligned in zip(aligned_uniprot[idx_start_match_alignment:], aligned_pdb[idx_start_match_alignment:]):
        if idx_alignment < idx_stop_match_alignment+1:
            if uniprot_res_aligned!="-" and pdb_res_aligned!="-": #because of mutations in pdb seq (add C-C bond) the residues can differ
                uniprot_pdb_dict[translate_MSA_seq[idx_alignment+1]+add_uniprot]=list(numbering_pdb.keys())[idx_start_uniprot_in_pdb]
            if pdb_res_aligned!="-":
                idx_start_uniprot_in_pdb +=1
            idx_alignment +=1

    result_file=json_folder+"/"+pdb_id+".json"
    json.dump(uniprot_pdb_dict,open(result_file,"w"),indent=1)
    return uniprot_pdb_dict

def access_mmCIF_numbering(pdb_id,uniprot,output_folder,type_gpcr):
    # Define the URL for the CIF file
    cif_url = f"https://www.ebi.ac.uk/pdbe/static/entry/{pdb_id.lower()}_updated.cif"

    # Fetch the CIF file content
    response = requests.get(cif_url)
    cif_content = response.text

    # Initialize a dictionary to store the UniProt to PDB residue mapping
    uniprot_to_pdb = {}
    found = False
    chain_found = ""
    first = True

    # Process the CIF file line by line
    for line in cif_content.splitlines():
        if line.startswith("ATOM"):
            columns = line.split()
            # if columns[6] == chain: #because sometimes you have 2 chain ids
            if type_gpcr =="chimera":
                if columns[-3] == uniprot.split('_')[0] or columns[-3] == uniprot.split('_')[1]: #chimeras
                    found = True
                    if first:
                        chain_found = columns[6]
                        first = False
            else:
                if uniprot in columns[-3]: #sometimes they add isoform number
                    found = True
                    if first:
                        chain_found = columns[6]
                        first = False
            if found and columns[6] == chain_found:
                # pdb_residue_number = int(columns[8])
                if "." in columns[16]:
                    print("Problem numbering with "+pdb_id, uniprot)
                    break 
                else:
                    pdb_residue_number = int(columns[16])
                uniprot_residue_number = int(columns[-2])
                uniprot_to_pdb[uniprot_residue_number] = pdb_residue_number
                found = False

    if len(uniprot_to_pdb)>0:
        # Save the mapping to a JSON file
        output_file = f"{pdb_id}.json"
        with open(output_folder+output_file, "w") as json_file:
            json.dump(uniprot_to_pdb, json_file, indent=4)
            return uniprot_to_pdb, True
    else:
        # print("Problem with "+pdb_id, uniprot)
        return uniprot_to_pdb, False

def map_PDB_uniprot(pdb_id,path_pdb,chain, uniprot,unaligned_seq, output_folder,type_gpcr):

    uniprot_pdb_dict, success=access_mmCIF_numbering(pdb_id,uniprot,output_folder,type_gpcr)
    
    if not success:
        uniprot_pdb_dict = match_uniprot_pdb_numbering(path_pdb,pdb_id,chain,uniprot,unaligned_seq,output_folder)

    return uniprot_pdb_dict