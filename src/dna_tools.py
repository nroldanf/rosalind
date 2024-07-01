import collections
from typing import Dict, List
from constants import *

def sequence_counter(sequence: str, vocabulary: List[str]) -> Dict[str, int]:    
    # Iterate and count the number of ocurrences of each token in
    # the vocabulary.
    sequence = sequence.upper()
    counts = {i:0 for i in vocabulary}
    for symbol in list(sequence):
        if not symbol in vocabulary:
            continue
        counts[symbol] += 1
    return counts

def sequence_counter_fast(sequence: str) -> Dict[str, int]:    
    # Iterate and count the number of ocurrences of each token in
    # the vocabulary.
    sequence = sequence.upper()
    return dict(collections.Counter(list(sequence)))

def dna_mrna_transcription(dna_sequence: str) -> str:
    # Transcribe DNA to mRNA
    # Note: the dna_sequence is consider to be the coding strand and not
    # the template strand.
    # same happens in biopython: https://biopython.org/wiki/Seq
    # In bioinformatics we normally assume the DNA is the coding strand (not the template strand) 
    # so this is a simple matter of replacing all the thymines with uracil:
    dna_sequence = dna_sequence.upper()
    return dna_sequence.replace("T", "U")

def dna_reverse_complement(dna_sequence: str) -> str:
    # Given the DNA strang 3'-5' (template)
    # return the reverse complement (the 5'-3' strand)
    # Reverse the sequence
    dna_sequence = dna_sequence[::-1]
    # Take the complement
    default_symbol = "*"
    return ''.join(dna_complements.get(symbol, default_symbol) for symbol in dna_sequence)

def reverse_transcription(mrna_sequence: str) -> str:
    # Given a mRNA sequence, use it as template to
    # generate a strand of DNA (cDNA)
    return "".join(reverse_transcription_mapping.get(symbol, "*") for symbol in mrna_sequence)

def translation(mrna_sequence: str, stop_symbol:str="*", to_stop: bool=False) -> str:
    # Translate mRNA into protein sequence
    # Heavily based on biopython's behavior
    current_index = 0
    protein = []
    # iterate the sequence in triplets
    while current_index < len(mrna_sequence) - 2:
        # Get the triplets
        codon = mrna_sequence[current_index:current_index+3]
        # Check if codong is a stop codon
        if codon in stop_codons:
            if to_stop:
                return "".join(protein)
            else:
                protein.append(stop_symbol)
                current_index += 3
                continue
        protein.append(std_codon_table.get(codon))
        current_index += 3
    return "".join(protein)

def get_number_rna_strings(protein_sequence: str) -> int:
    # Count the number of RNA strings that can be translated into the given protein sequence
    # This is a simple calculation based on the standard codon table
    # Note: This is not a real-world solution for large protein sequences
    # because the number of RNA strings grows exponentially with the length of the protein sequence.
    # However, for small protein sequences this should be sufficient.
    codons_of_interest_per_aa = 0
    for aa_in_protein in list(protein_sequence):
        # Filter the codons that codify for the given aa
        codons_of_interest = [codon for codon, aa_in_table in std_codon_table.items() if aa_in_protein == aa_in_table]
        # Multiple them for the previous count to obtain all the possible number of combinations
        # given the sequence
        if codons_of_interest_per_aa == 0:
            codons_of_interest_per_aa = len(codons_of_interest)
        else:
            codons_of_interest_per_aa *= len(codons_of_interest)
    # Finally, multiply by the number of possible stop codons that exists
    # to get the final number of possible combinations
    total_rna_strings = codons_of_interest_per_aa * len(stop_codons)
    # # Take the modulo 1M
    # total_rna_strings %= 10**6
    return total_rna_strings    

# https://biopython.org/docs/1.75/api/Bio.SeqUtils.html#Bio.SeqUtils.GC
def get_gc_content(dna_sequence: str) -> float:
    # Calculate GC content in the DNA sequence
    # i.e. porcentage of Guanine and Cytosine nucleobases
    dna_sequence = dna_sequence.upper()
    gc_count = sum(sequence_counter(dna_sequence, ["C", "G"]).values())
    total_count = len(dna_sequence)
    return (gc_count / total_count) * 100

# TODO: Switch to biopython
def read_fasta(filename: str) -> dict:
    # Read a FASTA file and return the sequences as a dictionary
    sequences = {}
    with open(filename, "r") as file:
        sequence_name = ""
        sequence_content = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if sequence_name:
                    # Assign the complete sequence to the sequence name
                    sequences[sequence_name] = sequence_content
                    # Clean up the sequence
                    sequence_content = ""
                # Assign the name at first iteration
                sequence_name = line[1:]
            else:
                # Accumulate the other part of the seuqences
                # note: usually fasta lines are no longer than 80 symbols
                sequence_content += line
    sequences[sequence_name] = sequence_content
    return sequences

# Mutations: types - deletion, duplication, inversion, insertion, translocation
# point mutations: change a base pair (changing one change the complement)

def hamming_distance(sequence_1: str, sequence_2: str) -> int:
    # Number of symbols that differ between two sequences
    # Note: Only considers sequences of equal length
    if len(sequence_1) == 0 or len(sequence_2) == 0:
        return  
    return len([symbol_1 for symbol_1, symbol_2 in zip(sequence_1, sequence_2) if symbol_1 != symbol_2])

if __name__ == "__main__":
    dna1 = "AGGACTATTAAGGAGCACCCGGGTGCAGTGGGCGTTTGTATCTTCCATCGCTAGGCTGCACGCCTGTTGTTTGCCGCTGGACGTTATACTAAGCTCACCGATTTATCCAGAAGCTTGCACGGGCCCATTGATCACCTTGGGCTACGAGTGGCCTCGACAAATAGTGATGCCCGACAAGGCACATATCCGGTGCCAGCGGAACTACCTCTACCGTGCAAGTAACCACAACGCAAGTTAGCGCAGGTATCTATGCTCGTCATAAGGGCAGGGTAACGGTAGCCGAGGCGGCAACTCTTTTGTAAGTGCGTGGAACGGCGATTAGAGCAGGCCCGCGTAACAAGCACTAGTCCGTGACGACACATAATGGCTGATAAGACTGCCTCGTGAAAGGACACAGGGGCCTTGAGTCCTCCGGGTGGTATCTATTCGTTTGCTATTAGGAGTCGGGCCCTGTCTCTCACTCGTACCGCTGCCTCCAGTCTTAGATTTCCCAGACGGAACGGGATAATGTCTCACGGGGAGATATAGTGTAAACTTAGGGGTCCCAGCACTTTCCTGCCTCACAAAATCGCTGTCGGCTGTCTACGTGTGTTTCTCCGGGAATACCTCGCTCACCGTATGAGAAATGTTGCTTTCAGTTTTGAGGTAGTGATTATGATGGGAAAGACAACTTTTAGTATAGTCTAGGCACTGTGATATCGAGAGGTGCAACCGAGGCGTCTAGGTTAGGAATCAGGGCCGCACGCCGGGTAACACTAGACCAGGAACGGTCACCTAGTCCCTCTTGCGAGTCGCACGTCGAAGCTTCAACTGTGCGAGAAGTTCATACATTAAAGAGATCCGTTAGAAGGTGCAACGTGCTGTCCGCTTAATATATCCGCAGGCGTCTGCGGTATGGCGAGCCGCTTCGGCTTTATTTATCA"
    dna2 = "CGAATTCTTTTCGAGTTTCTTGAAACAGGAAAATTTTAAGTCCATCGTCCCTCGCGTCCGTTTTGCTGGTTTGCCTCACTTCGTCAGACAAGCCTCGGCGATTTTGTCTCCAGCTGACCGAAGCGCAAGGATAACCGGAGGCAGGGCGGGACCACGGGCATCAGTCAAACTAGACAAGCTCCTTATCGCGAAAGAACTGGAAATACTCTTCCGTGCACTCAAGCAAATAGCAAATTAGCGCAGCTAACCAGTCTCTAGAAGTTGCACGTCTAACGCAAACCGTCACGAGCGCGCTTATACAGGTCTGCAATACGGCAGTGAGGGCGGCTGTGAACAGTTATCACTTGCTGCCGCGGAAGCACAACTGATAGTTAAACAACTACGTTGTAGCGGATCTGGGCCTTTATTACTTTCGGTGATTCCTCCTACTATTCCAACATGAGGCGGGGACTGGCGAATTATAGCAGGGCAGTTGCAGGACCTAGCTCTTCCATATTGGTATCTATGTGGACGCAAGGGGTAACAACGTCTAAACCTGTTTATGCTTTCGTTCGCCTGCAAAAGAAGATCTCCGCGTACCATCGCTCGCTTTTGGTCTGTGCTAACCGCCGTAACCGTATTAGAAATGTACCATCCAGAAATGCGCGCCATAGTCTGAAGTTACAGCTCACTGTTCGGATATGTTAGAGACTGTTCCGTCCACGCGTGCCCCAGATGCGCCGAGCCTCGTTATAATGGCAACACGTTGTGTATCCGACTAACGGGCATCGGTACACATTGACTGGCGAGACTCGCTCTTCGTTGCTTGAGCTCCGTACGTACGACAAGGATTAAAGATACTCAAGCCTAGCGGTAAAGGTTCTCCCTCTTAGTTGCGCCGCCGGAGACGGAGGATTGGTTAGCCGCTCTACTAGGTGTTGACA"
    print(hamming_distance(dna1, dna2))