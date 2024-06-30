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



if __name__ == "__main__":
    mrna = ""
    print(translation(mrna, to_stop=True))
    