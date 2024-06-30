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
    dna_sequence = dna_sequence.upper()
    return dna_sequence.replace("T", "U")

if __name__ == "__main__":
    sequence = "TTTTTGCAAAGTTAGCTGTCATAGGGCTTATGACCTCTACGACCGGGGACCAACTCTGCTTGAGCCCAATGCATAGGTTTTCCAGAATTGCATGTCCTCTATAAGTGGTTGTATTGTGGCGAAATCGCAGTTGTTTTATACGAGGATCTTCGTCTTTAATAGCTTCGGGATGTTTGTGGTGCCCATCGGAGTAGAACAGCGTTACCCGTGGATGACTGCCGTACGCCCTTGCCGACTGGGGTTCCCACTCTGTGTTCGTTCGTTCACCACAATCCGACTGGTATCGCTGGTCCGTGCAGTTATCAACGAAATAGCTGGTGTTGGCTGAGCAGACGATTTAACCAGTCGGGGTGGGCGCTAGTTCTCACCGGTATCTGGGAAGCCACTTGACGACGCGGTTGTTCCATGTAGACTCGTGGTTCTTATCCTGATGAGACGATCCGCTCTTTTTGGTGCCACCAAGACATACAACGCGGGGGTTATTACGGACCCCTTGCGAGGACAATAAGAGCCTGAGAGCGAAAAAAAGGACTGCCAGTTCTTTCTAACCCAATTAATTTTCTATTGGCTACAGTATAACGCGCCCAGAAGGGGCTGCTTCGGGGCCTGTGCTGATTTTCTTTCGTTACATTGATGGCCACTGATCGTATGCATTAATATTGAGGTCGGCCCCACCCCGCCCCCCACCGTGAAGCGCGTCCGCAAAGCCCTACACCTTCCGAAAGCTGTGCTCGCAGCCTGCCTATGGTTGTGCGGATGCACCAAAAAAGAACCAAAGACATCTTGGCAGTCCACGCCCCCGTGGAAGTTAACAAGAACGGGCTATATTCACCGTTTATGTTGCATGTCTTCGAGACACGGTAGATCCCTGAACCCATATAGTCGCAAAGTGGCGTGACGCGCCTAGGCTGACCGCAGAGTCCGTTCACTCTGACGCAGACCCGAGTGCAGGTTGTAATTCCAGGTGGAAA"
    print(dna_mrna_transcription(sequence))