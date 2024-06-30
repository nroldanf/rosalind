
# DNA vocab (nitrogenous bases)
dna_nucleobases = ["A", "C", "G", "T"]
# RNA vocab (nitrogenous bases)
rna_nucleobases = ["A", "C", "G", "U"]
# Complements: purines pairs with pyrimidines
dna_complements = {"A": "T", "T": "A", "C": "G", "G": "C"}
# reverse transcription mapping
reverse_transcription_mapping = {"U": "A", "A": "T", "G": "C", "C": "G"}

# codon tables https://github.com/biopython/biopython/blob/c2ba43db03a0e95593edfe42f1b335513b4889f3/Bio/Data/CodonTable.py#L603
# National Center for Biotechnology Information (NCBI) genetic codes: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
std_codon_table = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "UAU": "Y", "UAC": "Y",                          
    "UGU": "C", "UGC": "C",             "UGG": "W",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
stop_codons = [
    "UAA", "UAG", "UGA",
]