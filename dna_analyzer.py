"DNA Sequence Analyzer Module"

from Bio import Entrez
import matplotlib.pyplot as plt 


                    #CODON TABLE
codon_table = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

Entrez.email = "youremail@gmail.com"

def clean_sequence(seq):
    seq = seq.strip()
    seq = seq.upper()
    seq = seq.replace(" ", "")
    seq = seq.replace("\n", "")
    return seq

def is_valid_dna(seq):
    valid_bases = "ATGC"
    for base in seq:
        if base not in valid_bases:
            return False 
    return True

def count_nucleotides(seq):
    counts = {}
    counts["A"] = seq.count("A")
    counts["T"] = seq.count("T")
    counts["G"] = seq.count("G")
    counts["C"] = seq.count("C")
    return counts

def gc_content(seq):
    g_count = seq.count("G")
    c_count = seq.count("C")
    total = len(seq)
    gc = (g_count + c_count) / total * 100
    return gc

def reverse_complement(seq):
    complement_map = {"A": "T", "T": "A", "G": "C", "C": "G"}
    complement = ""
    for base in seq: 
        complement = complement + complement_map[base]
    reversed_complement = complement[::-1]
    return reversed_complement

def transcribe(seq):
    seq = seq.replace("T", "U")
    return seq

def translate(seq, frame=0):
    protein = ""
    for i in range(frame, len(seq), 3):
        codon = seq[i:i+3]
        if len(codon) < 3:
            break
        amino_acid = codon_table[codon]
        if amino_acid == "*":
            break
        else:
            protein = protein + amino_acid 
    return protein

def read_fasta(filename):
    file = open(filename, "r")
    header = ""
    sequence = ""

    for line in file:
        line = line.strip()
        if line.startswith(">"):
            header = line
        else: 
            sequence = sequence + line 

    file.close()
    return header, sequence

def fetch_from_ncbi(accession_number):
    handle = Entrez.efetch(db = "nucleotide", id = accession_number, rettype = "fasta", retmode = "text")
    record = handle.read()
    handle.close()
    return record 

def plot_nucleotide_counts(counts):
    plt.style.use("fivethirtyeight")
    plt.bar (counts.keys(), counts.values())
    plt.xlabel ("Nucleotide")
    plt.ylabel ("Count")
    plt.title ("NUCLEOTIDE FREQUENCY")
    plt.show ()




def analyze(seq):
    cleaned = clean_sequence(seq)

    if not is_valid_dna(cleaned):
        print("Error: sequence contains invalid characters (only A, T, G, C allowed)")

    counts = count_nucleotides(cleaned)
    gc = gc_content(cleaned)
    rev_comp = reverse_complement(cleaned)
    rna = transcribe(cleaned)

    print("DNA SEQUENCE ANALYSIS REPORT")
    print("Sequence:", cleaned)
    print("Length:", len(cleaned), "bp")
    print("Nucleotide Counts:", counts)
    print("GC Content: {:.2f}%".format(gc))
    print("Reverse Complement:", rev_comp)
    print("RNA(transcribed):", rna)
    print("Protein(translated):", translate(rna))






# =========================================
#           TEST / DEMO SECTION
# =========================================

if __name__ == "__main__":

    # --- Fetch a real gene from NCBI and run the full pipeline ---
    result = fetch_from_ncbi("NM_000207")

    output_file = open("insulin.fasta", "w")
    output_file.write(result)
    output_file.close()

    header, sequence = read_fasta("insulin.fasta")
    analyze(sequence)

    # --- Multiple reading frame translation demo ---
    print(translate("AUGGCCUAG", 0))
    print(translate("AUGGCCUAG", 1))
    print(translate("AUGGCCUAG", 2))

    # --- Visualization ---
    counts = count_nucleotides("ATGGCCATTGTAATGGGCCGCTGA")
    plot_nucleotide_counts(counts)

    # --- Older stage tests, kept for reference ---
    #   analyze("GGCCTAAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC")
    #   header, sequence = read_fasta("test_sequence.fasta")
    #   print(header)
    # analyze(sequence)