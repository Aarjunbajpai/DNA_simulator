# DNA SCRIPT SIMULATOR A
# Complete DNA codon → amino acid dictionary
DNA_codon_table = {
    'TTT': 'Phenylalanine', 'TTC': 'Phenylalanine',
    'TTA': 'Leucine', 'TTG': 'Leucine',
    'CTT': 'Leucine', 'CTC': 'Leucine', 'CTA': 'Leucine', 'CTG': 'Leucine',
    'ATT': 'Isoleucine', 'ATC': 'Isoleucine', 'ATA': 'Isoleucine',
    'ATG': 'Methionine',  # Start codon
    'GTT': 'Valine', 'GTC': 'Valine', 'GTA': 'Valine', 'GTG': 'Valine',
    'TCT': 'Serine', 'TCC': 'Serine', 'TCA': 'Serine', 'TCG': 'Serine',
    'CCT': 'Proline', 'CCC': 'Proline', 'CCA': 'Proline', 'CCG': 'Proline',
    'ACT': 'Threonine', 'ACC': 'Threonine', 'ACA': 'Threonine', 'ACG': 'Threonine',
    'GCT': 'Alanine', 'GCC': 'Alanine', 'GCA': 'Alanine', 'GCG': 'Alanine',
    'TAT': 'Tyrosine', 'TAC': 'Tyrosine',
    'TAA': 'Stop', 'TAG': 'Stop', 'TGA': 'Stop',  # Stop codons
    'CAT': 'Histidine', 'CAC': 'Histidine',
    'CAA': 'Glutamine', 'CAG': 'Glutamine',
    'AAT': 'Asparagine', 'AAC': 'Asparagine',
    'AAA': 'Lysine', 'AAG': 'Lysine',
    'GAT': 'Aspartic acid', 'GAC': 'Aspartic acid',
    'GAA': 'Glutamic acid', 'GAG': 'Glutamic acid',
    'TGT': 'Cysteine', 'TGC': 'Cysteine',
    'TGG': 'Tryptophan',
    'CGT': 'Arginine', 'CGC': 'Arginine', 'CGA': 'Arginine', 'CGG': 'Arginine',
    'AGT': 'Serine', 'AGC': 'Serine',
    'AGA': 'Arginine', 'AGG': 'Arginine',
    'GGT': 'Glycine', 'GGC': 'Glycine', 'GGA': 'Glycine', 'GGG': 'Glycine'
}

# Input DNA sequence
p = input("Enter the DNA sequence: ").upper()

# Validate sequence
for n in p:
    if n not in "ATGC":
        print("Invalid DNA sequence")
        break
else:
    # Transcribe to RNA
    RNA = p.replace("T", "U")
    print("\nRNA sequence:", RNA)

    # Reverse complement of DNA
    complement_DNA = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    rev_comp_DNA = ''.join(complement_DNA[b] for b in reversed(p))
    print("Reverse complement of DNA:", rev_comp_DNA)

    # Reverse complement of RNA
    complement_RNA = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
    rev_comp_RNA = ''.join(complement_RNA[b] for b in reversed(RNA))
    print("Reverse complement of RNA:", rev_comp_RNA)

    # Nucleotide counts
    a_count = p.count('A')
    t_count = p.count('T')
    g_count = p.count('G')
    c_count = p.count('C')
    total = len(p)
    print("\nNucleotide counts:")
    print("Adenine (A):", a_count)
    print("Thymine (T):", t_count)
    print("Guanine (G):", g_count)
    print("Cytosine (C):", c_count)
    print("Total nucleotides:", total)
    print("AT%:", round(((a_count + t_count) / total) * 100, 2))
    print("GC%:", round(((g_count + c_count) / total) * 100, 2))

    # Translate DNA to protein
    protein = []
    print("\nCodons and Amino Acids:")
    for x in range(0, len(p), 3):
        codon = p[x:x+3]
        if len(codon) < 3:
            continue  # skip incomplete codons
        amino_acid = DNA_codon_table.get(codon, "unknown")
        if amino_acid == "Stop":
            protein.append("[STOP]")
            continue
        print(codon, ":", amino_acid)
        protein.append(amino_acid)

    print("\nProtein chain:", "-".join(protein))

    # Count amino acids
    unique_amino_acids = set(protein)
    print("\nAmino acid counts:")
    for aa in unique_amino_acids:
        if aa == "[STOP]":
            continue
        print(f"{aa}: {protein.count(aa)}")
