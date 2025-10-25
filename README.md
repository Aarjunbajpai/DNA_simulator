DNA Script Simulator A

Author: Aarjun
Language: Python
Level: Beginner → Intermediate (10th grade project)

Overview

This is a Python-based DNA script simulator that converts a DNA sequence into its RNA sequence, calculates reverse complements, analyzes nucleotide composition (AT/GC content), and translates the DNA sequence into a protein chain. It also counts the occurrence of each amino acid in the resulting protein.

The tool is designed to mimic basic bioinformatics lab operations in a simple command-line interface, giving users immediate feedback on DNA sequences.

Features

Validate DNA sequences for correct nucleotides (A, T, G, C).

Generate RNA sequence by replacing T with U.

Compute reverse complement for both DNA and RNA sequences.

Count individual nucleotides (A, T, G, C) and calculate AT% and GC%.

Translate DNA into a protein chain using standard codons.

Count the occurrence of each amino acid in the protein.

Handles incomplete codons and stop codons gracefully.

Example Usage
Enter the DNA sequence: ATGCGTAAACCC
RNA sequence: AUGCGUAAACCC
Reverse complement of DNA: GGGTTTACGCAT
Reverse complement of RNA: GGGUUUACGCAU
Protein chain: Methionine-Arginine-Lysine-Proline
Amino acid counts:
Methionine: 1
Arginine: 1
Lysine: 1
Proline: 1
AT%: 50.0
GC%: 50.0

How to Run

Make sure Python 3 is installed.

Download dna_script_simulator.py.

Open a terminal and run:

python dna_script_simulator.py

Future Improvements

Add file handling to process multiple sequences.

Create a graphical interface for easier use.

Include sequence alignment tools for comparing DNA sequences.

Add GC skew and other nucleotide statistics.
