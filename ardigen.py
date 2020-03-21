#!/bin/python3

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# ---------------------------------- TASK 1 ----------------------------------
def open_seq(path):
    """
    :param path: Path to file with FASTA sequences
    :return: Content of a file
    """
    with open(path, 'r') as file:
        sequences = file.read()
    return sequences.split("\n")


def create_seq_dict(content):
    """
    :param content: Content of a opened file
    :return: Dictionary with sequence ID's and sequences
    """
    seq_dict = {}
    for element in content:
        if element.startswith(">"):
            tmp_key = element
            seq_dict[element] = ''
        else:
            seq_dict[tmp_key] += element
    return seq_dict


def get_correct_sequences(seq_dict):
    """
    :param seq_dict: Pre-check if sequences contain just allowed letters
    :return: Dictionary with sequences containing just allowed letters
    """
    correct_signs = "ATGC"
    correct_sequences = {}
    for seq_ID, sequence in seq_dict.items():
        correct = True
        for letter in sequence:
            if letter not in correct_signs:
                correct = False
        if correct:
            correct_sequences[seq_ID] = sequence
    return correct_sequences


def get_reverse_complementary(correct_sequences):
    """
    :param correct_sequences: Dictionary with correct DNA sequences
    :return: Dictionary with reverse-complementary sequences based on given dictionary with DNA sequences
    """
    complementarity = {"A": "T", "G": "C", "T": "A", "C": "G"}
    new_sequences = {}
    for sequence_ID in correct_sequences.keys():
        complementary_sequence = ''
        for letter in correct_sequences[sequence_ID]:
            complementary_sequence += complementarity[letter]
        new_sequences[sequence_ID] = complementary_sequence[::-1]
    return new_sequences


def get_mRNA(correct_sequences):
    """
    :param correct_sequences: Dictionary with correct DNA sequences
    :return: Dictionary with mRNA sequences based on given dictionary with DNA sequences
    """
    complementarity = {"A": "T", "G": "C", "T": "A", "C": "G"}
    new_sequences = {}
    for sequence_ID in correct_sequences.keys():
        mRNA = ''
        for letter in correct_sequences[sequence_ID]:
            mRNA += complementarity[letter]
        new_sequences[sequence_ID] = mRNA.replace("T", "U")
    return new_sequences


def translate(correct_sequences):
    """
    :param correct_sequences: Dictionary with correct DNA sequences
    :return: Dictionary with amino-acid sequences after translating given DNA/mRNA sequences
    """
    translated_dictionary = {}
    for sequence_ID in correct_sequences.keys():
        if "ATG" in correct_sequences[sequence_ID]:
            translated_seq = Seq(correct_sequences[sequence_ID], IUPAC.unambiguous_dna).translate()
            translated_dictionary[sequence_ID] = translated_seq
    return translated_dictionary


def save_to_fasta(complementary, mRNA, translated):
    """
    Saves output to three FASTA files

    :param complementary: Dictionary of complementary sequences with ID's
    :param mRNA: Dictionary of complementary sequences with ID's
    :param translated: Dictionary of complementary sequences with ID's
    """
    with open("complementary.fasta", "w") as fasta:
        for sequence_ID in complementary.keys():
            fasta.write(sequence_ID + "\n")
            fasta.write(complementary[sequence_ID])
            fasta.write("\n\n")

    with open("mRNA.fasta", "w") as fasta:
        for sequence_ID in mRNA.keys():
            fasta.write(sequence_ID + "\n")
            fasta.write(mRNA[sequence_ID])
            fasta.write("\n\n")

    with open("translated.fasta", "w") as fasta:
        for sequence_ID in translated.keys():
            fasta.write(sequence_ID + "\n")
            fasta.write(str(translated[sequence_ID]))
            fasta.write("\n\n")





if __name__ == "__main__":
    filename = 'sequences_zad2.fasta'
    seq = open_seq(filename)
    seq_dictionary = create_seq_dict(seq)
    correct_sequences_dictionary = get_correct_sequences(seq_dictionary)
    complementary_sequences_dictionary = get_reverse_complementary(correct_sequences_dictionary)
    mRNA_sequences_dictionary = get_mRNA(correct_sequences_dictionary)
    translated_sequences_dictionary = translate(correct_sequences_dictionary)
    save_to_fasta(complementary_sequences_dictionary, mRNA_sequences_dictionary, translated_sequences_dictionary)
