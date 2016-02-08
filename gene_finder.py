# -*- coding: utf-8 -*-
"""
Last updated: 1/31/16 9:00

This code file analyzes a DNA sequence and outputs snippets of DNA that are likely to be protein-coding genes

@author: Emma Price

"""

from load import load_seq
import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq

def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('T')
    'A'
    >>> get_complement('G')
    'C'
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    elif nucleotide  == 'T':
        return 'A'
    else:
        return None



def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    reversed_dna_string = ''
    dna_reverse = dna[::-1]
    for i in range(len(dna)):
        dna_nucleotide = dna_reverse[i]
        reversed_dna_string += get_complement(dna_nucleotide)
    return reversed_dna_string


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    end_dna = len(dna)
    for group in range(0, end_dna, 3):
        strand = dna[group:group+3]
        if strand in ['TGA', 'TAA', 'TAG']:
            return dna[0:group]
    return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    x = 3
    list_of_dna = []
    while (x < len(dna)):
        if dna[x-3:x] == 'ATG':
            current_dna = rest_of_ORF(dna[x-3:])
            list_of_dna.append(current_dna)
            x += len(current_dna)
        x += 3
    return list_of_dna


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    first_frame = find_all_ORFs_oneframe(dna[0:])
    second_frame = find_all_ORFs_oneframe(dna[1:])
    third_frame = find_all_ORFs_oneframe(dna[2:])
    list_of_frames = first_frame + second_frame + third_frame
    return list_of_frames
    

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    normal_dna = find_all_ORFs(dna)
    reverse_dna = find_all_ORFs(get_reverse_complement(dna))
    list_of_both_strands = normal_dna + reverse_dna
    return list_of_both_strands


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    i = 0
    both_ORFs = find_all_ORFs_both_strands(dna)
    longest = both_ORFs[0]
    while i + 1 < len(both_ORFs):
        if len(both_ORFs[i+1]) > len(longest):
            longest = both_ORFs[i+1]
        i+=1
    return longest

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        because this dna is shuffled, there is no way to run a doctest on it
        
        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    shuffled_dna = dna
    x = 0
    longest = 0
    while x < num_trials:
        shuffled_dna = shuffle_string(shuffled_dna)
        if len(str(longest_ORF(shuffled_dna))) > longest:
            longest = len(str(longest_ORF(shuffled_dna)))
        x +=1
    return longest

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents a protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    amino_acid = ''
    i = 0
    while i + 3 < len(dna) + 1:
        amino_acid += aa_table[dna[i:i+3]]
        i += 3
    return amino_acid


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna
        
        because this dna is shuffled, there is no way to run a doctest on it
        
        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 1500)
    list_of_aminos = []
    list_of_ORFs = find_all_ORFs_both_strands(dna)
    for an_ORF in list_of_ORFs:
        if an_ORF > threshold:
            amino_acid_ORF = coding_strand_to_AA(an_ORF)
            list_of_aminos.append(amino_acid_ORF)
    return list_of_aminos

if __name__ == "__main__":
    dna = load_seq("./data/X73525.fa")
    print gene_finder(dna)
#    import doctest
#    doctest.testmod()