#!/usr/bin/env python
# coding: utf-8

# In[1]:


from typing import List

def HammingDistance(s: str, t: str) -> int:
    """
    Calculate the Hamming distance between two DNA strings of equal length.
    
    Args:
        s (str): First DNA string
        t (str): Second DNA string
    
    Returns:
        int: Number of positions at which corresponding symbols differ
    
    Raises:
        ValueError: If strings are empty, unequal length, or contain invalid characters
    """
    if not s or not t:
        raise ValueError("Input strings cannot be empty")
    if len(s) != len(t):
        raise ValueError("Strings must be of equal length")
    
    valid_chars = set('ACGT')
    if not all(c in valid_chars for c in s) or not all(c in valid_chars for c in t):
        raise ValueError("Strings must contain only valid DNA nucleotides (A, C, G, T)")
        
    return sum(1 for i in range(len(s)) if s[i] != t[i])

def Translation(s: str) -> str:
    """
    Translate an RNA string into a protein string using the RNA codon table.
    
    Args:
        s (str): RNA string
    
    Returns:
        str: Protein string
    
    Raises:
        ValueError: If string is empty, invalid length, or contains invalid RNA nucleotides
    """
    if not s:
        raise ValueError("Input RNA string cannot be empty")
    
    if len(s) % 3 != 0:
        raise ValueError("RNA string length must be a multiple of 3")
    
    valid_chars = set('AUCG')
    if not all(c in valid_chars for c in s):
        raise ValueError("RNA string must contain only valid nucleotides (A, U, C, G)")
    
    CODON_TABLE = {
        'AUG': 'M', 'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'UCU': 'S',
        'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'UAU': 'Y', 'UAC': 'Y', 'UGU': 'C',
        'UGC': 'C', 'UGG': 'W', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAU': 'H', 'CAC': 'H',
        'CAA': 'Q', 'CAG': 'Q', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T',
        'ACG': 'T', 'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGU': 'S',
        'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V',
        'GUG': 'V', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAU': 'D',
        'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G',
        'GGG': 'G', 'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop'
    }
    
    protein = []
    start_codon_found = False
    
    for i in range(0, len(s), 3):
        codon = s[i:i+3]
        if len(codon) < 3:
            break
            
        if codon not in CODON_TABLE:
            raise ValueError(f"Invalid codon found: {codon}")
            
        amino_acid = CODON_TABLE[codon]
        
        if not start_codon_found:
            if codon == 'AUG':
                start_codon_found = True
            else:
                continue
                
        if amino_acid == 'Stop':
            break
        protein.append(amino_acid)
    
    if not start_codon_found:
        raise ValueError("No start codon (AUG) found in the sequence")
    
    return ''.join(protein)

def FindingMotif(s: str, t: str) -> List[int]:
    """
    Find all locations of a substring t in string s.
    
    Args:
        s (str): Main DNA string
        t (str): Motif to find
    
    Returns:
        List[int]: List of 1-based positions where t occurs in s
    
    Raises:
        ValueError: If strings are empty or contain invalid characters
    """
    if not s or not t:
        raise ValueError("Input strings cannot be empty")
    
    if len(t) > len(s):
        return []
        
    valid_chars = set('ACGT')
    if not all(c in valid_chars for c in s) or not all(c in valid_chars for c in t):
        raise ValueError("Strings must contain only valid DNA nucleotides (A, C, G, T)")
    
    positions = []
    for i in range(len(s) - len(t) + 1):
        if s[i:i+len(t)] == t:
            positions.append(i + 1)
    return positions

def RNASplicing(s: str, introns: List[str]) -> str:
    """
    Remove introns from DNA string and translate resulting RNA into protein.
    
    Args:
        s (str): DNA string
        introns (List[str]): List of intron sequences
    
    Returns:
        str: Resulting protein string
    
    Raises:
        ValueError: If input is invalid or contains invalid characters
    """
    if not s:
        raise ValueError("DNA string cannot be empty")
    
    if not isinstance(introns, list):
        raise ValueError("Introns must be provided as a list")
    
    valid_chars = set('ACGT')
    if not all(c in valid_chars for c in s):
        raise ValueError("DNA string must contain only valid nucleotides (A, C, G, T)")
    
    for intron in introns:
        if not intron or not all(c in valid_chars for c in intron):
            raise ValueError("All introns must be non-empty and contain valid DNA nucleotides")
        if len(intron) > len(s):
            raise ValueError("Intron cannot be longer than the DNA sequence")
    
    # Sort introns by length (descending) to handle nested introns correctly
    sorted_introns = sorted(introns, key=len, reverse=True)
    
    # Remove introns
    exon = s
    for intron in sorted_introns:
        exon = exon.replace(intron, '')
    
    if not exon:
        raise ValueError("After removing introns, no exon sequence remains")
    
    # Convert DNA to RNA
    rna = exon.replace('T', 'U')
    
    # Translate to protein
    return Translation(rna)

def LongestCommonSubstring(sequences: List[str]) -> str:
    """
    Find the longest common substring among DNA sequences.
    
    Args:
        sequences: List of DNA sequences
        
    Returns:
        str: The lexicographically smallest longest common substring
    
    Raises:
        ValueError: If input is invalid or contains invalid characters
    """
    if not sequences:
        raise ValueError("Input sequence list cannot be empty")
    
    if len(sequences) < 2:
        raise ValueError("At least two sequences are required")
    
    valid_chars = set('ACGT')
    for seq in sequences:
        if not seq:
            raise ValueError("Sequences cannot be empty")
        if not all(c in valid_chars for c in seq):
            raise ValueError("Sequences must contain only valid DNA nucleotides (A, C, G, T)")
    
    # Find the shortest string in the collection
    shortest_seq = min(sequences, key=len)
    n = len(shortest_seq)
    
    # Try substrings from longest to shortest
    for length in range(n, 0, -1):
        candidates = set()
        for start in range(n - length + 1):
            candidate = shortest_seq[start:start + length]
            if all(candidate in seq for seq in sequences):
                candidates.add(candidate)
        
        if candidates:
            return min(candidates)  # Return lexicographically smallest
            
    return ""

def FindingSubsequence(s: str, t: str) -> List[int]:
    """
    Find indices where string t appears as a subsequence in string s.
    Uses 1-based indexing.
    
    Args:
        s (str): The main DNA string to search in
        t (str): The subsequence to find
    
    Returns:
        List[int]: List of 1-based indices where characters of t appear in s
    
    Raises:
        ValueError: If input is invalid or contains invalid characters
    """
    if not s or not t:
        raise ValueError("Input strings cannot be empty")
    
    valid_chars = set('ACGT')
    if not all(c in valid_chars for c in s) or not all(c in valid_chars for c in t):
        raise ValueError("Strings must contain only valid DNA nucleotides (A, C, G, T)")
    
    if len(t) > len(s):
        return []
    
    result = []
    pos = 0
    
    for char in t:
        while pos < len(s) and s[pos] != char:
            pos += 1
            
        if pos >= len(s):
            return []
            
        result.append(pos + 1)
        pos += 1
        
    return result

def LongestCommonSubsequence(s: str, t: str) -> str:
    """
    Find the longest common subsequence of two DNA strings.
    
    Args:
        s (str): First DNA string
        t (str): Second DNA string
    
    Returns:
        str: Longest common subsequence
    
    Raises:
        ValueError: If input is invalid or contains invalid characters
    """
    if not s or not t:
        raise ValueError("Input strings cannot be empty")
    
    valid_chars = set('ACGT')
    if not all(c in valid_chars for c in s) or not all(c in valid_chars for c in t):
        raise ValueError("Strings must contain only valid DNA nucleotides (A, C, G, T)")
    
    m, n = len(s), len(t)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s[i - 1] == t[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
            else:
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1])
    
    # Reconstruct the LCS
    lcs = []
    i, j = m, n
    
    while i > 0 and j > 0:
        if s[i - 1] == t[j - 1]:
            lcs.append(s[i - 1])
            i -= 1
            j -= 1
        elif dp[i - 1][j] >= dp[i][j - 1]:
            i -= 1
        else:
            j -= 1
    
    return ''.join(reversed(lcs))

def ShortestCommonSupersequence(s: str, t: str) -> str:
    """
    Find a shortest common supersequence of two DNA strings.
    
    Args:
        s (str): First DNA string
        t (str): Second DNA string
    
    Returns:
        str: A shortest common supersequence
    
    Raises:
        ValueError: If input is invalid or contains invalid characters
    """
    if not s or not t:
        raise ValueError("Input strings cannot be empty")
    
    valid_chars = set('ACGT')
    if not all(c in valid_chars for c in s) or not all(c in valid_chars for c in t):
        raise ValueError("Strings must contain only valid DNA nucleotides (A, C, G, T)")
    
    m, n = len(s), len(t)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    backtrack = [[None] * (n + 1) for _ in range(m + 1)]
    
    # Fill dp table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s[i-1] == t[j-1]:
                dp[i][j] = dp[i-1][j-1] + 1
                backtrack[i][j] = 'match'
            else:
                if dp[i-1][j] >= dp[i][j-1]:
                    dp[i][j] = dp[i-1][j]
                    backtrack[i][j] = 'up'
                else:
                    dp[i][j] = dp[i][j-1]
                    backtrack[i][j] = 'left'
    
    # Build supersequence
    scs = []
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and backtrack[i][j] == 'match':
            scs.append(s[i-1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or backtrack[i][j] == 'up'):
            scs.append(s[i-1])
            i -= 1
        else:
            scs.append(t[j-1])
            j -= 1
    
    return ''.join(reversed(scs))

