#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import defaultdict
import sys
import numpy


# k-mer indexing minimizers/winnowing
# some functions are taken from the labs or the mapper.py found in example_file1

class MinimizerWinnowIndex:
    def __init__(self, text, k, w):
        self.k = k
        self.w = w
        self.kmers = defaultdict(list)      # k-mer -> positions in the reference
        for i in self.minimizers(text):
            self.kmers[text[i:i+k]].append(i)
 
    def minimizers(self, s):
        """
        pozycje minimizerow w s - w kazdym oknie zawierajacym w k-merow, 
        ten z najmniejsza wartoscia funkcji hashujacej jest wtbierany 
        (najbardziej z lewej jezeli wartosci sa rownae)
        """
        k, w = self.k, self.w
        h = [hash(s[i:i+k]) for i in range(len(s)-k+1)]
        pos = set()
        for j in range(len(h)-w+1):
            window = h[j:j+w]
            pos.add(j + window.index(min(window)))
        return sorted(pos)
 
    def query(self, read):
        """
        wszystkie pary (pozycja w odczycie, pozycja w sekwencji referencyjnej)
        gdzie minimizer z odczytu wystepuje w sekwencji referencyjnej
        """
        hits = []
        for i in self.minimizers(read):
            for p in self.kmers.get(read[i:i+self.k], []):
                hits.append((i, p))
        return hits

# choosing the clumps of kmers

# DP algorithm 

from Bio import SeqIO
from sys import argv

seq_rec_list=[seq_record for seq_record in SeqIO.parse(argv[1], "fasta")]
index = MinimizerWinnowIndex(str(seq_rec_list[0].seq), 20, 11)
del seq_rec_list

fout = open(argv[3], "w")
reads = SeqIO.parse(argv[2], "fasta")
for read in reads:
    hits = index.query(str(read.seq))

    fout.write("{}\t{}\n".format(read.id, len(hits)))
fout.close()



    
    
