#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import defaultdict
import sys
import numpy
import time
import psutil
import os
import threading
from Bio import SeqIO
from sys import argv


# k-mer indexing minimizers/winnowing
# some functions are taken from the labs or the mapper.py found in example_file1

class MinimizerWinnowIndex:
    def __init__(self, text, k, w):
        self.k = k
        self.w = w
        self.kmers = defaultdict(list)

        minimizer_positions = self.minimizers(text)

        for i in minimizer_positions:
            self.kmers[text[i:i+k]].append(i)


    def minimizers(self, s):
        """
        pozycje minimizerow w s - w kazdym oknie zawierajacym w k-merow, 
        ten z najmniejsza wartoscia funkcji hashujacej jest wtbierany 
        (najbardziej z lewej jezeli wartosci sa rownae)
        """
        k, w = self.k, self.w
        n = len(s) - k + 1
        if n < w:                      # sekwencja za krotka, nie ma zadnego okna
            return []
        h = numpy.fromiter((hash(s[i:i+k]) for i in range(n)), dtype=numpy.int64, count=n) #bierze kazdy kmer i haszuje go
        
        pos, last = [], -1
        for j in range(len(h)-w+1): #przechodzimy po kazdym oknie
            p = j + int(h[j:j+w].argmin()) #wyciagamy zahashowane kmery z okna

            if p != last:        # wybrana pozycja nigdy nie zmienia sie na lewo => nie potrzeba nam sortowania ani set'u
                pos.append(p)
                last = p

        return pos
 
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



class ReadMapper:
    def __init__(self, reference, index):
        self.reference = reference
        self.index = index

    def candidate_offsets(self, hits, tolerance=5):
        """
        Oblicza początek dopasowania danego reada w sekwencji referencyjnej i liczbę hitów w tej okolicy z dokładnością do zadanej tolerancji.
        Dla każdej takiej okolicy zwraca pierwsze i ostatnie dopasowanie oraz liczbę dopasowań.
        """

        if not hits:
            return []

        offsets = [ref_pos - read_pos for read_pos, ref_pos in hits]
        offsets.sort()


        common_offsets = [[offsets[0]]]

        for offset in offsets[1:]:
            if offset - common_offsets[-1][-1] <= tolerance:
                common_offsets[-1].append(offset)
            else:
                common_offsets.append([offset])


        candidates = []
        for common_offset in common_offsets:
            candidates.append((common_offset[0], common_offset[-1], len(common_offset)))

        return sorted(candidates, key=lambda x: x[2], reverse=True)

    def best_candidates(self, candidates, min_ratio=0.5):
        """
            wybieramy kandydata, dla którego jest najwięcej hitów oraz tych o zblizonej liczbie hitów, tzn. tych o ponad min_ratio*największa liczba hitów.
            max 5 kandydatów
        """
        if not candidates:
            return []

        best_support = candidates[0][2]

        selected_candidates = []

        for candidate in candidates:
            if candidate[2] >= min_ratio * best_support:
                selected_candidates.append(candidate)

            if len(selected_candidates) == 5:
                break


        return selected_candidates

    def candidate_ref_region(self, candidate, read_length, margin=100):

        offset_min, offset_max, support = candidate
        start = max(0, offset_min - margin)
        end = min(len(self.reference), offset_max + read_length + margin)

        return start, end


    def trace(self, D, read, reference_region, end):
        """Kod z labów: Odtwarza dopasowanie i zwraca jego początek oraz transcript."""
        i, j = len(read), end
        transcript = []

        while i > 0:
            diagonal = float("inf")
            vertical = float("inf")
            horizontal = float("inf")
            delta = None

            if i > 0 and j > 0:
                delta = 0 if read[i-1] == reference_region[j-1] else 1
                diagonal = D[i-1, j-1] + delta

            vertical = D[i-1, j] + 1

            if j > 0:
                horizontal = D[i, j-1] + 1

            if diagonal <= vertical and diagonal <= horizontal:
                transcript.append("R" if delta == 1 else "M")
                i -= 1
                j -= 1
            elif vertical <= horizontal:
                transcript.append("I")
                i -= 1
            else:
                transcript.append("D")
                j -= 1

        return j, "".join(transcript[::-1])

    def align(self, read, reference_region, lo, hi):
        """Znajduje najlepsze dopasowanie reada do fragmentu referencji."""
        n = len(read)
        m = len(reference_region)
        inf = 10**6

        rd = read.encode() #zmiana string na bajty, żeby łatwiej porównywa liczby w numpy
        ref = numpy.frombuffer(reference_region.encode(), dtype=numpy.uint8) #zmiana referencji na bajty
        neq = {c: (ref != c).astype(numpy.int32) for c in set(rd)} #dla kazdego znaku (klucza) z read tworzymy tabele z 1 jak ten sam znak w ref a 0 jak inny
        
        base = numpy.arange(m + 1, dtype=numpy.int32) #do przejsc poziomych??

        D = numpy.full((n + 1, m + 1), inf, dtype=numpy.int32)
        D[0, max(0, lo):min(m, hi) + 1] = 0 #0 tylko w rejonie, gdzie zakładamy początek dopasowania


        for i in range(1, n + 1):
            a = max(0, i + lo) #pas w ktorym obliczamy wartosci
            b = min(m, i + hi)

            if a > b:
                continue

            c = rd[i - 1]
            j0 = max(a, 1) #kolumna do ktorej przechodizmy, 0 oblsugiwane oddzielnie

            seg = numpy.minimum(D[i-1, j0-1:b] + neq[c][j0-1:b],   # diagonalnie - nowe wartosci w wycinku
                                D[i-1, j0:b+1] + 1)                 # pion

            row = D[i, a:b+1] #nasz wiersz

            if a == 0:
                row[0] = i
                row[1:] = seg
            else:
                row[:] = seg

            #przejscie poziome
            idx = base[:len(row)]
            row -= idx #odejmujemy wartosc indeksu od wartosci w wierszu?
            numpy.minimum.accumulate(row, out=row) #min z wartosci do tej pory od lewej, dla kazdego z zakresu
            row += idx

        # Znajdujemy minimum w ostatnim wierszu.
        best_end = int(numpy.argmin(D[n, :]))
        best_errors = int(D[n, best_end])

        # Odtwarzamy ścieżkę dopasowania.
        best_start, transcript = self.trace(
            D, read, reference_region, best_end
        )


        return best_errors, best_start, best_end, transcript


# DP algorithm 
seq_rec_list=[seq_record for seq_record in SeqIO.parse(argv[1], "fasta")]
reference = str(seq_rec_list[0].seq)


index = MinimizerWinnowIndex(reference, 20, 11)
mapper = ReadMapper(reference, index)
del seq_rec_list

fout = open(argv[3], "w")
reads = SeqIO.parse(argv[2], "fasta")

number_of_reads = 0


for read in reads:
    read_seq = str(read.seq)
    read_length = len(read_seq)
    number_of_reads += 1
    hits = index.query(read_seq)
    candidates = mapper.candidate_offsets(hits)
    best = mapper.best_candidates(candidates)

    best_alignment = None

    d = int(0.19 * read_length) + 1
    for candidate in best:
        start, end = mapper.candidate_ref_region(
            candidate, read_length, margin=d+5
        )

        reference_region = reference[start:end]

        lo = candidate[0] - start - d
        hi = candidate[1] - start + d

        errors, local_start, local_end, transcript = mapper.align(read_seq, reference_region, lo, hi)

        if best_alignment is None or errors < best_alignment[0]:
            best_alignment = (
                errors,
                start + local_start,
                start + local_end,
                transcript
            )


    if best_alignment is not None:
        errors, global_start, global_end, transcript = best_alignment

        if errors <= 0.19 * read_length:
            fout.write(
                f"{read.id}\t{global_start}\t{global_end}\n"
            )

fout.close()
