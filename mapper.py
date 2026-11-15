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
        """
        100 bo 10% bledow w dlugosci 1000 czyli np. 100 delcji/insercji
        """

        offset_min, offset_max, support = candidate
        start = max(0, offset_min - margin)
        end = min(len(self.reference), offset_max + read_length + margin)

        return start, end


    def trace(self, D, read, reference_region, end):
        """Odtwarza dopasowanie i zwraca jego początek oraz transcript."""
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

    def align(self, read, reference_region):
        """Znajduje najlepsze dopasowanie reada do fragmentu referencji."""
        n = len(read)
        m = len(reference_region)

        D = numpy.zeros((n + 1, m + 1), dtype=int)


        # Pierwszy wiersz wypełniony zerami.
        # Pierwsza kolumna zawiera kolejne koszty usunięcia znaków.
        D[1:, 0] = numpy.arange(1, n + 1)

        # Wypełnianie macierzy DP.
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                delta = 0 if read[i-1] == reference_region[j-1] else 1

                D[i, j] = min(
                    D[i-1, j-1] + delta,
                    D[i-1, j] + 1,
                    D[i, j-1] + 1
                )

        # Znajdujemy minimum w ostatnim wierszu.
        best_end = int(numpy.argmin(D[n, :]))
        best_errors = int(D[n, best_end])

        # Odtwarzamy ścieżkę dopasowania.
        best_start, transcript = self.trace(
            D, read, reference_region, best_end
        )


        return best_errors, best_start, best_end, transcript




from Bio import SeqIO
from sys import argv

seq_rec_list=[seq_record for seq_record in SeqIO.parse(argv[1], "fasta")]
reference = str(seq_rec_list[0].seq)

index = MinimizerWinnowIndex(reference, 20, 11)

mapper = ReadMapper(reference, index)
del seq_rec_list

fout = open(argv[3], "w")
reads = SeqIO.parse(argv[2], "fasta")

number_of_reads = 0


for read in reads:
    number_of_reads += 1
    hits = index.query(str(read.seq))
    candidates = mapper.candidate_offsets(hits)
    best = mapper.best_candidates(candidates)

    best_alignment = None

    for candidate in best:
        start, end = mapper.candidate_ref_region(
            candidate, len(str(read.seq)), margin=120
        )

        reference_region = reference[start:end]

        errors, local_start, local_end, transcript = mapper.align(
            str(read.seq), reference_region
        )

        if best_alignment is None or errors < best_alignment[0]:
            best_alignment = (
                errors,
                start + local_start,
                start + local_end,
                transcript
            )


    if best_alignment is not None:
        errors, global_start, global_end, transcript = best_alignment

        if errors <= 0.19 * len(str(read.seq)):
            fout.write(
                f"{read.id}\t{global_start}\t{global_end}\n"
            )

fout.close()