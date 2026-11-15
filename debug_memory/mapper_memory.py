#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import defaultdict
import sys
import numpy
import time
import psutil
import os
import threading

process = psutil.Process(os.getpid())


max_memory = 0
max_memory_label = ""
monitoring = True

def set_stage(label):
    global current_stage
    current_stage = label

def monitor_memory():
    global max_memory, max_memory_label

    while monitoring:
        memory = process.memory_info().rss

        if memory > max_memory:
            max_memory = memory
            max_memory_label = current_stage

        time.sleep(0.01)

current_stage = "Starting program"

def print_memory(label):
    process = psutil.Process(os.getpid())
    memory_mb = process.memory_info().rss / (1024 * 1024)
    print(f"{label}: {memory_mb:.2f} MB")

start_time = time.perf_counter()
monitor = threading.Thread(target=monitor_memory)
monitor.start()

# k-mer indexing minimizers/winnowing
# some functions are taken from the labs or the mapper.py found in example_file1

class MinimizerWinnowIndex:
    def __init__(self, text, k, w):
        self.k = k
        self.w = w
        self.kmers = defaultdict(list)

        print_memory("Index: before minimizers")

        minimizer_positions = self.minimizers(text)

        print_memory("Index: after minimizers")
        print(f"Number of minimizers: {len(minimizer_positions)}")

        for i in minimizer_positions:
            self.kmers[text[i:i+k]].append(i)

        print_memory("Index: after dictionary construction")
        print(f"Number of unique k-mers: {len(self.kmers)}")
    """
    def __init__(self, text, k, w):
        self.k = k
        self.w = w
        self.kmers = defaultdict(list)      # k-mer -> positions in the reference

        print_memory("Index: before minimizers")
        for i in self.minimizers(text):
            self.kmers[text[i:i+k]].append(i)
    """
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
        print_memory("Minimizers: after hashes")
        
        pos, last = [], -1
        for j in range(len(h)-w+1): #przechodzimy po kazdym oknie
            window = h[j:j+w]
            p = j + int(h[j:j+w].argmin()) #wyciagamy zahashowane kmery z okna
            if p != last:        # wybrana pozycja nigdy nie zmienia sie na lewo => nie potrzeba nam sortowania ani set'u
                pos.append(p)
                last = p
        print_memory("Minimizers: after positions")

        # print_memory("Minimizers: after sorting")
        return pos
 
    def query(self, read):
        """
        wszystkie pary (pozycja w odczycie, pozycja w sekwencji referencyjnej)
        gdzie minimizer z odczytu wystepuje w sekwencji referencyjnej
        """
        set_stage("Querying read")
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

        set_stage("Allocating DP matrix")
        D = numpy.zeros((n + 1, m + 1), dtype=int)


        # Pierwszy wiersz wypełniony zerami.
        # Pierwsza kolumna zawiera kolejne koszty usunięcia znaków.
        D[1:, 0] = numpy.arange(1, n + 1)

        set_stage("Computing DP")
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


# choosing the clumps of kmers



# DP algorithm 

from Bio import SeqIO
from sys import argv

seq_rec_list=[seq_record for seq_record in SeqIO.parse(argv[1], "fasta")]
reference = str(seq_rec_list[0].seq)


set_stage("Building index")

index = MinimizerWinnowIndex(reference, 20, 11)



set_stage("Mapping reads")
mapper = ReadMapper(reference, index)
del seq_rec_list

fout = open(argv[3], "w")
reads = SeqIO.parse(argv[2], "fasta")

number_of_reads = 0

set_stage("Mapping reads")
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
end_time = time.perf_counter()


monitoring = False
monitor.join()

elapsed = end_time - start_time

memory_mb = max_memory / (1024 * 1024)

limit = 100 + 3 * number_of_reads
print(f"Runtime: {elapsed:.2f} seconds, limit: {limit}")
print(f"Time requirement: {'PASS' if elapsed <= limit else 'FAIL'}")
print(f"Maximum memory: {memory_mb:.2f} MB")
print(f"Memory requirement: {'PASS' if memory_mb < 1024 else 'FAIL'}")
print(f"Peak occurred during: {max_memory_label}")
