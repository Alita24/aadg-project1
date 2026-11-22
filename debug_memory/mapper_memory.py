#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import defaultdict
import numpy
from Bio import SeqIO
from sys import argv


# Część funkcji pochodzi z laboratoriów i z pliku mapper.py z example_file1.

MAX_ERROR_RATE = 0.2 #Ponieważ "input data may contain unmappable reads that cannot be mapped to the reference sequence with error rate < 20%"

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
        Pozycje minimizerow w s - w każdym oknie zawierającym w k-merów, 
        ten z najmniejszą wartością funkcji hashującej jest wybierany 
        (najbardziej z lewej jeżeli wartości sa równe)
        """
        k, w = self.k, self.w
        n = len(s) - k + 1
        if n < w:                           #sekwencja za krotka, nie ma zadnego okna
            return []
        h = numpy.fromiter((hash(s[i:i+k]) for i in range(n)), dtype=numpy.int64, count=n) #bierze kazdy kmer i haszuje go
        
        pos, last = [], -1
        for j in range(len(h)-w+1):         #przechodzimy po kazdym oknie
            p = j + int(h[j:j+w].argmin())  #wyciagamy zahashowane kmery z okna

            if p != last:                   # wybrana pozycja nigdy nie zmienia sie na lewo => nie potrzeba nam sortowania ani set'u
                pos.append(p)
                last = p

        return pos
 
    def query(self, read):
        """
        Wszystkie pary (pozycja w odczycie, pozycja w sekwencji referencyjnej)
        gdzie minimizer z odczytu wystepuje w sekwencji referencyjnej
        """
        hits = []
        for i in self.minimizers(read):
            for p in self.kmers.get(read[i:i+self.k], []):
                hits.append((i, p))


        return hits



class ReadMapper:
    def __init__(self, reference):
        self.reference = reference

    def candidate_offsets(self, hits, tolerance=5):
        """
        Oblicza początek dopasowania danego reada w sekwencji referencyjnej i liczbę hitów w tej okolicy z dokładnością do zadanej tolerancji.
        Dla każdej takiej okolicy zwraca pozycję pierwszego i ostatniego dopasowania oraz liczbę hitów, posortowane malejąco po hitach.
        """

        if not hits:
            return []

        offsets = [ref_pos - read_pos for read_pos, ref_pos in hits]
        offsets.sort()


        common_offsets = [[offsets[0]]] #pogrupowane podobne przesunięcia

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
        Wybiera kandydata, dla którego jest najwięcej hitów oraz tych o zblizonej liczbie hitów, tzn. tych o ponad min_ratio*największa liczba hitów.
        Max 5 kandydatów.
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


    def trace(self, D, read, reference_region, end):
        """Zmodyfikowany kod z labów: odtwarza dopasowanie i zwraca jego początek"""
        i, j = len(read), end

        while i > 0:
            diagonal = float("inf")
            horizontal = float("inf")
            vertical = float("inf")
            delta = None

            if i>0 and j > 0:
                delta = 0 if read[i-1] == reference_region[j-1] else 1
                diagonal = D[i-1, j-1] + delta

            vertical = D[i-1, j] + 1

            if j > 0:
                horizontal = D[i, j-1] + 1

            if diagonal <= vertical and diagonal <= horizontal:
                i -= 1
                j -= 1
            elif vertical <= horizontal:
                i -= 1
            else:
                j -= 1

        return j

    def align(self, read, reference_region, lo, hi):
        """Znajduje najlepsze dopasowanie reada do fragmentu referencji (liczone tylko w pasie ograniczonym przez lo i hi)."""
        n = len(read)
        m = len(reference_region)
        inf = 10**6

        rd = read.encode() #zmiana string na bajty, żeby numpy szybciej porównywalo znaki
        ref = numpy.frombuffer(reference_region.encode(), dtype=numpy.uint8) #zmiana referencji na bajty
        neq = {c: (ref != c).astype(numpy.int32) for c in set(rd)} #dla każdego znaku (klucza) z read tworzymy tabele z 1 jak ten sam znak w ref a 0 jak inny
        
        base = numpy.arange(m + 1, dtype=numpy.int32) #do przejsc poziomych

        D = numpy.full((n + 1, m + 1), inf, dtype=numpy.int32)
        D[0, max(0, lo):min(m, hi) + 1] = 0 #0 tylko w rejonie, gdzie zakładamy początek dopasowania


        for i in range(1, n + 1):
            a = max(0, i + lo) #pas w ktorym obliczamy wartosci
            b = min(m, i + hi)

            if a > b:
                continue

            c = rd[i - 1]
            j0 = max(a, 1) #kolumna do ktorej przechodzimy, 0 obsługiwane oddzielnie

            seg = numpy.minimum(D[i-1, j0-1:b] + neq[c][j0-1:b],   # diagonalnie
                                D[i-1, j0:b+1] + 1)                # pion

            row = D[i, a:b+1] # nasz wycinek wiersza

            if a == 0:
                row[0] = i
                row[1:] = seg
            else:
                row[:] = seg

            #przejscie poziome: D[i, j] = min_{k<=j} (D[i, k] + j - k) = j + min_{k<=j} (D[i, k] - k)
            idx = base[:len(row)]
            row -= idx #odejmujemy wartosc indeksu od wartości w wierszu (j)
            numpy.minimum.accumulate(row, out=row) # min kumulatywne od lewej min_{k<=j} (D[i, k] - k)
            row += idx # dodajemy indeks z powrotem

        best_end = int(numpy.argmin(D[n, :]))
        best_errors = int(D[n, best_end])

        if best_errors >= inf: 
            return None

        best_start = self.trace(D, read, reference_region, best_end)


        return best_errors, best_start, best_end


seq_rec_list=[seq_record for seq_record in SeqIO.parse(argv[1], "fasta")]
reference = str(seq_rec_list[0].seq)


index = MinimizerWinnowIndex(reference, 20, 11)
mapper = ReadMapper(reference)
del seq_rec_list

fout = open(argv[3], "w")
reads = SeqIO.parse(argv[2], "fasta")


for read in reads:
    read_seq = str(read.seq)
    read_length = len(read_seq)
    hits = index.query(read_seq)
    candidates = mapper.candidate_offsets(hits)
    best = mapper.best_candidates(candidates)

    best_alignment = None

    d = int(MAX_ERROR_RATE * read_length)

    alignments = [] # dopasowania ze wszystkich kandydatow tego odczytu
    for candidate in best:
        start = max(0, candidate[0] - d)
        end = min(len(reference), candidate[1] + read_length + d)
        reference_region = reference[start:end] # wycinek referencji, z ktorym porownujemy odczyt
        lo = candidate[0] - start - d
        hi = candidate[1] - start + d

        result = mapper.align(read_seq, reference_region, lo, hi)
        if result is None: # brak dopasowania => pomijamy odczyt
            continue
        errors, local_start, local_end = result
        alignments.append((errors, start + local_start, start + local_end))

    if not alignments:  # zaden kandydat nie dal dopasowania => odczyt zostaje niezmapowany
        continue

    alignments.sort()  # sortuje po liczbie bledow, najlepsze na poczatku
    errors, global_start, global_end = alignments[0]

    tolerance = 20 # start dalej niz 20 baz = inne miejsce w referencji
    margin = 0.01 * read_length # ile bledow wiecej wciaz uznajemy za "prawie tak samo dobre"
    ambiguous = any(    # czy jest odczyt w innym miejscu o podobnej liczbie bledow
        abs(other_start - global_start) > tolerance and other_errors <= errors + margin
        for other_errors, other_start, _ in alignments[1:]
    )

    if errors < MAX_ERROR_RATE * read_length and not ambiguous: # zapis tylko gdy bledow <20% i brak innego alignmentu w innym miejscu
        fout.write(f"{read.id}\t{global_start}\t{global_end}\n")

fout.close()
