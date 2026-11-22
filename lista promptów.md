# Lista promptów AI
## Klasa MinimizerWinnowIndex: [wklejony i zmodyfikowany kod z labów - kmers]
- Prośba o zmodyfikowanie kodu by efektywniejsze były procesy → redukcja ilości zmiennych tymczasowych. 
- Pytanie o przedstawienie najprostszych funkcji hashujących.
- Zgłoszenie, że minimizery zajmują za dużo pamięci, i zapytanie o możliwości redukcji pamięci → propozycje: array/NumPy zamiast listy hashy i set() w minimizers, większe w, a ostatecznie indeks na posortowanych tablicach NumPy zamiast słownika.
## Funkcja align:  [wklejony zmodyfikowany kod z labów - kEditDp]
- Prośba o przerobienie align tak, żeby znając przybliżone położenie reada nie liczyć pełnych macierzy, z naciskiem na krótszy czas → wersja z pasmami.
- Pytanie, jak zapisać program z użyciem numpy, żeby był szybszy → dodanie m.in. np.minimum.accumulate, np.frombuffer.

## Całość:
- Prośba o przejrzenie kodu pod kątem zbędnych elementów → usunięcie kilku diagnostycznych printów i importów