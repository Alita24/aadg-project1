
# Wyniki Testów - Readme

Ten dokument zawiera wyniki porównań oraz informacje o ilości miejsca zajmowanego przez pierwszą wersję mappera. Dostarcza wglądu w wydajność i efektywność procesu mapowania - na serwerze students.

## Zawartość
- **Wyniki Porównań**: Szczegółowe rezultaty testów porównawczych.
- **Użycie Miejsca**: Informacje o wymaganiach dotyczących przestrzeni dyskowej oryginalnego mappera.

Wykorzystaj te dane do analizy i optymalizacji procesu mapowania.

## Pliki w folderze

W folderze znajdują się następujące pliki generowane podczas testów:

- **summary.log**: Krótki podsumowujący raport wyników wszystkich uruchomień.
- **[example_dir]_[reads_file]_to_[reference_file].out**: Wyniki mapowania dla danego zestawu danych.
- **[example_dir]_[reads_file]_to_[reference_file].console.log**: Logi konsoli zawierające informacje o użyciu pamięci i ewentualnych błędach.
- **[example_dir]_[reads_file]_to_[reference_file].compare.txt**: Szczegółowe porównanie wyników mapowania z odpowiedziami referencyjnymi.

Każdy plik wynikowy jest nazwany zgodnie z formatem `[example_dir]_[reads_file]_to_[reference_file]`, aby uniknąć konfliktów nazw i ułatwić analizę.


## podsumowanie
konstrukcja minimizerow przy wiekszych plikach przekraczala dozwolona granice <1gb pamieci wykorzystanej.

### proponowane wytlumaczenia gdzie jest duzo pamieci jest zuzywane:
- pythonowa lista przetrzymujaca zahashowane kmery wraz z set() 
- lista `h` ma po jednym obiekcie int na kazdy k-mer referencji (~20M), co w logu dalo skok z 86 MB do 415 MB ("after hashes")
- `set()` z pozycjami minimizerow, a potem `sorted(pos)`, podniosl zuzycie do ok. 736 MB ("after sorting")
- lista i set istnieja jednoczesnie, dlatego szczyt wystepuje wewnatrz funkcji `minimizers`
- slownik z unikatowymi stringami k-merow wraz z listami dla kazdego klucza
- ok. 3,1 mln unikatowych kluczy, kazdy to osobny string (20 znakow) i osobna lista pozycji

### kolejne kroki
- zamiast duzej pythonowej listy obiektow `int` wprowadzamy niskopoziomowa tablice numpy.array (w C)
- wyrzucamy `set()` i `sorted()`, bo przesuwajac okno w prawo, pozycja najmniejszego elementu nigdy nie cofa sie w lewo — wystarczy sprawdzac `if p != last:`, zeby nie dopisywac dwa razy tego samego miejsca