# Dokumentacja pipeline’u i danych

### A. Opis ogólny pipeline’u

  * Krótki opis projektu i celu pipeline’u (np. przygotowanie danych scRNA-seq z perturbacjami CRISPR do predykcji efektów gene knockdownów).
  * Schemat blokowy pipeline’u (np. w formie diagramu krok po kroku).
  * Wskazanie wersji oprogramowania i środowiska (np. Python 3.x, scanpy x.y, PyTorch z.z).
  * Informacje o źródłach danych (reference dataset, publiczne perturbation datasets).

### B. Szczegółowy opis każdego etapu

  * Opis działań wykonanych na danych (np. filtrowanie komórek, normalizacja, balansowanie).
  * Parametry użyte w każdym kroku (np. thresholdy filtrowania, target_sum do normalizacji).
  * Wykorzystywane narzędzia i funkcje (np. scanpy.pp.filter_cells(), scanpy.pp.normalize_total()).
  * Wyjaśnienie decyzji i kryteriów (dlaczego usuwamy komórki z MT > 15%, dlaczego wybrano dane publiczne itd.).
  * Efekty pośrednie i końcowe w postaci plików (np. train.h5ad, val.csv, public_pretrain.h5ad).

### C. Specyfikacja wejścia i wyjścia

  * Format i struktura plików wejściowych (np. AnnData .h5ad, CSV metadane).
  * Format plików wyjściowych po preprocessing (opis kolumn, typy danych, jakie metadane zawierają).
  * Założenia dotyczące struktury danych, np. rozmiar, liczba genów, skład metadanych.
  * Instrukcje dotyczące importerów danych do ML/AI (np. PyTorch DataLoadery).

### D. Instrukcja uruchamiania pipeline’u

  * Wskazówki, jak odtworzyć środowisko (np. plik conda environment.yml, Dockerfile).
  * Komendy uruchamiające poszczególne etapy pipeline’u.
  * Ścieżki do danych i plików konfiguracyjnych.
  * Informacje o seedach dla reproducibility.

### E. Opis danych i metadanych

  * Pełna lista kolumn metadanych w plikach obs (komórki) oraz var (geny).
  * Opis i format etykiet perturbacji (target_gene) i kontroli.
  * Informacje o dodatkowych cechach biologicznych w feature engineering.
  * Informacje o batchach, eksperymentach i innych zmiennych.

### F. Testy i sanity checks w pipeline

  * Opis testów jakościowych danych na każdym etapie.
  * Wyniki benchmarków bazowych z baseline models.
  * Raporty i wykresy opisujące efekty filtracji, normalizacji i integracji.

### G. Problemy i ograniczenia pipeline’u

  * Znane ograniczenia metody (np. trudności z imputacją brakujących genów).
  * Potencjalne pułapki i jak ich unikać (np. data leakage).
  * Propozycje dalszych usprawnień.
  
[powrót](../README.md)

