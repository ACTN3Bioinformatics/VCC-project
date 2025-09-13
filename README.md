# VCC-2025

Virtual Cell Challenge 2025 - Pipeline przygotowania danych


## Opis projektu

Pipeline służy do przetwarzania danych single-cell RNA-seq z perturbacjami CRISPR w komórkach H1-hESC.
Celem jest przygotowanie wysokiej jakości, spójnych i zbalansowanych zestawów danych do trenowania modeli ML/AI przewidujących efekty genowych perturbacji.

## Struktura repozytorium i plików

- `/scripts/` - skrypty pipeline, m.in. filtrowanie, balansowanie, agregacja, integracja, splitowanie, feature engineering, benchmarki
- `/config/` - pliki konfiguracyjne do sterowania przebiegiem pipeline’u
- `/docs/` - dokumentacja rozszerzona, raporty jakości, tutoriale
- `/data/` - przykładowe dane wejściowe i wyjściowe
- `/reports/` - wyniki, benchmarki i analizy końcowe

## Struktura pipeline’u

1. QC i filtracja danych ([filter_normalize.py](./scripts/filter_normalize.py))
2. Normalizacja i log-transformacja ([filter_normalize.py](./scripts/filter_normalize.py))
3. Balansowanie klas perturbacji ([balance.py](./scripts/balance.py))
4. Agregacja pseudobulkowa (pseudobulk.py)
5. Integracja publicznych datasetów i korekta batchów (integration.py)
6. Tworzenie splitów i cross-validation (split_data.py)
7. Feature engineering biologiczny (feature_engineering.py)
8. Budowa DataLoaderów dla ML (dataloader.py)
9. Sanity checks i baseline benchmarking (benchmark.py)

---

### Liczba i charakterystyka datasetów w VCC 2025

Typy datasetów:

| Zbiór          | Opis                                                                | Przeznaczenie                   | Uwagi dotyczące preprocessing |
|--------------|-------------------------------------------------------------|--------------------------------|----------------------------------------|
| Training       | ~300k scRNA-seq komórek H1-hESC z perturbacjami CRISPRi           | Trening modeli ML              | Pełny preprocessing z QC, balansowaniem, batch correction |
| Validation     | Zbiór podobny do treningowego, niezależny                           | Walidacja pośrednia            | Analogiczny preprocessing jak training |
| Final Test     | Zbiór z perturbacjami nieobecnymi w training i validation          | Test końcowy, unseen genes     | Minimalny preprocessing, brak balansowania, limitowana ingerencja |
| Public         | Zbiory publiczne i zewnętrzne (Perturb-seq, CROP-seq itd.)          | Rozszerzenie treningu/pretraining | Konieczna integracja, standaryzacja, różne batch correction |


### Różnicowanie przygotowania datasetów w pipeline

Aby zapewnić odpowiednie przygotowanie każdego datasetu, pipeline stosuje **konfigurowalny profil** per zbiór, definiowany w pliku YAML: [datasets.yaml](./config/datasets.yaml).


---

### Dobre praktyki

[Instrukcja do dokumentacji pipeline’u i danych](./docs/Pipeline_instr.md)

Środowisko conda - [plik environment.yml](./environment.yml)

Konfiguracja pipelinu - [plik config.yaml](./config/config.yaml)

Konfiguracja pipelinu i profilu datasetu - [datasets.yaml](./config/datasets.yaml)

---

## Instrukcja uruchamiania

```
### Environment setup

conda env create -f environment.yml
conda activate vcc2025

python scripts/filter_normalize.py --dataset_name training --config config/datasets.yaml

python scripts/balance.py --dataset_name training --config config/datasets.yaml
```

## Szczególowy opis pipeline'u

### 0. Pozyskanie i wstępne zapoznanie się z danymi

Pobranie danych: z oficjalnego repozytorium VCC 2025 (reference dataset, training set) w formacie AnnData (.h5ad) oraz dodatkowych public perturbation datasets (np. Perturb-seq).

Eksploracja danych:

  * Wczytanie danych w Pythonie za pomocą scanpy.read_h5ad().
  * Przegląd metadanych (obs = komórki, var = geny, X = macierz ekspresji).
  * Sprawdzenie wymiarów danych (np. liczba komórek, genów).
  * Przegląd warunków eksperymentalnych (perturbacje, kontrola).
  * Podstawowe podsumowania statystyczne dotyczące liczby UMIs, ekspresji i rozkładów.

Cel: uzyskanie świadomości co do jakości i struktury danych, co pozwoli zaplanować dalszy preprocessing.

### 1. Jakość i filtracja danych (Quality Control - QC)

Kontrola jakości komórek:

  * Liczba genów wykrytych w każdej komórce (n_genes_by_counts).
  * Suma UMIs na komórkę (total_counts).
  * Procent genów mitochondrialnych (np. geny rozpoczynające się na "MT-") — oznacza zdrowie/żywotność.
  * Detekcja dubletów

Filtracja:

  * Usuń komórki z niskim n_genes_by_counts — np. poniżej 200 genów (typowo usunięcie martwych/słabych komórek).
  * Usuń komórki o wysokim procencie ekspresji mitochondrialnej (np. >10-15%) — potencjalne artefakty.

Opcjonalnie: filtracja genów występujących w bardzo niewielu komórkach (np. mniej niż 3).

Narzędzia:

  * scanpy.pp.calculate_qc_metrics()
  * scanpy.pp.filter_cells(min_genes=200, max_genes=None, max_mt_pct=15)
  * Graficzne wizualizacje QC: histogramy liczby wykrytych genów, umi, procent MT — w scanpy.pl.violin lub scanpy.pl.scatter.

Cel: eliminacja błędnych lub słabych danych źródłowych, zwiększenie jakości downstream.

### 2. Normalizacja i transformacja danych

Normalizacja:

  * Skalowanie sumy ekspresji (UMI) per komórkę do stałej wartości (np. 10,000) — redukcja wpływu różnej głębokości sekwencjonowania.
  * Log-transformation: Zamiana surowych odczytów UMI na postać logarytmiczną (log1p) dla stabilizacji wariancji.

Skalowanie:

  * Skalowanie genów do średniej 0 i odchylenia 1.

Opcjonalnie: usuwanie efektu zmienności technicznej (regresja zmiennych takich jak suma UMI, procent MT).

Narzędzia:

  * scanpy.pp.normalize_total(target_sum=1e4)
  * scanpy.pp.log1p()
  * scanpy.pp.scale()
  * scanpy.pp.regress_out(variables=['total_counts', 'pct_counts_mt'])

Cel: standaryzowane dane do sensownej analizy i trenowania modeli ML.

### 3. Balansowanie danych względem perturbacji

Analiza rozkładu komórek:

  * Wyznaczenie liczebności komórek w każdej grupie perturbacji (target_gene).
  * Identyfikacja perturbacji z nadmierną lub zbyt małą liczbą komórek.

Metody balansowania:

  * Downsampling: zmniejszenie liczby komórek dla nadreprezentowanych perturbacji, aby wyrównać klasy.
  * Ewentualne oversampling (np. SMOTE) tylko w ostateczności i ostrożnie.

Narzędzia: Pandas do grupowania i losowego próbkowania (np. df.groupby('target_gene').sample(n=desired_count, random_state=42)).

Cel: zapobiegnięcie biasowi w modelach ML, które mogłyby nadmiernie faworyzować najliczniej reprezentowane perturbacje.

### 4. Agregacja pseudobulkowa (opcjonalna)

Procedura:

  * Grupowanie komórek po perturbacji (target_gene).
  * Agregacja ekspresji genów: średnia lub mediana wartości ekspresji dla każdej grupy.

Przydatność:

  * Redukcja szumu biologicznego i technicznego.
  * Zmniejszenie wymiaru danych – ułatwia klasyczne modele np. regresji liniowej.

Narzędzia: Pandas (groupby('target_gene').agg('mean')) lub scanpy przy agregacji na poziomie AnnData.

Cel: alternatywna reprezentacja danych dla bazowego modelu lub wyjściowego benchmarku.

### 5. Integracja zewnętrznych publicznych perturbation datasets

Pobranie publicznych danych:

  * Wybranie uznanych datasetów perturbacji CRISPR single-cell RNA-seq, np. Perturb-seq, CROP-seq, Replogle 2022, Norman 2019.
  * Pobranie danych w dostępnych formatach (h5ad, CSV).

Mapowanie genów:

  * Standaryzacja nazw genów na wspólny zestaw 18,080 genów VCC.
  * Konwersja nazw do Ensembl IDs lub HGNC symboli.
  * Obsługa problemów z nazwami: aliasy, synonimy, brakujące geny.

Unifikacja formatów:

  * Konwersja do spójnego formatu (AnnData lub pandas DataFrame).
  * Synchronizacja struktury metadanych (cell metadata, perturbation labels).

Obsługa brakujących danych:

  * Dla genów nieobecnych zero-fill (0 ekspresji).
  * Opcjonalnie imputacja za pomocą prostych algorytmów (np. KNN, MICE).

Korekta batchów:

  * Usunięcie efektów technicznych i różnic eksperymentalnych między zbiorami.
  * Narzędzia: Harmony, BBKNN (Python).

Efekt końcowy:

  * Jeden duży, spójny dataset gotowy do pre-trainingu i fine-tuningu modeli.
  * Zwiększenie różnorodności danych, poprawa generalizacji.

### 6. Tworzenie splitów danych i cross-validation

Zadania:

  * Podział danych na trening, walidację i test zgodnie z regułami VCC (train on 150 genów perturbowanych).
  * Cross-validation z pozostawieniem genów poza treningiem („leave-genes-out”) imitujący finalne unseen genes.

Rodzaje splitów:

  * Splity per gen – cała grupa komórek dla danego genu idzie do jednego z podzbiorów.
  * Krotna walidacja (np. 5-fold CV) z różnymi zestawami genów.
  * Alternatywnie splity per eksperyment/batch jeżeli wpływa na ocenę.

Cel:

  * Dokładna ocena zdolności generalizacji modelu na nowe, niewidziane perturbacje.
  * Uniknięcie przecieku danych między zbiorami.
  * Narzędzia: sklearn.model_selection.GroupKFold lub własne implementacje w Pythonie bazujące na metadanych AnnData.

### 7. Feature engineering z biologiczną wiedzą

Dodatkowe cechy biologiczne mogące zwiększyć moc predykcyjną:

  * Informacje o genach regulatorowych np. czy dany gen jest celem czynników transkrypcyjnych (TF).
  * Przynależność genów do kluczowych ścieżek sygnalizacyjnych (KEGG, Reactome).
  * Dane o interakcjach gen-regulator (Gene Regulatory Networks).
  * Współczynniki biologiczne opisujące ekspresję w stanach pluripotentnych H1-hESC (np. OCT4, SOX2, NANOG).

Źródła danych:

  * Bazy KEGG, Reactome, Enrichr.
  * Bioinformatyczne API takie jak biomaRt, Ensembl REST API.

Metodologia:

  * Stworzenie dodatkowej macierzy cech (feature_df) zawierającej powyższe informacje na poziomie genów.
  * Dołączenie tych cech do modelu ML jako priorytetowej warstwy informacji.

Cel: podniesienie biologicznej interpretowalności modelu oraz jego skuteczności w przewidywaniu efektów perturbacji.

### 8. Przygotowanie danych wejściowych do ML/AI

Przetwarzanie i zapis:

  * Konwersja znormalizowanych danych do formatu h5ad z podziałem na train/val/pretrain.
  * Uporządkowanie metadanych i cech w strukturze oczekiwanej przez modele ML.

Budowa DataLoaderów PyTorch:

  * Implementacja klas Dataset i DataLoader obsługujących dane AnnData z indeksacją i efektywnym dostępem do dużych zbiorów.
  * Wsparcie trybu „AnnData-backed” - nie ładowanie całych danych do pamięci, tylko strumieniowy dostęp.

Zabezpieczenia:

  * Walidacja spójności danych.
  * Zapewnienie reproducibility (seed, zapisywanie konfiguracji).

Cel:

  * Efektywne wejście dla trenowania AI/ML modeli bez przeciążania pamięci RAM.
  * Możliwość szybkiego eksperymentowania i zmiany hiperparametrów.

### 9. Sanity checks i benchmarki

Testy biologiczne sanity check:

  * Sprawdzenie, czy model przewiduje znane biologiczne efekty knockoutów np. KO OCT4 powodujące spadek ekspresji markerów pluripotencji.
  * Analiza, czy model nie uczy się "na pamięć" konkretnych genów czy artefaktów.

Baseline models:

  * Prostota modeli: średnia ekspresja per target gene (cell mean predictor).
  * Pseudobulk linear model na zagregowanych danych.
  * Zaawansowane modele transferowe: scVI, CPA, GPerturb itp.

Porównanie wyników:

  * Metryki: korelacja Pearsona/Spearmana, RMSE na walidacji.
  * Ocena, czy nowe modele poprawiają się względem baseline.

Cel:

  * Zapewnienie wiarygodnej bazy porównawczej.
  * Umożliwienie zespołowi ML szybkie zorientowanie się, czy ich pomysły rzeczywiście działają.



## Dodatkowe

### Usuwanie potencjalnych źródeł data leakage

Analiza:

Upewnienie się, że wejściowe cechy modelu nie zawierają jawnej informacji o tym, jaki gen jest perturbowany (uniknięcie łatwego "przecieku danych").

Działania:

  * Usunięcie kolumn metadanych zawierających target labels z cech treningowych.
  * Oddzielenie metadanych opisujących komórki/perturbacje od danych wejściowych.
  * Ustanowienie ścisłego podziału danych między treningiem, walidacją i testem, aby uniknąć przecieku informacji.

Cel: zapewnienie, że wyniki modelu są autentyczne i generalizują, a nie "uczą się na pamięć".

