Projekt, który powstaje w ramach przedmiotu "Wybrane zagadnienia matematyki dyskretnej" (1000-2M07MD) na MIMUW.
Składa się z:
- artykułu Hadamard-matrix.pdf, w którym zostały zdefiniowane macierze Hadamarda oraz zostały udowodnione poprawności algorytmicznych konstrukcji macierzy Paleya
- programu Hadamard-simple-cases.cpp napisanego w języku C++, który wyznacza macierze Hadamarda dla podanego rzędu (o ile rząd spełnia odpowiednie warunki)
  - n = potęga 2
  - n = q + 1, gdzie q = 3 mod 4 i q jest liczbą pierwszą
  - n = 2 * (q + 1), gdzie q = 1 mod 4 i q jest liczbą pierwszą
- programu IsHadamard.cpp sprawdzającego czy dana macierz jest macierzą Hadamarda
