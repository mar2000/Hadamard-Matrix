#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <unordered_set>

using namespace std;

typedef vector<int> Polynomial;

    // Funkcja sprawdzająca, czy liczba jest pierwsza.
    // Parametry:
    // - n: Liczba całkowita do sprawdzenia.
    // Zwraca true, jeśli liczba n jest pierwsza, w przeciwnym razie false.
bool is_prime(int n) {
    if (n <= 1) return false;       // Liczby mniejsze lub równe 1 nie są pierwsze.
    
    // Sprawdzanie dzielników od 2 do pierwiastka z n.
    for (int p = 2; p * p <= n; p++) {      
        if (n % p == 0) return false;       // Liczba nie jest pierwsza, jeśli ma dzielnik inny niż 1 i n.
    }
    
    return true;        // Liczba jest pierwsza.
}

    // Funkcja sprawdzająca, czy liczba jest potęgą liczby pierwszej.
    // Jeśli liczba q jest potęgą liczby pierwszej, funkcja zapisuje:
    // - odpowiednią liczbę pierwszą w zmiennej `p`,
    // - potęgę w zmiennej `m`.
    // Parametry:
    // - q: Liczba całkowita do sprawdzenia.
    // - p: Referencja do zmiennej, w której zostanie zapisana liczba pierwsza.
    // - m: Referencja do zmiennej, w której zostanie zapisana potęga.
    // Zwraca true, jeśli liczba q jest potęgą liczby pierwszej, w przeciwnym razie false.
bool is_power_of_prime(int q, int &p, int &m) {
    if (is_prime(q)) {      // Sprawdzanie, czy q jest liczbą pierwszą.
        p = q;              // Jeśli tak, ustawiamy p jako q.
        m = 1;              // Potęga wynosi 1, ponieważ q^1 = q.
        return true;
    }
    
    // Iterowanie przez możliwe liczby pierwsze do pierwiastka z q.
    for (p = 2; p <= sqrt(q); p++) {        
        if (is_prime(p)) {                  // Sprawdzanie, czy p jest liczbą pierwszą.
            int power = 1;                  // Inicjalizacja potęgi.
            int result = p;                 // Inicjalizacja wyniku.
            
            while (result < q) {            // Obliczanie kolejnych potęg liczby p.
                result *= p;
                power++;
            }
            
            if (result == q) {              // Jeśli wynik jest równy q, znaleziono liczbę pierwszą i jej potęgę.
                m = power;
                return true;
            }
        }
    }
    
    return false;       // Jeśli nie znaleziono, zwraca false.
}


    // Funkcja dzieląca wielomian w przez wielomian q i zwracająca resztę z dzielenia modulo p.
    // Parametry:
    // - w: Wielomian do podzielenia, reprezentowany jako wektor współczynników, gdzie w[i] oznacza współczynnik przy x^i.
    // - q: Wielomian dzielnik, reprezentowany jako wektor współczynników, gdzie q[i] oznacza współczynnik przy x^i.
    // - p: Liczba całkowita reprezentująca modulo, w którym przeprowadzane są operacje.
    // Zwraca:
    // - Wielomian reprezentujący resztę z dzielenia w przez q, znormalizowany w zakresie 0...p-1.
Polynomial divideWithRemainder(const Polynomial &w, const Polynomial &q, int p) {
    Polynomial remainder = w;       // Kopia wielomianu w, która będzie modyfikowana jako reszta.

    // Iteracja przez współczynniki wielomianu od najwyższego stopnia.
    for (int i = w.size() - 1; i >= int(q.size()) - 1; --i) {           
        if (remainder[i] != 0) {                                        // Jeśli współczynnik przy bieżącym stopniu jest różny od zera.
            int coef = remainder[i] / q.back();                         // Wyznaczanie współczynnika podziału.
            
            // Aktualizacja reszty przez odjęcie wielokrotności wielomianu q.
            for (int j = 0; j < q.size(); ++j) {
                remainder[i - j] -= coef * q[q.size() - 1 - j];         // Odjęcie współczynników.
                remainder[i - j] %= p;                                  // Redukcja modulo p.
                if (remainder[i - j] < 0) remainder[i - j] += p;        // Korekta ujemnych wartości.
            }
        }
    }

    // Usuwanie zer z końca wielomianu w celu znormalizowania stopnia reszty.
    while (!remainder.empty() && remainder.back() == 0) {
        remainder.pop_back();
    }

    return remainder;       // Zwracanie reszty jako wynik.
}


    // Funkcja do podniesienia wielomianu do kwadratu modulo p.
    // Parametry:
    // - s: Wielomian, który ma zostać podniesiony do kwadratu, reprezentowany jako wektor współczynników, 
    //      gdzie s[i] oznacza współczynnik przy x^i.
    // - p: Liczba całkowita reprezentująca modulo, w którym przeprowadzane są operacje.
    // Zwraca:
    // - Wielomian będący wynikiem podniesienia s do kwadratu, zredukowany modulo p.
Polynomial squarePolynomial(const Polynomial &s, int p) {
    int deg_s = s.size() - 1;                   // Stopień wejściowego wielomianu.
    Polynomial result(2 * deg_s + 1, 0);        // Inicjalizacja wynikowego wielomianu o odpowiednim rozmiarze.

    // Iteracja przez wszystkie pary współczynników wielomianu s.
    for (int i = 0; i <= deg_s; ++i) {
        for (int j = 0; j <= deg_s; ++j) {
            result[i + j] += s[i] * s[j];       // Dodanie iloczynu współczynników do odpowiedniego miejsca.
            result[i + j] %= p;                 // Redukcja wyniku modulo p.
        }
    }

    return result;      // Zwracanie wynikowego wielomianu.
}

    // Funkcja do generowania wszystkich elementów ciała GF(p^m).
    // Parametry:
    // - p: Liczba pierwsza, która jest podstawą ciała GF(p^m).
    // - m: Stopień rozszerzenia, określający liczbę elementów w wielomianach.
    // Zwraca:
    // - Wektor zawierający wszystkie elementy ciała GF(p^m) w postaci wielomianów,
    //   gdzie każdy wielomian jest reprezentowany jako wektor współczynników.
vector<Polynomial> generateFieldElements(int p, int m) {
    vector<Polynomial> elements;        // Wektor do przechowywania elementów ciała.
    int numElements = pow(p, m);        // Liczba elementów w ciele GF(p^m).

    // Generowanie wielomianów reprezentujących wszystkie elementy ciała.
    for (int i = 0; i < numElements; ++i) {
        Polynomial poly(m, 0);          // Inicjalizacja wielomianu o stopniu m-1 z zerowymi współczynnikami.
        int value = i;                  // Reprezentacja elementu jako liczby całkowitej.
        
        for (int j = 0; j < m; ++j) {
            poly[j] = value % p;        // Obliczanie współczynnika przy x^j jako reszty z dzielenia przez p.
            value /= p;                 // Przejście do kolejnego współczynnika.
        }
        elements.push_back(poly);       // Dodanie wielomianu do wektora elementów.
    }

    return elements;        // Zwracanie wszystkich elementów ciała GF(p^m).
}


    // Funkcja do znajdowania nieredukowalnego wielomianu dla danych wartości p i m.
    // Parametry:
    // - p: Liczba pierwsza, która określa podstawę ciała GF(p).
    // - m: Stopień rozszerzenia ciała GF(p^m).
    // Zwraca:
    // - Wektor współczynników nieredukowalnego wielomianu nad GF(p) stopnia m, gdzie współczynnik przy x^i
    //   znajduje się na pozycji i w wektorze (np. wielomian x^2 + 1 to {1, 0, 1}).
    // Wyjątek:
    // - Rzuca `runtime_error`, jeśli dla podanych wartości p i m nie można znaleźć nieredukowalnego wielomianu.
vector<int> find_irreducible_polynomial(int p, int m) {
    if (m == 1) {
        return {1, 1};
    } else if (m == 2) {
        if (p == 3) return {1, 0, 1};               // x^2 + 1
        if (p == 5) return {1, 0, 2};               // x^2 + 2
        if (p == 7) return {1, 0, 1};               // x^2 + 1
        if (p == 11) return {1, 0, 1};              // x^2 + 1
        if (p == 13) return {1, 0, 2};              // x^2 + 2
        if (p == 17) return {1, 1, 1};              // x^2 + x + 1
        if (p == 19) return {1, 0, 1};              // x^2 + 1
    } else if (m == 3) {
        if (p == 3) return {1, 0, 2, 1};            // x^3 + 2x + 1
        if (p == 5) return {1, 0, 2, 1};            // x^3 + 2x + 1
        if (p == 7) return {1, 0, 2, 1};            // x^3 + 2x + 1
    } else if (m == 4) {
        if (p == 3) return {1, 1, 1, 1, 1};         // x^4 + x^3 + x^2 + x + 1
    } else if (m == 5) {
        if (p == 3) return {1, 0, 0, 0, 2, 1};      // x^5 + 2x + 1
    }

    // Rzucanie wyjątku, jeśli dla podanych p i m nie znaleziono nieredukowalnego wielomianu.
    throw runtime_error("Nie znaleziono nieredukowalnego wielomianu dla podanych parametrów.");
}


    // Funkcja do obliczania wielomianu -r modulo p.
    // Parametry:
    // - poly: Wielomian wejściowy, reprezentowany jako wektor współczynników, gdzie poly[i] oznacza współczynnik przy x^i.
    // - p: Modulo, w którym przeprowadzane są operacje arytmetyczne.
    // Zwraca:
    // - Wielomian, w którym każdy współczynnik jest równy -(poly[i]) modulo p.
Polynomial negatePolynomial(const Polynomial &poly, int p) {
    Polynomial negated = poly;      // Kopia wejściowego wielomianu.
    
    // Iteracja przez współczynniki wielomianu.
    for (int &coeff : negated) {        
        coeff = (p - coeff) % p;        // Przemnożenie przez -1 modulo p.
        if (coeff < 0) coeff += p;      // Korekta dla ujemnych wyników, aby współczynniki były w zakresie 0..p-1.
    }
    
    return negated;         // Zwracanie wynikowego wielomianu.
}

    // Funkcja do sprawdzania, czy wielomian jest resztą kwadratową modulo p.
    // Parametry:
    // - diff: Wielomian do sprawdzenia, reprezentowany jako wektor współczynników.
    // - remainders: Wektor wielomianów, które mogą być resztami kwadratowymi.
    // - n: Liczba całkowita, określająca parzystość w ramach (mod 4).
    // - p: Modulo, w którym przeprowadzane są operacje.
    // Zwraca:
    // - true, jeśli diff jest resztą kwadratową modulo p w ciele GF(p), false w przeciwnym razie.
bool isQuadraticResidue(const Polynomial &diff, const vector<Polynomial> &remainders, int n, int p) {
    // Iteracja przez potencjalne reszty kwadratowe.
    for (const Polynomial &r : remainders) {                    
        if (diff == r && n % 4 == 3) {                          // Sprawdzenie dla przypadku n ≡ 3 (mod 4).
            return true;
        }
        
        if (n % 4 == 1) {                                       // Sprawdzenie dla przypadku n ≡ 1 (mod 4).
            Polynomial negated_r = negatePolynomial(r, p);      // Obliczenie -r modulo p.
            if (diff == r || diff == negated_r) {               // Sprawdzenie, czy diff jest r lub -r.
                return true;
            }
        }
    }
    return false;       // Zwracanie false, jeśli diff nie jest resztą kwadratową.
}


    // Funkcja do generowania macierzy Jacobsthala.
    // Parametry:
    // - fieldElements: Wektor wielomianów reprezentujących elementy ciała GF(p^m).
    // - remainders: Wektor wielomianów będących potencjalnymi resztami kwadratowymi.
    // - p: Liczba pierwsza określająca modulo operacji arytmetycznych.
    // Zwraca:
    // - Macierz Jacobsthala (n x n), gdzie n to liczba elementów w ciele GF(p^m). Elementy macierzy to:
    //   - 1, jeśli różnica między dwoma elementami jest resztą kwadratową modulo p.
    //   - -1, jeśli różnica między dwoma elementami nie jest resztą kwadratową modulo p.
    //   - 0, na przekątnej (dla elementów porównywanych ze sobą).
vector<vector<int>> generateJacobianMatrix(const vector<Polynomial> &fieldElements, const vector<Polynomial> &remainders, int p) {
    int n = fieldElements.size();                           // Liczba elementów w ciele GF(p^m).
    vector<vector<int>> matrix(n, vector<int>(n, 0));       // Inicjalizacja macierzy n x n wypełnionej zerami.

    // Iteracja przez wszystkie pary elementów ciała GF(p^m).
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                matrix[i][j] = 0; // Na przekątnej macierzy wartość wynosi 0.
            } else {
                Polynomial diff = fieldElements[i];             // Kopia elementu i.
                
                // Obliczenie różnicy między elementami fieldElements[i] i fieldElements[j].
                for (int k = 0; k < diff.size(); ++k) {
                    if (k < fieldElements[j].size()) {
                        diff[k] -= fieldElements[j][k];         // Odjęcie współczynnika elementu j.
                    }
                    diff[k] %= p;                               // Redukcja modulo p.
                    if (diff[k] < 0) diff[k] += p;              // Korekta dla ujemnych wyników.
                }
                
                // Usunięcie wiodących zer w wielomianie różnicy.
                while (!diff.empty() && diff.back() == 0) {
                    diff.pop_back();
                }
                
                // Sprawdzenie, czy różnica jest resztą kwadratową.
                matrix[i][j] = isQuadraticResidue(diff, remainders, n, p) ? 1 : -1;
            }
        }
    }

    return matrix;      // Zwracanie wygenerowanej macierzy Jacobsthala.
}


    // Konstrukcja macierzy Hadamarda dla potęgi 2.
    // Parametry:
    // - n: Rozmiar macierzy (n musi być potęgą liczby 2, np. 2, 4, 8, ...).
    // Zwraca:
    // - Kwadratową macierz Hadamarda o rozmiarze n x n, gdzie elementy to 1 lub -1.
    // Macierz Hadamarda o rozmiarze potęgi 2 jest generowana rekurencyjnie za pomocą relacji:
    // H(2n) = [ H(n)  H(n) ]
    //         [ H(n) -H(n) ]
vector<vector<int>> construct_power_of_two_hadamard(int n) {
    // Inicjalizacja macierzy Hadamarda H(1) = [[1]].
    vector<vector<int>> H = {{1}};
    
    // Iteracyjne generowanie macierzy Hadamarda do żądanego rozmiaru n.
    while (H.size() < n) {
        int m = H.size();                                           // Aktualny rozmiar macierzy H(n).
        vector<vector<int>> new_H(2 * m, vector<int>(2 * m));       // Nowa macierz H(2n).

        // Wypełnianie nowej macierzy Hadamarda.
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                new_H[i][j] = H[i][j];                  // Górny lewy kwadrat: H(n).
                new_H[i][j + m] = H[i][j];              // Górny prawy kwadrat: H(n).
                new_H[i + m][j] = H[i][j];              // Dolny lewy kwadrat: H(n).
                new_H[i + m][j + m] = -H[i][j];         // Dolny prawy kwadrat: -H(n).
            }
        }
        
        H = move(new_H);        // Przypisanie wygenerowanej macierzy do H.
    }
    
    return H;       // Zwracanie wygenerowanej macierzy Hadamarda.
}


    // Konstrukcja macierzy Hadamarda dla q + 1, gdzie q = 4k + 3.
    // Parametry:
    // - n: Rozmiar macierzy (q + 1, gdzie q jest liczbą spełniającą 4k + 3).
    // - p: Liczba pierwsza definiująca ciało GF(p^m).
    // - m: Stopień rozszerzenia ciała GF(p^m).
    // Zwraca:
    // - Kwadratową macierz Hadamarda o rozmiarze n x n, gdzie n = q + 1.
    // Proces oparty jest na generowaniu ciała GF(p^m), obliczaniu reszt kwadratowych oraz generowaniu macierzy Jacobsthala.
vector<vector<int>> construct_q_plus_one_hadamard(int n, int p, int m) {
    // Generowanie elementów ciała GF(p^m).
    vector<Polynomial> fieldElements = generateFieldElements(p, m);

    // Znalezienie nieredukowalnego wielomianu q definiującego ciało GF(p^m).
    Polynomial q = find_irreducible_polynomial(p, m);

    // Obliczanie reszt kwadratowych dla każdego elementu ciała GF(p^m).
    vector<Polynomial> remainders;
    for (const Polynomial &s : fieldElements) {
        Polynomial w = squarePolynomial(s, p);                      // Obliczanie kwadratu wielomianu s.
        Polynomial remainder = divideWithRemainder(w, q, p);        // Reszta z dzielenia przez nieredukowalny wielomian q.
        remainders.push_back(remainder);                            // Dodanie reszty do wektora reszt.
    }

    // Generowanie macierzy Jacobsthala na podstawie elementów ciała i reszt kwadratowych.
    vector<vector<int>> Q = generateJacobianMatrix(fieldElements, remainders, p);

    // Konstrukcja macierzy Hadamarda o rozmiarze n x n.
    vector<vector<int>> H(n, vector<int>(n, 1));        // Inicjalizacja macierzy z 1 na wszystkich pozycjach.
    for (int i = 1; i < n; i++) {
        H[0][i] = 1;                                    // Ustawienie wartości w pierwszym wierszu.
        H[i][0] = -1;                                   // Ustawienie wartości w pierwszej kolumnie.
        for (int j = 1; j < n; j++) {
            H[i][j] = Q[i - 1][j - 1];                  // Wypełnienie pozostałych pozycji wartościami z macierzy Jacobsthala.
        }
        H[i][i] = 1;                                    // Ustawienie wartości na przekątnej.
    }

    return H;       // Zwracanie wygenerowanej macierzy Hadamarda.
}


    // Konstrukcja macierzy Hadamarda dla 2(q + 1), gdzie q = 4k + 1.
    // Parametry:
    // - n: Rozmiar macierzy (2 * (q + 1), gdzie q = 4k + 1).
    // - p: Liczba pierwsza definiująca ciało GF(p^m).
    // - m: Stopień rozszerzenia ciała GF(p^m).
    // Zwraca:
    // - Kwadratową macierz Hadamarda o rozmiarze 2(q + 1) x 2(q + 1).
    // Proces oparty jest na generowaniu ciała GF(p^m), obliczaniu reszt kwadratowych, generowaniu macierzy Jacobsthala oraz łączeniu jej z macierzą jednostkową.
vector<vector<int>> construct_two_q_plus_one_hadamard(int n, int p, int m) {
    // Generowanie elementów ciała GF(p^m).
    vector<Polynomial> fieldElements = generateFieldElements(p, m);

    // Znalezienie nieredukowalnego wielomianu q definiującego ciało GF(p^m).
    Polynomial q = find_irreducible_polynomial(p, m);

    // Obliczanie reszt kwadratowych dla każdego elementu ciała GF(p^m).
    vector<Polynomial> remainders;
    for (const Polynomial &s : fieldElements) {
        Polynomial w = squarePolynomial(s, p);                      // Obliczanie kwadratu wielomianu s.
        Polynomial remainder = divideWithRemainder(w, q, p);        // Reszta z dzielenia przez nieredukowalny wielomian q.
        remainders.push_back(remainder);                            // Dodanie reszty do wektora reszt.
    }

    // Generowanie macierzy Jacobsthala na podstawie elementów ciała i reszt kwadratowych.
    vector<vector<int>> Q = generateJacobianMatrix(fieldElements, remainders, p);
    int sizer = pow(p, m) + 1;      // Rozmiar macierzy Jacobsthala to q + 1.

    // Konstrukcja macierzy C na podstawie macierzy Jacobsthala.
    vector<vector<int>> C(sizer, vector<int>(sizer));
    for (int i = 1; i <= sizer - 1; i++) {         // Poza pierwszą kolumną i pierwszym rzędem jest macierz Q.
        for (int j = 1; j <= sizer - 1; j++) {
            C[i][j] = Q[i - 1][j - 1];
        }
    }
    for (int i = 1; i <= sizer - 1; i++) {
        C[0][i] = C[i][0] = 1;                     // Wartości w pierwszym wierszu i kolumnie.
    }
    C[0][0] = 0;                                   // Element w lewym górnym rogu.

    // Konstrukcja macierzy jednostkowej I.
    vector<vector<int>> I(sizer, vector<int>(sizer, 0));
    for (int i = 0; i < sizer; i++) {
        I[i][i] = 1;        // Jedynki na przekątnej.
    }

    // Konstrukcja końcowej macierzy Hadamarda H.
    vector<vector<int>> H(2 * sizer, vector<int>(2 * sizer));
    for (int i = 0; i < sizer; i++) {
        for (int j = 0; j < sizer; j++) {
            H[i][j] = C[i][j] + I[i][j];                            // Lewy górny kwadrat: C + I.
            H[i][j + sizer] = C[i][j] - I[i][j];                    // Prawy górny kwadrat: C - I.
            H[i + sizer][j] = C[i][j] - I[i][j];                    // Lewy dolny kwadrat: C - I.
            H[i + sizer][j + sizer] = -(C[i][j] + I[i][j]);         // Prawy dolny kwadrat: -(C + I).
        }
    }

    return H;       // Zwracanie wygenerowanej macierzy Hadamarda.
}


    // Funkcja obliczająca iloczyn Kroneckera macierzy A i B.
    // Iloczyn Kroneckera dwóch macierzy A (o wymiarach rowsA x colsA) i B (o wymiarach rowsB x colsB)
    // to macierz wynikowa o wymiarach (rowsA * rowsB) x (colsA * colsB), gdzie każdy element macierzy A
    // zostaje pomnożony przez całą macierz B.
    // Parametry:
    // - A: Pierwsza macierz wejściowa o wymiarach rowsA x colsA.
    // - B: Druga macierz wejściowa o wymiarach rowsB x colsB.
    // Zwraca:
    // - Macierz wynikowa będąca iloczynem Kroneckera o wymiarach (rowsA * rowsB) x (colsA * colsB).
vector<vector<int>> kroneckerProduct(const vector<vector<int>> &A, const vector<vector<int>> &B) {
    int rowsA = A.size();           // Liczba wierszy w macierzy A.
    int colsA = A[0].size();        // Liczba kolumn w macierzy A.
    int rowsB = B.size();           // Liczba wierszy w macierzy B.
    int colsB = B[0].size();        // Liczba kolumn w macierzy B.

    // Wymiary macierzy wynikowej.
    int rowsResult = rowsA * rowsB;         // Liczba wierszy w macierzy wynikowej.
    int colsResult = colsA * colsB;         // Liczba kolumn w macierzy wynikowej.

    // Inicjalizacja macierzy wynikowej wypełnionej zerami.
    vector<vector<int>> result(rowsResult, vector<int>(colsResult, 0));

    // Obliczanie iloczynu Kroneckera.
    for (int i = 0; i < rowsA; ++i) {                   // Iteracja po wierszach macierzy A.
        for (int j = 0; j < colsA; ++j) {               // Iteracja po kolumnach macierzy A.
            for (int k = 0; k < rowsB; ++k) {           // Iteracja po wierszach macierzy B.
                for (int l = 0; l < colsB; ++l) {       // Iteracja po kolumnach macierzy B.
                    // Obliczenie i przypisanie odpowiedniego elementu w macierzy wynikowej.
                    result[i * rowsB + k][j * colsB + l] = A[i][j] * B[k][l];
                }
            }
        }
    }

    return result;      // Zwrócenie macierzy wynikowej.
}


    // Główna funkcja do konstrukcji macierzy Hadamarda.
    // Funkcja ta generuje odpowiednią macierz Hadamarda dla podanego rozmiaru n.
    // Jeżeli n jest potęgą liczby 2, używa funkcji do konstrukcji macierzy Hadamarda dla potęgi 2.
    // W przeciwnym razie, dla innych wartości n, sprawdza, czy n spełnia odpowiednie warunki 
    // (w szczególności, czy q = n-1 lub q = n/2-1 jest potęgą liczby pierwszej oraz spełnia warunki modulo 4).
    // Na podstawie tego, wybiera odpowiednią metodę konstrukcji macierzy Hadamarda.
    // Parametry:
    // - n: Rozmiar macierzy Hadamarda, musi być liczbą dodatnią.
    // Zwraca:
    // - Macierz Hadamarda odpowiedniego rozmiaru.
vector<vector<int>> construct_hadamard(int n) {
    // Sprawdzenie, czy n jest potęgą liczby 2.
    if ((n & (n - 1)) == 0) { 
        return construct_power_of_two_hadamard(n);      // Jeśli n jest potęgą 2, wywołaj funkcję do konstrukcji Hadamarda.
    }

    int q, p, m;
    int power_of_two = 0;

    // Pętla, w której dzielimy n przez 2, aż uzyskamy n = 2^k.
    while (n % 2 == 0) {
        
        // Sprawdzenie, czy q = n-1 jest potęgą liczby pierwszej i q % 4 == 3.
        q = n - 1;
        if (is_power_of_prime(q, p, m) && (q % 4 == 3)) {
            // Konstrukcja macierzy Hadamarda dla q + 1 (gdzie q = 4k+3) oraz dla macierzy Hadamarda potęgi 2.
            vector<vector<int>> H_paleya = construct_q_plus_one_hadamard(n, p, m);
            vector<vector<int>> H2_power = construct_power_of_two_hadamard(pow(2, power_of_two));
            return kroneckerProduct(H2_power, H_paleya);        // Zwrócenie iloczynu Kroneckera.

        }
        
        // Sprawdzenie, czy q = n/2 - 1 jest potęgą liczby pierwszej i q % 4 == 1.
        q = n / 2 - 1;
        if (is_power_of_prime(q, p, m) && (q % 4 == 1)) {
            // Konstrukcja macierzy Hadamarda dla 2(q+1) (gdzie q = 4k+1) oraz dla macierzy Hadamarda potęgi 2.
            vector<vector<int>> H_paleya = construct_two_q_plus_one_hadamard(n, p, m);
            vector<vector<int>> H2_power = construct_power_of_two_hadamard(pow(2, power_of_two));
            return kroneckerProduct(H2_power, H_paleya);        // Zwrócenie iloczynu Kroneckera.
        }

        n /= 2;                 // Dzielimy n przez 2, aby uzyskać kolejne potencjalne wartości n.
        power_of_two++;         // Zwiększamy potęgę 2.
    }

    throw runtime_error("Invalid Hadamard order"); // Rzucenie wyjątku, jeśli n nie spełnia żadnych warunków.
}

    // Główna funkcja programu, która wykonuje konstrukcję macierzy Hadamarda dla zadanego rozmiaru n.
    // Program najpierw sprawdza, czy liczba n jest podzielna przez 4, ponieważ macierz Hadamarda
    // może być skonstruowana tylko dla rozmiarów będących wielokrotnościami 4. Następnie wywołuje funkcję 
    // konstrukcji odpowiedniej macierzy Hadamarda i wyświetla wynik na ekranie.
    // W przypadku wystąpienia błędów (np. niepoprawny rozmiar n) zwrócony zostanie komunikat o błędzie.
    // Parametry:
    // - n: Rozmiar macierzy Hadamarda (dodatnia liczba całkowita).
    // Zwraca:
    // - 0: Jeśli operacja zakończy się powodzeniem.
    // - 1: Jeśli wystąpił błąd (np. n nie jest podzielne przez 4 lub nie wyznaczono macierzy).
int main() {
    int n;
    cin >> n;       // Wczytanie wartości n (rozmiaru macierzy Hadamarda).

    // Sprawdzenie, czy n jest podzielne przez 4 i większe od 0 (macierz Hadamarda może być skonstruowana tylko dla n % 4 == 0, n > 0).
    if (n % 4 != 0 && n <= 0) {
        cout << "The Hadamard matrix does not exist";       // Komunikat o błędzie, jeśli n nie jest podzielne przez 4 lub n <= 0.
        return 1;                                           // Zakończenie programu z kodem błędu 1.
    }

    try {                       
        // Wywołanie funkcji do konstrukcji macierzy Hadamarda i przypisanie wyniku do zmiennej H.
        vector<vector<int>> H = construct_hadamard(n);

        // Wyświetlanie macierzy Hadamarda.
        for (const auto &row : H) {     // Iteracja po wierszach macierzy H.
           for (int val : row) {        // Iteracja po elementach wiersza.
               cout << val << " ";      // Wyświetlanie wartości wiersza.
           }
           cout << endl;                // Nowa linia po każdym wierszu.
        }
        
    } catch (const runtime_error &e) {                  // Obsługa wyjątków, jeśli wystąpi błąd w konstrukcji macierzy Hadamarda.
        cout << "Invalid Hadamard order" << endl;       // Komunikat o błędzie, jeśli n nie spełnia odpowiednich warunków.
        return 1;                                       // Zakończenie programu z kodem błędu 1.
    }

    return 0;       // Zakończenie programu z kodem sukcesu 0.
}
