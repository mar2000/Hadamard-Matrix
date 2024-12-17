// Program działa poprawnie dla n, które jest potęgą liczby 2 lub jest postaci (q+1) dla q = 3 mod 4 lub 2(q+1) dla q = 1 mod 4, gdzie q=p^m i p jest liczbą pierwszą.

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <unordered_set>

using namespace std; 

typedef vector<int> Polynomial;

// Funkcja sprawdzająca, czy liczba jest pierwsza
bool is_prime(int n) {
    if (n <= 1) return false;
    for (int p = 2; p * p <= n; p++) {
        if (n % p == 0) return false;
    }
    return true;
}

// Funkcja sprawdzająca, czy liczba jest potęgą liczby pierwszej
bool is_power_of_prime(int q, int &p, int &m) {
    if (is_prime(q)) {
        p = q;
        m = 1;
        return true;
    }
    for (p = 2; p <= sqrt(q); p++) {
        if (is_prime(p)) {
            int power = 1;
            int result = p;
            while (result < q) {
                result *= p;
                power++;
            }
            if (result == q) {
                m = power;
                return true;
            }
        }
    }
    
    return false;
}

// Konstrukcja macierzy Hadamarda dla potęgi 2
vector<vector<int>> construct_power_of_two_hadamard(int n) {
    vector<vector<int>> H = {{1}};
    while (H.size() < n) {
        int m = H.size();
        vector<vector<int>> new_H(2 * m, vector<int>(2 * m));
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                new_H[i][j] = H[i][j];
                new_H[i][j + m] = H[i][j];
                new_H[i + m][j] = H[i][j];
                new_H[i + m][j + m] = -H[i][j];
            }
        }
        H = move(new_H);
    }
    return H;
}



// Funkcja do dzielenia wielomianu w przez q, zwracająca resztę
Polynomial divideWithRemainder(const Polynomial& w, const Polynomial& q, int p) {
    Polynomial remainder = w;

    for (int i = w.size() - 1; i >= int(q.size()) - 1; --i) {
        if (remainder[i] != 0) { // Jeżeli współczynnik jest różny od 0
            int coef = remainder[i] / q.back();

            for (int j = 0; j < q.size(); ++j) {
                remainder[i - j] -= coef * q[q.size() - 1 - j];
                remainder[i - j] %= p; // Redukcja modulo p
                if (remainder[i - j] < 0) remainder[i - j] += p; // Poprawka do zakresu 0...(p-1)
            }
        }
    }

    // Obcięcie nadmiarowych współczynników
    while (!remainder.empty() && remainder.back() == 0) {
        remainder.pop_back();
    }

    return remainder;
}

// Funkcja do podniesienia wielomianu s do kwadratu
Polynomial squarePolynomial(const Polynomial& s, int p) {
    int deg_s = s.size() - 1;
    Polynomial result(2 * deg_s + 1, 0);

    for (int i = 0; i <= deg_s; ++i) {
        for (int j = 0; j <= deg_s; ++j) {
            result[i + j] += s[i] * s[j];
            result[i + j] %= p;
        }
    }

    return result;
}

// Funkcja do generowania wszystkich elementów ciala GF(p^m)
vector<Polynomial> generateFieldElements(int p, int m) {
    vector<Polynomial> elements;
    int numElements = pow(p, m);

    for (int i = 0; i < numElements; ++i) {
        Polynomial poly(m, 0);
        int value = i;
        for (int j = 0; j < m; ++j) {
            poly[j] = value % p;
            value /= p;
        }
        elements.push_back(poly);
    }

    return elements;
}

// Funkcja do sprawdzania, czy wielomian jest resztą kwadratową
// bool isQuadraticResidue(const Polynomial& diff, const vector<Polynomial>& remainders) {
//     for (const Polynomial& r : remainders) {
//         if (diff == r) {
//             return true;
//         }
//     }
//     return false;
// }

// Funkcja do obliczania wielomianu -r modulo p
Polynomial negatePolynomial(const Polynomial& poly, int p) {
    Polynomial negated = poly;
    for (int& coeff : negated) {
        coeff = (p - coeff) % p; // Przemnożenie przez -1 modulo p
        if (coeff < 0) coeff += p; // Poprawka, aby współczynniki były dodatnie
    }
    return negated;
}

// Funkcja do sprawdzania, czy wielomian jest resztą kwadratową
bool isQuadraticResidue(const Polynomial& diff, const vector<Polynomial>& remainders, int n, int p) {
    for (const Polynomial& r : remainders) {
        if (diff == r && n % 4 == 3) {
            return true; // Dla n ≡ 3 (mod 4)
        }
        if (n % 4 == 1) { // Dla n ≡ 1 (mod 4)
            Polynomial negated_r = negatePolynomial(r, p); // Obliczamy -r
            if (diff == r || diff == negated_r) {
                return true;
            }
        }
    }
    return false;
}


// Funkcja do generowania macierzy Jacobsthala
vector<vector<int>> generateJacobianMatrix(const vector<Polynomial>& fieldElements, const vector<Polynomial>& remainders, int p) {
    int n = fieldElements.size();
    vector<vector<int>> matrix(n, vector<int>(n, 0));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                matrix[i][j] = 0;
            } else {
                // Obliczenie różnicy modulo p
                Polynomial diff = fieldElements[i];
                for (int k = 0; k < diff.size(); ++k) {
                    if (k < fieldElements[j].size()) {
                        diff[k] -= fieldElements[j][k];
                    }
                    diff[k] %= p;
                    if (diff[k] < 0) diff[k] += p;
                }

                // Usunięcie zer z końca wielomianu
                while (!diff.empty() && diff.back() == 0) {
                    diff.pop_back();
                }

                // Sprawdzenie, czy różnica jest resztą kwadratową
                matrix[i][j] = isQuadraticResidue(diff, remainders, n, p) ? 1 : -1;
            }
        }
    }

    return matrix;
}



vector<int> find_irreducible_polynomial(int p, int m) {
    // Nieredukowalne wielomiany dla różnych przypadków
    if (m == 1) {
        return {1, 1};
    }
    if (m == 2) {
        if (p == 2) return {1, 1, 1};       // x^2 + x + 1
        if (p == 3) return {1, 0, 1};       // x^2 + 1
        if (p == 5) return {1, 0, 2};       // x^2 + 2
        if (p == 7) return {1, 0, 3};       // x^2 + 3
        if (p == 11) return {1, 0, 2};      // x^2 + 2
        if (p == 13) return {1, 0, 2};      // x^2 + 2
        if (p == 17) return {1, 0, 3};      // x^2 + 3
    } else if (m == 3) {
        if (p == 3) return {1, 0, 2, 1};    // x^3 + 2x + 1
        if (p == 5) return {1, 0, 2, 1};    // x^3 + 2x + 1
        if (p == 7) return {1, 0, 2, 1};    // x^3 + 2x + 1
    }
    
    throw runtime_error("Nieredukowalny wielomian nie został zdefiniowany dla podanych p i m");
}


// Konstrukcja macierzy Hadamarda dla q+1, gdzie q = 4k+3
vector<vector<int>> construct_q_plus_one_hadamard(int n, int p, int m) {
    // Generowanie elementów ciala GF(p^m)
    vector<Polynomial> fieldElements = generateFieldElements(p, m);

    Polynomial q = find_irreducible_polynomial(p, m);
    
    // Obliczanie reszt kwadratowych dla każdego elementu
    vector<Polynomial> remainders;
    for (const Polynomial& s : fieldElements) {
        Polynomial w = squarePolynomial(s, p);
        Polynomial remainder = divideWithRemainder(w, q, p);
        remainders.push_back(remainder);
    }

    // Generowanie macierzy Jacobsthala
    vector<vector<int>> Q = generateJacobianMatrix(fieldElements, remainders, p);

    vector<vector<int>> H(n, vector<int>(n, 1));
    for (int i = 1; i < n; i++) {
        H[0][i] = 1;
        H[i][0] = -1;
        for (int j = 1; j < n; j++) {
            H[i][j] = Q[i - 1][j - 1];
        }
        H[i][i] = 1;
    }
    return H;
}

// Konstrukcja macierzy Hadamarda dla 2(q+1), gdzie q = 4k+1
vector<vector<int>> construct_two_q_plus_one_hadamard(int n, int p, int m) {
    // Generowanie elementów ciala GF(p^m)
    vector<Polynomial> fieldElements = generateFieldElements(p, m);

    Polynomial q = find_irreducible_polynomial(p, m);
    
    // Obliczanie reszt kwadratowych dla każdego elementu
    vector<Polynomial> remainders;
    for (const Polynomial& s : fieldElements) {
        Polynomial w = squarePolynomial(s, p);
        Polynomial remainder = divideWithRemainder(w, q, p);
        remainders.push_back(remainder);
    }

    // Generowanie macierzy Jacobsthala
    vector<vector<int>> Q = generateJacobianMatrix(fieldElements, remainders, p);
    int sizer = pow(p, m) + 1;
        
    // Construct C matrix
    vector<vector<int>> C(sizer, vector<int>(sizer));
    for (int i = 1; i <= sizer - 1; i++) {
        for (int j = 1; j <= sizer - 1; j++) {
            C[i][j] = Q[i - 1][j - 1];
        }
    }
    for (int i = 1; i <= sizer - 1; i++) {
        C[0][i] = C[i][0] = 1;
    }
    C[0][0] = 0; 

    // Identity matrix I
    vector<vector<int>> I(sizer, vector<int>(sizer, 0));
    for (int i = 0; i < sizer; i++) {
        I[i][i] = 1;
    }

    // Construct Hadamard matrix H
    vector<vector<int>> H(2 * (sizer), vector<int>(2 * (sizer)));
    for (int i = 0; i < sizer; i++) {
        for (int j = 0; j < sizer; j++) {
            H[i][j] = C[i][j] + I[i][j];
            H[i][j + sizer] = C[i][j] - I[i][j];
            H[i + sizer][j] = C[i][j] - I[i][j];
            H[i + sizer][j + sizer] = -(C[i][j] + I[i][j]);
        }
    }

    return H;
}






// Główna funkcja do konstrukcji macierzy Hadamarda
vector<vector<int>> construct_hadamard(int n) {
    assert(n > 0);

    if ((n & (n - 1)) == 0) { // Potęga 2
        return construct_power_of_two_hadamard(n);
    }

    int q, p, m;
    q = n - 1;
    if (is_power_of_prime(q, p, m) && (q % 4 == 3)) { // q+1, q = 4k+3
        return construct_q_plus_one_hadamard(n, p, m);
    }

    q = n / 2 - 1;
    if (is_power_of_prime(q, p, m) && (q % 4 == 1)) { // 2(q+1), q = 4k+1
         return construct_two_q_plus_one_hadamard(n, p, m);
    }

    throw runtime_error("Invalid Hadamard Order");
}

int main() {
    int n;
    cin >> n;

    if (n % 4 != 0) {
        cout << "The number must be divisible by 4.";
        return 1;
    }
    
    try {
        vector<vector<int>> H = construct_hadamard(n);
        for (const auto &row : H) {
            for (int val : row) {
                cout << val << " ";
            }
            cout << endl;
        }
    } catch (const runtime_error &e) {
        cout << "Invalid Hadamard order" << endl;
    }

    return 0;
}


