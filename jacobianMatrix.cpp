#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

typedef vector<int> Polynomial;

// Funkcja do wypisania wielomianu
void printPolynomial(const Polynomial& p) {
    bool first = true;
    for (int i = p.size() - 1; i >= 0; --i) {
        if (p[i] != 0) { // Ignoruj współczynniki równe 0
            if (!first && p[i] > 0) cout << "+";
            cout << p[i];
            if (i > 0) cout << "x";
            if (i > 1) cout << "^" << i;
            first = false;
        }
    }
    if (first) cout << "0";
    cout << endl;
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
bool isQuadraticResidue(const Polynomial& diff, const vector<Polynomial>& remainders) {
    for (const Polynomial& r : remainders) {
        if (diff == r) {
            return true;
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
                matrix[i][j] = isQuadraticResidue(diff, remainders) ? 1 : -1;
            }
        }
    }

    return matrix;
}

int main() {
    int p, m;
    Polynomial q;

    // Wczytanie liczby p i m
    cout << "Podaj p (modulo): ";
    cin >> p;
    cout << "Podaj m (stopień wielomianu - 1): ";
    cin >> m;

    // Wczytanie wielomianu q
    cout << "Podaj wspolczynniki wielomianu q od najwyzszej potegi (stopien = " << m << "): \n";
    q.resize(m + 1);
    for (int i = m; i >= 0; --i) {
        cin >> q[i];
        q[i] %= p; // Redukcja modulo p
        if (q[i] < 0) q[i] += p;
    }

    // Generowanie elementów ciala GF(p^m)
    vector<Polynomial> fieldElements = generateFieldElements(p, m);

    // Obliczanie reszt kwadratowych dla każdego elementu
    vector<Polynomial> remainders;
    for (const Polynomial& s : fieldElements) {
        Polynomial w = squarePolynomial(s, p);
        Polynomial remainder = divideWithRemainder(w, q, p);
        remainders.push_back(remainder);
    }

    // Generowanie macierzy Jacobsthala
    vector<vector<int>> jacobianMatrix = generateJacobianMatrix(fieldElements, remainders, p);

    // Wyświetlenie macierzy Jacobsthala
    cout << "Macierz Jacobsthala:\n";
    for (const auto& row : jacobianMatrix) {
        for (int val : row) {
            cout << val << " ";
        }
        cout << endl;
    }

    return 0;
}
