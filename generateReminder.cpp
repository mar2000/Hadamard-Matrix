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
            int coef = remainder[i] / q[q.size() - 1];

            for (int j = 0; j < q.size(); ++j) {
                remainder[i - j] -= coef * q[q.size() - 1 - j];
                remainder[i - j] %= p; // Redukcja modulo p
                if (remainder[i - j] < 0) remainder[i - j] += p; // Poprawka do zakresu 0...(p-1)
            }
        }
    }

    // Redukcja modulo p dla całej reszty
    for (int& coef : remainder) {
        coef %= p;
        if (coef < 0) coef += p;
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
            if (result[i + j] < 0) result[i + j] += p;
        }
    }

    return result;
}

// Funkcja do generowania wszystkich elementów ciała GF(p^m)
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
    q.resize(m+1);
    for (int i = m; i >= 0; --i) {
        cin >> q[i];
        q[i] %= p; // Redukcja modulo p
        if (q[i] < 0) q[i] += p;
    }

    // Generowanie elementów ciała GF(p^m)
    vector<Polynomial> fieldElements = generateFieldElements(p, m);

    // Dla każdego elementu oblicz jego kwadrat i resztę z dzielenia przez q
    cout << "Elementy ciała i odpowiadające im reszty z dzielenia:\n";
    int i=0;
    for (const Polynomial& s : fieldElements) {
        Polynomial w = squarePolynomial(s, p);
        Polynomial remainder = divideWithRemainder(w, q, p);
        i++;
        cout << "Element: ";
        printPolynomial(s);
        cout << "Kwadrat: ";
        printPolynomial(w);
        cout << "Reszta: ";
        printPolynomial(remainder);
        cout << "-----------------------------\n";
    }
    cout << i;
    return 0;
}
