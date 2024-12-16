#include <iostream>
#include <vector>
#include <iomanip>
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

int main() {
    Polynomial w, q, remainder;
    int deg_w, deg_q, p;

    // Wczytanie liczby p
    cout << "Podaj p (modulo): ";
    cin >> p;

    // Wczytanie wielomianu w
    cout << "Podaj stopien wielomianu w: ";
    cin >> deg_w;
    w.resize(deg_w + 1);
    cout << "Podaj wspolczynniki wielomianu w od najwyzszej potegi: \n";
    for (int i = deg_w; i >= 0; --i) {
        cin >> w[i];
        w[i] %= p; // Redukcja modulo p
        if (w[i] < 0) w[i] += p;
    }

    // Wczytanie wielomianu q
    cout << "Podaj stopien wielomianu q: ";
    cin >> deg_q;
    q.resize(deg_q + 1);
    cout << "Podaj wspolczynniki wielomianu q od najwyzszej potegi: \n";
    for (int i = deg_q; i >= 0; --i) {
        cin >> q[i];
        q[i] %= p; // Redukcja modulo p
        if (q[i] < 0) q[i] += p;
    }

    // Dzielenie wielomianów
    remainder = divideWithRemainder(w, q, p);

    // Wynik
    cout << "Reszta: ";
    printPolynomial(remainder);

    return 0;
}
