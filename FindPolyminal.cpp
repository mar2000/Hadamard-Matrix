#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>

using namespace std;

// Funkcja pomocnicza do obliczania modulo
int mod(int a, int p) {
    return ((a % p) + p) % p;
}

// Funkcja do mnożenia wielomianów w Z_p[x]
vector<int> multiplyPolynomials(const vector<int>& a, const vector<int>& b, int p) {
    vector<int> result(a.size() + b.size() - 1, 0);
    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < b.size(); ++j) {
            result[i + j] = mod(result[i + j] + a[i] * b[j], p);
        }
    }
    return result;
}

// Funkcja do dzielenia wielomianów z resztą w Z_p[x]
vector<int> dividePolynomials(const vector<int>& dividend, const vector<int>& divisor, int p) {
    vector<int> result = dividend;
    while (result.size() >= divisor.size()) {
        int degreeDiff = result.size() - divisor.size();
        int leadingCoeff = mod(result.back() * divisor.back(), p);

        vector<int> subtractor(degreeDiff + 1, 0);
        subtractor.back() = leadingCoeff;
        subtractor = multiplyPolynomials(subtractor, divisor, p);

        for (size_t i = 0; i < subtractor.size(); ++i) {
            result[result.size() - subtractor.size() + i] = mod(result[result.size() - subtractor.size() + i] - subtractor[i], p);
        }

        while (!result.empty() && result.back() == 0) {
            result.pop_back();
        }
    }
    return result;
}

// Funkcja sprawdzająca, czy wielomian jest nierozkładalny
bool isIrreducible(const vector<int>& poly, int p) {
    if (poly.size() <= 1) return false;

    vector<int> x = {0, 1}; // Wielomian x
    vector<int> power = x;

    for (size_t d = 1; d <= poly.size() / 2; ++d) {
        power = multiplyPolynomials(power, x, p);
        power = dividePolynomials(power, poly, p);

        if (dividePolynomials(power, poly, p).empty()) {
            return false;
        }
    }

    return true;
}

// Funkcja generująca losowy wielomian stopnia m w Z_p[x]
vector<int> generateRandomPolynomial(int m, int p) {
    vector<int> poly(m + 1);
    for (int i = 0; i <= m; ++i) {
        poly[i] = rand() % p;
    }
    poly[m] = rand() % (p - 1) + 1; // Współczynnik najwyższego stopnia nie może być zerem
    return poly;
}

// Funkcja generująca wielomian nierozkładalny stopnia m w Z_p[x]
vector<int> generateIrreduciblePolynomial(int p, int m) {
    srand(time(nullptr));
    while (true) {
        vector<int> candidate = generateRandomPolynomial(m, p);
        if (isIrreducible(candidate, p)) {
            return candidate;
        }
    }
}

int main() {
    int p, m;
    cout << "Podaj liczbę pierwszą p: ";
    cin >> p;
    cout << "Podaj stopień wielomianu m: ";
    cin >> m;

    vector<int> irreduciblePolynomial = generateIrreduciblePolynomial(p, m);

    cout << "Współczynniki wielomianu nierozkładalnego: ";
    for (int coeff : irreduciblePolynomial) {
        cout << coeff << " ";
    }
    cout << endl;

    return 0;
}

