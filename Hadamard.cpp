// Program działa poprawnie dla n, które jest potęgą liczby 2 lub jest postaci (p+1) dla p = 3 mod 4 lub 2(q+1) dla p = 1 mod 4, gdzie p jest liczbą pierwszą.

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <unordered_set>

using namespace std;

bool is_prime(int n) {
    if (n <= 1) return false;
    for (int p = 2; p * p <= n; p++) {
        if (n % p == 0) return false;
    }
    return true;
}

bool is_valid_hadamard_order(int n) {
    if (n <= 1 || n >= 500) return false;
    if ((n & (n - 1)) == 0) return true; // Power of 2

    int q;
    if (is_prime(n - 1) && ((n - 1) % 4 == 3)) return true; // q+1, q = 4k+3

    if (n % 2 == 0 && is_prime(n / 2 - 1) && ((n / 2 - 1) % 4 == 1)) return true; // 2(q+1), q = 4k+1

    return false;
}

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

vector<vector<int>> construct_q_plus_one_hadamard(int n, int q) {
    vector<int> chi(q, 0);
    for (int x = 1; x < q; x++) {
        chi[(x * x) % q] = 1;
    }
    vector<vector<int>> Q(q, vector<int>(q));
    for (int i = 0; i < q; i++) {
        for (int j = 0; j < q; j++) {
            Q[i][j] = (i == j) ? 0 : (chi[(i - j + q) % q] ? 1 : -1);
        }
    }

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

vector<vector<int>> construct_two_q_plus_one_hadamard(int n, int q) {
    int m = q + 1;

    // Construct Q matrix for chi values
    vector<int> chi(q, 0);
    for (int x = 1; x < q; x++) {
        chi[(x * x) % q] = 1;
    }
    vector<vector<int>> Q(q, vector<int>(q));
    for (int i = 0; i < q; i++) {
        for (int j = 0; j < q; j++) {
            Q[i][j] = (i == j) ? 0 : (chi[(i - j + q) % q] ? 1 : -1);
        }
    }

    // Construct C matrix
    vector<vector<int>> C(m, vector<int>(m));
    for (int i = 1; i <= q; i++) {
        for (int j = 1; j <= q; j++) {
            C[i][j] = Q[i - 1][j - 1];
        }
    }
    for (int i = 1; i <= q; i++) {
        C[0][i] = C[i][0] = 1;
    }
    C[0][0] = 0; // Top-left corner

    // Identity matrix I
    vector<vector<int>> I(m, vector<int>(m, 0));
    for (int i = 0; i < m; i++) {
        I[i][i] = 1;
    }

    // Construct Hadamard matrix H
    vector<vector<int>> H(2 * m, vector<int>(2 * m));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            H[i][j] = C[i][j] + I[i][j];
            H[i][j + m] = C[i][j] - I[i][j];
            H[i + m][j] = C[i][j] - I[i][j];
            H[i + m][j + m] = -(C[i][j] + I[i][j]);
        }
    }

    return H;
}


vector<vector<int>> construct_hadamard(int n) {
    assert(is_valid_hadamard_order(n));

    if ((n & (n - 1)) == 0) { // Power of 2
        return construct_power_of_two_hadamard(n);
    }

    int q;
    if (is_prime(n - 1) && ((n - 1) % 4 == 3)) { // q+1 case
        return construct_q_plus_one_hadamard(n, n - 1);
    }

    if (n % 2 == 0 && is_prime(n / 2 - 1) && ((n / 2 - 1) % 4 == 1)) { // 2(q+1) case
        return construct_two_q_plus_one_hadamard(n, n / 2 - 1);
    }

    throw runtime_error("Invalid Hadamard order");
}

int main() {
    int n;
    cin >> n;

    if (!is_valid_hadamard_order(n)) {
        cout << "Invalid Hadamard order" << endl;
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

// Działające rozmiary n in {2, 4, 8, 12, 16, 20, 24, 28, 32, 36, 44, 48, 60, 64, 68, 72, 76, 80, 84}
// Niedziałąjące rozmiary n in {28, 36, 40, 52, 56, 76, 88, 92, 96, 100}
