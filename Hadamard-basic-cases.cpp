// Pierwsza próba implementowania rozwiązania dla q = p^m

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <unordered_set>

using namespace std;

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


// Funkcja generująca elementy ciała GF(p^m)
vector<vector<int>> generateField(int p, int m) {
    int fieldSize = pow(p, m);
    vector<vector<int>> field(fieldSize, vector<int>(m, 0));

    // Generowanie elementów w reprezentacji wektorowej
    for (int i = 0; i < fieldSize; ++i) {
        int value = i;
        for (int j = 0; j < m; ++j) {
            field[i][j] = value % p;
            value /= p;
        }
    }

    return field;
}

// Funkcja obliczająca charakter kwadratowy w GF(p^m)
int quadraticCharacter(int a, int p) {
    if (a == 0) return 0; // Jeśli a = 0, charakter = 0
    int power = (p - 1) / 2; // Dla ciała GF(p), to jest (p-1)/2
    int result = 1;
    int base = a % p;

    // Iteracyjne potęgowanie modularne
    while (power > 0) {
        if (power % 2 == 1) result = (result * base) % p;
        base = (base * base) % p;
        power /= 2;
    }

    return result == 1 ? 1 : -1; // Zwraca 1, jeśli kwadrat, -1, jeśli niekwadrat
}

// Funkcja tworząca macierz Jacobsthala
vector<vector<int>> generateJacobianMatrix(int p, int m) {
    int fieldSize = pow(p, m);
    vector<vector<int>> matrix(fieldSize, vector<int>(fieldSize, 0));

    for (int i = 0; i < fieldSize; ++i) {
        for (int j = 0; j < fieldSize; ++j) {
            int difference = (i - j + fieldSize) % fieldSize; // (a - b) mod p^m
            matrix[i][j] = quadraticCharacter(difference, p);
        }
    }

    return matrix;
}


// Konstrukcja macierzy Hadamarda dla q+1, gdzie q = 4k+3
vector<vector<int>> construct_q_plus_one_hadamard(int n, int p, int m) {
    vector<vector<int>> Q = generateJacobianMatrix(p, m);

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
    vector<vector<int>> Q = generateJacobianMatrix(p, m);
    int q = pow(p, m);
        
    // Construct C matrix
    vector<vector<int>> C(q+1, vector<int>(q+1));
    for (int i = 1; i <= q; i++) {
        for (int j = 1; j <= q; j++) {
            C[i][j] = Q[i - 1][j - 1];
        }
    }
    for (int i = 1; i <= q; i++) {
        C[0][i] = C[i][0] = 1;
    }
    C[0][0] = 0; 

    // Identity matrix I
    vector<vector<int>> I(q+1, vector<int>(q+1, 0));
    for (int i = 0; i < q+1; i++) {
        I[i][i] = 1;
    }

    // Construct Hadamard matrix H
    vector<vector<int>> H(2 * (q+1), vector<int>(2 * (q+1)));
    for (int i = 0; i < q+1; i++) {
        for (int j = 0; j < q+1; j++) {
            H[i][j] = C[i][j] + I[i][j];
            H[i][j + q+1] = C[i][j] - I[i][j];
            H[i + q+1][j] = C[i][j] - I[i][j];
            H[i + q+1][j + q+1] = -(C[i][j] + I[i][j]);
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
