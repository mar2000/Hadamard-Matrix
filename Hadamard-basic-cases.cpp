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


// Funkcja generująca ciało GF(p)
vector<int> generate_field(int p) {
    vector<int> field(p);
    for (int i = 0; i < p; ++i) field[i] = i;
    return field;
}

// Funkcja do obliczenia znaku kwadratowego w GF(q)
int quadratic_character(int a, int q) {
    if (a == 0) return 0;
    for (int b = 1; b < q; b++) {
        if ((b * b) % q == a) return 1;
    }
    return -1;
}

// Funkcja generująca macierz Jakobstala
vector<vector<int>> generate_jacobsthal_matrix(int q) {
    vector<int> field = generate_field(q);
    vector<vector<int>> jacobsthal_matrix(q, vector<int>(q, 0));

    for (int i = 0; i < q; ++i) {
        for (int j = 0; j < q; ++j) {
            int diff = (field[i] - field[j] + q) % q;
            jacobsthal_matrix[i][j] = quadratic_character(diff, q);
        }
    }

    return jacobsthal_matrix;
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































// Pierwsza próba implementowania rozwiązania dla q = p^m

// #include <iostream>
// #include <vector>
// #include <cmath>
// #include <cassert>
// #include <unordered_set>

// using namespace std;

// // Funkcja sprawdzająca, czy liczba jest pierwsza
// bool is_prime(int n) {
//     if (n <= 1) return false;
//     for (int p = 2; p * p <= n; p++) {
//         if (n % p == 0) return false;
//     }
//     return true;
// }

// // Funkcja sprawdzająca, czy liczba jest potęgą liczby pierwszej
// bool is_power_of_prime(int q, int &p, int &m) {
//     if (is_prime(q)) {
//         p = q;
//         m = 1;
//         return true;
//     }
//     for (p = 2; p <= sqrt(q); p++) {
//         if (is_prime(p)) {
//             int power = 1;
//             int result = p;
//             while (result < q) {
//                 result *= p;
//                 power++;
//             }
//             if (result == q) {
//                 m = power;
//                 return true;
//             }
//         }
//     }
    
//     return false;
// }

// // Konstrukcja macierzy Hadamarda dla potęgi 2
// vector<vector<int>> construct_power_of_two_hadamard(int n) {
//     vector<vector<int>> H = {{1}};
//     while (H.size() < n) {
//         int m = H.size();
//         vector<vector<int>> new_H(2 * m, vector<int>(2 * m));
//         for (int i = 0; i < m; i++) {
//             for (int j = 0; j < m; j++) {
//                 new_H[i][j] = H[i][j];
//                 new_H[i][j + m] = H[i][j];
//                 new_H[i + m][j] = H[i][j];
//                 new_H[i + m][j + m] = -H[i][j];
//             }
//         }
//         H = move(new_H);
//     }
//     return H;
// }


// vector<int> find_irreducible_polynomial(int p, int m) {
//     // Nieredukowalne wielomiany dla różnych przypadków
//     if (m == 2) {
//         if (p == 2) return {1, 1, 1};       // x^2 + x + 1
//         if (p == 5) return {1, 0, 2};       // x^2 + 2
//         if (p == 7) return {1, 0, 3};       // x^2 + 3
//         if (p == 11) return {1, 0, 2};      // x^2 + 2
//         if (p == 13) return {1, 0, 2};      // x^2 + 2
//         if (p == 17) return {1, 0, 3};      // x^2 + 3
//     } else if (m == 3) {
//         if (p == 7) return {1, 1, 0, 3};    // x^3 + x + 3
//     }
    
//     // Rzucenie wyjątku, jeśli przypadek nie jest obsłużony
//     throw runtime_error("Nieredukowalny wielomian nie został zdefiniowany dla podanych p i m");
// }


// // Dodawanie wielomianów w GF(p)
// vector<int> add_polynomials(const vector<int>& a, const vector<int>& b, int p) {
//     vector<int> result(max(a.size(), b.size()), 0);
//     for (size_t i = 0; i < a.size(); ++i) result[i] = a[i];
//     for (size_t i = 0; i < b.size(); ++i) result[i] = (result[i] + b[i]) % p;
//     return result;
// }

// // Mnożenie wielomianów w GF(p)
// vector<int> multiply_polynomials(const vector<int>& a, const vector<int>& b, int p) {
//     vector<int> product(a.size() + b.size() - 1, 0);
//     for (size_t i = 0; i < a.size(); ++i) {
//         for (size_t j = 0; j < b.size(); ++j) {
//             product[i + j] = (product[i + j] + a[i] * b[j]) % p;
//         }
//     }
//     return product;
// }

// // Redukcja wielomianu modulo inny wielomian
// vector<int> reduce_modulo(const vector<int>& poly, const vector<int>& mod, int p) {
//     vector<int> result = poly;
//     while (result.size() >= mod.size()) {
//         int coeff = result.back();
//         int degree_diff = result.size() - mod.size();
//         for (size_t i = 0; i < mod.size(); ++i) {
//             result[degree_diff + i] = (result[degree_diff + i] - coeff * mod[i]) % p;
//             if (result[degree_diff + i] < 0) result[degree_diff + i] += p;
//         }
//         result.pop_back();
//     }
//     return result;
// }

// // Potęgowanie modularne w GF(p^m)
// vector<int> power_mod(const vector<int>& base, int exp, const vector<int>& mod, int p) {
//     vector<int> result = {1}; // Element neutralny
//     vector<int> current = base;

//     while (exp > 0) {
//         if (exp % 2 == 1) {
//             result = reduce_modulo(multiply_polynomials(result, current, p), mod, p);
//         }
//         current = reduce_modulo(multiply_polynomials(current, current, p), mod, p);
//         exp /= 2;
//     }

//     return result;
// }

// // Generowanie elementów ciała GF(p^m)
// vector<vector<int>> generate_galois_field(int p, int m, const vector<int>& irreducible) {
//     int field_size = pow(p, m);
//     vector<vector<int>> field(field_size, vector<int>(m, 0));

//     // Każdy element GF(p^m) reprezentowany jako wielomian
//     for (int i = 0; i < field_size; ++i) {
//         int value = i;
//         for (int j = 0; j < m; ++j) {
//             field[i][j] = value % p;
//             value /= p;
//         }
//     }

//     return field;
// }

// // Obliczanie znaku kwadratowego
// int quadratic_character(const vector<int>& element, const vector<vector<int>>& field, const vector<int>& irreducible, int p) {
//     if (element == vector<int>(element.size(), 0)) return 0; // Znak dla 0
//     int field_size = field.size();
//     vector<int> result = power_mod(element, (field_size - 1) / 2, irreducible, p); // Element^(q-1)/2
//     return (result == vector<int>(result.size(), 1)) ? 1 : -1; // 1 dla kwadratów, -1 dla niekwadratów
// }

// // Generowanie macierzy Jakobstala
// vector<vector<int>> generateJacobianMatrix(int p, int m) {
//     vector<int> irreducible = find_irreducible_polynomial(p, m);
//     vector<vector<int>> field = generate_galois_field(p, m, irreducible);

//     int field_size = field.size();
//     vector<vector<int>> jacobian_matrix(field_size, vector<int>(field_size, 0));

//     for (int i = 0; i < field_size; ++i) {
//         for (int j = 0; j < field_size; ++j) {
//             vector<int> difference = add_polynomials(field[i], field[j], p); // Obliczanie (a - b)
//             jacobian_matrix[i][j] = quadratic_character(difference, field, irreducible, p);
//         }
//     }

//     return jacobian_matrix;
// }





// // Konstrukcja macierzy Hadamarda dla q+1, gdzie q = 4k+3
// vector<vector<int>> construct_q_plus_one_hadamard(int n, int p, int m) {
//     vector<vector<int>> Q = generateJacobianMatrix(p, m);

//     vector<vector<int>> H(n, vector<int>(n, 1));
//     for (int i = 1; i < n; i++) {
//         H[0][i] = 1;
//         H[i][0] = -1;
//         for (int j = 1; j < n; j++) {
//             H[i][j] = Q[i - 1][j - 1];
//         }
//         H[i][i] = 1;
//     }
//     return H;
// }

// // Konstrukcja macierzy Hadamarda dla 2(q+1), gdzie q = 4k+1
// vector<vector<int>> construct_two_q_plus_one_hadamard(int n, int p, int m) {
//     vector<vector<int>> Q = generateJacobianMatrix(p, m);
//     int q = pow(p, m);
        
//     // Construct C matrix
//     vector<vector<int>> C(q+1, vector<int>(q+1));
//     for (int i = 1; i <= q; i++) {
//         for (int j = 1; j <= q; j++) {
//             C[i][j] = Q[i - 1][j - 1];
//         }
//     }
//     for (int i = 1; i <= q; i++) {
//         C[0][i] = C[i][0] = 1;
//     }
//     C[0][0] = 0; 

//     // Identity matrix I
//     vector<vector<int>> I(q+1, vector<int>(q+1, 0));
//     for (int i = 0; i < q+1; i++) {
//         I[i][i] = 1;
//     }

//     // Construct Hadamard matrix H
//     vector<vector<int>> H(2 * (q+1), vector<int>(2 * (q+1)));
//     for (int i = 0; i < q+1; i++) {
//         for (int j = 0; j < q+1; j++) {
//             H[i][j] = C[i][j] + I[i][j];
//             H[i][j + q+1] = C[i][j] - I[i][j];
//             H[i + q+1][j] = C[i][j] - I[i][j];
//             H[i + q+1][j + q+1] = -(C[i][j] + I[i][j]);
//         }
//     }

//     return H;
// }






// // Główna funkcja do konstrukcji macierzy Hadamarda
// vector<vector<int>> construct_hadamard(int n) {
//     assert(n > 0);

//     if ((n & (n - 1)) == 0) { // Potęga 2
//         return construct_power_of_two_hadamard(n);
//     }

//     int q, p, m;
//     q = n - 1;
//     if (is_power_of_prime(q, p, m) && (q % 4 == 3)) { // q+1, q = 4k+3
//         return construct_q_plus_one_hadamard(n, p, m);
//     }

//     q = n / 2 - 1;
//     if (is_power_of_prime(q, p, m) && (q % 4 == 1)) { // 2(q+1), q = 4k+1
//          return construct_two_q_plus_one_hadamard(n, p, m);
//     }

//     throw runtime_error("Invalid Hadamard Order");
// }

// int main() {
//     int n;
//     cin >> n;

//     if (n % 4 != 0) {
//         cout << "The number must be divisible by 4.";
//         return 1;
//     }
    
//     try {
//         vector<vector<int>> H = construct_hadamard(n);
//         for (const auto &row : H) {
//             for (int val : row) {
//                 cout << val << " ";
//             }
//             cout << endl;
//         }
//     } catch (const runtime_error &e) {
//         cout << "Invalid Hadamard order" << endl;
//     }

//     return 0;
// }
