#include <iostream>
#include <vector>

using namespace std;

int main() {
    // Wczytywanie macierzy M
    int n, rows, cols;
    cin >> n;
    cols = n;
    rows = n;

    vector<vector<int>> M(rows, vector<int>(cols));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            cin >> M[i][j];
        }
    }

    // Utworzenie macierzy wynikowej M^T * M
    vector<vector<int>> result(cols, vector<int>(cols, 0));

    // Obliczanie M^T * M
    for (int i = 0; i < cols; ++i) {
        for (int j = 0; j < cols; ++j) {
            for (int k = 0; k < rows; ++k) {
                result[i][j] += M[k][i] * M[k][j];
            }
        }
    }

    // Wyświetlenie wyniku
    for (const auto& row : result) {
        for (int value : row) {
            cout << value << " ";
        }
        cout << endl;
    }

    return 0;
}
