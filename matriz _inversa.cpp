#include <iostream>
#include <vector>

using namespace std;

class matrizInversa{
public:
vector<vector<double>> getcofactor(const vector<vector<double>>& m, int i, int j) {
    vector<vector<double>> result;
    for (int row = 0; row < m.size(); ++row) {
        if (row != i) {
            vector<double> new_row;
            for (int col = 0; col < m[row].size(); ++col) {
                if (col != j) {
                    new_row.push_back(m[row][col]);
                }
            }
            result.push_back(new_row);
        }
    }
    return result;
}

double determinantOfMatrix(const vector<vector<double>>& mat) {
    if (mat.size() == 1) {
        return mat[0][0];
    }
    if (mat.size() == 2) {
        return mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];
    }
    double Sum = 0;
    for (int current_column = 0; current_column < mat.size(); ++current_column) {
        int sign = (current_column % 2 == 0) ? 1 : -1;
        vector<vector<double>> sub_matrix = getcofactor(mat, 0, current_column);
        double sub_det = determinantOfMatrix(sub_matrix);
        Sum += (sign * mat[0][current_column] * sub_det);
    }
    return Sum;
}

void pivot(vector<vector<double>>& mat, int col) {
    int tam = mat.size();
    if (mat[col][col] == 0) {
        int sel = col;
        for (int i = col + 1; i < tam; ++i) {
            if (mat[i][col] != 0) {
                sel = i;
                break;
            }
        }
        if (sel != col) {
            swap(mat[col], mat[sel]);
        }
    }
}

vector<vector<double>> aumenta(const vector<vector<double>>& mat, int tam) {
    vector<vector<double>> result = mat;
    for (int i = 0; i < tam; ++i) {
        vector<double> identity_row;
        for (int j = 0; j < tam; ++j) {
            if (i == j) {
                identity_row.push_back(1);
            } else {
                identity_row.push_back(0);
            }
        }
        result[i].insert(result[i].end(), identity_row.begin(), identity_row.end());
    }
    for (int i = 0; i < tam; ++i) {
        vector<double> zeros(tam, 0);
        result.push_back(zeros);
    }
    return result;
}

vector<vector<double>> escalonamento(vector<vector<double>>& mat, int tam) {
    for (int i = 0; i < tam; ++i) {
        if (mat[i][i] != 0) {
            pivot(mat, i);
        }
        for (int j = 0; j < tam; ++j) {
            if (i == j) {
                continue;
            }
            double mul = mat[j][i] / mat[i][i];
            for (int k = 0; k < tam * 2; ++k) {
                mat[j][k] = mat[j][k] - (mat[i][k] * mul);
            }
        }
    }
    return mat;
}

vector<vector<double>> arruma(vector<vector<double>>& mat, int tam) {
    for (int i = 0; i < tam; ++i) {
        double aux = mat[i][i];
        for (int j = 0; j < tam * 2; ++j) {
            mat[i][j] = mat[i][j] / aux;
        }
    }
    return mat;
}

vector<vector<double>> inversa(vector<vector<double>>& mat) {
    if (determinantOfMatrix(mat) == 0) {
        vector<vector<double>> result;
        result.push_back(vector<double>{0.0});
        return result;
    }
    int tam = mat.size();
    mat = aumenta(mat, tam);
    mat = escalonamento(mat, tam);
    mat = arruma(mat, tam);

    return mat;
}    
};

int main() {
    int tam;
    cout << "Tamanho da matriz? ";
    cin >> tam;

    vector<vector<double>> mat(tam, vector<double>(tam));

    for (int i = 0; i < tam; ++i) {
        cout << "Linha " << i + 1 << ": ";
        for (int j = 0; j < tam; ++j) {
            cin >> mat[i][j];
        }
    }
    matrizInversa matrizInv;
    vector<vector<double>> result = matrizInv.inversa(mat);

    if (result[0][0] == 0) {
        cout << "A matriz nao eh inversivel" << endl;
    } else {
        mat = result;
        for (int i = 0; i < tam; ++i) {
            for (int j = tam; j < tam * 2; ++j) {
                if (mat[i][j] == 0) {
                    mat[i][j] = 0.0;
                }
                cout << mat[i][j] << " \t";
            }
            cout << endl;
        }
    }

    return 0;
}