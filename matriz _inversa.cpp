// Alunos: Arthur Hayashi e Vitor Siqueira
// Trabalho Computacional 1 Matriz Inversa

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// Função para imprimir a matriz
void imprimirMatriz(const vector<vector<double>>& matriz) {
    int linhas = matriz.size();
    int colunas = matriz[0].size();

    for (int i = 0; i < linhas; i++) {
        for (int j = 0; j < colunas; j++) {
            cout << matriz[i][j] << "\t";
        }
        cout << endl;
    }
}

// Função para permutar linhas da matriz
void permutarLinhas(vector<vector<double>>& matriz, int i, int j) {
    int colunas = matriz[0].size();

    for (int k = 0; k < colunas; k++) {
        double temp = matriz[i][k];
        matriz[i][k] = matriz[j][k];
        matriz[j][k] = temp;
    }
}

// Função para calcular a matriz inversa
vector<vector<double>> calcularMatrizInversa(vector<vector<double>>& matriz) {
    int n = matriz.size();

    // Criar uma matriz expandida com a matriz original e uma matriz identidade
    vector<vector<double>> matrizExpandida(n, vector<double>(2 * n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrizExpandida[i][j] = matriz[i][j];
        }
        matrizExpandida[i][i + n] = 1.0;
    }

    // Aplicar a eliminação de Gauss
    for (int i = 0; i < n; i++) {
        if (matrizExpandida[i][i] == 0.0) {
            // Se o pivô é zero, permuta a linha com uma linha não nula abaixo
            bool encontrouLinhaNaoNula = false;
            for (int j = i + 1; j < n; j++) {
                if (matrizExpandida[j][i] != 0.0) {
                    permutarLinhas(matrizExpandida, i, j);
                    encontrouLinhaNaoNula = true;
                    break;
                }
            }
            if (!encontrouLinhaNaoNula) {
                cout << "A matriz não é invertível." << endl;
                return matrizExpandida;
            }
        }

        for (int j = 0; j < n; j++) {
            if (i != j) {
                double ratio = matrizExpandida[j][i] / matrizExpandida[i][i];

                for (int k = 0; k < 2 * n; k++) {
                    matrizExpandida[j][k] -= ratio * matrizExpandida[i][k];
                }
            }
        }
    }

    // Normalizar as linhas da matriz expandida
    for (int i = 0; i < n; i++) {
        double divisor = matrizExpandida[i][i];
        for (int j = 0; j < 2 * n; j++) {
            matrizExpandida[i][j] /= divisor;
        }
    }

    // Extrair a matriz inversa da matriz expandida
    vector<vector<double>> matrizInversa(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrizInversa[i][j] = matrizExpandida[i][j + n];
        }
    }

    return matrizInversa;
}

int main() {
    int n;
    cout << "Digite a ordem da matriz: ";
    cin >> n;

    vector<vector<double>> matriz(n, vector<double>(n, 0.0));

    cout << "Digite os elementos da matriz: " << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> matriz[i][j];
        }
    }

    cout << "Matriz original:" << endl;
    imprimirMatriz(matriz);

    vector<vector<double>> matrizInversa = calcularMatrizInversa(matriz);

    cout << "Matriz inversa:" << endl;
    imprimirMatriz(matrizInversa);

    return 0;
}