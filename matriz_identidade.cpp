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

vector<vector<double>> calcularMatrizIdentidade(vector<vector<double>>& matriz){
    int n = matriz.size();
    vector<vector<double>> matrizIdentidade(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if(i == j){
                matrizIdentidade[i][j] = 1;
            }
        }
    }
    return matrizIdentidade;
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

    vector<vector<double>> matrizIdentidade = calcularMatrizIdentidade(matriz);

    cout << "Matriz identidade:" << endl;
    imprimirMatriz(matrizIdentidade);

    return 0;
}