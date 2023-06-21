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

int main() {
    int l, c, maxmin;
    cout << "Digite 0 se o problema for de maximizar e 1 se for para minimizar: ";
    cin >> maxmin;

    cout << "Digite a quantidade de linhas e colunas da matriz considerando a seguinte base: ";
    cout << "col:  x   y   z   w    Result"  << endl; 
    cin >> l >> c;

    vector<vector<double>> matriz(l, vector<double>(c, 0.0));

    cout << "Digite os elementos da matriz como o exemplo a seguir: " << endl;
    if(maxmin == 0){
        cout << "max: -x + 2y -z -3w"  << endl; 
    }else{
        cout << "min: -x + 2y -z -3w"  << endl; 
    }
    cout << "col:  x   y   z   w    Result"  << endl; 
    // cout << " x -y   0  4w  >=  50"  << endl;
    // cout << "-x  0   z   w   =  50"   << endl;
    // cout << " 0  y  -z  -w  <=  50"  << endl;
    for (int i = 0; i < l; i++) {
        for (int j = 0; j < c; j++) {
            cin >> matriz[i][j];
        }
    }
    if(maxmin == 0){
        for (int j = 0; j < c; j++) {
            matriz[1][j] *= -1; 
        }
    }
    cout << " x  y   z   w    Result"  << endl; 
    cout << "Matriz inicial:" << endl;

    imprimirMatriz(matriz);

    return 0;
}