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
        for (int j = 0; j < n; j++) {
            if (i != j) {
                double ratio = matrizExpandida[j][i] / matrizExpandida[i][i];

                for (int k = 0; k < 2 * n; k++) {
                    matrizExpandida[j][k] -= ratio * matrizExpandida[i][k];
                }
            }
        }
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

    // Extrair a matriz inversa da matriz expandida
    vector<vector<double>> matrizInversa(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrizInversa[i][j] = matrizExpandida[i][j + n];
        }
    }

    return matrizInversa;
}
void multiplicaMatriz(vector<vector<double>>& b,vector<vector<double>>& B) {
	int linhas1 = b.size();
    int colunas1 = b.size();
	int linhas2 = B.size();
    int colunas2 = B.size();
	
	vector<vector<double>>product ={{0, 0, 0}}; 

    for (int row = 0; row < linhas1; row++) {
        for (int col = 0; col < colunas2; col++) {
            // Multiply the row of A by the column of B to get the row, column of product.
            for (int inner = 0; inner < linhas2; inner++) {
                product[row][col] += b[row][inner] * B[inner][col];
				// cout << b[row][inner] << " " << B[inner][col] << endl;
            }
            std::cout << product[row][col] << "  ";
        }
        std::cout << "\n";
    }
}


int main() {
	// N = aN1 aN2 = a3 a4
    vector<vector<double>> N = {{1,0},
								{0,1},
								{0,0}};
	//B = aB1 aB2 aB3 = a1 a2 a5
	vector<vector<double>> B = {{1,1,0},
								{1,0,1},
								{0,1,1}};
	imprimirMatriz(B);
	vector<vector<double>> bT = {{-2,-1, 0}};
	// {-2,-1,0} * {1,1,0} 
	// 			   {1,0,1} Linha por Coluna.
	// 			   {0,1,1}
	vector<vector<double>> matrizInversa = calcularMatrizInversa(B);
	imprimirMatriz(bT);	
	imprimirMatriz(matrizInversa);	
	multiplicaMatriz(bT,matrizInversa);
	vector<vector<double>> A;

	vector<vector<double>> X, Y;
	
    return 0;
}