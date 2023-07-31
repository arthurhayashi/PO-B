#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <limits>
using namespace std;

class Aux {
public:
    static void print(string txt, vector<vector<double>>& matriz) {
        cout << txt << endl;
        for (vector<double>& linha : matriz) {
            print("$", linha);
        }
        cout << endl;
    }

    static void print(string txt, vector<double>& vetor) {
        if (txt != "$") {
            cout << txt << endl;
        }

        cout << "[ ";
        for (double valor : vetor) {
            if (valor >= 0) {
                cout << " ";
            }

            if (valor == 0) {
                cout << "0.0 ";
            } else {
                cout << valor << " ";
            }
        }
        cout << " ]" << endl;

        if (txt != "$") {
            cout << endl;
        }
    }

    static int fat(int num) {
        int res = num;
        if (num > 1) {
            return res * fat(num - 1);
        }
        return res;
    }

    static double min(vector<double>& vetor) {
        double res = 0;
        for (double valor : vetor) {
            if (valor <= res && valor > 0) {
                res = valor;
            }
        }
        return res;
    }

    static vector<double> concatenarVetores(const vector<double>& vetor1, const vector<double>& vetor2) {
        vector<double> res;
        res.reserve(vetor1.size() + vetor2.size());

        for (double valor : vetor1) {
            res.push_back(valor);
        }

        for (double valor : vetor2) {
            res.push_back(valor);
        }
        return res;
    }


    static bool compararMatrizes(vector<vector<double>>& matriz1, vector<vector<double>>& matriz2) {
        if (matriz1.size() != matriz2.size()) {
            throw runtime_error("Tamanho diferente entre matrizes!");
        }

        for (int i = 0; i < matriz1.size(); i++) {
            if (matriz1[i].size() != matriz2[i].size()) {
                throw runtime_error("Tamanho diferente entre matrizes!");
            }

            for (int j = 0; j < matriz1[i].size(); j++) {
                if (matriz1[i][j] != matriz2[i][j]) {
                    return false;
                }
            }
        }

        return true;
    }

    static double multiplicar(vector<double>& vetor1, vector<double>& vetor2) {
        if (vetor1.size() != vetor2.size()) {
            throw runtime_error("Tamanho diferente entre vetores!");
        }

        double res = 0;

        for (int i = 0; i < vetor1.size(); i++) {
            res += vetor1[i] * vetor2[i];
        }

        return res;
    }

    static vector<double> multiplicar(vector<double>& vetor, vector<vector<double>>& matriz) {
        if (vetor.size() != matriz.size()) {
            throw runtime_error("Tamanho diferente entre vetor e matriz!");
        }

        vector<double> res(matriz[0].size(), 0.0);

        for (int i = 0; i < matriz[0].size(); i++) {
            for (int j = 0; j < vetor.size(); j++) {
                res[i] += vetor[j] * matriz[j][i];
            }
        }

        return res;
    }

    static vector<vector<double>> multiplicar(vector<vector<double>>& matriz1, vector<vector<double>>& matriz2) {
        if (matriz1.size() != matriz2[0].size()) {
            throw runtime_error("Tamanho diferente entre matrizes!");
        }

        vector<vector<double>> res(matriz1.size(), vector<double>(matriz1.size(), 0.0));

        for (int i = 0; i < matriz1.size(); i++) {
            for (int k = 0; k < matriz1.size(); k++) {
                for (int j = 0; j < matriz1.size(); j++) {
                    res[i][j] += matriz1[i][k] * matriz2[k][j];
                }
            }
        }

        return res;
    }

    static vector<vector<double>> verifPivoNulo(vector<vector<double>>& matriz, int it) {
        for (int i = 0; i < matriz.size(); i++) {
            if (matriz[i][i] == 0) {
                matriz = trocaLinhasMatriz(matriz, i, it++);
            }
        }

        return matriz;
    }

    static vector<vector<double>> trocaLinhasMatriz(vector<vector<double>>& matriz, int linha, int it) {
        if (it > fat(matriz.size()) || it <= 0) {
            throw runtime_error("Matriz inválida");
        }

        int troca = 0;

        for (int i = 0; i < matriz.size(); i++) {
            if (matriz[i][linha] != 0 && matriz[linha][i] != 0 && i != linha) {
                if (i != linha) {
                    continue;
                }

                double aux;
                troca = i;

                for (int j = 0; j < matriz[0].size(); j++) {
                    aux = matriz[i][j];
                    matriz[i][j] = matriz[linha][j];
                    matriz[linha][j] = aux;
                }
            }
        }

        print("L" + to_string(linha + 1) + " <-> L" + to_string(troca + 1), matriz);

        // Confere se precisa de mais trocas
        verifPivoNulo(matriz, it);

        return matriz;
    }

    static vector<vector<double>> trocaCol(vector<vector<double>>& matriz, int colInicial, int colFinal) {
        for (int i = 0; i < matriz.size(); i++) {
            double aux = matriz[i][colInicial];
            matriz[i][colInicial] = matriz[i][colFinal];
            matriz[i][colFinal] = aux;
        }

        return matriz;
    }
};

class DecompLU {
private:
    vector<vector<double>> l;
    vector<vector<double>> u;

public:
    DecompLU(vector<vector<double>>& matriz) {
        int n = matriz.size();
        l = vector<vector<double>>(n, vector<double>(n, 0.0));
        u = matriz;

        for (int i = 0; i < n; i++) {
            l[i][i] = 1;
        }

        for (int i = 0; i < n; i++) {
            for (int k = i; k < n; k++) {
                double soma = 0;

                for (int j = 0; j < i; j++) {
                    soma += l[i][j] * u[j][k];
                }

                u[i][k] = matriz[i][k] - soma;

                if (u[i][i] == 0) {
                    matriz = Aux::trocaLinhasMatriz(matriz, i, 1);
                    i = 0;
                }
            }

            for (int k = i; k < n; k++) {
                if (i != k) {
                    double soma = 0;
                    for (int j = 0; j < i; j++) {
                        soma += (l[k][j] * u[j][i]);
                    }

                    l[k][i] = (matriz[k][i] - soma) / u[i][i];
                }
            }
        }
    }

    vector<double> resolveSistema(vector<double>& b) {
        int n = b.size();
        vector<double> y(n, 0.0);
        vector<double> x(n, 0.0);

        // Ly = b
        for (int i = 0; i < n; i++) {
            double soma = 0;

            for (int j = 0; j <= i - 1; j++) {
                soma += l[i][j] * y[j];
            }

            y[i] = b[i] - soma;
        }

        // Ux = y
        for (int i = n - 1; i >= 0; i--) {
            double soma = 0;

            for (int j = i + 1; j < n; j++) {
                soma += u[i][j] * x[j];
            }

            x[i] = (y[i] - soma) / u[i][i];
        }

        return x;
    }
};

class Matriz {
private:
    double det(vector<vector<double>>& matriz) {
        double determinante = 0;

        if (matriz.size() == 1) {
            determinante = matriz[0][0];
            return determinante;
        }

        if (matriz.size() == 2) {
            return matriz[0][0] * matriz[1][1] - matriz[0][1] * matriz[1][0];
        }

        for (int i = 0; i < matriz.size(); i++) {
            vector<vector<double>> mat_temp(matriz.size() - 1, vector<double>(matriz[0].size() - 1, 0.0));
            for (int j = 1; j < matriz[0].size(); j++) {
                for (int k = 0; k < matriz.size(); k++) {
                    if (k < i) {
                        mat_temp[j - 1][k] = matriz[j][k];
                    } else if (k > i) {
                        mat_temp[j - 1][k - 1] = matriz[j][k];
                    }
                }
            }

            determinante += matriz[0][i] * pow(-1, i) * det(mat_temp);
        }

        return determinante;
    }

    vector<vector<double>> identidade(vector<vector<double>>& matriz) {
        vector<vector<double>> identidade(matriz.size(), vector<double>(matriz[0].size(), 0.0));

        for (int i = 0; i < matriz.size(); i++) {
            identidade[i][i] = 1;
        }

        return identidade;
    }

    vector<vector<double>> transposta(vector<vector<double>>& matriz) {
        vector<vector<double>> transposta(matriz[0].size(), vector<double>(matriz.size(), 0.0));

        for (int i = 0; i < matriz.size(); i++) {
            for (int j = 0; j < matriz[0].size(); j++) {
                transposta[j][i] = matriz[i][j];
            }
        }
        return transposta;
    }

public:
    vector<vector<double>> inversa(vector<vector<double>>& matriz) {
        // Tratamento para matriz singular
        if (det(matriz) == 0) {
            throw runtime_error("Matriz singular!");
        }

        vector<vector<double>> inversa(matriz.size(), vector<double>(matriz[0].size(), 0.0));
        vector<vector<double>> identityMat = identidade(matriz); // Rename the variable

        DecompLU dLU(matriz);

        // Resolve LUx = b
        for (int i = 0; i < identityMat[0].size(); i++) {
            vector<double> b(identityMat.size(), 0.0);

            for (int j = 0; j < identityMat.size(); j++) {
                b[j] = identityMat[j][i];
            }

            vector<double> aux = dLU.resolveSistema(b);

            for (int j = 0; j < identityMat.size(); j++) {
                inversa[j][i] = aux[j];
            }
        }

        return inversa;
    }

};

class Simplex {
private:
    bool max;
    int tamN;
    vector<double> cB, cR, b, n, xChapeu;
    vector<vector<double>> B, N, BInversa;

public:
    Simplex(bool max, vector<double>& cB, vector<vector<double>>& matriz, vector<int>& simbolos) {
        this->cB = cB;

        if (max) {
            for (int i = 0; i < cB.size(); i++) {
                cB[i] = cB[i] * (-1);
            }
        }

        tamN = 0;
        for (int i = 0; i < simbolos.size(); i++) {
            if (simbolos[i] != 0) {
                tamN++;
            }
        }

        B = vector<vector<double>>(matriz.size(), vector<double>(matriz[0].size()));
        N = vector<vector<double>>(matriz.size(), vector<double>(matriz[0].size() - 1));

        b = vector<double>(matriz.size());
        n = vector<double>(matriz[0].size() - 1);

        for (int i = 0; i < matriz.size(); i++) {
            for (int j = 0; j < matriz[0].size(); j++) {
                if (j < matriz.size() - 1) {
                    B[i][j] = matriz[i][j];
                } else {
                    b[i] = matriz[i][j];
                }
            }
        }

        B[matriz.size() - 1][matriz.size() - 1] = simbolos[simbolos.size() - 1]; // DEF [x1 x2 xn]

        // DEF [x3 xi (i < n)]
        for (int i = 0; i < N[0].size(); i++) {
            N[i][i] = simbolos[i];
        }

        // Call resolver function
        cR = resolver(simbolos, 1);
    }

private:
    void trocaColunas(int col1, int col2) {
        vector<vector<double>> BN(B.size(), vector<double>(B[0].size() + N[0].size()));

        for (int i = 0; i < BN.size(); i++) {
            for (int j = 0; j < BN[0].size(); j++) {
                if (j < B.size()) {
                    BN[i][j] = B[i][j];
                } else {
                    BN[i][j] = N[i][j - B.size()];
                }
            }
        }

        BN = Aux::trocaCol(BN, col1, col2);

        for (int i = 0; i < BN.size(); i++) {
            for (int j = 0; j < BN[0].size(); j++) {
                if (j < B.size()) {
                    B[i][j] = BN[i][j];
                } else {
                    N[i][j - B.size()] = BN[i][j];
                }
            }
        }
    }

    vector<double> resolver(const vector<int>& simbolos, int it) {
        cout << "teste " << it << endl;

        Matriz m;
        BInversa = m.inversa(B);

        Aux::print("B:", B);
        Aux::print("N:", N);
        Aux::print("b:", b);
        Aux::print("n:", n);
        Aux::print("cB", cB);
        Aux::print("B⁻¹:", BInversa);

        xChapeu = Aux::concatenarVetores(Aux::multiplicar(b, BInversa), n); // b_chapeu e n_chapeu


        Aux::print("xChapeu", xChapeu);

        for (int i = 0; i < B.size(); i++) {
            if (xChapeu[i] < 0) {
                //throw runtime_error("Sistema Ax = b fere a condição de não-negatividade");
            }
        }

        vector<double> lambda_t = Aux::multiplicar(cB, BInversa);
        Aux::print("lambda", lambda_t);

        vector<double> res = Aux::multiplicar(lambda_t, N);
        vector<double> y(N[0].size());

        int kMenor = -1;
        double menor = numeric_limits<double>::infinity();

        for (int i = 0; i < y.size(); i++) {
            y[i] = xChapeu[i + tamN] - res[i];

            if (y[i] < menor) {
                menor = y[i];
                kMenor = i;
            }
        }

        Aux::print("y: ", y);

        vector<double> direcao(B.size());

        for (int i = 0; i < B.size(); i++) {
            direcao[i] = N[i][kMenor];
        }

        direcao = Aux::multiplicar(direcao, BInversa);

        Aux::print("direcao: ", direcao);

        if (menor < 0) {
            if (it > Aux::fat(B.size()) || it <= 0) {
                throw runtime_error("Sistema sem solução");
            }

            it++;

            vector<double> tamanho_passo(B.size());

            int kMin = -1;
            double min = numeric_limits<double>::infinity();

            for (int i = 0; i < tamanho_passo.size(); i++) {
                if (B[i][kMenor] > 0 && direcao[i] != 0 && (xChapeu[i] / direcao[i]) >= 0) {
                    tamanho_passo[i] = xChapeu[i] / direcao[i];
                    if (tamanho_passo[i] < min) {
                        min = tamanho_passo[i];
                        kMin = i;
                    }
                }
            }

            trocaColunas(kMenor, kMin);

            y = resolver(simbolos, it);
        }

        return y;
    }
};

int main() {
    // Example usage:
    vector<double> cB = { 3, 1 };
    vector<vector<double>> matriz = { { 1, 1, 1 }, { 1, -1, 1 } };
    vector<int> simbolos = { 1, -1, -1, -1 }; // 1 for <=, -1 for >=

    Simplex simplex(false, cB, matriz, simbolos);

    return 0;
}
