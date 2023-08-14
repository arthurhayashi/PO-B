#include <iostream>
#include <vector>

using namespace std;

const double EPSILON = 1e-10;

// Função para encontrar a coluna pivô
int encontrarColunaPivo(const vector<vector<double>>& tableau) {
    int m = tableau.size() - 1;
    int n = tableau[0].size() - 1;

    int colunaPivo = 1;
    for (int i = 2; i <= n; ++i) {
        if (tableau[m][i] < tableau[m][colunaPivo]) {
            colunaPivo = i;
        }
    }

    if (tableau[m][colunaPivo] >= -EPSILON) {
        return -1;  // Nenhuma coluna negativa, a solução é otimizada
    }

    return colunaPivo;
}

// Função para encontrar a linha pivô
int encontrarLinhaPivo(const vector<vector<double>>& tableau, int colunaPivo) {
    int m = tableau.size() - 1;

    int linhaPivo = -1;
    for (int i = 1; i <= m; ++i) {
        if (tableau[i][colunaPivo] > EPSILON) {
            if (linhaPivo == -1) {
                linhaPivo = i;
            } else {
                double ratio = tableau[i][0] / tableau[i][colunaPivo];
                double minRatio = tableau[linhaPivo][0] / tableau[linhaPivo][colunaPivo];
                if (ratio < minRatio) {
                    linhaPivo = i;
                }
            }
        }
    }

    return linhaPivo;
}

// Função Simplex
void simplex(vector<vector<double>>& tableau) {
    while (true) {
        int colunaPivo = encontrarColunaPivo(tableau);
        if (colunaPivo == -1) {
            break;  // Solução otimizada alcançada
        }

        int linhaPivo = encontrarLinhaPivo(tableau, colunaPivo);
        if (linhaPivo == -1) {
            cout << "Problema ilimitado." << endl;
            return;
        }

        // Fazer o pivô
        double pivo = tableau[linhaPivo][colunaPivo];
        for (int i = 0; i < tableau[0].size(); ++i) {
            tableau[linhaPivo][i] /= pivo;
        }
        for (int i = 0; i < tableau.size(); ++i) {
            if (i != linhaPivo) {
                double ratio = tableau[i][colunaPivo];
                for (int j = 0; j < tableau[0].size(); ++j) {
                    tableau[i][j] -= ratio * tableau[linhaPivo][j];
                }
            }
        }
    }
}

int main() {
    int m, n;
    cout << "Digite o número de restrições (linhas): ";
    cin >> m;
    cout << "Digite o número de variáveis de decisão (colunas): ";
    cin >> n;


    vector<vector<double>> tableau(m + 1, vector<double>(n + m + 1));
    
    cout << "Digite os coeficientes da função objetivo (linha 0): ";
    for (int j = 0; j <= n-1; ++j) {
        cin >> tableau[0][j];
        cout << j << endl;
    }

    cout << "Digite os coeficientes das restrições (linhas 1 a " << m << "):" << endl;
    for (int i = 1; i <= m; ++i) {
        cout << "Restrição " << i << ": ";
        for (int j = 0; j <= n-1; ++j) {
            cin >> tableau[i][j];
            
        }
        cout << "Lado direito " << i << ": ";
        cin >> tableau[i][n + m];
    }

    // Preencha as variáveis de folga e de excesso
    for (int i = 1; i <= m; ++i) {
        for (int j = n + 1; j <= n + m; ++j) {
            if (j - n == i) {
                tableau[i][j] = 1.0;
            } else {
                tableau[i][j] = 0.0;
            }
        }
    }

    simplex(tableau);

    cout << "Solução ótima encontrada:" << endl;
    cout << "Valor da função objetivo: " << -tableau[0][0] << endl;
    cout << "Valores das variáveis de decisão:" << endl;
    for (int j = 1; j <= n; ++j) {
        bool isBasic = false;
        int basicRow = -1;
        for (int i = 1; i <= m; ++i) {
            if (tableau[i][j] == 1.0 && basicRow == -1) {
                isBasic = true;
                basicRow = i;
            } else if (tableau[i][j] != 0.0) {
                isBasic = false;
                break;
            }
        }
        if (isBasic) {
            cout << "x" << j << " = " << tableau[basicRow][n + m] << endl;
        } else {
            cout << "x" << j << " = 0 (não básica)" << endl;
        }
    }

    return 0;
}
