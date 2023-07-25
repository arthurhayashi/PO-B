#include <iostream>
#include <vector>

#include "matriz.cpp"
#include "decompLU.cpp"
#include "simplex.cpp"
#include ".\Utils\auxiliar.cpp"

using namespace std;

int main()
{
    /*
        // Testes manuais
    vector<vector<double>> m1 = {
        {1, 2, 3},
        {0, 1, 0},
        {1, 0, 2}
    };

    vector<vector<double>> m2 = {
        {2, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    };

    vector<vector<double>> m3 = {
        {0, 2, 4},
        {0, 1, 0},
        {1, 0, 4}
    };

    vector<vector<double>> m4 = {
        {2, 2, 3},
        {4, 4, 2},
        {0, 5, 5}
    };

    vector<vector<double>> m6 = {
        {1, 1, 4},
        {1, 0, 3},
        {0, 1, 3.5}
    };
        */
    vector<vector<double>> A;
    vector<int> simbolos;
    vector<double> fx;

    int tam;
    cout << "Digite o tamanho da matriz: >> " << endl;
    cin >> tam;
    cout << endl;

    A = vector<vector<double>>(tam, vector<double>(tam));
    simbolos = vector<int>(tam);
    fx = vector<double>(tam);

    bool max = false; // Definindo max e min

    /*
    cout << "(min/max) \n>> ";
    string tipoFuncao;
    cin >> tipoFuncao;
    cout << endl;

    if (tipoFuncao == "max") {
        max = true;
    }
    */

    for (int i = 0; i < A.size(); i++)
    {
        cout << "Valor de x" << (i + 1) << "\n>> ";
        cin >> fx[i];
    }

    cout << endl;

    for (int i = 0; i < tam; i++)
    {
        for (int j = 0; j < tam - 1; j++)
        {
            cout << "Valor de a" << (i + 1) << (j + 1) << "\n>> ";
            cin >> A[i][j];
        }

        char aux;
        cout << "\nSinal do sistema \n>> ";
        cin >> aux;

        if (aux == '<')
        {
            simbolos[i]++;
        }
        else if (aux == '>')
        {
            simbolos[i]--;
        }

        cout << "\nValor de a" << (i + 1) << tam << "\n>> ";
        cin >> A[i][tam - 1];
    }

    cout << endl;
    cout << "det(A) = " << det(A) << "\n\n";
    Print2double("A:", A);

    Simplex(max, fx, A, simbolos);
    return 0;
}