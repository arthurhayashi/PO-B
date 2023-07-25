#include <iostream>
#include <vector>
#include <cmath>

#include ".\Utils\auxiliar.cpp"
#include "decompLU.cpp"

using namespace std;

double det(const vector<vector<double>> &matriz)
{
    vector<vector<double>> mat_temp;
    double determinante = 0;

    if (matriz.size() == 1)
    {
        determinante = matriz[0][0];
        return determinante;
    }

    if (matriz.size() == 2)
    {
        return matriz[0][0] * matriz[1][1] - matriz[0][1] * matriz[1][0];
    }

    for (size_t i = 0; i < matriz.size(); i++)
    {
        mat_temp = vector<vector<double>>(matriz.size() - 1, vector<double>(matriz[0].size() - 1));
        for (size_t j = 1; j < matriz[0].size(); j++)
        {
            for (size_t k = 0; k < matriz.size(); k++)
            {
                if (k < i)
                {
                    mat_temp[j - 1][k] = matriz[j][k];
                }
                else if (k > i)
                {
                    mat_temp[j - 1][k - 1] = matriz[j][k];
                }
            }
        }

        determinante += matriz[0][i] * pow(-1, i) * det(mat_temp);
    }

    return determinante;
}

vector<vector<double>> transposta(const vector<vector<double>> &matriz)
{
    vector<vector<double>> transposta(matriz[0].size(), vector<double>(matriz.size(), 0));

    for (size_t i = 0; i < matriz.size(); i++)
    {
        for (size_t j = 0; j < matriz[0].size(); j++)
        {
            transposta[j][i] = matriz[i][j];
        }
    }
    return transposta;
}

vector<vector<double>> inversa(const vector<vector<double>> &matriz)
{
    // Tratamento para matriz singular
    if (det(matriz) == 0)
    {
        throw runtime_error("Matriz singular!");
    }

    vector<vector<double>> inversa(matriz.size(), vector<double>(matriz[0].size(), 0));
    vector<vector<double>> id = identidade(matriz);

    // Resolve LUx = b
    for (size_t i = 0; i < id.size(); i++)
    {
        vector<double> b(id[0].size(), 0);

        for (size_t j = 0; j < id[0].size(); j++)
        {
            b[j] = id[j][i];
        }

        vector<double> aux = resolveSistema(b);

        for (size_t j = 0; j < id[0].size(); j++)
        {
            inversa[j][i] = aux[j];
        }
    }

    return inversa;
}

vector<vector<double>> identidade(const vector<vector<double>> &matriz)
{
    vector<vector<double>> identidade(matriz.size(), vector<double>(matriz[0].size(), 0));

    for (size_t i = 0; i < matriz.size(); i++)
    {
        identidade[i][i] = 1;
    }

    return identidade;
}