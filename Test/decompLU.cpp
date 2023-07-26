#include <iostream>
#include <vector>

#include ".\Utils\auxiliar.h"

using namespace std;

vector<vector<double>> l;
vector<vector<double>> u;

int DecompLU(const vector<vector<double>> &matriz)
{
    l = vector<vector<double>>(matriz.size(), vector<double>(matriz.size(), 0));
    u = vector<vector<double>>(matriz.size(), vector<double>(matriz.size(), 0));

    for (size_t i = 0; i < matriz.size(); i++)
    {
        l[i][i] = 1;
        for (size_t j = 0; j < matriz.size(); j++)
        {
            u[i][j] = matriz[i][j];
        }
    }

    for (size_t i = 0; i < matriz.size(); i++)
    {
        // U
        for (size_t k = i; k < matriz.size(); k++)
        {
            double soma = 0;

            for (size_t j = 0; j < i; j++)
            {
                soma += l[i][j] * u[j][k];
            }

            u[i][k] = matriz[i][k] - soma;

            if (u[i][i] == 0)
            {
                matriz = trocaLinhas(matriz, i, 1);
                i = 0;
            }
        }

        // L
        for (size_t k = i; k < matriz.size(); k++)
        {
            if (i != k)
            {
                double soma = 0;
                for (size_t j = 0; j < i; j++)
                {
                    soma += (l[k][j] * u[j][i]);
                }

                l[k][i] = (matriz[k][i] - soma) / u[i][i];
            }
        }
    }
    return 0;
}

vector<double> resolveSistema(const vector<double> &b)
{
    vector<double> y(b.size());
    vector<double> x(b.size());

    int i, j;

    // Ly = b
    for (i = 0; i < static_cast<int>(b.size()); i++)
    {
        double soma = 0;

        for (j = 0; j <= i - 1; j++)
        {
            soma += l[i][j] * y[j];
        }

        y[j] = b[i] - soma;
    }

    // Ux = y
    for (i = static_cast<int>(b.size()) - 1; i >= 0; i--)
    {
        double soma = 0;

        for (j = i + 1; j < static_cast<int>(b.size()); j++)
        {
            soma += u[i][j] * x[j];
        }

        x[i] = (y[i] - soma) / u[i][i];
    }

    return x;
}