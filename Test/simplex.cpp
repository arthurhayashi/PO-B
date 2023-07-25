#include <iostream>
#include <vector>
#include <cmath>

#include ".\Utils\auxiliar.cpp"
#include "decompLU.cpp"
#include "matriz.cpp"

using namespace std;

bool max;
int tamN;
vector<double> cB, cR, b, n, xChapeu;
vector<vector<double>> B, N, BInversa;

int Simplex(bool max, vector<double> CB, vector<vector<double>> matriz, vector<int> simbolos)
{
    cB = CB;

    if (max)
    {
        for (size_t i = 0; i < cB.size(); i++)
        {
            cB[i] = cB[i] * (-1);
        }
    }

    tamN = 0;

    for (size_t i = 0; i < simbolos.size(); i++)
    {
        if (simbolos[i] != 0)
        {
            tamN++;
        }
    }

    B = vector<vector<double>>(matriz.size(), vector<double>(matriz[0].size(), 0));
    N = vector<vector<double>>(matriz.size(), vector<double>(matriz[0].size() - 1, 0));

    b = vector<double>(matriz.size(), 0);
    n = vector<double>(matriz[0].size() - 1, 0);

    for (size_t i = 0; i < matriz.size(); i++)
    {
        for (size_t j = 0; j < matriz[0].size(); j++)
        {
            if (j < matriz.size() - 1)
            {
                B[i][j] = matriz[i][j];
            }
            else
            {
                b[i] = matriz[i][j];
            }
        }
    }

    B[matriz.size() - 1][matriz.size() - 1] = simbolos[simbolos.size() - 1]; // DEF [x1 x2 xn]

    // DEF [x3 xi (i < n)]
    for (size_t i = 0; i < N[0].size(); i++)
    {
        N[i][i] = simbolos[i];
    }

    cR = resolver(simbolos, 1);
    return 0;
}

void trocaColunas(int col1, int col2)
{
    vector<vector<double>> BN(B.size(), vector<double>(B[0].size() + N[0].size(), 0));

    for (size_t i = 0; i < BN.size(); i++)
    {
        for (size_t j = 0; j < BN[0].size(); j++)
        {
            if (j < B.size())
            {
                BN[i][j] = B[i][j];
            }
            else
            {
                BN[i][j] = N[i][j - B.size()];
            }
        }
    }

    BN = trocaColuna(BN, col1, col2);

    for (size_t i = 0; i < BN.size(); i++)
    {
        for (size_t j = 0; j < BN[0].size(); j++)
        {
            if (j < B.size())
            {
                B[i][j] = BN[i][j];
            }
            else
            {
                N[i][j - B.size()] = BN[i][j];
            }
        }
    }
}

vector<double> resolver(const vector<int> &simbolos, int it)
{
    cout << "teste " << it << endl;

    BInversa = inversa(B);

    Print2double("B:", B);
    Print2double("N:", N);
    Print1double("b:", b);
    Print1double("n:", n);
    Print1double("cB", cB);
    Print2double("B⁻¹:", BInversa);

    xChapeu = concatena(mult(b, BInversa), n); // b_chapeu e n_chapeu

    Print1double("xChapeu", xChapeu);

    for (size_t i = 0; i < B.size(); i++)
    {
        if (xChapeu[i] < 0)
        {
            // throw runtime_error("Sistema Ax = b fere a condição de não-negatividade");
        }
    }

    vector<double> lambda_t = mult(cB, BInversa);
    Print1double("lambda", lambda_t);

    vector<double> res = mult(lambda_t, N);
    vector<double> y(N[0].size(), 0);

    int kMenor = -1;
    double menor = numeric_limits<double>::infinity();

    for (size_t i = 0; i < y.size(); i++)
    {
        y[i] = xChapeu[i + tamN] - res[i];

        if (y[i] < menor)
        {
            menor = y[i];
            kMenor = i;
        }
    }

    Print1double("y: ", y);

    vector<double> direcao(B.size(), 0);

    for (size_t i = 0; i < B.size(); i++)
    {
        direcao[i] = N[i][kMenor];
    }

    direcao = mult(direcao, BInversa);

    Print1double("direcao: ", direcao);

    if (menor < 0)
    {
        if (it > fatorial(B.size()) || it <= 0)
        {
            throw runtime_error("Sistema sem solução");
        }

        it++;

        vector<double> tamanho_passo(B.size(), 0);

        int kMin = -1;
        double min = numeric_limits<double>::infinity();

        for (size_t i = 0; i < tamanho_passo.size(); i++)
        {
            if (B[i][kMenor] > 0 && direcao[i] != 0 && (xChapeu[i] / direcao[i]) >= 0)
            {
                tamanho_passo[i] = xChapeu[i] / direcao[i];
                if (tamanho_passo[i] < min)
                {
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
