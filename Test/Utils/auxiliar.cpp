#include <iostream>
#include <vector>
#include <cmath>
#include <string>

using namespace std;

void Print2double(string txt, vector<vector<double>> &matriz)
{
    cout << txt << endl;

    for (auto &linha : matriz)
    {
        Print1double("$", linha);
    }

    cout << endl;
}

void Print1double(string txt, vector<double> &vetor)
{
    if (txt != "$")
    {
        cout << txt << endl;
    }

    cout << "[ ";

    for (double valor : vetor)
    {
        if (valor >= 0)
        {
            cout << " ";
        }

        if (valor == 0)
        {
            cout << "0.0 ";
        }
        else
        {
            cout << valor << " ";
        }
    }

    cout << " ]" << endl;

    if (txt != "$")
    {
        cout << endl;
    }
}

int fatorial(int num)
{
    int res = num;
    if (num > 1)
    {
        return res * fatorial(num - 1);
    }
    return res;
}

double min(vector<double> &vetor)
{
    double res = 0;

    for (double valor : vetor)
    {
        if (valor <= res && valor > 0)
        {
            res = valor;
        }
    }

    return res;
}

vector<double> concatena(vector<double> &vetor1, vector<double> &vetor2)
{
    vector<double> res(vetor1.size() + vetor2.size());

    for (size_t i = 0; i < res.size(); i++)
    {
        if (i < vetor1.size())
        {
            res[i] = vetor1[i];
        }
        else
        {
            res[i] = vetor2[i - vetor1.size()];
        }
    }

    return res;
}

bool comparaMatriz(vector<vector<double>> &matriz1, vector<vector<double>> &matriz2)
{
    if (matriz1.size() != matriz2.size())
    {
        throw runtime_error("Tamanho diferente entre matrizes!");
    }

    for (size_t i = 0; i < matriz1.size(); i++)
    {
        for (size_t j = 0; j < matriz1[i].size(); j++)
        {
            if (matriz1[i][j] != matriz2[i][j])
            {
                return false;
            }
        }
    }

    return true;
}

double mult(vector<double> &vetor1, vector<double> &vetor2)
{
    if (vetor1.size() != vetor2.size())
    {
        throw runtime_error("Tamanho diferente entre vetores!");
    }

    double res = 0;

    for (size_t i = 0; i < vetor1.size(); i++)
    {
        res += vetor1[i] * vetor2[i];
    }

    return res;
}

vector<double> mult(vector<double> &vetor, vector<vector<double>> &matriz)
{
    if (vetor.size() != matriz.size())
    {
        throw runtime_error("Tamanho diferente entre vetor e matriz!");
    }

    vector<double> res(vetor.size(), 0);

    for (size_t i = 0; i < matriz[0].size(); i++)
    {
        for (size_t j = 0; j < vetor.size(); j++)
        {
            res[i] += vetor[j] * matriz[j][i];
        }
    }

    return res;
}

vector<vector<double>> mult(vector<vector<double>> &matriz1, vector<vector<double>> &matriz2)
{
    if (matriz1.size() != matriz2[0].size())
    {
        throw runtime_error("Tamanho diferente entre matrizes!");
    }

    vector<vector<double>> res(matriz1.size(), vector<double>(matriz1.size(), 0));

    for (size_t i = 0; i < matriz1.size(); i++)
    {
        for (size_t k = 0; k < matriz1.size(); k++)
        {
            for (size_t j = 0; j < matriz1.size(); j++)
            {
                res[i][j] += matriz1[i][k] * matriz2[k][j];
            }
        }
    }

    return res;
}

vector<vector<double>> pivoNulo(vector<vector<double>> matriz, int it)
{
    for (size_t i = 0; i < matriz.size(); i++)
    {
        if (matriz[i][i] == 0)
        {
            matriz = trocaLinhas(matriz, i, it++);
        }
    }

    return matriz;
}

vector<vector<double>> trocaLinhas(vector<vector<double>> matriz, int linha, int it)
{
    if (it > fatorial(matriz.size()) || it <= 0)
    {
        throw runtime_error("Matriz inválida");
    }

    int troca = 0;

    for (size_t i = 0; i < matriz.size(); i++)
    {
        if (matriz[i][linha] != 0 && matriz[linha][i] != 0 && i != linha)
        {
            if (i != linha)
            {
                continue;
            }

            double aux;
            troca = i;

            for (size_t j = 0; j < matriz[0].size(); j++)
            {
                aux = matriz[i][j];
                matriz[i][j] = matriz[linha][j];
                matriz[linha][j] = aux;
            }
        }
    }
    Print2double("L" + to_string(linha + 1) + " <-> L" + to_string(troca + 1), matriz);
    // Confere se precisa de mais trocas
    pivoNulo(matriz, it);

    return matriz;
}

vector<vector<double>> trocaColuna(vector<vector<double>> matriz, int colInicial, int colFinal)
{
    for (size_t i = 0; i < matriz.size(); i++)
    {
        double aux = matriz[i][colInicial];
        matriz[i][colInicial] = matriz[i][colFinal];
        matriz[i][colFinal] = aux;
    }

    return matriz;
}