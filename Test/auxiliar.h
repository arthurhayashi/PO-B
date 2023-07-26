#include <iostream>
#include <vector>
#include <cmath>
#include <string>

using namespace std;

bool max;
int tamN;
vector<double> cB, cR, b, n, xChapeu;
vector<vector<double> > B, N, BInversa, l, u;

vector<double> resolver(const vector<int> &simbolos, int it);
vector<double> resolveSistema(const vector<double> &b);
vector<vector<double> > trocaLinhas(vector<vector<double> > matriz, int linha, int it);
vector<vector<double> > identidade(const vector<vector<double> > &matriz);
vector<vector<double> > inversa(const vector<vector<double> > &matriz);

void Print(string txt, vector<vector<double> > &matriz)
{
    cout << txt << endl;

    for (auto &linha : matriz)
    {
        Print("$", linha);
    }

    cout << endl;
}

void Print(string txt, vector<double> &vetor)
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

bool comparaMatriz(vector<vector<double> > &matriz1, vector<vector<double> > &matriz2)
{
    if (matriz1.size() != matriz2.size())
    {
        cout <<("Tamanho diferente entre matrizes!");
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

double multiplica(vector<double> &vetor1, vector<double> &vetor2)
{
    if (vetor1.size() != vetor2.size())
    {
        cout <<("Tamanho diferente entre vetores!");
    }

    double res = 0;

    for (size_t i = 0; i < vetor1.size(); i++)
    {
        res += vetor1[i] * vetor2[i];
    }

    return res;
}

vector<double> multiplica(vector<double> &vetor, vector<vector<double> > &matriz)
{
    if (vetor.size() != matriz.size())
    {
        cout <<("Tamanho diferente entre vetor e matriz!");
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

vector<vector<double> > multiplica(vector<vector<double> > &matriz1, vector<vector<double> > &matriz2)
{
    if (matriz1.size() != matriz2[0].size())
    {
        cout <<("Tamanho diferente entre matrizes!");
    }

    vector<vector<double> > res(matriz1.size(), vector<double>(matriz1.size(), 0));

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

vector<vector<double> > pivoNulo(vector<vector<double> > matriz, int it)
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

vector<vector<double> > trocaLinhas(vector<vector<double> > matriz, int linha, int it)
{
    if (it > fatorial(matriz.size()) || it <= 0)
    {
        cout <<("Matriz inválida");
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
    Print("L" + to_string(linha + 1) + " <-> L" + to_string(troca + 1), matriz);
    // Confere se precisa de mais trocas
    pivoNulo(matriz, it);

    return matriz;
}

vector<vector<double> > trocaColuna(vector<vector<double> > matriz, int colInicial, int colFinal)
{
    for (size_t i = 0; i < matriz.size(); i++)
    {
        double aux = matriz[i][colInicial];
        matriz[i][colInicial] = matriz[i][colFinal];
        matriz[i][colFinal] = aux;
    }

    return matriz;
}

void Simplex(bool max, vector<double> CB, vector<vector<double> > matriz, vector<int> simbolos)
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

    B = vector<vector<double> >(matriz.size(), vector<double>(matriz[0].size(), 0));
    N = vector<vector<double> >(matriz.size(), vector<double>(matriz[0].size() - 1, 0));

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
    Print("cR: ", cR);
}

void trocaColunas(int col1, int col2)
{
    vector<vector<double> > BN(B.size(), vector<double>(B[0].size() + N[0].size(), 0));

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

    Print("B:", B);
    Print("N:", N);
    Print("b:", b);
    Print("n:", n);
    Print("cB", cB);
    Print("B⁻¹:", BInversa);

    xChapeu = concatena(multiplica(b, BInversa), n); // b_chapeu e n_chapeu

    Print("xChapeu", xChapeu);

    for (size_t i = 0; i < B.size(); i++)
    {
        if (xChapeu[i] < 0)
        {
            // cout <<("Sistema Ax = b fere a condição de não-negatividade");
        }
    }

    vector<double> lambda_t = multiplica(cB, BInversa);
    Print("lambda", lambda_t);

    vector<double> res = multiplica(lambda_t, N);
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

    Print("y: ", y);

    vector<double> direcao(B.size(), 0);

    for (size_t i = 0; i < B.size(); i++)
    {
        direcao[i] = N[i][kMenor];
    }

    direcao = multiplica(direcao, BInversa);

    Print("direcao: ", direcao);

    if (menor < 0)
    {
        if (it > fatorial(B.size()) || it <= 0)
        {
            cout <<("Sistema sem solução");
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

double det(const vector<vector<double> > &matriz)
{
    vector<vector<double> > mat_temp;
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
        mat_temp = vector<vector<double> >(matriz.size() - 1, vector<double>(matriz[0].size() - 1));
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

vector<vector<double> > transposta(const vector<vector<double> > &matriz)
{
    vector<vector<double> > transposta(matriz[0].size(), vector<double>(matriz.size(), 0));

    for (size_t i = 0; i < matriz.size(); i++)
    {
        for (size_t j = 0; j < matriz[0].size(); j++)
        {
            transposta[j][i] = matriz[i][j];
        }
    }
    return transposta;
}

vector<vector<double> > inversa(const vector<vector<double> > &matriz)
{
    // Tratamento para matriz singular
    if (det(matriz) == 0)
    {
        cout <<("Matriz singular!");
    }

    vector<vector<double> > inversa(matriz.size(), vector<double>(matriz[0].size(), 0));
    vector<vector<double> > id = identidade(matriz);

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

vector<vector<double> > identidade(const vector<vector<double> > &matriz)
{
    vector<vector<double> > identidade(matriz.size(), vector<double>(matriz[0].size(), 0));

    for (size_t i = 0; i < matriz.size(); i++)
    {
        identidade[i][i] = 1;
    }

    return identidade;
}

void DecompLU(const vector<vector<double> > &matriz)
{
    l = vector<vector<double> >(matriz.size(), vector<double>(matriz.size(), 0));
    u = vector<vector<double> >(matriz.size(), vector<double>(matriz.size(), 0));

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