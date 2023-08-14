#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

class funcaoEx
{
public:
    vector<vector<double>> Multiplicacao_matrizes(const vector<vector<double>> &matrizA, const vector<vector<double>> &matrizB)
    {
        int rowsA = matrizA.size();
        int colsA = matrizA[0].size();
        int colsB = matrizB[0].size();

        vector<vector<double>> matrizC(rowsA, vector<double>(colsB, 0.0));

        for (int i = 0; i < rowsA; ++i)
        {
            for (int j = 0; j < colsB; ++j)
            {
                for (int k = 0; k < colsA; ++k)
                {
                    matrizC[i][j] += matrizA[i][k] * matrizB[k][j];
                }
            }
        }

        return matrizC;
    }

    double Multiplicacao_vetores(const vector<double> &vetorA, const vector<double> &vetorB)
    {
        double resultado = 0.0;
        for (size_t i = 0; i < vetorA.size(); ++i)
        {
            resultado += vetorA[i] * vetorB[i];
        }
        return resultado;
    }

    vector<vector<double>> Transposta(const vector<vector<double>> &matriz)
    {
        int rows = matriz.size();
        int cols = matriz[0].size();

        vector<vector<double>> matrizTransposta(cols, vector<double>(rows, 0.0));

        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                matrizTransposta[j][i] = matriz[i][j];
            }
        }

        return matrizTransposta;
    }

    double Determinante(const vector<vector<double>> &m)
    {
        int n = m.size();
        if (n == 1)
        {
            return m[0][0];
        }
        if (n == 2)
        {
            return m[0][0] * m[1][1] - m[1][0] * m[0][1];
        }
        double det = 0.0;
        for (int i = 0; i < n; ++i)
        {
            vector<vector<double>> submatrix;
            for (int j = 1; j < n; ++j)
            {
                vector<double> row;
                for (int k = 0; k < n; ++k)
                {
                    if (k != i)
                    {
                        row.push_back(m[j][k]);
                    }
                }
                submatrix.push_back(row);
            }
            det += pow(-1.0, i) * m[0][i] * Determinante(submatrix);
        }
        return det;
    }
};

class matrizInversa
{
public:
    vector<vector<double>> Inversa(const vector<vector<double>> &matrizA, const vector<double> &independentes)
    {
        int n = matrizA.size();
        funcaoEx funcao;
        // Check if determinant is not equal to 0
        if (funcao.Determinante(matrizA) == 0)
        {
            vector<vector<double>> result;
            result.push_back(vector<double>{0.0});
            return result;
        }

        vector<vector<double>> inv(n, vector<double>(n, 0.0));
        vector<vector<double>> mat(matrizA);

        for (int i = 0; i < n; ++i)
        {
            inv[i][i] = 1.0;
        }

        // Performing elementary operations
        for (int i = 0; i < n; ++i)
        {
            if (mat[i][i] == 0)
            {
                int c = 1;
                while ((i + c) < n && mat[i + c][i] == 0)
                {
                    c += 1;
                }
                if ((i + c) == n)
                {
                    break;
                }

                int j = i;

                for (int k = 0; k < n; ++k)
                {
                    swap(mat[j][k], mat[j + c][k]);
                    swap(inv[j][k], inv[j + c][k]);
                }
            }

            for (int j = 0; j < n; ++j)
            {
                if (i != j)
                {
                    double p = mat[j][i] / mat[i][i];
                    for (int k = 0; k < n; ++k)
                    {
                        mat[j][k] -= mat[i][k] * p;
                        inv[j][k] -= inv[i][k] * p;
                    }
                }
            }
        }

        for (int i = 0; i < n; ++i)
        {
            if (mat[i][i] != 1)
            {
                double p = 1.0 / mat[i][i];
                for (int j = 0; j < n; ++j)
                {
                    mat[i][j] *= p;
                    inv[i][j] *= p;
                }
            }
        }

        return inv;
    }

    vector<vector<double>> Inversa_(const vector<vector<double>> &matriz)
    {
        return Inversa(matriz, vector<double>{}); // Empty vector for consistency
    }
};

class submatriz
{
public:
    vector<vector<double>> Cria_submatriz(const vector<vector<double>> &matrizA, const vector<int> &vetorX)
    {
        int numRows = matrizA.size();
        int numCols = vetorX.size();

        vector<vector<double>> submatriz(numRows, vector<double>(numCols, 0.0));

        for (int j = 0; j < numRows; ++j)
        {
            for (int i = 0; i < numCols; ++i)
            {
                submatriz[j][i] = matrizA[j][vetorX[i]];
            }
        }

        return submatriz;
    }
};

class separaMatriz
{
public:
    struct SeparacaoResult
    {
        vector<vector<double>> matrizA;
        vector<int> basicas;
        vector<int> naoBasicas;
        vector<double> independentes;
    };

    SeparacaoResult Separacao_da_matriz(const vector<double> &funcaoZ, const vector<vector<double>> &funcoes)
    {
        SeparacaoResult result;

        vector<double> modifiedFuncaoZ = funcaoZ;

        vector<double> inequacoes;
        int tam = funcoes.size();

        for (int i = 0; i < tam; ++i)
        {
            char condicao = funcoes[i][funcoes[i].size() - 2];
            if (condicao != '=')
            {
                modifiedFuncaoZ.push_back(0.0);
                if (condicao == '<')
                {
                    inequacoes.push_back(1.0);
                }
                else
                {
                    inequacoes.push_back(-1.0);
                }
            }
            else
            {
                inequacoes.push_back(0.0);
            }
        }

        int tam2 = inequacoes.size();
        int tam3 = modifiedFuncaoZ.size();

        result.matrizA.resize(tam, vector<double>(tam2, 0.0));
        for (int i = 0; i < tam; ++i)
        {
            for (int j = 0; j < tam2; ++j)
            {
                if (i == j)
                {
                    result.matrizA[i][j] = inequacoes[j];
                }
                else
                {
                    result.matrizA[i][j] = 0.0;
                }
            }
        }

        result.basicas.resize(tam3, 0);
        result.naoBasicas.resize(tam3 - tam, 0);

        for (int i = 0; i < tam3; ++i)
        {
            if (i < tam3 - tam)
            {
                result.naoBasicas[i] = i;
            }
            else
            {
                result.basicas[i - (tam3 - tam)] = i;
            }
        }

        int tam4 = funcoes.size();
        result.independentes.resize(tam4, 0.0);
        for (int i = 0; i < tam4; ++i)
        {
            result.independentes[i] = funcoes[i][funcoes[i].size() - 1];
        }

        return result;
    }
};

class calculoRelativo
{
public:
    vector<double> Calculo_x_relativo(const vector<vector<double>> &BInversa, const vector<double> &b)
    {
        funcaoEx func;
        vector<vector<double>> bTransposta = func.Transposta(vector<vector<double>>({b}));
        vector<vector<double>> xRelativoBasico = func.Multiplicacao_matrizes(BInversa, bTransposta);

        // Extract the first row of xRelativoBasico to a vector
        vector<double> xRelativo = xRelativoBasico[0];

        return xRelativo;
    }
};

class fcusto
{
public:
    vector<double> Custo(const std::vector<double> &funcaoZ, const std::vector<int> &variaveis)
    {
        std::vector<double> custoBasico(variaveis.size(), 0.0);
        for (size_t i = 0; i < custoBasico.size(); ++i)
        {
            custoBasico[i] = funcaoZ[variaveis[i]];
        }
        return custoBasico;
    }
};

class calculaLamb
{
public:
    vector<double> Calcula_lambda(const std::vector<double> &custoBasico, const std::vector<std::vector<double>> &basicaInversa)
    {
        int rows = basicaInversa.size();
        int cols = basicaInversa[0].size();

        std::vector<double> lambdaSimplex(cols, 0.0);

        for (int j = 0; j < cols; ++j)
        {
            for (int k = 0; k < rows; ++k)
            {
                lambdaSimplex[j] += custoBasico[k] * basicaInversa[k][j];
            }
        }

        return lambdaSimplex;
    }
};


int main()
{
    // Test the functions here if needed
    return 0;
}
