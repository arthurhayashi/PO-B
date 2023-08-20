#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

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
    vector<double> Custo(const vector<double> &funcaoZ, const vector<int> &variaveis)
    {
        vector<double> custoBasico(variaveis.size(), 0.0);
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
    vector<double> Calcula_lambda(const vector<double> &custoBasico, const vector<vector<double>> &basicaInversa)
    {
        int rows = basicaInversa.size();
        int cols = basicaInversa[0].size();

        vector<double> lambdaSimplex(cols, 0.0);

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

class custosRel
{
public:
    vector<double> Custos_Relativos(const vector<double> &lambdaSimplex, const vector<double> &custoNaoBasico, const vector<vector<double>> &matrizNaoBasica)
    {
        int rows = matrizNaoBasica.size();
        int cols = matrizNaoBasica[0].size();

        vector<double> custoRelativoNaoBasico = custoNaoBasico;

        for (int j = 0; j < cols; ++j)
        {
            for (int k = 0; k < rows; ++k)
            {
                custoRelativoNaoBasico[j] -= lambdaSimplex[k] * matrizNaoBasica[k][j];
            }
        }

        return custoRelativoNaoBasico;
    }
};

class calculk
{
public:
    int Calcula_k(const vector<double> &custoRelativoNaoBasico)
    {
        auto minElementIt = min_element(custoRelativoNaoBasico.begin(), custoRelativoNaoBasico.end());
        int k = distance(custoRelativoNaoBasico.begin(), minElementIt);
        return k;
    }
};

class otimil
{
public:
    bool Otimalidade(const vector<double> &custoRelativoNaoBasico, int k)
    {
        if (custoRelativoNaoBasico[k] >= 0)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
};

class direcao
{
public:
    vector<double> Direcao_simplex(const vector<vector<double>> &BasicaInversa, const vector<vector<double>> &matrizA, int k, const vector<int> &naoBasicas)
    {
        int numRows = matrizA.size();
        funcaoEx fun;
        vector<double> colunaK(numRows);

        for (int i = 0; i < numRows; ++i)
        {
            colunaK[i] = matrizA[i][naoBasicas[k]];
        }

        vector<vector<double>> colunaKMat(1, colunaK); // Convert colunaK to a matrix

        vector<double> y = fun.Multiplicacao_matrizes(BasicaInversa, fun.Transposta(colunaKMat))[0]; // Extract the first row of the result

        return y;
    }
};

class CalculaL
{
public:
    int Calcula_l(const vector<double> &y, const vector<double> &xRelativoBasico)
    {
        const double MAXINT = numeric_limits<double>::max(); // Define the maximum value for comparison

        // Check if y is not greater than 0
        bool seguro = false;
        for (size_t i = 0; i < y.size(); ++i)
        {
            if (y[i] > 0)
            {
                seguro = true;
                break;
            }
        }
        if (!seguro)
        {
            return -1; // Returning -1 to indicate "false"
        }

        vector<double> razoes;
        for (size_t i = 0; i < xRelativoBasico.size(); ++i)
        {
            if (y[i] <= 0)
            {
                razoes.push_back(MAXINT);
            }
            else
            {
                razoes.push_back(xRelativoBasico[i] / y[i]);
            }
        }

        double passo = *min_element(razoes.begin(), razoes.end());
        int l = distance(razoes.begin(), find(razoes.begin(), razoes.end(), passo));

        return l;
    }
};

class TrocaKAndL
{
public:
    pair<vector<int>, vector<int>> Troca_k_l(const vector<int> &basicas, const vector<int> &naoBasicas, int k, int l)
    {
        vector<int> newBasicas = basicas;
        vector<int> newNaoBasicas = naoBasicas;

        int aux = newBasicas[l];
        newBasicas[l] = newNaoBasicas[k];
        newNaoBasicas[k] = aux;

        return make_pair(newBasicas, newNaoBasicas);
    }
};

class ValorFuncao
{
public:
    double Valor_funcao(const vector<double> &funcaoZ, const vector<double> &xRelativoBasico, const vector<int> &basicas)
    {
        double resultado = 0.0;
        for (size_t i = 0; i < xRelativoBasico.size(); ++i)
        {
            resultado += funcaoZ[basicas[i]] * xRelativoBasico[i];
        }
        return resultado;
    }
};

class leitura
{
public:
    struct SeparacaoLeituras
    {
        vector<double> funcaoZ;
        string minMax;
        vector<vector<double>> funcoes;
    };

    SeparacaoLeituras Leituras()
    {
        SeparacaoLeituras leitura;
        string entrada;
        int numeroFuncoes;

        cout << "digite a funcao z separada por espacos (2 -4 3): ";
        getline(cin, entrada);
        int pos = 0;
        while ((pos = entrada.find(' ')) != string::npos)
        {
            leitura.funcaoZ.push_back(stod(entrada.substr(0, pos)));
            entrada.erase(0, pos + 1);
        }
        leitura.funcaoZ.push_back(stod(entrada));
        cout << endl;

        cout << "digite min para minimizar e max para maximizar: ";
        cin >> leitura.minMax;
        cout << endl;

        cout << "qual o numero de funcoes? ";
        cin >> numeroFuncoes;
        cout << endl;
        cout << "insira as funcoes separadas por enter:" << endl;

        for (int i = 0; i < numeroFuncoes; ++i)
        {
            string novaFuncao;
            getline(cin, novaFuncao);
            vector<double> funcaoAtual;
            int pos = 0;

            while ((pos = novaFuncao.find(' ')) != string::npos)
            {
                if (novaFuncao.substr(0, pos) != novaFuncao.substr(novaFuncao.size() - 2))
                {
                    funcaoAtual.push_back(stod(novaFuncao.substr(0, pos)));
                }
                novaFuncao.erase(0, pos + 1);
            }
            leitura.funcoes.push_back(funcaoAtual);
        }

        return leitura;
    }
};
/// @brief
/// @return
int main()
{
    leitura::SeparacaoLeituras leitura = leitura::Leituras();                                                   // funcaoZ, funcoes, minMax = Leitura();
    separaMatriz::SeparacaoResult matriz = separaMatriz::Separacao_da_matriz(leitura.funcaoZ, leitura.funcoes); // matrizA, basicas, naoBasicas, independentes = Separacao_da_matriz(funcaoZ, funcoes)

    // indFixo = deepcopy(independentes)
    separaMatriz::SeparacaoResult indFixo;
    indFixo.independentes = matriz.independentes;

    // funcaoFin = deepcopy(funcaoZ)
    leitura::SeparacaoLeituras funcaoFin;
    funcaoFin.funcaoZ = leitura.funcaoZ;

    // tam = len(funcaoZ)
    int tam = leitura.funcaoZ.size();

    // if(minMax == 'max'):
    //     for i in range(tam):
    //         funcaoZ[i] *= -1
    if (leitura.minMax == 'max')
    {
        for (int i = 0; i < tam; i++)
        {
            leitura.funcaoZ[i] *= -1;
        }
    }

    // it = 1
    int it = 1;

    // maxit = 10
    int maxit = 10;

    // solucaoOtima = []
    vector<double> solucaoOtima = [];

    // deu = True
    bool funciona = true;

    // while(it < maxit):
    while (it < maxit)
    {
        //     print()
        //     independentes = indFixo
        //     print('it: ', it)
        //     matrizBasica = Cria_submatriz(matrizA, basicas)
        //     matrizNaoBasica = Cria_submatriz(matrizA, naoBasicas)
        //     //print('basica: ',matrizBasica)
        cout << endl;
        matriz.independentes = indFixo.independentes;
        cout << "it: " << it << endl;
        vector<double> matrizBasicas = submatriz::Cria_submatriz(matriz.matrizA, matriz.basicas);
        vector<double> matrizNaoBasicas = submatriz::Cria_submatriz(matriz.matrizA, matriz.naoBasicas);
        //     matrizBasicaInversa, independentes = Inversa(
        //         matrizBasica, independentes)
        //     if(matrizBasicaInversa == False):
        //         deu = False
        //         break
        vector<vector<double>> matrizBasicaInversa = matrizInversa::Inversa(matriz.matrizA, matriz.independentes);
        if (matrizBasicaInversa == false)
        {
            funciona = false;
            break;
        }
        //     xRelativo = Calculo_x_relativo(matrizBasicaInversa, independentes)
        vector<double> xRelativo = calculoRelativo::Calculo_x_relativo(matrizBasicaInversa, matriz.independentes);
        //     custoBasico = Custo(funcaoZ, basicas)
        //     lambdaTransposto = Calcula_lambda(custoBasico, matrizBasicaInversa)
        vector<double> custoBasico = fcusto::Custo(leitura.funcaoZ, matriz.basicas);
        vector<double> lambdaTransposto = calculaLamb::Calcula_lambda(custoBasico, matrizBasicaInversa);
        //     custoNaoBasico = Custo(funcaoZ, naoBasicas)
        //     custoRelativoNaoBasico = Custos_Relativos(lambdaTransposto, custoNaoBasico, matrizNaoBasica)
        vector<double> custoNaoBasico = fcusto::Custo(leitura.funcaoZ, matriz.naoBasicas);
        vector<double> custoRelativoNaoBasico = custosRel::Custos_Relativos(lambdaTransposto, custoNaoBasico, matrizNaoBasicas);
        //     k = Calcula_k(custoRelativoNaoBasico)
        int k = calculk::Calcula_k(custoRelativoNaoBasico);
            //     if(Otimalidade(custoRelativoNaoBasico, k)):
            //         print("Otimo!")
            //         solucaoOtima = xRelativo
            //         deu = True
            //         break
            //     print("Nao otimo!")
        if (otimil::Otimalidade(custoRelativoNaoBasico, k))
        {
            cout << "Otimo!" << endl;
            solucaoOtima = xRelativo;
            funciona = true;
            break;
        }
        cout << "Não otimo!" << endl;
        //     y = Direcao_simplex(matrizBasicaInversa, matrizA, k, naoBasicas)
        vector<double> y = direcao::Direcao_simplex(matrizBasicaInversa, matriz.matrizA, k, matriz.naoBasicas);
        //     l = Calcula_l(y, xRelativo)
        //     if(isinstance(l, bool) and l == False):
        //         deu = False
        //         break
        int l = CalculaL::Calcula_l(y, xRelativo);
        if (l == false)
        {
            funciona = false;
            break;
        }
        //     basicas, naoBasicas = Troca_k_l(basicas, naoBasicas, k, l)
        pair<vector<int>, vector<int>> TrocaAux = TrocaKAndL::Troca_k_l(matriz.basicas, matriz.naoBasicas, k, l);
        matriz.basicas = TrocaAux.first;
        matriz.naoBasicas = TrocaAux.second;
        //     it += 1
        it += 1;
    }
    // # fim do laco de repeticao simplex

    // if(deu):
    //     print("A solucao factivel otima eh:")
    //     tam = len(solucaoOtima)
    //     for i in range(tam):
    //         print(f'x{basicas[i]} = {solucaoOtima[i]}, ', end=' ')
    //     print(f'z = {Valor_funcao(funcaoFin, solucaoOtima, basicas)}')
    // else:
    //     print('Em algum momento nao deu para fazer a inversa porque o determinante deu 0\nou a direcao simplex deu <= 0')
    if (funciona)
    {
        cout << "A solucao factivel otima é:" << endl;
        tam = solucaoOtima.size();
        for (int i = 0; i < tam; i++)
        {
            cout << "" << endl;
        }
        cout << "" << endl;
    }
    else
    {
        cout << "Em algum momento nao foi possivel fazer a inversa, o determinante deu 0 ou a direcao simplex deu <= 0" << endl;
    }
}
