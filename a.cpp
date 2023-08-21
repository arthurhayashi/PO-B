#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

class funcaoEx
{
public:
    vector<vector<double> > Multiplicacao_matrizes(const vector<vector<double> > &matrizA, const vector<vector<double> > &matrizB)
    {
        int rowsA = matrizA.size();
        int colsA = matrizA[0].size();
        int colsB = matrizB[0].size();

        vector<vector<double> > matrizC(rowsA, vector<double>(colsB, 0.0));

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

    vector<vector<double> > Transposta(const vector<vector<double> > &matriz)
    {
        int rows = matriz.size();
        int cols = matriz[0].size();

        vector<vector<double> > matrizTransposta(cols, vector<double>(rows, 0.0));

        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                matrizTransposta[j][i] = matriz[i][j];
            }
        }

        return matrizTransposta;
    }

    double Determinante(const vector<vector<double> > &m)
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
            vector<vector<double> > submatrix;
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
    vector<vector<double> > Inversa(const vector<vector<double> > &matrizA, const vector<double> &independentes)
    {
        int n = matrizA.size();
        funcaoEx funcao;
        // Check if determinant is not equal to 0
        if (funcao.Determinante(matrizA) == 0)
        {
            vector<vector<double> > result;
            result.push_back(vector<double>{0.0});
            return result;
        }

        vector<vector<double> > inv(n, vector<double>(n, 0.0));
        vector<vector<double> > mat(matrizA);

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

    vector<vector<double> > Inversa_(const vector<vector<double> > &matriz)
    {
        return Inversa(matriz, vector<double>{}); // Empty vector for consistency
    }
};

class submatriz
{
public:
    vector<vector<double> > Cria_submatriz(const vector<vector<double> > &matrizA, const vector<int> &vetorX)
    {
        int numRows = matrizA.size();
        int numCols = vetorX.size();

        vector<vector<double> > submatriz(numRows, vector<double>(numCols, 0.0));

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
        vector<vector<double> > matrizA;
        vector<int> basicas;
        vector<int> naoBasicas;
        vector<double> independentes;
    };

    SeparacaoResult Separacao_da_matriz(const vector<double> &funcaoZ, const vector<vector<double> > &funcoes)
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
    vector<double> Calculo_x_relativo(const vector<vector<double> > &BInversa, const vector<double> &b)
    {
        funcaoEx func;
        vector<vector<double> > bTransposta = func.Transposta(vector<vector<double> >({b}));
        vector<vector<double> > xRelativoBasico = func.Multiplicacao_matrizes(BInversa, bTransposta);

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
    vector<double> Calcula_lambda(const vector<double> &custoBasico, const vector<vector<double> > &basicaInversa)
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
    vector<double> Custos_Relativos(const vector<double> &lambdaSimplex, const vector<double> &custoNaoBasico, const vector<vector<double> > &matrizNaoBasica)
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
        float menor = *std::min_element(custoRelativoNaoBasico.begin(), custoRelativoNaoBasico.end());
        int index = std::distance(custoRelativoNaoBasico.begin(), std::find(custoRelativoNaoBasico.begin(), custoRelativoNaoBasico.end(), menor));
        return index;
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
    vector<double> Direcao_simplex(const vector<vector<double> > &BasicaInversa, const vector<vector<double> > &matrizA, int k, const vector<int> &naoBasicas)
    {
        int numRows = matrizA.size();
        funcaoEx fun;
        vector<double> colunaK(numRows);

        for (int i = 0; i < numRows; ++i)
        {
            colunaK[i] = matrizA[i][naoBasicas[k]];
        }

        vector<vector<double> > colunaKMat(1, colunaK); // Convert colunaK to a matrix

        vector<double> y = fun.Multiplicacao_matrizes(BasicaInversa, fun.Transposta(colunaKMat))[0]; // Extract the first row of the result

        return y;
    }
};

class CalculaL
{
public:
    int Calcula_l(const std::vector<double>& y, const std::vector<double>& xRelativoBasico) {
    const double MAXINT = std::numeric_limits<double>::max(); // Defina MAXINT como o valor máximo de float
    bool seguro = false;
    for (size_t i = 0; i < y.size(); ++i) {
        if (y[i] > 0) {
            seguro = true;
            break;
        }
    }
    if (!seguro) {
        return -1; // Retorna -1 para indicar que não há solução ótima finita
    }
    std::vector<float> razoes;
    for (size_t i = 0; i < xRelativoBasico.size(); ++i) {
        if (y[i] <= 0) {
            razoes.push_back(MAXINT);
        } else {
            razoes.push_back(xRelativoBasico[i] / y[i]);
        }
    }
    float passo = *std::min_element(razoes.begin(), razoes.end());
    int l = std::distance(razoes.begin(), std::find(razoes.begin(), razoes.end(), passo));

    return l;
}
};

class TrocaKAndL
{
public:
    pair<vector<int>, vector<int> > Troca_k_l(const vector<int> &basicas, const vector<int> &naoBasicas, int k, int l)
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

class leitura {
public:
    struct SeparacaoLeituras {
        vector<double> funcaoZ;
        string minMax;
        vector<vector<double>> funcoes;
    };

    SeparacaoLeituras Leituras() {
        SeparacaoLeituras leitura;
        string inputString;
        int numeroFuncoes;

        cout << "Digite a função z separada por espaços (2 -4 3): ";
        getline(cin, inputString);
        parseInputToVector(inputString, leitura.funcaoZ, true);

        cout << "Digite 'min' para minimizar ou 'max' para maximizar: ";
        cin >> leitura.minMax;

        cout << "Digite o número de funções: ";
        cin >> numeroFuncoes;

        cout << "Insira as funções separadas por enter:" << endl;
        cin.ignore(); // Clear newline character
        for (int i = 0; i < numeroFuncoes; ++i) {
            cout << "Digite a função " << i + 1 << " separada por espaços (2 -4 3 <= 5): ";
            getline(cin, inputString);
            vector<double> novaFuncao;
            parseInputToVector(inputString, novaFuncao, false);
            leitura.funcoes.push_back(novaFuncao);
        }

        return leitura;
    }

    void parseInputToVector(const string& inputString, vector<double>& outputVector, bool allowSymbols) {
        istringstream iss(inputString);
        string token;
        while (iss >> token) {
            if (allowSymbols) {
                outputVector.push_back(parseToken(token));
            } else {
                try {
                    double value = parseToken(token);
                    outputVector.push_back(value);
                } catch (const std::exception&) {
                    cout << "Entrada inválida: " << token << endl;
                    outputVector.clear(); // Clear the vector if parsing fails
                    break;
                }
            }
        }
    }

    double parseToken(const string& token) {
        try {
            return stod(token);
        } catch (const std::exception&) {
            return 0.0; // Treat symbols as 0.0
        }
    }
};

/// @brief
/// @return
int main()
{
    //Criação de objetos auxiliares
    submatriz SubMatrizAux;
    matrizInversa inversaAux;
    calculoRelativo RelativoAux;
    fcusto CustoAux;
    calculaLamb LambidaAux;
    custosRel CustoRelativoAux;
    calculk CalculaAux;
    otimil OtimilidadeAux;
    direcao DirecaoAux;
    CalculaL CalculaLAux;
    TrocaKAndL TrocaKAux;
    ValorFuncao ValorAux;
    leitura LeituraAux;
    separaMatriz MatrizAux;

    leitura::SeparacaoLeituras leitura = LeituraAux.Leituras();                                             
    cout <<"Passou por aqui Leitura"<< endl;
    separaMatriz::SeparacaoResult matriz = MatrizAux.Separacao_da_matriz(leitura.funcaoZ, leitura.funcoes);
    cout <<"Passou por aqui Separação"<< endl;
    
    separaMatriz::SeparacaoResult indFixo;
    indFixo.independentes = matriz.independentes;
    cout <<"Passou por aqui Indice Fixo"<< endl;
    
    leitura::SeparacaoLeituras funcaoFin;
    funcaoFin.funcaoZ = leitura.funcaoZ; 
    cout <<"Passou por aqui Func Fin"<< endl;

    int tam = leitura.funcaoZ.size(); 
    cout <<"Passou por aqui Tam"<< endl;

    if (leitura.minMax == "max")
    {
        for (int i = 0; i < tam; i++)
        {
            leitura.funcaoZ[i] *= -1;
        }
    }

    int it = 1;
    int maxit = 10;
    vector<double> solucaoOtima = {0};
    bool funciona = true;

    while (it < maxit)
    {
        cout << endl;
        matriz.independentes = indFixo.independentes;
        cout << "it: " << it << endl;
        vector<vector<double>> matrizBasicas = SubMatrizAux.Cria_submatriz(matriz.matrizA, matriz.basicas);
        vector<vector<double>> matrizNaoBasicas = SubMatrizAux.Cria_submatriz(matriz.matrizA, matriz.naoBasicas);
        vector<vector<double> > matrizBasicaInversa = inversaAux.Inversa(matriz.matrizA, matriz.independentes);
        if (matrizBasicaInversa != matrizBasicaInversa)
        {
            funciona = false;
            break;
        }
        vector<double> xRelativo = RelativoAux.Calculo_x_relativo(matrizBasicaInversa, matriz.independentes);
        vector<double> custoBasico = CustoAux.Custo(leitura.funcaoZ, matriz.basicas);
        vector<double> lambdaTransposto = LambidaAux.Calcula_lambda(custoBasico, matrizBasicaInversa);
        vector<double> custoNaoBasico = CustoAux.Custo(leitura.funcaoZ, matriz.naoBasicas);
        vector<double> custoRelativoNaoBasico = CustoRelativoAux.Custos_Relativos(lambdaTransposto, custoNaoBasico, matrizNaoBasicas);
        int k = CalculaAux. Calcula_k(custoRelativoNaoBasico);

        if (OtimilidadeAux.Otimalidade(custoRelativoNaoBasico, k))
        {
            cout << "Otimo!" << endl;
            solucaoOtima = xRelativo;
            funciona = true;
            break;
        }
        cout << "Não otimo!" << endl;
        
        vector<double> y = DirecaoAux.Direcao_simplex(matrizBasicaInversa, matriz.matrizA, k, matriz.naoBasicas);
       
        int l = CalculaLAux.Calcula_l(y, xRelativo);
        if (l == false)
        {
            funciona = false;
            break;
        }

        pair<vector<int>, vector<int> > TrocaAux = TrocaKAux.Troca_k_l(matriz.basicas, matriz.naoBasicas, k, l);
        matriz.basicas = TrocaAux.first;
        matriz.naoBasicas = TrocaAux.second;
        it += 1;
    }
    // fim do laco de repeticao simplex

    if (funciona)
    {
        cout << "A solucao factivel otima é:" << endl;
        tam = solucaoOtima.size();
        for (int i = 0; i < tam; i++)
            cout << "x"<< matriz.basicas[i] <<" = "<<solucaoOtima[i] << "end= ";
        cout << "z = "<< ValorAux.Valor_funcao(funcaoFin.funcaoZ, solucaoOtima, matriz.basicas) << endl;
    }
    else
    {
        cout << "Em algum momento nao foi possivel fazer a inversa, o determinante deu 0 ou a direcao simplex deu <= 0" << endl;
    }
}
