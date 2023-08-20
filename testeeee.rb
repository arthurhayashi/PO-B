require 'matrix'

def Multiplicacao_matrizes(matrizA, matrizB)
  matrizC = matrizA * matrizB
  return matrizC
end

def Multiplicacao_vetores(vetorA, vetorB)
  resultado = Vector.elements(vetorA) * Vector.elements(vetorB)
  return resultado
end

def Transposta(matriz)
  return matriz.transpose
end

def Determinante(m)
  return Matrix[*m].determinant
end

def Inversa(matrizA, independentes)
  matrizA = Matrix[*matrizA]
  n = matrizA.row_count

  if Determinante(matrizA.to_a) == 0
    return false
  end

  inv = Matrix.build(n, n) { |i, j| i == j ? 1.0 : 0.0 }

  n.times do |i|
    if matrizA[i, i] == 0
      c = 1
      while (i + c) < n && matrizA[i + c, i] == 0
        c += 1
      end

      if (i + c) == n
        break
      end

      j = i

      n.times do |k|
        temp = matrizA[j, k]
        temp2 = inv[j, k]

        matrizA[j, k] = matrizA[j + c, k]
        inv[j, k] = inv[j + c, k]

        matrizA[j + c, k] = temp
        inv[j + c, k] = temp2
      end
    end

    n.times do |j|
      next if i == j

      p = matrizA[j, i] / matrizA[i, i]

      n.times do |k|
        matrizA[j, k] -= matrizA[i, k] * p
        inv[j, k] -= inv[i, k] * p
      end
    end

    if matrizA[i, i] != 1
      p = 1.0 / matrizA[i, i]
      n.times do |j|
        matrizA[i, j] *= p
        inv[i, j] *= p
      end
    end
  end

  return inv.to_a, independentes
end

def Cria_submatriz(matrizA, vetorX)
  submatriz = []

  matrizA.each do |linha|
    sublinha = []
    vetorX.each do |i|
      sublinha << linha[i]
    end
    submatriz << sublinha
  end

  return submatriz
end

def Separacao_da_matriz(funcaoZ, funcoes)
  inequacoes = []
  tam = funcoes.size

  funcaoZ_copy = funcaoZ.dup

  for i in 0...tam
    condicao = funcoes[i][-2]
    if condicao != '='
      funcaoZ_copy.push(0.0)
      if condicao == '<='
        inequacoes.push(1.0)
      else
        inequacoes.push(-1.0)
      end
    else
      inequacoes.push(0.0)
    end
  end

  matrizA = []
  for i in 0...tam
    linha = funcoes[i][0...-2]
    tam2 = inequacoes.size
    for j in 0...tam2
      if i == j
        linha.push(inequacoes[j])
      else
        linha.push(0.0)
      end
    end
    matrizA.push(linha)
  end

  basicas = []
  naoBasicas = []
  tam3 = funcaoZ_copy.size
  for i in 0...tam3
    if i < tam3 - tam
      naoBasicas.push(i)
    else
      basicas.push(i)
    end
  end
  basicas.sort!
  naoBasicas.sort!

  independentes = []
  for i in 0...tam
    independentes.push(funcoes[i][-1])
  end

  return matrizA, basicas, naoBasicas, independentes
end

def Calculo_x_relativo(bInversa, b)
  xRelativoBasico = Multiplicacao_matrizes(bInversa, Transposta(Matrix.rows([b])))
  return xRelativoBasico
end

def Custo(funcaoZ, variaveis)
  custoBasico = []
  variaveis.each { |index| custoBasico.push(funcaoZ[index]) }
  return custoBasico
end

def Calcula_lambda(custoBasico, basicaInversa)
  lambdaSimplex = Multiplicacao_matrizes(custoBasico, basicaInversa)
  return lambdaSimplex
end

def Custos_Relativos(lambdaSimplex, custoNaoBasico, matrizNaoBasica)
  naoBasicaTransposta = Transposta(Matrix.rows(matrizNaoBasica))
  for i in 0...custoNaoBasico.size
    custoNaoBasico[i] -= Multiplicacao_vetores(lambdaSimplex, naoBasicaTransposta.column(i).to_a)
  end
  return custoNaoBasico
end

def Calcula_k(custoRelativoNaoBasico)
  menor = custoRelativoNaoBasico.min
  return custoRelativoNaoBasico.index(menor)
end

def Otimalidade(custoRelativoNaoBasico, k)
  return custoRelativoNaoBasico[k] >= 0
end

def Direcao_simplex(basicaInversa, matrizA, k, naoBasicas)
  colunaK = []
  matrizA.each { |row| colunaK.push(row[naoBasicas[k]]) }
  colunaK = Transposta(colunaK)
  y = Multiplicacao_matrizes(basicaInversa, colunaK)
  return y
end

def Calcula_l(y, xRelativoBasico)
  seguro = false
  y.each do |value|
    if value > 0
      seguro = true
      break
    end
  end

  if !seguro
    return false
  end

  razoes = []
  xRelativoBasico.each_with_index do |value, i|
    if y[i] <= 0
      razoes.push(MAXINT)
    else
      razoes.push(value / y[i])
    end
  end

  passo = razoes.min
  l = razoes.index(passo)

  return l
end

def Troca_k_l(basicas, naoBasicas, k, l)
  aux = basicas[l]
  basicas[l] = naoBasicas[k]
  naoBasicas[k] = aux
  return [basicas, naoBasicas]
end

def Valor_funcao(funcaoZ, xRelativoBasico, basicas)
  resultado = 0
  xRelativoBasico.each_with_index do |value, i|
    resultado += funcaoZ[basicas[i]] * value
  end
  return resultado
end

def Leitura
  print "digite a funcao z separada por espacos (2 -4 3): "
  funcaoZ = gets.chomp
  funcaoZ = funcaoZ.split(' ').map(&:to_f)
  print 'digite min para minimizar e max para maximizar: '
  minMax = gets.chomp
  print 'qual o numero de funcoes? '
  numeroFuncoes = gets.chomp.to_i
  puts 'insira as funcoes separadas por enter'
  funcoes = []

  numeroFuncoes.times do |i|
    print "digite a funcao #{i + 1} separada por espacos (2 -4 3 <= 5): "
    novaFuncao = gets.chomp.split.map(&:to_f)
    
    if novaFuncao.any? { |valor| valor != novaFuncao[-2] }
      novaFuncao[-2] = novaFuncao[-2].to_f
    end
    
    funcoes << novaFuncao
  end

  return funcaoZ, funcoes, minMax
end
  
def main
  funcaoZ, funcoes, minMax = Leitura()
  puts
  matrizA, basicas, naoBasicas, independentes = Separacao_da_matriz(
    funcaoZ, funcoes)
  indFixo = independentes.dup
  funcaoFin = funcaoZ.dup
  tam = funcaoZ.length
  if minMax == 'max'
    tam.times { |i| funcaoZ[i] *= -1 }
  end

  it = 1
  maxit = 10
  solucaoOtima = []
  deu = true

  while it < maxit
    puts
    independentes = indFixo.dup
    puts "it: #{it}"
    matrizBasica = Cria_submatriz(matrizA, basicas)
    matrizNaoBasica = Cria_submatriz(matrizA, naoBasicas)
    puts "basica: ",matrizBasica
    matrizBasicaInversa, independentes = Inversa(
      matrizBasica, independentes)

    if matrizBasicaInversa == false
      deu = false
      break
    end

    xRelativo = Calculo_x_relativo(matrizBasicaInversa, independentes)

    custoBasico = Custo(funcaoZ, basicas)
    lambdaTransposto = Calcula_lambda(
      custoBasico, matrizBasicaInversa)

    custoNaoBasico = Custo(funcaoZ, naoBasicas)
    custoRelativoNaoBasico = Custos_Relativos(
      lambdaTransposto, custoNaoBasico, matrizNaoBasica)

    k = Calcula_k(custoRelativoNaoBasico)

    if Otimalidade(custoRelativoNaoBasico, k)
      puts "Otimo!"
      solucaoOtima = xRelativo
      deu = true
      break
    end
    puts "Nao otimo!"

    y = Direcao_simplex(matrizBasicaInversa, matrizA, k, naoBasicas)

    l = Calcula_l(y, xRelativo)
    if l.is_a?(FalseClass) && l == false
      deu = false
      break
    end

    basicas, naoBasicas = Troca_k_l(basicas, naoBasicas, k, l)
    it += 1
  end
  # fim do laco de repeticao simplex

  if deu
    puts "A solucao factivel otima eh:"
    tam = solucaoOtima.length
    tam.times { |i| print "x#{basicas[i]} = #{solucaoOtima[i]}, " }
    puts "z = #{Valor_funcao(funcaoFin, solucaoOtima, basicas)}"
  else
    puts 'Em algum momento nao deu para fazer a inversa porque o determinante deu 0'
    puts 'ou a direcao simplex deu <= 0'
  end
end

main()