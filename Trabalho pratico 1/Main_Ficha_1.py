from Fuctions_Ficha_1 import getInfo, getOcorrencia, histograma, entropia, mediaHuffman
from Fuctions_Ficha_1 import getEntropiaAgrupado, varianciaHuffman, percorreTarget, graficoInformacaoMutua
import numpy as np

imagem1 = "MRI.bmp"
imagem2 = "MRIbin.bmp"
imagem3 = "landscape.bmp"
som1 = "soundMono.wav"
texto1 = "lyrics.txt"

filename = texto1

#EXERCICIO_1
[fonte, alfabeto] = getInfo(filename)
arrayOcorrencias = getOcorrencia(fonte, alfabeto)

#EXERCICIO_2
entropy = entropia(arrayOcorrencias, fonte)
print("Entropia ===>", entropy, " bits/simbolo")

#EXERCICIO_3
histograma(arrayOcorrencias, alfabeto)

#EXERCICIO_4
[mediaComprimento, lengths] = mediaHuffman(fonte, arrayOcorrencias)
print("Media de comprimento ===>", mediaComprimento)
variancia = varianciaHuffman(mediaComprimento, fonte, arrayOcorrencias, lengths)
print("Variancia ==> ", variancia)

#EXERCICIO_5
entropiaAgrupado = getEntropiaAgrupado(fonte, alfabeto, filename)
print("Entropia de simbolos agrupados ===>", entropiaAgrupado, " bits/simbolo")

#Exercicio_6A
querySimulacao = np.asarray([2, 6, 4, 10, 5, 9, 5, 8, 0, 8])
targetSimulacao = np.asarray([6, 8, 9, 7, 2, 4, 9, 9, 4, 9, 1, 4, 8, 0, 1, 2, 2, 6, 3, 2, 0, 7, 4, 9, 5, 4, 8, 5, 2, 7, 
          8, 0, 7, 4, 8, 5, 7, 4, 3, 2, 2, 7, 3, 5, 2, 7, 4, 9, 9, 6])
alfabetoSimulacao = np.arange(0, 11)
passoSimulacao = 1
[arraySimulacao, eixoXX] = percorreTarget(passoSimulacao, targetSimulacao, querySimulacao, alfabetoSimulacao)
print("Informacao mutua simulacao ==> ", arraySimulacao)

#Exercicio_6B
target01 = "target01 - repeat.wav"
target02 = "target02 - repeatNoise.wav"
queryName = "saxriff.wav"
[fonte, alfabetoQuery] = getInfo(queryName)
[fonteTarget, useless] = getInfo(target01)
query = fonte
target = fonteTarget
passo = int(1/4 * len(query))
#TARGET01
[arrayInfoMutua, eixoXX] = percorreTarget(passo, target, query, alfabetoQuery)
graficoInformacaoMutua(arrayInfoMutua, eixoXX)
print("saxriff.wav ==> target01 - repeat.wav:\n ", arrayInfoMutua)
#TARGET02
[fonteTarget2, useless] = getInfo(target02)
target = fonteTarget2
[arrayInfoMutua2, eixoXX] = percorreTarget(passo, target, query, alfabetoQuery)
graficoInformacaoMutua(arrayInfoMutua2, eixoXX)
print("saxriff.wav ==> target02 – repeatNoise.wav:\n ", arrayInfoMutua2)

#EXERCICIO_6C
nomes = ["Song01.wav", "Song02.wav", "Song03.wav", "Song04.wav", "Song05.wav", "Song06.wav", "Song07.wav"]
informacaoMaxima = list()
[query, alfabeto] = getInfo("saxriff.wav") 
passo = int (0.25 * len(query))
for i in range(len(nomes)):
    [target, useless] = getInfo(nomes[i])
    [arrayInfoMutua, eixoXX] = percorreTarget(passo, target, query, alfabeto)
    informacaoMaxima.append(max(arrayInfoMutua))
#associamos cada informacao maxima a cada nome
associacaoInfo = dict(zip(informacaoMaxima, nomes))
#ordena-se os valores
informacaoMaxima = np.sort(informacaoMaxima)
for i in range(len(informacaoMaxima) - 1, -1, -1):
    print("Nome ==> ",associacaoInfo[informacaoMaxima[i]], " Informacao Mutua Maxima ==> ", informacaoMaxima[i])




