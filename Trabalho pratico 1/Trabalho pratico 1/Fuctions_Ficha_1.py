import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
import matplotlib.image as mpimg
from huffmancodec import HuffmanCodec


def getInfo(fileName):
    if ".txt" in fileName:
        # Guarda numa lista o conteudo do ficheiro de texto
        with open(fileName, 'r') as f:
            text = list(f.read())
            charText = list()
        # Iteramos sobre text para filtrar apenas simbolos regulares do alfabeto
        # e guardamos em charText
        for i in range(len(text)):
            toAscii = ord(text[i])
            if ((toAscii > 64 and toAscii < 91) or (toAscii > 96 and toAscii < 123)):
                charText.append(toAscii)
        fonte = np.asarray(charText)
        # Criamos um alfabeto de ord('A') -> ord('Z') e de ord('a') -> ord('z'
        alfabeto = list()
        for i in range(ord('A'), ord('Z') + 1):
            alfabeto.append(i)
        for i in range(ord('a'), ord('z') + 1):
            alfabeto.append(i)
        alfabeto = np.asarray(alfabeto)
        return [fonte, alfabeto]

    elif ".wav" in fileName:
        [fs, fonte] = wavfile.read(fileName)
        # No caso de a fonte ter mais que 1 canal, utilizamos apenas 1 deles
        if (len(fonte.shape) > 1):
            fonte = fonte[:, ]
        fonte = fonte.flatten()
        # Criamos um alfabeto de 0 -> 255
        alfabeto = np.arange(0, 256)
        return [fonte, alfabeto]

    elif ".bmp" in fileName:
        # Criamos um alfabeto de 0 -> 255
        alfabeto = np.arange(0, 256)
        fonte = mpimg.imread(fileName)
        fonte = fonte.flatten()
        return [fonte, alfabeto]


def getOcorrencia(fonte, alfabeto):
    listaOcorrencias = list()
    # np.where devolve um tuplo com as posicoes onde alfabeto[i] ocorre na fonte
    # Length desse tuplo sera o numero de ocorrencias
    for i in range(len(alfabeto)):
        aux = np.where(alfabeto[i] == fonte)
        numeroOcorrencias = len(aux[0])
        listaOcorrencias.append(numeroOcorrencias)
    finalArray = np.asarray(listaOcorrencias)
    return finalArray


def histograma(ocorrencias, alfabeto):
    plt.xlabel("alfabeto")
    plt.ylabel("ocorrencias")
    plt.bar(alfabeto, ocorrencias)
    plt.show()


def entropia(ocorrencias, fonte):
    lenFonte = np.sum(ocorrencias)
    arrayProbabilidade = ocorrencias / lenFonte
    # Elimina todas as posicoes onde probabilidade = 0, pois log2(0) nao existe
    arrayProbabilidade = arrayProbabilidade[arrayProbabilidade > 0]
    entropy = np.sum(arrayProbabilidade * np.log2(1/arrayProbabilidade))
    return entropy


def mediaHuffman(fonte, ocorrencias):
    codec = HuffmanCodec.from_data(fonte)
    symbols, lengths = codec.get_code_len()
    # Array com os comprimentos respetivos de cada um dos simbolos da fonte
    lengths = np.asarray(lengths)
    # Cria um novo array com todas as ocorrencias, exceto com valor 0
    arrayAux = ocorrencias[ocorrencias > 0]
    media = np.sum((arrayAux/len(fonte)) * lengths)
    return media, lengths


def varianciaHuffman(mediaComprimento, fonte, arrayOcorrencias, lengths):
    variancia = 0
    # Cria um novo array com todas as ocorrencias, exceto com valor 0
    auxList = arrayOcorrencias[arrayOcorrencias > 0]
    variancia = np.sum((auxList/len(fonte)) * (lengths - mediaComprimento)**2)
    return variancia


def getEntropiaAgrupado(fonte, alfabeto, filename):
    if ".txt" in filename:
        auxList = np.zeros((122, 122), np.float32)
    else:
        auxList = np.zeros((len(alfabeto), len(alfabeto)),np.float32)
    lenFonte = len(fonte)
    if (lenFonte % 2 == 1):
        lenFonte = lenFonte - 1
    for i in range(0, lenFonte, 2):
        auxList[fonte[i], fonte[i+1]] += 1
    auxList = auxList/len(fonte)
    auxList = auxList[auxList > 0]
    entropia = np.sum(auxList * np.log2(1/auxList))
    return entropia


def targetQueryAgrupado(query, target, inicio):
    auxList = list()    
    for i in range(len(query)):
            auxList.append((query[i], target[i + inicio]))
    agrupado = auxList
    return agrupado

def ocorrenciasNormal(fonte, alfabeto):
    arrayOcorrencias = np.zeros(len(alfabeto))  
    for i in range (len(fonte)):
        arrayOcorrencias[fonte[i]] += 1
    return arrayOcorrencias 


def entropyConjunta(query, alfabeto, target):
    probabilidadeConjunta = np.zeros((len(alfabeto), len(alfabeto)),np.float32)
    #percorre a query e o target ao mesmo tempo e na matriz probabilidadeConjunta incrementa-se 1 na posicao query[i]target[i]
    for i in range(len(query)):
        probabilidadeConjunta[query[i]][target[i]] += 1
    probabilidadeConjunta = probabilidadeConjunta / len(query)
    probabilidadeConjunta = probabilidadeConjunta[probabilidadeConjunta>0]
    entropiaConjunta = - np.sum(probabilidadeConjunta * np.log2(probabilidadeConjunta))
    return entropiaConjunta


def getInformacaoMutua(query, target, alfabeto, agrupado):
    ocorrenciasQuery = ocorrenciasNormal(query, alfabeto)
    ocorrenciasTarget = ocorrenciasNormal(target, alfabeto)
    #calcula-se a entropia da query e do target
    entropiaQuery = entropia(ocorrenciasQuery, alfabeto)
    entropiaTarget = entropia(ocorrenciasTarget, alfabeto)
    #calcula-se a entropia conjugada da query e do target
    entropiaConjunta = entropyConjunta(query, alfabeto, target)
    #calcula-se a informacao mutua
    informacaoMutua = entropiaQuery + entropiaTarget - entropiaConjunta
    return informacaoMutua


def percorreTarget(passo, target, query, alfabeto):
    informacaoMutua = list()
    inicio = 0
    fim = len(query)
    tempo = 0
    for i in range(len(target / passo) + 1):
        tempo += 1
        #cria-se um array com os elementos da query e do target agrupados
        tqAgrupado = targetQueryAgrupado(query, target, inicio)
        #adiciona-se os valores da informacao mutua a um vetor
        informacaoMutua.append(getInformacaoMutua(query, target[inicio:fim], alfabeto, tqAgrupado))
        inicio += passo
        fim += passo
        if (fim > len(target)):
            break
    eixoXX = np.arange(tempo)
    informacaoMutua = np.asarray(informacaoMutua)
    return informacaoMutua, eixoXX   

def graficoInformacaoMutua(informacaoMutua, tempo):
    plt.xlabel("TEMPO")
    plt.ylabel("INFORMACAO MUTUA")
    plt.plot(tempo, informacaoMutua)
    plt.show()
    

    

