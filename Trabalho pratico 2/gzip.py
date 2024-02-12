import sys
from huffmantree import HuffmanTree
import numpy as np

class GZIPHeader:
    ''' class for reading and storing GZIP header fields '''

    ID1 = ID2 = CM = FLG = XFL = OS = 0
    MTIME = []
    lenMTIME = 4
    mTime = 0

    # bits 0, 1, 2, 3 and 4, respectively (remaining 3 bits: reserved)
    FLG_FTEXT = FLG_FHCRC = FLG_FEXTRA = FLG_FNAME = FLG_FCOMMENT = 0   
    
    # FLG_FTEXT --> ignored (usually 0)
    # if FLG_FEXTRA == 1
    XLEN, extraField = [], []
    lenXLEN = 2
    
    # if FLG_FNAME == 1
    fName = ''  # ends when a byte with value 0 is read
    
    # if FLG_FCOMMENT == 1
    fComment = ''   # ends when a byte with value 0 is read
        
    # if FLG_HCRC == 1
    HCRC = []
        
        
    
    def read(self, f):
        ''' reads and processes the Huffman header from file. Returns 0 if no error, -1 otherwise '''

        # ID 1 and 2: fixed values
        self.ID1 = f.read(1)[0]  
        if self.ID1 != 0x1f: return -1 # error in the header
            
        self.ID2 = f.read(1)[0]
        if self.ID2 != 0x8b: return -1 # error in the header
        
        # CM - Compression Method: must be the value 8 for deflate
        self.CM = f.read(1)[0]
        if self.CM != 0x08: return -1 # error in the header
                    
        # Flags
        self.FLG = f.read(1)[0]
        
        # MTIME
        self.MTIME = [0]*self.lenMTIME
        self.mTime = 0
        for i in range(self.lenMTIME):
            self.MTIME[i] = f.read(1)[0]
            self.mTime += self.MTIME[i] << (8 * i)                 
                        
        # XFL (not processed...)
        self.XFL = f.read(1)[0]
        
        # OS (not processed...)
        self.OS = f.read(1)[0]
        
        # --- Check Flags
        self.FLG_FTEXT = self.FLG & 0x01
        self.FLG_FHCRC = (self.FLG & 0x02) >> 1
        self.FLG_FEXTRA = (self.FLG & 0x04) >> 2
        self.FLG_FNAME = (self.FLG & 0x08) >> 3
        self.FLG_FCOMMENT = (self.FLG & 0x10) >> 4
                    
        # FLG_EXTRA
        if self.FLG_FEXTRA == 1:
            # read 2 bytes XLEN + XLEN bytes de extra field
            # 1st byte: LSB, 2nd: MSB
            self.XLEN = [0]*self.lenXLEN
            self.XLEN[0] = f.read(1)[0]
            self.XLEN[1] = f.read(1)[0]
            self.xlen = self.XLEN[1] << 8 + self.XLEN[0]
            
            # read extraField and ignore its values
            self.extraField = f.read(self.xlen)
        
        def read_str_until_0(f):
            s = ''
            while True:
                c = f.read(1)[0]
                if c == 0: 
                    return s
                s += chr(c)
        
        # FLG_FNAME
        if self.FLG_FNAME == 1:
            self.fName = read_str_until_0(f)
        
        # FLG_FCOMMENT
        if self.FLG_FCOMMENT == 1:
            self.fComment = read_str_until_0(f)
        
        # FLG_FHCRC (not processed...)
        if self.FLG_FHCRC == 1:
            self.HCRC = f.read(2)
            
        return 0
            



class GZIP:
    ''' class for GZIP decompressing file (if compressed with deflate) '''

    gzh = None
    gzFile = ''
    fileSize = origFileSize = -1
    numBlocks = 0
    f = None
    

    bits_buffer = 0
    available_bits = 0        

    
    def __init__(self, filename):
        self.gzFile = filename
        self.f = open(filename, 'rb')
        self.f.seek(0,2)
        self.fileSize = self.f.tell()
        self.f.seek(0)

        
    #Guardamos no array Clens os comprimentos dos codigos do "alfabeto de comprimento
    #de codigos"
    def comprimentosCodigosHCLEN(self, HCLEN):
        CLens = np.zeros(19, np.int16)
        posicao = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15]
        
        for i in range(HCLEN + 4):
            CLens[posicao[i]] = self.readBits(3)
        #print(CLens)
        return CLens
    
    
    #Guardamos no array HLIT_lens os comprimentos dos codigos referente ao
    #alfabeto de comprimentos
    def arrayHLIT(self, HLIT, tree):
        pos = 0
        HLIT_lens = np.zeros(HLIT + 257, np.int16)      
        while (pos < HLIT + 257):
            bit = str(self.readBits(1))
            value = tree.nextNode(bit)
            while(value < 0):
                bit = str(self.readBits(1))
                value = tree.nextNode(bit)
            if value < 16:
                HLIT_lens[pos] = value
                pos += 1
                
            elif value == 16:
                number = HLIT_lens[pos - 1]
                for i in range(3 + self.readBits(2)):
                    HLIT_lens[pos] = number
                    pos += 1
                    
            elif value == 17:
                for i in range(3 + self.readBits(3)):
                    HLIT_lens[pos] = 0
                    pos += 1
                    
            elif value == 18:
                for i in range(11 + self.readBits(7)):
                    HLIT_lens[pos] = 0
                    pos += 1
            tree.resetCurNode()   
        return HLIT_lens
    
    
    #Guardamos no array HDIST_lens os comprimentos dos codigos referentes ao 
    #alfabeto de comprimentos, tal como, na função arrayHLIT
    def arrayHDIST(self, HDIST, tree):
        pos = 0
        HDIST_lens = np.zeros(HDIST + 1, np.int16)      
        while (pos < HDIST + 1):
            bit = str(self.readBits(1))
            value = tree.nextNode(bit)
            while(value < 0):
                bit = str(self.readBits(1))
                value = tree.nextNode(bit)
            if value < 16:
                HDIST_lens[pos] = value
                pos += 1
                
            elif value == 16:
                number = HDIST_lens[pos - 1]
                for i in range(3 + self.readBits(2)):
                    HDIST_lens[pos] = number
                    pos += 1
                    
            elif value == 17:
                for i in range(3 + self.readBits(3)):
                    HDIST_lens[pos] = 0
                    pos += 1
                    
            elif value == 18:
                for i in range(11 + self.readBits(7)):
                    HDIST_lens[pos] = 0
                    pos += 1
            tree.resetCurNode()
        return HDIST_lens
    
    #algoritmo slide 34
    #obtemos o bl_count ==> array com a frequencia dos comprimentos dos codigos
    #na posicao 1 do bl_count temos quantas vezes ocorre o comprimento '1'
    #o comprimento maximo é 8,  no caso da HCLEN, pois foi lido em 3 bits (0->7, há 8 numeros)
    #o comprimento maximo é 16, no caso da HLIT e HDIST, pois foi lido em 4 bits(0->15, há 16 numeros)
    #obtemos os Values, ou seja, os primeiros valores de cada comprimento
    def conversaoCodigosHuffman(self, Lens, Max_Bits):
        bl_count = np.zeros(Max_Bits, np.int16)
        for i in range(len(Lens)):
            bl_count[Lens[i]] += 1
        bl_count[0] = 0
        
        code = 0
        next_code = np.zeros(Max_Bits, np.int16)
        for i in range(1, Max_Bits):
            code = (code + bl_count[i-1]) << 1
            next_code[i] = code
        
        Values = np.zeros(len(Lens), np.int16)
        for i in range(len(Lens)):
            length = Lens[i]
            if(length != 0):
                Values[i] = next_code[length]
                next_code[length] += 1;
        #print(Values)
        return Values
    
    #Com os comprimentos(cLens) e os valores(cValues) criamos os codigos binarios dos simbolos
    #caso o comprimento binario seja inferior ao respetivo valor de cLens, adicionamos
    #'0's no final até que os dois comprimentos sejam iguais
    def decimalToBinario(self, Lens, Values):
        string = list()
        
        for i in range(len(Lens)):
            binario = list()
            if Lens[i] != 0:
                value = Values[i]
                while(value > 0):
                    binario.insert(0,value%2)
                    value = value//2
                while (len(binario) < Lens[i]):
                    binario.insert(0, 0)
                string.append(''.join(map(str, binario)))
            else:
                string.append('-1')
        #print(string)
        return string
    
    #Criamos uma nova folha na arvore de Huffman com o simbolo i
    #Com o codigo de Huffman string[i]
    def fillTree(self, tree, string):
        for i in range(len(string)):
            if string[i] != '-1':
                tree.addNode(string[i], i)
                
    
    #Descompactação
    def descompactacao(self, HLIT_tree, HDIST_tree):
        output_stream = list()
        while(1):
            value = -1
            value1 = -1
            #percorre a arvore HLIT ate chegar a uma folha
            while(value < 0):
                bit = str(self.readBits(1))
                value = HLIT_tree.nextNode(bit)
            HLIT_tree.resetCurNode()
            #caso o valor seja menor que 256 add a output_stream o valor
            if value < 256:
                output_stream.append(value)
            #caso seja igual a 256 significa que chegamos ao fim
            elif value == 256:
                break
            #caso seja maior que 256, obtemos um length com o auxilio das tabelas
            #logica usada: value - (code - primeiro length de cada extra bit) + (value - primeiro code de cada length)
            # * (diferenca entre dois primeiros length consecutivos do mesmo extra bit - 1)
            #Desta forma arranjamos uma solução para encontrar o length de cada extra bit de diferentes codes
            elif value > 256 and value < 265:
                length = value - 254
            elif value >  264 and value < 269:
                length = self.readBits(1) + (value - 254 + value - 265)
            elif value > 268 and value < 273:
                length = self.readBits(2) + (value - 250 + (value - 269) * 3)
            elif value > 272 and value < 277:
                length = self.readBits(3) + (value - 238 + (value - 273) * 7)
            elif value > 276 and value < 281:
                length = self.readBits(4) + (value - 210 + (value - 277) * 15)
            elif value > 280 and value < 285:
                length = self.readBits(5) + (value - 150 + (value - 281) * 31)
            elif value == 285:
                length = 258
        
            #caso value > 256 precisamos de obter a distancia a recuar atraves da HDIST_tree
            if value > 256:
                #percorremos a arvore HDIST ate chegar a uma folha
                while(value1 < 0):
                    bit1 = str(self.readBits(1))
                    value1 = HDIST_tree.nextNode(bit1)
                HDIST_tree.resetCurNode()
                #Obtemos o valor da distancia a recuar atraves das tabelas
                #Logica usada: value1 + (diferenca entre o primeiro length com o primeiro code de cada
                #extra bit) + (value1 - primeiro code de cada extra bit) * (2**extrabits - 1)
                #Desta forma arranjamos uma solução para encontrar a dist de cada extra bit de diferentes codes
                if value1 >= 0 and value1 < 4:
                    dist = value1 + 1
                elif value1 > 3 and value1 < 6:
                    dist = self.readBits(1) + value1 + 1 + value1 - 4
                elif value1 > 5  and value1 < 8:
                    dist = self.readBits(2) + value1 + 3 + (value1 - 6) * 3
                elif value1 > 7 and value1 < 10:
                    dist = self.readBits(3) + value1 + 9 + (value1 - 8) * 7
                elif value1 > 9 and value1 < 12:
                    dist = self.readBits(4) + value1 + 23 + (value1 - 10) * 15
                elif value1 > 11 and value1 < 14:
                    dist = self.readBits(5) + value1 + 53 + (value1 - 12) * 31
                elif value1 > 13 and value1 < 16:
                    dist = self.readBits(6) + value1 + 115 + (value1 - 14) * 63
                elif value1 > 15 and value1 < 18:
                    dist = self.readBits(7) + value1 + 241 + (value1 - 16) * 127
                elif value1 > 17 and value1 < 20:
                    dist = self.readBits(8) + value1 + 495 + (value1 - 18) * 255
                elif value1 > 19 and value1 < 22:
                    dist = self.readBits(9) + value1 + 1005 + (value1 - 20) * 511
                elif value1 > 21 and value1 < 24:
                    dist = self.readBits(10) + value1 + 2027 + (value1 - 22) * 1023
                elif value1 > 23 and value1 < 26:
                    dist = self.readBits(11) + value1 + 4073 + (value1 - 24) * 2047
                elif value1 > 25 and value1 < 28:
                    dist = self.readBits(12) + value1 + 8167 + (value1 - 26) * 4095
                elif value1 > 27 and value1 < 30:
                    dist = self.readBits(13) + value1 + 16357 + (value1 - 28) * 8191
                
                #Ao chegarmos aqui já temos a distancia a recuar e o numero de vezes a copiar os numeros
                # da stream (por exemplo: dist = 2 && length = 3 ==> recuamos 2 posições no array e copiamos 
                #os 3 valores seguintes)
                pos = len(output_stream) - dist
                for i in range(length):
                    output_stream.append(output_stream[pos])
                    pos += 1
        return output_stream
    
    def gravarDados(self, fname ,output_stream):
        with open(fname,'w') as fich:
            for dados in output_stream:
                fich.write('%s' %chr(dados))
            fich.close()
        
    def decompress(self):
        ''' main function for decompressing the gzip file with deflate algorithm '''
        
        numBlocks = 0

        # get original file size: size of file before compression
        origFileSize = self.getOrigFileSize()
        print(origFileSize)
        
        # read GZIP header
        error = self.getHeader()
        if error != 0:
            print('Formato invalido!')
            return
        
        # show filename read from GZIP header
        print(self.gzh.fName)
        
        
        # MAIN LOOP - decode block by block
        BFINAL = 0    
        while not BFINAL == 1:    
            
            BFINAL = self.readBits(1)
            
            BTYPE = self.readBits(2)                    
            if BTYPE != 2:
                print('Error: Block %d not coded with Huffman Dynamic coding' % (numBlocks+1))
                return
            
            HLIT = self.readBits(5)
            HDIST = self.readBits(5)
            HCLEN = self.readBits(4)
            
            #Obtenção da arvore HCLEN
            tree_HCLEN = HuffmanTree()
            HCLEN_lens = self.comprimentosCodigosHCLEN(HCLEN)
            HCLEN_values = self.conversaoCodigosHuffman(HCLEN_lens, 8)
            strings_HCLEN = self.decimalToBinario(HCLEN_lens, HCLEN_values)
            self.fillTree(tree_HCLEN, strings_HCLEN)
            
            #Obtenção da arvore HCLIT
            tree_HCLIT = HuffmanTree()
            HLIT_lens = self.arrayHLIT(HLIT, tree_HCLEN)
            HLIT_values = self.conversaoCodigosHuffman(HLIT_lens, 16)
            strings_HLIT = self.decimalToBinario(HLIT_lens, HLIT_values)
            self.fillTree(tree_HCLIT, strings_HLIT)
            
            #Obtenção da arvore HDIST
            tree_HDIST = HuffmanTree()
            HDIST_lens = self.arrayHDIST(HDIST, tree_HCLEN)
            HDIST_values = self.conversaoCodigosHuffman(HDIST_lens, 16)
            strings_HDIST = self.decimalToBinario(HDIST_lens, HDIST_values)
            self.fillTree(tree_HDIST, strings_HDIST)
            
            output_stream=self.descompactacao(tree_HCLIT,tree_HDIST)
        
            self.gravarDados(self.gzh.fName,output_stream)
                                                                                                                                             
            # update number of blocks read
            numBlocks += 1

        
        # close file            
        
        self.f.close()    
        print("End: %d block(s) analyzed." % numBlocks)
    
    
    def getOrigFileSize(self):
        ''' reads file size of original file (before compression) - ISIZE '''
        
        # saves current position of file pointer
        fp = self.f.tell()
        
        # jumps to end-4 position
        self.f.seek(self.fileSize-4)
        
        # reads the last 4 bytes (LITTLE ENDIAN)
        sz = 0
        for i in range(4): 
            sz += self.f.read(1)[0] << (8*i)
        
        # restores file pointer to its original position
        self.f.seek(fp)
        
        return sz        
    

    
    def getHeader(self):  
        ''' reads GZIP header'''

        self.gzh = GZIPHeader()
        header_error = self.gzh.read(self.f)
        return header_error
        

    def readBits(self, n, keep=False):
        ''' reads n bits from bits_buffer. if keep = True, leaves bits in the buffer for future accesses '''

        while n > self.available_bits:
            self.bits_buffer = self.f.read(1)[0] << self.available_bits | self.bits_buffer
            self.available_bits += 8
        
        mask = (2**n)-1
        value = self.bits_buffer & mask

        if not keep:
            self.bits_buffer >>= n
            self.available_bits -= n

        return value

    


if __name__ == '__main__':

    # gets filename from command line if provided
    fileName = "FAQ.txt.gz"
    if len(sys.argv) > 1:
        fileName = sys.argv[1]            

    # decompress file
    gz = GZIP(fileName)
    gz.decompress()
