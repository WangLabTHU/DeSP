import logging
import random
import sys
from math import sqrt
from reedsolo import RSCodec
import numpy as np

#-----------------segmentation----------------#
def file_to_indexed_dnas(file_name, chunk_size, index_length = None):
    data = preprocess(file_name, chunk_size)
    index_l = index_length
    if index_l == None:
        index_l = index_len(len(data))
    return data_to_dnas(data,index_l), index_l

def preprocess(file_name, chunk_size, is_text = False):
    data = data_from_file(file_name,is_text)
    return segments(data,chunk_size,is_text)

def data_from_file(file_name, is_text = False):
    try:
      if is_text == False:
        f = open(file_name, 'rb')
      else:
        f = open(file_name,'r')
    except: 
      logging.error("%s file not found", file_name)
      sys.exit(0)
    data = f.read()
    f.close()    
    return data

def segments(data, chunk_size,is_text = False):

    pad = -len(data) % chunk_size
    if pad > 0:
      logging.debug("Padded the file with %d zero to have a round number of blocks of data", pad)    
    if is_text == False:
        data += b'\0' * pad #zero padding.
    else:
        data += ' ' * pad 
    size = len(data)

    chunk_num = int(size/chunk_size)
    data_array = [None]*chunk_num
    for num in range(chunk_num):
        start = chunk_size * num
        end = chunk_size * (num+1)
        chunk_binary = data[start:end]
        data_array[num] = chunk_binary

    return data_array, pad

def lines_from_file(file_name):
    lines = []
    with open(file_name,'r') as f:
        while True:
            l = f.readline().split('\n')[0]
            if(l == ''):
                 break
            lines.append(l)
        f.close()
        return lines


def parse_int(f):
    line = f.readline()
    return int(line.split('\n')[0].split(' ')[1])

def load_dna(file_name):
    with open(file_name) as f:
        dnas = f.readlines()
    in_dnas = [dna.split('\n')[0] for dna in dnas]
    return in_dnas


#-----------------scanner-------------------#S
class Scanner:
    def __init__(self,max_repeat = 3, gc_interval = [0.45,0.55]):
        self.max_repeat = max_repeat
        self.gc_interval = gc_interval
        
    def scan_repeats(self,dna,record_position = False):
        repeats = []
        prv = dna[0]
        r_num = 1
        for i,c in enumerate(dna[1:]):
            if prv == c:
                r_num += 1
            else:
                if(r_num > self.max_repeat):
                    if(record_position):
                        repeats.append([prv,r_num,i-r_num+1])
                    else:
                        repeats.append([prv,r_num])
                r_num = 1
                prv = c

        if(r_num > self.max_repeat): 
                    if(record_position):
                        repeats.append([prv,r_num,i-r_num+1])
                    else:
                        repeats.append([prv,r_num])
        return repeats
    
    def max_repeats(self,dna):
        rs = self.scan_repeats(dna)
        if rs == []: return 0
        else: return max([r[1] for r in rs])
        
    def repeats_point(self,dna):
        rs = self.scan_repeats(dna)
        if rs == []: return 0
        else:
            return sum([r[1] / (self.max_repeat + 1) for r in rs])

    def Gc(self,dna):
        gc = dna.count('G') + dna.count('C')
        l = len(dna)
        return float(gc) / l
    
    def gc_pass(self,dna):
        if self.gc_interval[0]  < self.Gc(dna) < self.gc_interval[1]:
            return True
        else:
            return False
    
    def Pass(self,dna,with_primer = False):
        if self.gc_pass(dna) and self.repeats_point(dna)  == 0:
            return True
        else:
            return False
    
    def ave_gc(self,dnas):
        return sum([self.Gc(dna) for dna in dnas])/len(dnas)

    def rp_total(self,dnas):
        return(sum([self.repeats_point(dna) for dna in dnas]))

    def select_best(self,dnas):
        min_rp = 10000
        best_dna = dnas[0]
        for dna in dnas:
            if(self.gc_pass(dna)):
                if self.repeats_point(dna) < min_rp:
                    min_rp = self.repeats_point(dna)
                    best_dna = dna
        return best_dna,min_rp

    def analyze(self,dnas):
        gcs = [self.Gc(dna) for dna in dnas]
        gc_out_num = sum([not self.gc_pass(dna) for dna in dnas])
        ave_Gc = self.ave_gc(dnas)
        
        rps = [self.repeats_point(dna) for dna in dnas]
        rp_too_long = sum([rp > 0 for rp in rps])
        dic = {
            'gc_list': gcs,
            'gc_out': gc_out_num,
            'homo_list': rps,
            'average_gc': ave_Gc,
            'homo_too_long': rp_too_long
        }
        return dic


#------------------xor---------------------#

xor_map = {
    'A':{
        'A': 'A',
        'C': 'C',
        'G': 'G',
        'T': 'T'
    },
    'C':{
        'A': 'C',
        'C': 'A',
        'G': 'T',
        'T': 'G'
    },
    'G':{
        'A': 'G',
        'C': 'T',
        'G': 'A',
        'T': 'C'
    },
    'T':{
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }
}

def xor_dna(d1, d2):
    if(len(d1) != len(d1)):
        logging.error("length not equal")
    dr = ''.join([xor_map[c1][c2] for c1,c2 in zip(d1,d2)])
    return dr

def xor_ord(ord_array1, ord_array2):
       return [ord1 ^ ord2 for (ord1,ord2) in zip(ord_array1,ord_array2)]

def xor(byte_array1, byte_array2):
    return bytes([b1 ^ b2 for (b1,b2) in zip(byte_array1,byte_array2)])

#------------------random------------------#

def happen(prob):
    i = random.uniform(0,1)
    if i<=prob:
        return 1
    else:
        return 0
    
def random_base():
    r = random.random()
    if(r < 0.25): 
        return 'A'
    elif(r < 0.5):
        return 'C'
    elif(r < 0.75):
        return 'G'
    else:
        return 'T'
    
def random_dna(num):
    return ''.join([random_base() for i in range(num)])

#----------------transfromation functions---------------------#
BASE = ['A','C','G','T']
QUANT = {'A': 0, 'C':1, 'G':2, 'T':3}

def dna_to_int_array(dna_str):
    #convert a string like ACTCA to an array of ints like [10, 2, 4]
    s = ''.join('{0:02b}'.format(QUANT[dna_str[t]]) for t in range(0, len(dna_str),1))
    return [int(s[t:t+8],2) for t in range(0,len(s), 8)]


#transform a number to qua_len bases
#example: num_to_dna(6,3) returns 'ACC'
def num_to_dna(num, qua_len):
    arr = []
    while True:
        lef = num % 4
        arr.append(lef)
        num = int(num / 4)
        if 0 == num:
            break

    outString = ''
    for n in arr:
        outString = BASE[n] + outString
    
    dl = qua_len - len(arr)
    if(dl < 0):
        logging.error('space not enough for encoding num')
        return -1
    elif(dl == 0):
        return outString
    else:
        return 'A' * dl + outString
    return outString

# dna <-> bytes
def byte_to_dna(s):
    #convert byte data (\x01 \x02) to DNA data: ACTC
    bin_data = ''.join('{0:08b}'.format(s[t]) for t in range(0,len(s)))
    return bin_to_dna(bin_data)

def dna_to_byte(dna):
    #convert a string like ACTCA to a string of bytes like \x01 \x02
    num = [QUANT[b] for b in dna]
    s = ''.join('{0:02b}'.format(num[t]) for t in range(0, len(num),1))
    data = b''.join(bytes([int(s[t:t+8],2)]) for t in range(0, len(s), 8))
    return data

def bin_to_dna(bin_str):
    s = ''.join(BASE[int(bin_str[t:t+2],2)] for t in range(0, len(bin_str),2)) 
    return s

def dna_to_num(dna):
    return sum([num * 4**i for i,num in enumerate([QUANT[b] for b in dna][::-1])])


#-----------------------indexing dna chunks------------------------#
def data_to_dnas(data,index_length = 8):
    dnas = []
    for i,d in enumerate(data):
        d = num_to_dna(i,index_length) + byte_to_dna(d)
        dnas.append(d)
    return dnas

def dnas_to_data(dnas,chunk_num,index_length = 8):
    data_chunks = [b'' for b in range(chunk_num)]
    for dna in dnas:
        index = dna_to_num(dna[:index_length])
        payload = dna_to_byte(dna[index_length:])
        data_chunks[index] = payload
    return b''.join(data_chunks)

def index_len(chunk_num):
    return int(sqrt(sqrt(chunk_num))) + 1

#------------------------RS------------------------------------------#
def rs_decode(data, rs_obj = None, rs = None, max_hamming = 1):
    if not rs_obj:
        if rs: rs_obj = RSCodec(rs)
        else: 
            print('rs decoder not assigned.')
            return
    try:
        data_corrected = list(rs_obj.decode(data)[0])
    except:
        logging.debug('can not correct ori data')
        return -1, None 
    #we will encode the data again to evaluate the correctness of the decoding
    data_again = list(rs_obj.encode(data_corrected))
    if np.count_nonzero(data != list(data_again)) > max_hamming: #measuring hamming distance between raw input and expected raw input
        #too many errors to correct in decoding
        logging.debug('too many errors!')
        return -1, None
    return 0, data_corrected