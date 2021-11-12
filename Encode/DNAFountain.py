import struct
import random
import os
import logging
import sys
import operator
import math
import numpy as np
from collections import defaultdict

from reedsolo import RSCodec
from Encode.Helper_Functions import *
from Encode.RPNG import * 


#----------------------------------------------------Droplet-------------------------------------------------#  
class Droplet:
    def __init__(self, data, seed, num_chunks = None, rs = 0, rs_obj = None, degree = None):
        #num_chunks is a list of the orignal packets numbers used to xor 
        #rs is the number of Reed Solomon symbols to add to the message

        self.data = data
        self.seed = seed
        self.num_chunks = set(num_chunks)
        self.rs = rs
        self.rs_obj = rs_obj
        self.degree = degree

        self.DNA = None
    
    def toDNA(self, flag = None):
        #this function wraps the seed, data payload, and Reed Solomon.
        if self.DNA is not None:
            return self.DNA
        self.DNA = byte_to_dna(self._package())
        return self.DNA
    
    def chunkStr(self):
        num = 0
        s = ''
        for i in self.num_chunks:
            if(6 == num):
                s += '...'
                break
            s+= str(i) + ' '
        num += 1
        return s
        
    def _package(self):
        #this function converts the seed to a list of 4bytes HARD CODED!!!
        #adds the seed to the data (list of integers)
        #computes a reed solomon on the seed+data.
        #returns everything.

        seed_ord = self.seed.to_bytes(4, byteorder = 'big')
        # print(type(self.data[0]))
        message = seed_ord + bytes(self.data)
        # message = seed_ord + bytearray(self.data)

        
        if self.rs > 0:
            message = self.rs_obj.encode(message) #adding RS symbols to the message

        return message
    
#----------------------------------------------------Fountain-------------------------------------------------#      
class DNAFountain:

    def __init__(self, 
                file_in,
                alpha, 
                stop = None,
                rs = 0, 
                c_dist = 0.1, 
                delta = 0.5, 
                scanner = None
                ):

        #alpha is the redundency level
        #stop is whether we have a limit on the number of oligos
        #chunk_size and file_size are in bytes
        #rs is the number of bytes for reed-solomon error correcting code over gf(2^8).
        #c_dist is a parameter of the degree distribution
        #delta is a parameter of the degree distribution
        #np: should we use numpy random number generator? Faster, but incompatible for previous versions
        #max_homopolymer: the largest homopolymer allowed
        #gc: the allowable range of gc +- 50%

        #data:
        self.file_in = file_in
        self.chunk_size = len(file_in[0])
        self.num_chunks = len(file_in)
        
        #reduancy:
        self.alpha = alpha
        self.stop = stop
        self.final = self.calc_stop()

        #random mnumber generator
        self.lfsr = lfsr(lfsr32s(), lfsr32p()) #starting an lfsr with a certain state and a polynomial for 32bits.
        self.lfsr_l = len('{0:b}'.format( lfsr32p())) - 1 #calculate the length of lsfr in bits 
        self.seed = self.lfsr.__next__()

        self.PRNG = PRNG(K = self.num_chunks, delta = delta, c = c_dist, np = False) #creating the solition distribution object
        self.PRNG.set_seed(self.seed)

        #error correcting code:
        self.rs = rs #the number of symbols (bytes) to add
        self.rs_obj = RSCodec(self.rs) #initalizing an reed solomon object

        #biological screens:
        self.scanner = scanner
        if self.scanner == None:
            self.scanner = Scanner()

        self.tries = 0 #number of times we tried to create a droplet
        self.good = 0 #droplets that were screened successfully.
        
        self.oligo_l = self.calc_oligo_length()
        # store the generated droplets
        self.dna_df = None
        self.dna_dl = []
    
    
    def calc_oligo_length(self):
        #return the number of nucleotides in an oligo:
        bits = self.chunk_size * 8 + self.lfsr_l + self.rs * 8
        return bits/4


    def calc_stop(self):
        if self.stop is not None:
            return self.stop
        stop = int(self.num_chunks*(1+self.alpha))+1
        return stop

    def droplet(self):
        #creating a droplet.
        data = None

        d, num_chunks = self.rand_chunk_nums() #creating a random list of segments.

        # print(num_chunks)
        for num in num_chunks: #iterating over each segment
            if data is None: #first round. data payload is empty.
                data = self.chunk(num) #just copy the segment to the payload.
            else: #more rounds. Xor the new segments with the payload.
                data = xor(data,self.chunk(num))  #map(operator.xor, data, self.chunk(num))
        
        # print(data)
        self.tries +=  1 

        #we have a droplet:
        return Droplet(data = data, 
                       seed = self.seed, 
                       rs = self.rs,
                       rs_obj = self.rs_obj,
                       num_chunks = num_chunks,
                       degree = d)

    def chunk(self, num):
        #return the num-th segment from the file
        return self.file_in[num]

    #-------------------generate random chunk numebers----------------# 
    def updateSeed(self):
        #This function creates a fresh seed for the droplet and primes the solition inverse cdf sampler
        self.seed = self.lfsr.__next__() #deploy one round of lfsr, and read the register.
        self.PRNG.set_seed(self.seed) #update the seed with the register

    def rand_chunk_nums(self):
        #This funcation returns a subset of segments based on the solition distribution.
        #It updates the lfsr to generates a new seed.
        self.updateSeed() #get a fresh seed and prime the solition inverse cdf sampler.
        blockseed, d, ix_samples = self.PRNG.get_src_blocks_wrap()
        return d, ix_samples #return a list of segments.

    #----------------screen generated droplets----------------------#
    def screen(self, droplet):
        if self.scanner.Pass(droplet.toDNA()):
            self.good += 1
            dna = droplet.toDNA()
            degree = droplet.degree
            chunk_str = droplet.chunkStr()
            seed = droplet.seed
            self.dna_dl.append([dna,seed,degree,chunk_str])
            return 1
        return 0

    def save(self,file_name = 'out.dna'):
        with open(file_name, 'w') as f:
            # f.write('Fountain code\n')
            # f.write('CN: ' + str(self.num_chunks) +'\n')
            # f.write('CL: ' + str(self.chunk_size) + '\n')
            # f.write('RS: ' + str(self.rs) + '\n')
            f.writelines('\n'.join([d[0] for d in self.dna_dl]))
            f.close()

    def encode(self):
        self.dl = []
        self.tries = 0
        self.good = 0
        while self.good < self.final:
            self.screen(self.droplet())
            if self.tries%2000 == 0:
                logging.info("generate %d chunks after %d tries",self.good, self.tries)
                # print("generate %d chunks after %d tries" % (self.good, self.tries))
                
        logging.info("Finish generating %d chunks after %d tries", self.good,self.tries)
        # print("Finish generating %d chunks after %d tries"% (self.good,self.tries))
        return self.good, self.tries    

#----------------------------------------------------Glass-------------------------------------------------#        
class Glass:
    def __init__(self, in_file_name, chunk_num, header_size = 4, 
                 rs = 0, c_dist = 0.1, delta = 0.5, 
                flag_correct = True, gc = 0.05, max_homopolymer = 3, 
                max_hamming = 100, chunk_size = 32, exDNA = False, np = False, truth = None):
        
        self.entries = []
        self.droplets = set()
        self.num_chunks = chunk_num
        self.chunks = [None] * self.num_chunks
        self.header_size = header_size
        self.chunk_size = chunk_size
        self.exDNA = exDNA
        self.np = np
        self.chunk_to_droplets = defaultdict(set)
        self.done_segments = set()
        self.truth = truth
        self.in_file_name = in_file_name
        self.max_hamming = max_hamming

        self.PRNG = PRNG(K = self.num_chunks, delta = delta, c = c_dist, np = np)

        self.rs = rs
        self.RSCodec = None
        self.correct = flag_correct
        self.seen_seeds = set()
        
        if self.rs > 0:
            self.RSCodec = RSCodec(rs)
    
    def add_dna(self, dna_string):
        # transfer data to int
        data = dna_to_int_array(dna_string)
         
        # try error correcting, if rs code is added and we want to correct error
        if self.rs > 0:
            if self.correct: 
                flag, data_corrected = rs_decode(data, self.RSCodec)
                if flag == -1:
                    return -1, None
            else: #if we don't want to evaluate the error correcting code, just delete the rs code
                data_corrected  = data[0:len(data) - self.rs] 
        else:
            data_corrected = data
        
        # split seed and payload
        seed_array = data_corrected[:self.header_size]
        seed = sum([   int(x)*256**i        for i, x in enumerate(seed_array[::-1])   ])
        payload = data_corrected[self.header_size:]
        self.add_seed(seed)
        
        # decode
        self.PRNG.set_seed(seed)
        blockseed, d, ix_samples = self.PRNG.get_src_blocks_wrap() #reconstruct the linear combination
        # print(blockseed,d,ix_samples,payload)
        d = Droplet(payload, seed, ix_samples) #reconstruct droplet
        self.addDroplet(d) 
        return seed, data

    def addDroplet(self, droplet):
        self.droplets.add(droplet)
        for chunk_num in droplet.num_chunks:
            self.chunk_to_droplets[chunk_num].add(droplet) #we document for each chunk all connected droplets        
        self.updateEntry(droplet) #one round of message passing

        
    def updateEntry(self, droplet):
        #removing solved segments from droplets
        for chunk_num in (droplet.num_chunks & self.done_segments):
            droplet.data = xor(droplet.data,self.chunks[chunk_num]) #subtract (ie. xor) the value of the solved segment from the droplet.
            droplet.num_chunks.remove(chunk_num)
            self.chunk_to_droplets[chunk_num].discard(droplet) #remove the edge between droplet and input segment.

        #solving segments when the droplet have exactly 1 segment
        if len(droplet.num_chunks) == 1: #the droplet has only one input segment
            lone_chunk = droplet.num_chunks.pop()
            self.chunks[lone_chunk] = droplet.data #assign the droplet value to the input segment (=entry[0][0])
            self.done_segments.add(lone_chunk) #add the lone_chunk to a data structure of done segments.
            self.droplets.discard(droplet) #remove the edge between the droplet and input segment
            self.chunk_to_droplets[lone_chunk].discard(droplet) #remove the edge between the input segment and the droplet
            #update other linked droplets
            for other_droplet in self.chunk_to_droplets[lone_chunk].copy():
                self.updateEntry(other_droplet)

    def add_seed(self, seed):
        self.seen_seeds.add(seed)

    def len_seen_seed(self):
        return len(self.seen_seeds)

    def isDone(self):
        if self.num_chunks - len(self.done_segments) > 0:
            return None 
        return True

    def chunksDone(self):
        return len(self.done_segments)

    def String(self):
        res = ''
        for x in self.chunks:
            res += ''.join(map(chr, x))
        return res
    
    def StringNoPadding(self):
        return self.String().rstrip('\0')

    def removePadding(self,pad):
        if pad != -1:
            self.chunks[-1] = self.chunks[-1][:-pad]

        crp = []
        for b in self.chunks[-1]:
            if 0 == b:
                break 
            crp.append(b)
        self.chunks[-1] = crp
        return crp
    
    def save(self,file_name, pad = -1):
        self.removePadding(pad)
        with open(file_name,'wb') as f:
            for c in self.chunks:
                f.write(bytes(c))
#             logging.info('saved')
            print('saved')
            f.close()
        
    def binString(self):
        bs = b''
        for c in self.chunks:
            bs += bytes(c)
        return bs
    
    def bchunks(self):
        chunks = []
        for c in self.chunks:
            chunks.append(bytes(c))
        return chunks
    
    def print_chunks(self):
        print(self.chunks)
        
    def display_chunks(self):
        i = 0
        not_none = []
        for x in self.chunks:
            print(i,''.join(map(chr, x)))
            i+=1
            if x!= None:
                not_none.append(i)
        return not_none
    
    def decode(self):
        f = open(self.in_file_name,'r')
        line = 0
        errors = 0
        solve_num = []
        while True:
            #read line
            try:     
                dna = f.readline().rstrip('\n')
            except:
                logging.info("After reading %d lines, %d chunks are done. So far: %d rejections (%f) %d barcodes", line, self.chunksDone(), errors, errors/(line+0.0), self.len_seen_seed())
                logging.info("Finished reading input file!")
                # print("After reading %d lines, %d chunks are done. So far: %d rejections (%f) %d barcodes" % (line, self.chunksDone(), errors, errors/(line+0.0), self.len_seen_seed()))
                # print('Finished reading input file!')
                return -1, solve_num, line, self.chunksDone(), errors
            if len(dna) == 0:
                logging.info("After reading %d lines, %d chunks are done. So far: %d rejections (%f) %d barcodes", line, self.chunksDone(), errors, errors/(line+0.0), self.len_seen_seed())
                logging.info("Finished reading input file!")
                # print("After reading %d lines, %d chunks are done. So far: %d rejections (%f) %d barcodes" % (line, self.chunksDone(), errors, errors/(line+0.0), self.len_seen_seed()))
                # print("Finished reading input file. Failed to decode!")
                return -1, solve_num, line, self.chunksDone(), errors
            line += 1
            
            seed, data = self.add_dna(dna)
            if seed == -1:
                errors += 1
            #logging
            if line % 200 == 0:
                logging.info("After reading %d lines, %d chunks are done. So far: %d rejections (%f) %d barcodes", line, self.chunksDone(), errors, errors/(line+0.0), self.len_seen_seed())
                # print("After reading %d lines, %d chunks are done. So far: %d rejections (%f) %d barcodes" % (line, self.chunksDone(), errors, errors/(line+0.0), self.len_seen_seed()))
                pass
            solve_num.append(self.chunksDone())

            if self.isDone():
                logging.info("After reading %d lines, %d chunks are done. So far: %d rejections (%f) %d barcodes", line, self.chunksDone(), errors, errors/(line+0.0), self.len_seen_seed())
                logging.info("Done!")
                # print("After reading %d lines, %d chunks are done. So far: %d rejections (%f) %d barcodes" % (line, self.chunksDone(), errors, errors/(line+0.0), self.len_seen_seed()))
                # print('done!')
                f.close()
                return 0, solve_num, line, self.chunksDone(), errors


