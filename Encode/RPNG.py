from math import log, floor, sqrt
import numpy as np
from numpy.random import RandomState
import random
import logging

    
def LCG(seed,mi,ma,n,key=[2**32,1103,12345]):
    rd = []
    m,a,c=key
    for i in range(n):
        while(1):
            seed = (a * seed + c) % m
            tmp=int((ma-mi) * seed /float(m-1) + mi)

            if(tmp not in rd):
                break
        rd.append(tmp)
    return rd

def gen_tau(S, K, delta):
    """The Robust part of the RSD, we precompute an
    array for speed
    """
    pivot = int(floor(K/S))

    val1 =  [S/K * 1/d for d in range(1, pivot)] 
    val2 =  [S/K * log(S/delta)] 
    val3 =  [0 for d in range(pivot, K)] 
 
    return val1 + val2 + val3

def gen_rho(K):
    """The Ideal Soliton Distribution, we precompute
    an array for speed
    """
    return [1.0/K] + [1.0/(d*(d-1)) for d in range(2, K+1)]


def gen_mu(K, S, delta):
    """The Robust Soliton Distribution on the degree of 
    transmitted blocks
    """

    tau = gen_tau(S, K, delta)
    rho = gen_rho(K)

    Z = sum(rho) + sum(tau)
    mu = [(rho[d] + tau[d])/Z for d in range(K)]
    return (mu, Z)

def gen_rsd_cdf(K, S, delta):
    """The CDF of the RSD on block degree, precomputed for
    sampling speed"""

    mu, Z = gen_mu(K, S, delta)
    cdf = np.cumsum(mu)
    return cdf, Z

class PRNG(object):
    """A Pseudorandom Number Generator that yields samples
    from the set of source blocks using the RSD degree
    distribution described above.
    """

    def __init__(self, K, delta, c, np = None, enc_num=1,enc_key=[2**32,1103,12345]):
        """Provide RSD parameters on construction
        # K is the number of segments
        # delta and c are parameters that determine the distribution
        #np is to use numpy random number generator which is faster
        """

        self.K = float(K)
        self.K_int = int(K)
        self.delta = delta
        self.c = c

        S = self.calc_S()
        cdf, Z = gen_rsd_cdf(K, S, delta)
        self.cdf = cdf
        self.Z = Z

        #self.inter = inter.interp1d(np.concatenate(([0], cdf)), range(0,K+1))
        self.np_rand = RandomState(1)
        self.np = np

        self.state = 1
        
        self.enc_num=enc_num
        self.enc_key=enc_key

    def calc_S(self):
        """ A helper function to calculate S, the expected number of degree=1 nodes
        """
  
        K = self.K
        S = self.c * log(self.K/self.delta) * sqrt(self.K) 
        self.S = S
        return S


    def get_S(self):
        return self.S


    def set_seed(self, seed):
        """Reset the state of the PRNG to the 
        given seed
        """

        self.state = seed
    
    def get_state(self):
        """Returns current state of the linear PRNG
        """

        return self.state


    def get_src_blocks_wrap(self, seed=None):
        #a wrapper function to get source blocks.
        #if np flag is on, it will use a numpy-based method.
        #otherwise, it will use the native python random function.
        #np is faster but in compatible with python random which implemented in previous versions.
        if self.enc_num:
            return self.get_src_blocks_enc(seed)
        elif self.np:
            return self.get_src_blocks_np(seed)
        else:
            return self.get_src_blocks(seed)
        
    def get_src_blocks_enc(self,seed=None):
        if seed:
            self.state = seed

        blockseed = self.state
        random.seed(self.state)
        
        d = self._sample_d()

        nums = LCG(blockseed,0,self.K_int,d,self.enc_key)
        return blockseed, d, nums
        
    def get_src_blocks(self, seed=None):
        """Returns the indices of a set of `d` source blocks
        sampled from indices i = 1, ..., K-1 uniformly, where
        `d` is sampled from the RSD described above.
        """

        if seed:
            self.state = seed

        blockseed = self.state
        random.seed(self.state)
        
        d = self._sample_d()

        nums = random.sample(range(self.K_int), d)
        return blockseed, d, nums


    def get_src_blocks_np(self, seed=None):
        """Returns the indices of a set of `d` source blocks
        sampled from indices i = 1, ..., K-1 uniformly, where
        `d` is sampled from the RSD described above.
        Uses numpy for speed.
        """

        if seed:
            self.state = seed


        blockseed = self.state
        self.np_rand.seed(self.state)
        
        d = self._sample_d_np()
        nums = self.np_rand.randint(0, self.K_int, d)
        return blockseed, d, nums

    def _sample_d_np(self):
        """Samples degree given the precomputed
        distributions above. Uses numpy for speed"""

        p = self.np_rand.rand()
        for ix, v in enumerate(self.cdf):
            if v > p:
                return ix + 1
        return ix + 1        


    def _sample_d_inter(self):
        """Samples degree given the precomputed
        distributions above using interpolation
        """

        p = random.random()
        return int(self.inter(p))+1 #faster than math.ceil albeit can return the wrong value...

    # Samples from the CDF of mu
    def _sample_d(self):
        """Samples degree given the precomputed
        distributions above"""

        p = random.random()

        for ix, v in enumerate(self.cdf):
            if v > p:
                return ix + 1
        
        return ix + 1
        
    
def lfsr(state, mask):
    #Galois lfsr:
    result = state
    nbits = mask.bit_length()-1
    while True:
        result = (result << 1)
        xor = result >> nbits
        if xor != 0:
            result ^= mask

        yield result

def lfsr32p():
    #this function returns a hard coded polynomial (0b100000000000000000000000011000101).
    #The polynomial corresponds to 1 + x^25 + x^26 + x^30 + x^32, which is known 
    #to repeat only after 32^2-1 tries. Don't change unless you know what you are doing.
    return 0b100000000000000000000000011000101

def lfsr32s():
    #this function returns a hard coded state for the lfsr (0b001010101)
    #this state is the inital position in the register. You can change it without a major implication.
    return 0b001010101

def test():
    #run the test to see a stream of seed by the polynomial
    for pattern in lfsr(0b001, 0b100000000000000000000000011000101):
        print(pattern)
        