import numpy as np
from math import sqrt,log
import copy 
import Model.config as config
import time

BASE = np.array(['A','C','G','T'])
QUANT = {'A': 0, 'C':1, 'G':2, 'T':3}
qua2str = lambda qua: ''.join(BASE[qua])
str2qua = lambda dna: np.array([QUANT[base] for base in dna],dtype = 'uint8')

class DNA_Channel_Model:
    def __init__(self, Modules, arg = config.DEFAULT_PASSER):
        if Modules:
            self.Modules = Modules
        else:
            self.Modules = [
                ('synthesizing',Synthesizer(arg)),
                # ('decaying',Decayer(arg)),
                ('pcring',PCRer(arg = arg)),
                ('sampling',Sampler(arg = arg)),
                ('sequencing',Sequencer(arg))
            ]

    def __call__(self, dnas, inspectFunction = None, print_state = True):
        if print_state: print('Model running...  ', end = '\r')
        ast = time.time()
        for stage_name, module in self.Modules:
            st = time.time()
            if print_state: print(stage_name + '.... ', end = '')
            dnas = module(dnas)
            if print_state: print(' Done. Time spent: ' + str(round(time.time() - st,3)) +'s')
            if inspectFunction: inspectFunction(dnas)
        if print_state: print(f'Simulation done. Time spent: {round(time.time()-ast,3)}s')
        return dnas
    
class Synthesizer:
    def __init__(self, arg):
        self.Yield = arg.syn_yield
        self.N = arg.syn_number
        self.pcrc = arg.syn_pcrc
        self.pcrp = arg.syn_pcrp
        self.performPCR = arg.syn_performPCR
        self.probS = arg.syn_sub_prob
        self.probD = arg.syn_del_prob
        self.probI = arg.syn_ins_prob

        self.syn = Syn_D(self.Yield, self.N)
        self.err = ErrorAdder(self.probS, self.probD, self.probI)
        self.pcr = PCRer(self.pcrc,self.pcrp)
    
    def __call__(self, dnas):
        dnas = self.syn(dnas)
        dnas = self.err(dnas)
        if self.performPCR: dnas = self.pcr(dnas)
        return dnas

class Decayer:
    def __init__(self, arg):
        self.error_rate = arg.decay_er
        self.loss_rate = arg.decay_loss_rate
        self.constructTm()
        
        self.sam = Sampler(1-self.loss_rate)
        self.err = ErrorAdder(probS = 0, probD = 0, probI = 0, TM = self.TM)

    def constructTm(self):
        Tm = [[0 for i in range(4)] for i in range(4)]
        Tm[1][3] = Tm[2][1] = self.error_rate
        Tm[1][1] = Tm[2][2] = 1-self.error_rate # C2T, G2A
        Tm[0][0] = Tm[3][3] = 1
        self.TM = Tm
    
    def __call__(self, dnas):
        dnas = self.sam(dnas)
        dnas = self.err(dnas)
        return dnas

class Sequencer:
    def __init__(self,arg):
        self.copies_required = arg.seq_copies
        self.pcrp = arg.seq_prcp
        self.performPCR = arg.seq_performPCR
        
        self.seq_depth = arg.seq_depth
        
        self.TM = arg.seq_TM
    
    def __call__(self, dnas):
        if self.performPCR:
            dnas = self.pcr(dnas)
        dnas = self.sample(dnas)
        self.E = ErrorAdder(probI = 0.00001, probD= 0.00001, TM = self.TM)
        dnas = self.E(dnas)
        return dnas
    
    def pcr(self,dnas):
        rNs = [dna['num'] for dna in dnas]
        average_copies = sum(rNs) / len(rNs)
        amplify_ratio = self.copies_required / average_copies
        self.pcrc = int(log(amplify_ratio) / log(self.pcrp+1))
        dnas = PCRer(self.pcrc, self.pcrp)(dnas)
        return dnas        

    def sample(self, dnas):
        rNs = [dna['num'] for dna in dnas]
        average_copies = sum(rNs) / len(rNs)
        self.sample_ratio = self.seq_depth / average_copies
        dnas = Sampler(self.sample_ratio)(dnas)
        return dnas
        
class Syn_D:
    def __init__(self, Yield = 0.99, N = 30):
        self.Yield = Yield
        self.N = N
        
    def distribution(self):
        return np.random.binomial(self.N, self.p)
    
    def __call__(self,dnas):
        self.L = len(dnas[0])
        self.p = self.Yield ** self.L

        out = []
        for dna in dnas:
            n = self.distribution()
            out.append({'ori':dna, 'num':n,'re':[[n,[]]]})
        return out
    
class Sampler:
    def __init__(self, p=0.001, sam_to_number = False, arg = None):
        if arg: 
            self.p = arg.sam_ratio
            self.sam_to_number = arg.sam_to_number
        else: 
            self.p = p
            self.sam_to_number = sam_to_number
    
    def distribution(self,N):
        return np.random.binomial(N,self.p)

    def run(self,re_dnas):
        markers = []
        for i,dna in enumerate(re_dnas):
            dna[0] = self.distribution(dna[0])
            if dna[0] > 0: markers.append(i)
        re_dnas = [re_dnas[i] for i in markers]
        return re_dnas
    
    def __call__(self,dnas, in_place = False):
        if not in_place:
            out_dnas = copy.deepcopy(dnas)
        else:
            out_dnas = dnas
        
        if self.sam_to_number:
            rNs = [dna['num'] for dna in dnas]
            average_copies = sum(rNs) / len(rNs)
            self.p = self.sam_to_number / average_copies
            
        for dna in out_dnas:
            dna['re'] = self.run(dna['re'])
            dna['num'] = sum([tp[0] for tp in dna['re']])
        return out_dnas

# class PCRer:
#     def __init__(self,N = 16, p = 0.7, arg = None):
#         if arg: 
#             p = arg.pcrp
#             N = arg.pcrc
#         self.p = p
#         self.N = N
#         self.u0 = (1+p)**N
#         self.sigma0 = np.sqrt((1-p) / (1+p) * ((1+p)**(2*N) - (1+p)**N))
    
#     def distribution(self,ori):
#         assert ori >= 0
#         return max(int(np.random.normal(self.u0 * ori, self.sigma0 * sqrt(ori))),0)

#     def run(self,re_dnas):
#         out = []
#         for dna in re_dnas:
#             dna[0] = self.distribution(dna[0])
#             if dna[0] > 0:
#                 out.append(dna)
#         return out
    
#     def __call__(self,dnas,in_place = False):
#         if not in_place:
#             out_dnas = copy.deepcopy(dnas)
#         else:
#             out_dnas = dnas
            
#         for dna in out_dnas:
#             dna['re'] = self.run(dna['re'])
#             dna['num'] = sum([tp[0] for tp in dna['re']])
#         return out_dnas
    
class PCRer:
    def __init__(self,N = 16, p = 0.7, pBias = 0.05, arg = None):
        if arg: 
            p = arg.pcrp
            N = arg.pcrc
            pBias = arg.pcrBias
        self.p = p
        self.N = N
        self.pBias = pBias
        
        self.u0 = (1+p)**N
        self.sigma0 = np.sqrt((1-p) / (1+p) * ((1+p)**(2*N) - (1+p)**N))
    
    def distribution(self,ori):
        assert ori >= 0
        p = np.random.uniform(self.p - self.pBias, self.p + self.pBias)
        N = self.N
        u0 = (1+p)**N
        sigma0 = np.sqrt((1-p) / (1+p) * ((1+p)**(2*N) - (1+p)**N))
        return max(int(np.random.normal(u0 * ori, sigma0 * sqrt(ori))),0)

    def run(self,re_dnas):
        out = []
        for dna in re_dnas:
            dna[0] = self.distribution(dna[0])
            if dna[0] > 0:
                out.append(dna)
        return out
    
    def __call__(self,dnas,in_place = False):
        if not in_place:
            out_dnas = copy.deepcopy(dnas)
        else:
            out_dnas = dnas
            
        for dna in out_dnas:
            dna['re'] = self.run(dna['re'])
            dna['num'] = sum([tp[0] for tp in dna['re']])
        return out_dnas

class ErrorAdder:
    def __init__(self,probS = 0.001, probD = 0.0005, probI = 0.0005, TM = None):
        if TM != None:
            self.TM = TM
            self.all_equal = 0
        else:
            self.TM = genTm(probS)
            self.all_equal = 1
        
        self.probD = probD
        self.probI = probI

    def genNewError(self,dna):
        Errors = []
        for i,base in enumerate(['A','C','G','T']):
            Pi = np.where(dna==base)[0]
            subi = np.random.choice(['A','C','G','T'],size = Pi.size, p = self.TM[i])
            subPi = np.where(subi != base)[0]
            for pos in subPi:
                Errors.append((Pi[pos],'s', subi[pos]))
        delP = np.where(np.random.choice([False,True],size = len(dna), p = [1-self.probD,self.probD]))[0]
        insP = np.where(np.random.choice([False,True],size = len(dna), p = [1-self.probI,self.probI]))[0]
        Errors += ([(pos,'-',dna[pos]) for pos in delP] + [(pos,'+',np.random.choice(['A','T','C','G'])) for pos in insP])
        return Errors

    def run(self,ori_dna,re_dnas):
        ori_dna = np.array(list(ori_dna))
        new_types = []
        for re_dna in re_dnas:
            for i in range(re_dna[0]):
                new_error = self.genNewError(ori_dna)
                if len(new_error) > 0:
                    new_types.append([1, re_dna[1] + new_error])
                    re_dna[0] -= 1
        return re_dnas + new_types

    def __call__(self,dnas,in_place = False, apply = True):
        if not in_place:
            out_dnas = copy.deepcopy(dnas)
        else:
            out_dnas = dnas
            
        for dna in out_dnas:
            dna['re'] = self.run(dna['ori'], dna['re'])

        if apply:
            out_dnas = self.apply_batch(out_dnas)
        return out_dnas
    
    # apply errors to dnas
    def apply(self, ori_dna, errors):
        dna = list(ori_dna)
        errors.sort(key = lambda x: x[0])
        # substitutions
        for error in errors:
            pos, tp, base = error
            if tp == 's':
                dna[pos] = base
        # del / insertions
        for error in errors:
            bias = 0
            pos, tp, base = error
            if tp == '-':
                try:
                    dna.pop(pos + bias)
                except:
                    # print('pop index error:', pos + bias)
                    break
                bias -= 1
            elif tp == '+':
                dna.insert(pos,base)
                bias += 1
        dna = ''.join(dna)
        return dna
    
    def apply_batch(self, dnas):
        for dna in dnas:
            ori_dna = dna['ori']
            re = []
            for re_dna in dna['re']:
                if re_dna[0] == 0: pass
                re.append([re_dna[0],re_dna[1],self.apply(ori_dna, re_dna[1])])
            dna['re'] = re
        return dnas

def genTm(prob):
    tm = []
    for i in range(4):
        row = [prob for i in range(4)]
        row[i] = 1 - 3* prob 
        tm.append(row)
    return tm