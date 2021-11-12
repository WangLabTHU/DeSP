from Analysis.html_printer import html_table, color_bold,background_color_bold
from prettytable import PrettyTable

import numpy as np
from copy import deepcopy

import plotly.graph_objects as go 
import matplotlib.pyplot as plt
import seaborn as sns

from difflib import ndiff

def inspect_distribution(dnas, num_th = 1, show = True):
    plt.figure(figsize = (7,3))
    plt.subplot(1,2,1)
    zN = plot_oligo_number_distribution(dnas)

    plt.subplot(1,2,2)
    esn = plot_error_distribution(dnas)
    
    if show: plt.show()
    else: return zN, esn

def inspect_number_only(dnas, num_th = 1):
    rNs = [dna['num'] for dna in dnas]
    lost_num = sum([rN == 0 for rN in rNs])
    error_nums = error_distribution(dnas)
    error_num = sum([n >= num_th for n in error_nums])
    # print(f'{lost_num} lost. {error_num} strands have {num_th} errors or more. sum: {error_num  + lost_num}')
    return lost_num, error_num

def plot_oligo_number_distribution(dnas):
    N = len(dnas)
    rNs = [dna['num'] for dna in dnas]
    zN = np.sum(np.array(rNs)==0)
    label = f'{zN}({round(zN/N*100,3)}%) lost' + f'\n{round(sum(rNs)/len(rNs),2)} copies'
    
    ax = sns.distplot(rNs, hist=False, color="r", kde_kws={"shade": True})
    plt.xlabel('oligo number')
    plt.ylabel('frequency')

    plt.text(0.95,0.95,f'{round(sum(rNs)/len(rNs),2)}',fontsize = 12,  transform = ax.transAxes, va = 'top', ha = 'right')
    plt.text(0.05,0.95,f'{zN}({round(zN/N*100,3)}%)',fontsize = 12,color = 'r',  transform = ax.transAxes, va = 'top', ha = 'left')
    return zN

def plot_error_distribution(dnas, th = 1):
    error_nums = error_distribution(dnas)
    esn = sum([n >= th for n in error_nums])
    label = f'{esn}({round(esn/len(error_nums)*100,2)}%) with more than {th} errors'
    ax = sns.distplot(error_nums, hist= False, color="r", kde_kws={"shade": True, 'bw': 0.1})

    plt.text(0.95,0.95,f'{esn}({round(esn/len(error_nums)*100,2)}%)',fontsize = 12, color = 'r',transform = ax.transAxes, va = 'top', ha = 'right')

    plt.xlabel('error number')
    plt.ylabel('frequency')
    return esn

def examine_strand(dnas, index = 0):
    dna = dnas[index]
    dc = dna_chunk(dna)
    dc.plot_re_dnas()
    return dc.plot_voting_result()

def error_distribution(dnas):
    error_nums = []
    for dna in dnas:
        if dna['num'] == 0: continue
        error_num = dna_chunk(dna).voting_error()
        error_nums.append(error_num)
    return error_nums

def save_simu_result(dnas, file_name = 'simu_res.dna',ignore_index = None):
    with open(file_name, 'w') as f:
        for i,dna in enumerate(dnas):
            if dna['num'] == 0: continue
            if ignore_index and (i in ignore_index): continue
            re_dna = dna_chunk(dna).voting_result()
            f.write(re_dna + '\n')
    
class dna_chunk:
    def __init__(self,dna_in,env = 'jupyter'):
        self.ori_dna = dna_in['ori']
        self.re_dnas = deepcopy(dna_in['re'])
        self.rN = dna_in['num']
        
        self.env = env
        self.BASE = ['A','T','C','G']
        self.CMAP_PY = {'s':'33',
        '-':'31',
        '+':'36'
       }
        self.CMAP = {'s':'red',
            '-':'#9900FF',
            '+':'#FFFF00'
        }
        self.BASE_COLOR = {
            'A':'green',
            'T':'goldenrod',
            'C':'blue',
            'G':'red',
        }
        
        self.R = None

    def plot_re_dnas(self):
        if self.env == 'jupyter':
            return self.plot_re_dnas_jupyter()
        else:
            return self.plot_re_dnas_html()

    def plot_re_dnas_html(self):
        table = []
        for re_dna in self.re_dnas:
            num,error,dna = re_dna
            dna = self.plot_error_dna_html(dna,error)
            table.append([(dna,'style = "color:#cccccc"'), num])
        return html_table().print(table, ['DNA','Num'])

    def plot_error_dna_html(self,dna, error):
        out = ''
        prv = cur = 0
        for e in error:
            cur,tp,base = e
            c = self.CMAP[tp]
            out += dna[prv:cur] + color_bold(base,c)
            prv = cur+1
        out += dna[cur+1:]
        return out

    def plot_re_dnas_jupyter(self, compress = True):
        table = PrettyTable(['read num','re_dna'])
        for i in range(len(self.re_dnas)):
            num, error,dna = self.re_dnas[i]
            table.add_row([num,self.plot_error_dna_jupyter(error,compress)])
        print(table)
    
    def plot_error_dna_jupyter(self, error,compress = True):
        dna = ''
        prv = 0
        color = lambda c,s: '\033[1;'+ c + 'm' + s + '\033[0m'
        for e in error:
            pos, typ, base = e
            if compress:
                pos = int(e[0]/2)
            else: 
                pos = e[0]
            gap = pos - prv
            prv = pos
            c = self.CMAP_PY[typ]
            base_change = color(c,base)
            dna = dna + ' '*gap + base_change
        return dna
    
    def vote(self):
        if self.R: return self.R
        if self.rN == 0: return None
        R = [{'A':0,'T':0,'G':0,'C':0} for i in range(len(self.ori_dna))]
        # compute subs
        for re_dna in self.re_dnas:
            read_num, error, _ = re_dna
            for e in error:
                pos,tp,base = e
                if tp == 's':
                    R[min(len(self.ori_dna)-1,pos)][base]+= read_num
        # original base
        for i in range(len(self.ori_dna)):
            base = self.ori_dna[i]
            r = R[i]
            r[base] = self.rN - sum([r[b] for b in self.BASE])
        # to prob
        for r in R:
            for base in self.BASE:
                r[base] = r[base] / self.rN * 100
        self.R = R
        return R
    
    def voting_error(self):
        self.vote()
        error_num = 0 
        for i in range(len(self.ori_dna)):
            comp = self.R[i]
            ori_base = self.ori_dna[i]
            for base in ['A','T','C','G']:
                if base is not ori_base:
                    if comp[base] >= comp[ori_base]:
                        error_num += 1
                        break
        return error_num
    
    def voting_result(self):
        if not self.vote(): return None
        re_dna = ''
        for comp in self.R:
            maxP = -1
            next_base = 'N'
            for base in ['A','T','C','G']:
                if comp[base] > maxP:
                    next_base = base
                    maxP = comp[base]
            re_dna = re_dna + next_base
        return re_dna
                    
    def plot_voting_result(self):
        if not self.vote():
            print('NOOOOOOO! Strand is lost.')
            return None
        data = []
        y = [self.R[i][self.ori_dna[i]] for i in range(len(self.ori_dna))]
        GT = go.Scatter(y=y,marker_color='gray',mode='lines',line_width = 1)
        data.append(GT)

        for b in self.BASE:
            y = [p[b] for p in self.R]
            prop = go.Scatter(y=y,marker_color = self.BASE_COLOR[b],mode='markers',marker_size = 4,hovertext = [str(p) for p in self.R])
            data.append(prop)
        
        fig = go.Figure(data) 
        fig.update_layout(height=300,width = 900, showlegend = False,title="Voting Result", xaxis_title="Position",yaxis_title="Frequency",)
        return fig
    
def diff(a,b):
    return sp_to_final(diff_sp(a,b))

def diff_sp(a,b):
    diff_arr = []
    for i,s in enumerate(ndiff(a,b)):
        if(s[0] != ' '):
            diff_arr.append([i,s[0],s[2]])
            
    diff_arr_sp = []
    sp = []
    for diff in diff_arr:
        if(sp == []):
            sp.append(diff)
        else:
            if(diff[0] == sp[-1][0]+1):
                sp.append(diff)
            else:
                diff_arr_sp.append(sp)
                sp = [diff]
    if(sp!=[]):
        diff_arr_sp.append(sp)
    return diff_arr_sp

def sp_to_final(diff_arr_sp):
    #d a to s
    pos_d = 0
    diff_arr = []
    for sp in diff_arr_sp:
        diff_type = [d[1] for d in sp]
        d = diff_type.count('-')
        a = diff_type.count('+')
        
        if( d-a <= 0):
            pos_d += d
            for df in sp:
                df[0] -= pos_d
            
            for df in sp[d:d*2]:
                df[1] = 's'
            
            diff_arr += sp[d:]
        else:
            for df in sp[:d-a]:
                df[0] -= pos_d
            
            pos_d += a
            for df in sp[d:]:
                df[0] -= pos_d
                df[1] = 's'
              
            diff_arr += sp[:d-a] + sp[d:]
    return diff_arr
    
    
    
    
    
    

