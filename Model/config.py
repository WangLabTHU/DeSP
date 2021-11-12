def genTm(prob):
    tm = []
    for i in range(4):
        row = [prob for i in range(4)]
        row[i] = 1 - 3* prob 
        tm.append(row)
    return tm

'ACGT'
TM_NGS = [
    [0.9993,0.00015,0.0004,0.00015],
    [0.0004,0.9994,0.00005,0.00015],
    [0.0003,0.00045,0.999,0.00025],
    [0.0001,0.00025,0.00035,0.9993],
]

TM_NNP = [
    [0.788,0.02,0.144,0.048],
    [0,0.992,0.008,0],
    [0.02,0.01,0.97,0],
    [0,0.015,0,0.985],
]

DEFAULT_DIC = {
    'syn_yield': 0.99,
    'syn_number': 30,
    'syn_pcrc': 9,
    'syn_pcrp': 0.8,
    'syn_performPCR': False,
    'syn_sub_prob': 0.001, # total: 0.004
    'syn_ins_prob': 0.002,
    'syn_del_prob': 0.001,
    'decay_er': 0.001,
    'decay_loss_rate': 0.3,
    'pcrc': 12,
    'pcrp':0.8,
    'pcrBias': 0.05,
    'sam_ratio': 0.01,
    'sam_to_number': 25,
    'seq_copies': 5000,
    'seq_prcp': 0.8,
    'seq_performPCR': False,
    'seq_depth': 10,
    'seq_TM': genTm(0.0015)
}

''' 
----------------------------------------------------------------------------------------
Parameters provided in arg passer format.
----------------------------------------------------------------------------------------
''' 
class ArgumentPasser:
    """Simple Class for passing arguments in arg object.
        Init all arttributes from a dictionary.
    """
    def __init__(self, dic):
        self.__dict__.update(dic)

DEFAULT_PASSER = ArgumentPasser(DEFAULT_DIC)