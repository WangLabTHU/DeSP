import streamlit as st
import plotly.io as pio
import matplotlib.pyplot as plt
from Model.Model import * 
import plotly.express as px
from Analysis.Analysis import dna_chunk, plot_oligo_number_distribution, plot_error_distribution, save_simu_result
from Analysis.Fountain_analyzer import FT_Analyzer_Simplified
from Encode.DNAFountain import *
pio.templates.default = "plotly_white"
import Model.config as config

# ---------------------------- Choosing parameters ---------------------------- # 

def inspect(dnas, num_th = 1, inspect_index = 20):
    fig = plt.figure(figsize = (12,6))
    plt.subplot(1,2,1)
    plot_oligo_number_distribution(dnas)
    plt.subplot(1,2,2)
    plot_error_distribution(dnas,th = num_th)
    st.write(fig)
    dc = dna_chunk(dnas[inspect_index],'html')
    table = dc.plot_re_dnas()
    st.markdown(table, unsafe_allow_html = True)
    st.write(dc.plot_voting_result())

# --------------------------- assigning parameters ------------------------ # 
arg = config.DEFAULT_PASSER
st.sidebar.subheader('Parameters of DNA data storage channel')
arg.syn_number = st.sidebar.slider('Syn number', min_value = 10, max_value = 50, value = 30)
arg.syn_sub_prob = st.sidebar.number_input('Syn Error rate', min_value = 0.0, max_value = 0.1, value = 0.01) / 3 # 3 kinds of substitutions
arg.syn_yield = st.sidebar.slider('Syn Yield', min_value = 0.98, max_value = 0.995, value = 0.99)

arg.pcrc = st.sidebar.slider('PCR cycle',min_value = 0, max_value =20,value = 12)
arg.pcrp = st.sidebar.number_input('PCR prob',min_value = 0.5, max_value = 1.0,value = 0.8)

arg.sam_ratio = st.sidebar.number_input('Sampling ratio',min_value = 0.0, max_value =1.0,value = 0.005)
arg.seq_depth = st.sidebar.slider('Seq Depth', min_value = 1, max_value = 100, value = 10)
seq_platform = st.sidebar.selectbox('Sequencing Platform',['Illumina Sequencing','Nanopore'])

index = st.sidebar.slider('inspect index', max_value = 600, value = 0)

st.sidebar.subheader('Parameters of Fountain code')
alpha = st.sidebar.slider('Alpha', min_value = 0.25, max_value = 1.0, value = 0.5)
rs = st.sidebar.slider('RS', min_value = 4, max_value = 10, value = 4)
num_th = int(rs / 2)

if seq_platform == 'Illumina Sequencing':
    arg.seq_TM = config.TM_NGS
else:
    arg.seq_TM = config.TM_NNP

# -------------------------- encoding ------------------------------------ #
st.header('DNA-D2S: a systematic error simulation Model for DNA Data Storage channel')
st.markdown('In this demonstrative web app, you will encode lena.jpg to DNA with DNA fountain code,\
     run the simulation process to see how errors are generated and passed through different stages,\
      and optimize encoding parameters according to the noise structures of the channel.')
st.markdown('You can assign parameters of DNA data storage channel and fountain code in the sidebar.\
    Play with different combinations of parameters to see how the noise structures and optimal encoding designs are influenced.') 

file_name = 'lena.jpg'
file_name, suffix = file_name.split('.')
file_name = 'files/' + file_name

in_file_name = file_name + '.' + suffix
in_dna_name = file_name + '.dna'
out_dna_name = file_name + '_simu.dna'
out_file_name = file_name + '_re.' + suffix

print(in_file_name)
st.header('Encoding the file into DNA')
if suffix in ['jpg','png']:
    st.image(in_file_name, width = 300)
data,pad = preprocess(in_file_name,20)
in_file_name, ' loaded and split into ', len(data), ' data chunks.'
f = DNAFountain(data,alpha,rs = rs)
good, tries = f.encode()
'Data encoded into ' ,good, ' DNA strands after ', tries, ' tries.'
'Saved to ', in_dna_name
f.save(in_dna_name)


# --------------------------- error simulation ---------------------------------- #

st.header('Error simulation of the DNA data storage channel')

st.subheader('Load Data')
with open(in_dna_name) as f:
    dnas = f.readlines()
in_dnas = [dna.split('\n')[0] for dna in dnas]

in_dna_name, ' loaded: ', len(in_dnas), ' strands of length ', len(in_dnas[0])
'Sequence ', index, ' will be inspected in detail to show how errors are formed in one sequence.\
     You can choose another sequence to inspect by altering the **inspect index** in the sidebar.' 
'Three figures will be depicted for each stage:'
'1. Oligo copy number distribution and voting error number distribution.'
'2. Error types of sequence ', index
'3. Voting results of current copies of sequence ', index
'The simulation process now begins:'

st.subheader('Synthesis')
SYN = Synthesizer(arg)
dnas_syn = SYN(in_dnas)
inspect(dnas_syn,inspect_index = index)

st.subheader('Decay')
DEC = Decayer(arg)
dnas_dec = DEC(dnas_syn)
inspect(dnas_dec,inspect_index = index)

st.subheader('PCR')
PCR = PCRer(arg = arg)
dnas_pcr = PCR(dnas_dec)
inspect(dnas_pcr,inspect_index = index)

st.subheader('Sampling')
SAM = Sampler(arg = arg)
dnas_sam = SAM(dnas_pcr)
inspect(dnas_sam,inspect_index = index)

st.subheader('Sequencing')
SEQ = Sequencer(arg)
dnas_seq = SEQ(dnas_sam)
inspect(dnas_seq,inspect_index = index)
save_simu_result(dnas_seq,out_dna_name)
'Simulation results saved to ', out_dna_name

# --------------------------- decoding ---------------------------- #
st.header('Decoding')

def plot_solve_num(solve_num, ret):
    label = f'{len(solve_num)} reads -> {solve_num[-1]} solved'
    if ret == 0 :label += ';Succeed.'
    else: label += ';Fail'
    fig = plt.figure(figsize = (6,4), dpi = 200)
    plt.plot(solve_num,color = 'r', label = label)
    plt.xlabel('read nums')
    plt.ylabel('solve nums')
    plt.legend()
    st.write(fig)

'Trying to decode from sequencing readouts.'
g = Glass(out_dna_name, len(data), rs = rs)
ret, solve_num, lineRead, chunksDone, errors = g.decode()
plot_solve_num(solve_num, ret)
if ret == 0: 
    st.write('Decoding succeeded!')
    g.save(out_file_name,pad)
    st.markdown(f'**{out_file_name}**:')
    st.image(in_file_name, width = 300)
else:
    st.write('Decoding Failed-.-')

# ------------------------ optimizing ----------------------------- #
st.header('Encoding Optimization')
'The encoding-simulation-decoding process should have been run several times to estimate the distributions.'
'To save computation here, we only run the process once and use the observed values and pre-computed values to perform the estimation. Some bias might be introduced because of this.'

loc = 500
scale = 7.5
N = len(data)
Ld = len(data[0])
FA = FT_Analyzer_Simplified(N,Ld,alpha,loc,scale,dnas_seq)

st.subheader('Choosing RS length: ')
fig, rs = FA.choose_rs()
st.write(fig)
'According to D(k), rs can be set to ', rs

st.subheader('Choosing Alpha ')
fig = FA.choose_alpha()
st.write(fig)
'Alpha can be selected from the second graph to meet a specific success possibility requirement.'