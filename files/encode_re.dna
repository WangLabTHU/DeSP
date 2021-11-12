import streamlit as st
import plotly.io as pio
import matplotlib.pyplot as plt
from Model.Model import * 
import plotly.express as px
from Analysis.Analysis import dna_chunk, plot_oligo_number_distribution, plot_error_distribution,save_simu_result
from Analysis.Helper_Functions import preprocess
from Encode.DNAFountain import *
pio.templates.default = "plotly_white"
import Model.config as config 

def inspect(dnas, num_th = 1):
    fig = plt.figure(figsize = (12,6))
    plt.subplot(1,2,1)
    plot_oligo_number_distribution(dnas, show = False)
    plt.subplot(1,2,2)
    plot_error_distribution(dnas)
    st.write(fig)

# ------------ assigning parameters --------------- #
arg = config.DEFAULT_PASSER

arg.syn_number = st.sidebar.slider('Syn number', min_value = 10, max_value = 50, value = 30)
arg.syn_sub_prob = st.sidebar.number_input('Syn Error rate', min_value = 0.0, max_value = 0.1, value = 0.001) / 3 # 3 kinds of substitutions
arg.syn_yield = st.sidebar.slider('Syn Yield', min_value = 0.98, max_value = 0.995, value = 0.99)

arg.sam_ratio = st.sidebar.number_input('Sampling ratio',min_value = 0.0, max_value =1.0,value = 0.005)
arg.seq_depth = st.sidebar.slider('Seq Depth', min_value = 1, max_value = 100, value = 10)

alpha = st.sidebar.slider('Alpha', min_value = 0.1, max_value = 0.5, value = 0.15)
rs = st.sidebar.slider('RS', min_value = 0, max_value = 8, value = 2)
num_th = int(rs / 2)

st.header('Choosing proper reduancy with the model')
# Encode
st.header('Encode')
data = preprocess('lena.jpg',20)
st.image('lena.jpg',width = 300)
'Lena.jpg split into ', len(data), ' data chunks.'
f = DNAFountain(data,alpha,rs = rs)
good, tries = f.encode()
'Data encoded into ' ,good, ' DNA strands after ', tries, ' tries.'
'Saved to lena.dna.'
f.save('lena.dna')
#Model
st.header('Channel')
'running..'
with open('lena.dna') as f:
    dnas = f.readlines()
in_dnas = [dna.split('\n')[0] for dna in dnas]

Model = DNA_Channel_Model(None,arg)
out_dnas = Model(in_dnas)
inspect(out_dnas)
save_simu_result(out_dnas,'lena_simu.dna')

# ---------------- Decoding --------------- #
st.header('Decode')

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
    

g = Glass('files/lena_simu.dna', len(data), rs = rs)
solve_num, ret = g.decode()
plot_solve_num(solve_num, ret)
if ret == 0: 
    st.write('Decoding succeeded!')
    g.save('relena.jpg')
else:
    st.write('Decoding Failed-.-')
    st.write('Change the parameter, and run the model again.')



