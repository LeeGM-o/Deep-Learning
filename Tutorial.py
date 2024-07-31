#1. Initialization
import torch
import os
from pyfasta import Fasta
import numpy as np
import pandas as pd
os.environ['CUDA_VISIBLE_DEVICES']='0'
from model import EpiGePT
from model.utils import *
from model.config import *

#2. Load pretrained model
model = EpiGePT.EpiGePT(WORD_NUM,TF_DIM,BATCH_SIZE)
model.load_state_dict(torch.load('pretrainModel/model.ckpt',map_location='cuda:0')['state_dict'])
genome = Fasta('hg38.fa')
model.eval()
model = model.cuda()

#3. Predict
acgt2num = {'A': 0,
            'C': 1,
            'G': 2,
            'T': 3}

def seq2mat(seq):
    seq = seq.upper()
    h = 4
    w = len(seq)
    mat = np.zeros((h, w), dtype=bool)  # True or false in mat
    for i in range(w):
        if seq[i] != 'N':
            mat[acgt2num[seq[i]], i] = 1.
    return mat

def model_predict(model,chrom,start,end,tf_feature):
    seq = genome[chrom][start:end]
    seq_embeds = seq2mat(seq)
    seq_embeds = np.array(seq_embeds,dtype='float16')#(50,)
    seq_embeds = np.expand_dims(seq_embeds,axis=0)
    seq_embeds = torch.from_numpy(seq_embeds)
    tf_feature = np.pad(tf_feature,((0, 0), (0, 1)),'constant',constant_values = (0,0))
    tf_feature = np.expand_dims(tf_feature,axis=0)
    tf_feature = torch.from_numpy(tf_feature)
    seq_embeds = seq_embeds.type(torch.FloatTensor)
    tf_feature = tf_feature.type(torch.FloatTensor)
    seq_embeds = seq_embeds.to('cuda')
    tf_feature = tf_feature.to('cuda')
    signals = model(seq_embeds,tf_feature)
    np_signals = signals.cpu().detach().numpy()
    return np_signals

#3-1 Prepare input: motif score for transcription factor binding and gene expression
#1) Prepare tf expression
# get information about 711 TFs
tf_list, motif2tf, isContain = get_tf_motif_match()
print(tf_list[:10])

ref_tf_expression = pd.read_csv('reference_TPM_tf_expr.csv',header=0,sep='\t',index_col=[0])

def quantile_norm(df,ref_exp):
    """quantile normalization
    Arg:
        df: Pandas DataFrame with TFs x cells
    Return:
        normalized Pandas DataFrame 
    """
    rank_mean = ref_exp.stack().groupby(ref_exp.rank(method='first').stack().astype(int)).mean()
    return df.rank(method='min').stack().astype(int).map(rank_mean).unstack()

cell_id = 0 # index for the cell type
N = 64 # number for cell types
tf_TPM_expression = pd.read_csv('EXAMPLE_TPM.csv',header=0,sep='\t',index_col=[0]) # TPM value of 711 TFs
tf_expression = quantile_norm(tf_TPM_expression,ref_tf_expression)
tf_expression = np.array(tf_expression.iloc[:,cell_id]) # quntile transformed TPM values

tf_TPM_expression.head()

#2) Prepare tf motif score
# scan the motif binding score using Homer tool
! findMotifsGenome.pl bins.bed hg38.fa ./ -find all_motif_rmdup.motif -p 32 -cache 320000 > homer_motifscore.txt

"""parse the motifscan results by homer tool.
    Args:
        bin_file: path to bed file that was used for motif scanning
        motifscan_file: path to the motifscan file from Homer tool (e.g. homer_motifscore.txt)
"""
tf_motifscore = get_motif_score(tf_list, motif2tf, isContain, bin_file, motifscan_file, save = True)

chrom, start, end = 'chr1', 768000, 896000
region = [chrom, start, end]
tf_motifscore = get_motifscore(region,tf_list,motif2tf) # motif binding status from homer tools

# TF feature
tf_feature = tf_motifscore * np.log(tf_expression + 1)
print(tf_feature.shape)

#3.2 Prepare input: genomic region with 128kbp length
chrom, start, end = 'chr1', 768000, 896000
region = [chrom, start, end]
bins = [chrom+':'+str(start+i*128)+'-'+str(start+(i+1)*128) for i in range(1000)]

#3.3 Model prediction
# bin level prediction
predict = model_predict(model,region[0],region[1],region[2],tf_feature)
pd_predict = pd.DataFrame(predict[0,:,:],index=bins,columns = ['DNase','CTCF','H3K27ac','H3K4me3','H3K36me3','H3K27me3','H3K9me3','H3K4me1'])
pd_predict
