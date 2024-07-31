# EpiGePT Project
![EpiGePT overview](https://github.com/user-attachments/assets/6b2cf334-89a6-4478-9276-1e7e88685b2c)
EpiGePT is a pretrained transformer-based model for predicting context-specific epigenomic signals and chromatin contacts by taking long DNA sequence and transcription factor profile as inputs. By taking the context-specific activities of transcription factors (TFs) and 3D genome interactions into consideration, EpiGePT offers wider applicability and deeper biological insights than models trained on DNA sequence only. Given the expression profile of hundreds of TFs from a cellular context, EpiGePT is able to predict the genome-wide chromatin states.

## Component
#### 1) Sequence module   
* size: (128000(128 kb), 4)   
* one-hot matrix(A = [0, 0, 0, 1], C = [0, 1, 0, 0], G = [0, 0, 1, 0], T = [0, 0, 0, 1])
#### 2) TF module
* HOCOMOCO database: Potential binding sites for a set of 711 human transcription factors were scanned in the input bin using the Homer tool. We combined two vectors (motif features & expression features) made as follows and concatenated the results to the output of the sequence module.
> 1. Selected the maximum score of reported binding status for each transcription factor to obtain a vector of 711 dimensions as the motif feature for each DNA bin.
> 2. For gene expression, focused on log-transformed TPM values of the 711 transcription factors and obtained a vector of 711 dimensions after quantile normalization as the expression feature.
#### 3) Transformer module
* The input word embedding (X) of the transformer encoder = (Sequence length, embedding dim) -> (1000, 968).
> Input genomic bin sequence has a length of 1000   
> Each genomic bin has an embedded representation that combines the sequence information with cell-type-specific features with dimension of 968.
* For position embedding, we employed absolute position embedding to represent the positional information of the 1000 genomic bins in the input 128kbp DNA sequence, with dimensions of (1000, 968).
* Each Transformer encoder includes a multi-head self-attention mechanism and a feed-forward neural network. For self-attention in each head, the calculation is based on the matrix operation.
![attention function](https://github.com/user-attachments/assets/2f2472cb-bd16-423e-845b-d19b97bfd415)
> 
#### 4) Multi-task prediction module

## Requirements
* Pytorch-lightning==1.2.10
* Pytorch==1.7.1
* Python==3.6.9
* Pyfasta==0.5.2
