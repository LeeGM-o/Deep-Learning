# EpiGePT Project
![EpiGePT overview](https://github.com/user-attachments/assets/6b2cf334-89a6-4478-9276-1e7e88685b2c)
EpiGePT is a pretrained transformer-based model for predicting context-specific epigenomic signals and chromatin contacts by taking long DNA sequence and transcription factor profile as inputs. By taking the context-specific activities of transcription factors (TFs) and 3D genome interactions into consideration, EpiGePT offers wider applicability and deeper biological insights than models trained on DNA sequence only. Given the expression profile of hundreds of TFs from a cellular context, EpiGePT is able to predict the genome-wide chromatin states.

## Component
#### 1) Sequence module   
* size: (128000(128 kb), 4)   
* one-hot matrix(A = [0, 0, 0, 1], C = [0, 1, 0, 0], G = [0, 0, 1, 0], T = [0, 0, 0, 1])
#### 2) TF module
#### 3) Transformer module
#### 4) Multi-task module

## Requirements
* Pytorch-lightning==1.2.10
* Pytorch==1.7.1
* Python==3.6.9
* Pyfasta==0.5.2
