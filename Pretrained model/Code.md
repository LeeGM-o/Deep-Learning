### Data downloading

We provided `Download_raw_data.sh` for download RNA-seq data (.tsv), DNase-seq data (.bam) and ChIP-seq data (.bam) from the ENCODE project
We pre-defined cell type ID from 1-37. After downloading the meta data from ENCODE website (`head -n 1 files.txt|xargs -L 1 curl -O -L`), one can run the following script:

```python
bash Download_raw_data.bash  -c <CELL_ID> -r -c -d
-c  CELLID: pre-defined cell ID (from 1 to 55)
-r  download RNA-seq data (.tsv)
-d  download chromatin accessible readscount from DNase-seq data (.bam)
-c  download readscount from ChIP-seq data (.bam)
```
one can also run ```bash Download_raw_data.bash  -h``` to show the script instructions. Note that `.bam` files downloading may take time. After downloading the raw data, the raw data folder will be organized by `cell-assay-experiment-file` order. Note that each experiment may contain multiple replicates. See an example of the folder tree:

```
data/
    |-- raw_data/
    |   |-- 1/
    |   |   |-- dseq/
    |   |   |   |-- ENCSR000EIE/
    |   |   |   |   |-- ENCFF953HEA.bed.gz
    |   |   |   |   |-- ENCFF983PML.bam
    |   |   |   |-- ENCSR000ELW/
    |   |   |   |   |...
    |   |   |-- rseq/
    |   |   |   |-- ENCSR000BXY/
    |   |   |   |   |-- ENCFF110IED.tsv
    |   |   |   |   |-- ENCFF219FVQ.tsv
    |   |   |   |-- ENCSR000BYH/
    |   |   |   |   |...
```

### Data preprocessing

After downloading the raw DNase-seq, RNA-seq and ChIP-seq data, you can align BAM files to the reference sequence to obtain read counts. Users can use the following command for alignment:

```shell
bash runme.sh
```

Then we merge multiple replicate of RNA-seq data by taking the average expression of each gene across replicates in a cell type. As for DNase-seq data and ChIP-seq data, we only keep bins that appear in more than half of the replicates with respect to a cell type. One can run the following scripts to merge relicates of both DNase-seq, ChIP-seq and RNA-seq data, and generate TF gene expression matrix (N x C where N is the number of TFs and C is the number of cell types).

Additionally, users can also obtain the motif binding score of motifscan results by homer tool.

```shell
python preprocess.py
```

As for the motif binding score, we conducted scanning of 711 transcription factors across the entire genome and integrated them into the [EpiGePT-online](http://health.tsinghua.edu.cn/epigept/) platform, enabling users to directly perform predictions through the web interface. The file of genomic regions and samll demo data for motif score and target data can be found [here](data).


### Model training


The main python script `train.py` is used for implementing EpiGePT for predicting 8 epigenomic profiles. Model architecture for EpiGePT can be find in `./model/EpiGePT.py`. Data loader or data sampler can be find in `./model/dataset.py`.

One can run the following commond to train a EpiGePT model, the preprocessed data of TF expression value can be downloaded from the `Supplementary Materials` of EpiGePT, the motif score file and target data for training used in the paper can be download from the Download page of [EpiGePT-online](http://health.tsinghua.edu.cn/epigept/download.php).

```shell
CUDA_VISIBLE_DEVICES=0 python train.py --train True  --num_train_region 10000 
[train]  --  whethre use train mode
[num_train_region]  --  number of genomic regions used for training
[cell_idxs_path] -- indexes of cell types used for training
```




 Next, users can also train a EpiGePT model specially for the prediction of the chromatin accessibility (DNase-seq) by running the following commond.


```shell
CUDA_VISIBLE_DEVICES=0 python train_DNase.py --train True --cell_idxs_path data
[train]  --  whethre use train mode
[cell_idxs_path] -- indexes of cell types used for training
```

Train cell type index file and test cell type index file can be found as [data/train_idxs_5_fold.npy](./data/train_idxs_5_fold.npy) and [data/test_idxs_5_fold.npy](./data/test_idxs_5_fold.npy). When training the DNase-specific pre-trained model with a batch size of 128, it takes approximately 2-3 hours to complete one epoch of training.


### Model testing

For model prediction, we provide three different approaches. Given that EpiGePT allows prediction of epigenomic signal values in any specified genomic region and cell type, we employed three testing approaches to evaluate model performance within the collected 28 cell types: `cross-cell-type`, `cross-region`, and `cross-both`. one can run the following command to conduct cross cell type prediction

```shell
CUDA_VISIBLE_DEVICES=0 python train.py --train False --pred_method 0 --num_train_region 10000 --cell_idxs_path train_cell_type_idxs.npy --pretrained_model_path checkpoint/best_model.ckpt
[train]  --  whethre use train mode
[pred_method] -- Method for testing the model, 0 for cross cell type prediction, 1 for cross genomic region prediction and 2 for cross both prediction
[num_train_region]  --  number of genomic regions used for training
[cell_idxs_path] -- indexes of cell types used for training
[pretrained_model_path] -- path for the pretrained EpiGePT model
```

and can run the following command to conduct cross genomic region prediction
```shell
CUDA_VISIBLE_DEVICES=0 python train.py --train False --pred_method 1 --num_train_region 10000 --cell_idxs_path test_cell_type_idxs.npy --pretrained_model_path checkpoint/best_model.ckpt
```


Similarly, we provide separate scripts for cross-cell type prediction of DNase, allowing users to predict DNase-seq signals across the entire genome in specified test cell types among the 129 available cell types using the following script.

```shell
CUDA_VISIBLE_DEVICES=0 python train_DNase.py --train False --cell_idxs_path test_cell_type_idxs.npy --pretrained_model_path checkpoint/best_model.ckpt
[train]  --  whethre use train mode
[cell_idxs_path] -- indexes of cell types used for training
[pretrained_model_path] -- path for the pretrained EpiGePT model
```

### Pretrain Models

We provide various of pretrain models for a quick implementation of EpiGePT. First, one needs to download the pretrain models from the [EpiGePT-online](http://health.tsinghua.edu.cn/epigept/download.php) website. Then uncompress it under `EpiGePT` folder. For the above models that use `predict.py` for model prediction. For an example, one can run 

```python
python predict.py  --pretrained_model_path PATH/model.ckpt
```
