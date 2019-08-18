# GWAS_Flow

## Introduction 

GWAS_Flow is an open source python based software provding a GPU-accelerated framework for performing genome association studies, published under the MIT-License. 
GWAS is a set  of  major algorithms in quantitative genetics  to find associations between  phenotypes and their respective genotypes. With a broad range of applications ranging from plant breeding to medicine. 
In recent years the data sets used for those studies increased rapidly in size, and accordingly the time necessary too perform these on conventional CPU-powered machines increased exponentially. Here we used TensorFlow a framework that is commonly used for machine learning applications to utilize graphical processing units (GPU) for GWAS. 

## Requirements


## Required Software
- [python](https://www.python.org/ "Python programming language") (v.3.7.3)

## Required python packages
- [tensorflow](https://www.tensorflow.org/ "tensorflow") (v.1.14.0)
- [numpy](https://numpy.org/ "numerical python") (v.1.16.4)
- [pandas](https://pandas.pydata.org/ "import and manipulate data frames")(v.24.2)
- [scipy](https://www.scipy.org/ "scientific python") (v.1.3.0)
- [h5py](https://www.h5py.org/ "import and manipulated hdf files") (v.2.9.0)


## Docker and Singularity
- [Docker](https://www.docker.com/) (v.19.03.1)
- [Singularity](https://singularity.lbl.gov/) (v.2.5.2)

## Installation 

### git and anaconda 

clone the directory directly from git 

```shell
git clone https://github.com/Joyvalley/GWAS_Flow
``` 

create and anaconda environment and install the necassary packages
```shell
conda create -n gwas_flow python=3.7.3
conda activate gwas_flow
conda install tensorflow==1.14 # conda install tensorflow-gpu==1.4 for gpu usage
conda install scipy pandas numpy h5py
pip install limix

```

### docker 

```shell 
git clone https://github.com/Joyvalley/GWAS_Flow.git 

docker build  -t gwas_flow  docker

```

### singularity

```shell 
git clone https://github.com/Joyvalley/gwas_tf.git 

docker build  -t gwas_flow docker

docker run -v /var/run/docker.sock:/var/run/docker.sock -v /path/to/git/gwastf:/output --privileged -t singularityware/docker2singularity:1.11 ggwas:latest
```
## Execution 
To run the gwas with default settings and the sample data provided in gwas_sample_data/ 
Make sure to have all the required packages installed 

```shell
python gwas.py -x gwas_sample_data/AT_geno.hdf5 -y gwas_sample_data/phenotype.csv -k gwas_sample_data/kinship_ibs_binary_mac5.h5py

```

Flgas and options are 	
```shell
-x , --genotype : file containing marker information in csv or hdf5 format of size
-y , --phenotype : file container phenotype information in csv format
-k , --kinship : file containing kinship matrix of size k X k in csv or hdf5 format
-m : integer specifying the column of phenotype file to use. Default -m 0
-a , --mac_min : integer specifying the minimum minor allele count necessary for a marker to be included. Default a = 1
-bs, --batch-size : integer specifying the number of markers processed at once. Default -bs 500000
-p , --perm : 
-o , --out : name of output file. Default -o results.csv  
-h , --help : prints help and command line options

```

use `python gwas.py -h` to see the command line options


## Performance Benchmarking and Recommendations
![benchmark](https://user-images.githubusercontent.com/26280192/63228473-bf2b2400-c1f3-11e9-86c2-081ca86127bd.png)

The image displays the average time of 10 runs with 10000 markers each and varying number of phenotypes for GWAS_Flow on GPU and CPUs and a standard [R-Script](https://github.com/arthurkorte/GWAS "GWAS") for GWAS.  The computational time growths exponentially with increasing number of phenotypes. With lower numbers of phenotypes (< 800), the CPU version is faster than the GPU Version. This gets more and more lopsided the more phenotypes are included. 
All calculations have been performed on 16 i9 vCPUS and a NVIDIA Tesla P100 graphic card.
