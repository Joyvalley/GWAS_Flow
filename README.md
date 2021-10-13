# GWAS_Flow

## Citing
`GWAS-Flow` was written and published in the hope that you might find it useful. If you do and use it for your research please cite the paper published alongside the software, which is currently publicly accessible on the BiorXiv preprint server. https://www.biorxiv.org/content/10.1101/783100v1  doi: 10.1101/783100 



## Introduction 

`GWAS_Flow` is an open source python based software provding a GPU-accelerated framework for performing genome-wide association studies (GWAS), published under the MIT-License. 
GWAS is a set  of  major algorithms in quantitative genetics to find associations between phenotypes and their respective genotypes.
With a broad range of applications ranging from plant breeding to medicine. 
In recent years the data sets used for those studies increased rapidly in size, and accordingly the time necessary to perform these on conventional CPU-powered machines increased exponentially.
Here we use TensorFlow a framework that is commonly used for machine learning applications to utilize graphical processing units (GPU) for GWAS. 

## Requirements

### Required Software
- [python](https://www.python.org/ "Python programming language") (v.3.9)
- [anaconda](https://www.anaconda.com/ "Anaconda virtual environments")
- [git](https://git-scm.com/)

### Required python packages
- [tensorflow](https://www.tensorflow.org/ "tensorflow") (v.2.6.0)
- [numpy](https://numpy.org/ "numerical python")
- [pandas](https://pandas.pydata.org/ "import and manipulate data frames")(v.1.0.5)
- [scipy](https://www.scipy.org/ "scientific python")
- [h5py](https://www.h5py.org/ "import and manipulated hdf files") 
- [matplotlib](https://matplotlib.org/ "plot library for python")

### Docker and Singularity
- [Docker](https://www.docker.com/) (v.19.03.1)
- [Singularity](https://singularity.lbl.gov/) (v.2.5.2)

## Installation 

### git and anaconda 
This has been tested on multiple linux systems with anconda versions > 4.7 

clone the repository directly with git 

```shell
git clone https://github.com/Joyvalley/GWAS_Flow
``` 

create an anaconda environment and install the necessary packages using the gwas_flow_env.yaml configuration file

```shell
###  optional: 
conda create -n gwas_flow python==3.9 pip
conda activate gwas_flow
### set up environment with pip 
pip install -r requirements.txt

```

### docker 
For the installation with docker the only required software is docker itself.

```shell 
git clone https://github.com/Joyvalley/GWAS_Flow.git 
cd GWAS_Flow
docker build  -t gwas_flow .
```

Then you can run GWAS_Flow using your user id and files in your current working directory like this:
```
docker run -u $UID:$GID -v $PWD:/data --rm gwas_flow -x gwas_sample_data/G_sample.csv -y gwas_sample_data/Y_sample.csv -k gwas_sample_data/K_sample.csv -o docker_out.csv
```


### singularity

```shell 
git clone https://github.com/Joyvalley/GWAS_Flow.git 

docker build  -t gwas_flow .

!! make sure to change /PATH/TO/FOLDER
docker run -v /var/run/docker.sock:/var/run/docker.sock -v /PATH/TO/FOLDER:/output --privileged -t singularityware/docker2singularity:1.11 gwas_flow:latest
change the name of e.g. gwas_flow_latest-2019-08-19-8c98f492dd54.img to gwas_flow_sing.img
```

## Execution with anaconda installation
### Input data 
GWAS_Flow is designed to work with several different input data formats. For all of them there is are sample data avaialble in the folder gwas_sample_data/
The minimal requirement is to provide a genotype and a phenotype file if no kinship matrix is provided a kinship matrix according to van Raden ist caluculated from the provided marker information. Depending on the size of the marker matrix this can take a while.

#### hdf5 input

```shell
python gwas.py -x gwas_sample_data/AT_geno.hdf5 -y gwas_sample_data/phenotype.csv -k gwas_sample_data/kinship_ibs_binary_mac5.h5py

```
#### csv input 

```shell
python gwas.py -x gwas_sample_data/G_sample.csv -y gwas_sample_data/Y_sample.csv -k gwas_sample_data/K_sample.csv

```
#### plink input
To use PLINK data format add a bed bim and fam file with the same prefix to the folder. You can tell GWAS-Flow to use those files by using prefix.plink as the option for the genotype file

```shell
python gwas.py -x gwas_sample_data/my_plink.plink -y gwas_sample_data/pheno2.csv -k gwas_sample_data/kinship_ibs_binary_mac5.h5py
```


Flgas and options are 	
```shell
-x , --genotype : file containing marker information in csv or hdf5 format of size
-y , --phenotype : file container phenotype information in csv format
-k , --kinship : file containing kinship matrix of size k X k in csv or hdf5 format
-m : name of column to be used in phenotype file. Default m='phenotype_value' 
--cof: file with cofactor information (only one co-factor as of now)
-a , --mac_min : integer specifying the minimum minor allele count necessary for a marker to be included. Default a = 1
-bs, --batch-size : integer specifying the number of markers processed at once. Default -bs 500000
-p , --perm : perform n permutations
--out_perm : output individual resulst of the permuation. Default False, enable with arbitary string (e.g. --out_perm yo)
--plot : create manhattanplot 
-o , --out : name of output file. Default -o results.csv  
-h , --help : prints help and command line options
```

use `python gwas.py -h` to see the command line options

### Execution with docker and singularity 

Execute the docker container with the sample data
```shell
docker run --rm -u $UID:$GID -v $PWD:/data gwas_flow:latest  -x gwas_sample_data/AT_geno.hdf5 -y gwas_sample_data/phenotype.csv -k gwas_sample_data/kinship_ibs_binary_mac5.h5py
```

On Windows you can use something like this after activating the file sharing for the drive the repo is stored on:

```cmd
cd c:\PATH\TO\REPO\GWAS_Flow
docker run -v c:/PATH/TO/REPO/GWAS_Flow:/data gwas_flow:latest -x gwas_sample_data/AT_geno.hdf5 -y gwas_sample_data/phenotype.csv -k gwas_sample_data/kinship_ibs_binary_mac5.h5py
```

!! The GPU versions of docker and singularity are still under development and might or might not work properly with your setup.
To run the GWAS-Flow on GPUs as of now we recommand the usage of anaconda environments

Execute the singularity image with the sample data
```shell
singularity run  gwas_flow_sing.img -x gwas_sample_data/AT_geno.hdf5 -y gwas_sample_data/phenotype.csv -k gwas_sample_data/kinship_ibs_binary_mac5.h5py
```

### further options 
#### Co-factor
So far GWAS-Flow is capable of using on co-factor the co-factor is added to the analysis with the flag `--cof FILENAME` 
e.g 
```shell
 python gwas.py -x gwas_sample_data/G_sample.csv -y gwas_sample_data/Y_sample.csv -k gwas_sample_data/K_sample.csv --cof gwas_sample_data/cof.csv 
```


#### Permutation
add the flag `--perm 100` to calculate a significance threshold based on 100 permutations. Change 100 to any integer larger 2 to perform n permutations

#### Manhattan plot 
By default there is no plot generated if you add `--plot True` a manhattan plot is generated

![manhattan](https://user-images.githubusercontent.com/26280192/71103427-6a57e400-21ba-11ea-8eab-aa40abcc46ec.png)

The dash-dotted line is the bonferroni threshold of significance and the dashed line the permutation based threshold
The latter is only calculated if the flag `--perm n` was used with n > 2.



## Performance Benchmarking and Recommendations
![time_plot](https://user-images.githubusercontent.com/26280192/71102215-61661300-21b8-11ea-90cc-7d4690e645be.png)

The image displays the average time of 10 runs with 10000 markers each and varying number of phenotypes for `GWAS_Flow` on GPU and CPUs and a standard [R-Script](https://github.com/arthurkorte/GWAS "GWAS") for GWAS.
The computational time growths exponentially with increasing number of phenotypes.
With lower numbers of phenotypes (< 800), the CPU version is faster than the GPU Version.
This gets more and more lopsided the more phenotypes are included. 
All calculations have been performed on 16 i9 vCPUS and a NVIDIA Tesla P100 graphic card.


## Unit tests 

The unit tests can be run one the console with:

```shell 
python -m unittest tests/test.py
```

All the necassary test data is stored in test_data
