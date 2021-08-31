FROM ubuntu:20.04

LABEL maintainer="jan.freudenthal@uni-wuerzburg.de"

ENV DEBIAN_FRONTEND noninteractive

RUN apt update 

RUN apt install --yes && \
    apt install --yes \
    	locales
	
RUN locale-gen en_US.UTF-8

RUN apt install -y emacs git 
RUN apt install -y wget bzip2 sudo 
    		 
RUN adduser --disabled-password --gecos '' ubuntu
RUN adduser ubuntu sudo
RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
USER ubuntu
WORKDIR /home/ubuntu/
RUN chmod a+rwx /home/ubuntu/

#### change anaconda to miniconda 

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b 
RUN rm Miniconda3-latest-Linux-x86_64.sh 
ENV PATH /home/ubuntu/miniconda3/bin:$PATH


# Updating Anaconda packages
RUN conda update conda

# Install python packages
ADD requirements.txt .
RUN conda install python=3.7 pip
RUN pip install -r requirements.txt

# Add scripts
ADD gwas.py .
ADD gwas_flow ./gwas_flow

WORKDIR /data
ENTRYPOINT ["python","-u","/home/ubuntu/gwas.py"]

# docker build  -t gwas_flow  docker

### Run docker container
## docker run --rm -u $UID:$GID -v $PWD:/data gwas_flow:latest  -x gwas_sample_data/AT_geno.hdf5 -y gwas_sample_data/phenotype.csv -k gwas_sample_data/kinship_ibs_binary_mac5.h5py

### Build singulartiy container from docker container locally 

## docker run -v /var/run/docker.sock:/var/run/docker.sock -v ../GWAS_Flow:/output --privileged -t singularityware/docker2singularity:1.11 tf_image:latest
