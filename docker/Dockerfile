FROM ubuntu:18.04


LABEL maintainer="jan.freudenthal@uni-wuerzburg.de"

RUN apt update 

RUN apt install --yes && \
    apt install --yes \
    	locales
	
RUN locale-gen en_GB.UTF-8


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
RUN conda install tensorflow==1.14
RUN conda install scipy
RUN conda install numpy
RUN conda install pandas
RUN pip install limix
RUN conda install h5py


ADD gwas.py .
ADD herit.py .
ADD main.py .

WORKDIR /data
ENTRYPOINT ["python","-u","/home/ubuntu/gwas.py"]

# docker build  -t tf_image  docker

### Run docker container
## docker run --rm -u $UID:$GID -v $PWD:/data gwas_flow:latest  -x gwas_sample_data/AT_geno.hdf5 -y gwas_sample_data/phenotype.csv -k gwas_sample_data/kinship_ibs_binary_mac5.h5py

### Build singulartiy container from docker container locally 

## docker run -v /var/run/docker.sock:/var/run/docker.sock -v /home/joyvalley/PHD/Projects/gwastf:/output --privileged -t singularityware/docker2singularity:1.11 tf_image:latest
