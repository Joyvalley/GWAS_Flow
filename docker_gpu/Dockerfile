FROM ubuntu:18.04
LABEL maintainer="jan.freudenthal@uni-wuerzburg.de"

RUN apt update 
RUN apt install -y locales
RUN export DEBIAN_FRONTEND=noninteractive	
RUN locale-gen en_GB.UTF-8

RUN apt install -y vim git 
RUN apt install -y wget bzip2 sudo 
RUN apt install -y aptitude
RUN aptitude install -y nvidia-driver-418

RUN adduser --disabled-password --gecos '' ubuntu
RUN adduser ubuntu sudo
RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
USER ubuntu
WORKDIR /home/ubuntu/
RUN chmod a+rwx /home/ubuntu/


#### change anaconda to miniconda 
#RUN wget https://repo.continuum.io/miniconda/Miniconda3-4.3.14-Linux-x86_64.sh
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda*.sh -b -p /home/ubuntu/con
RUN rm Miniconda*
ENV PATH /home/ubuntu/con/bin:$PATH


# Updating Anaconda packages
RUN conda update conda
RUN conda install scipy 
RUN conda install numpy
RUN conda install pandas
RUN pip install limix
RUN conda install h5py
RUN conda install tensorflow-gpu

ADD gwas.py .
ADD herit.py .
ADD main.py .


WORKDIR /data
RUN which python
ENTRYPOINT ["python","-u","/home/ubuntu/gwas.py"]

# docker build  -t gwas_flow_gpu  docker

### Run docker container
## docker run --rm -u $UID:$GID -v $PWD:/data tf_image:latest  -x gwas_sample_data/KE_geno.csv -y gwas_sample_data/KE_pheno.csv -k gwas_sample_data/K.csv -m 2 -perm No

### Build singulartiy container from docker container locally 

## docker run -v /var/run/docker.sock:/var/run/docker.sock -v /home/joyvalley/PHD/Projects/gwastf:/output --privileged -t singularityware/docker2singularity tf_image:latest
