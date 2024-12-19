# Use a base image with the desired Linux distribution (e.g., Ubuntu)
FROM ubuntu:latest

# Set non-interactive mode for package installations
ENV DEBIAN_FRONTEND=noninteractive

# Install necessary build tools
RUN apt-get update && apt-get install -y build-essential wget nano less curl  && apt-get install -y wget

RUN apt-get install -y make && apt install -y gcc && apt-get install -y libz-dev && apt install -y bzip2 && apt-get install -y git

RUN apt-get install -y findutils

RUN apt-get install -y libncurses-dev && apt-get install -y liblzma-dev && apt install -y libbz2-dev && apt-get install -y parallel

RUN apt-get install -y libcurl4-openssl-dev && apt install -y unzip && apt-get install -y openjdk-17-jdk

RUN apt-get install libffi-dev

RUN wget https://www.python.org/ftp/python/3.12.8/Python-3.12.8.tgz && tar zxf Python-3.12.8.tgz && cd Python-3.12.8 && ./configure --prefix=`pwd` && make && make install
RUN mv /Python-3.12.8/bin/python3 /Python-3.12.8/bin/python
RUN mv /Python-3.12.8/bin/pip3 /Python-3.12.8/bin/pip

# ADD python to ENV #
ENV PATH=$PATH:/Python-3.12.8/bin/

RUN pip install whatshap -i http://mirrors.aliyun.com/pypi/simple/ --trusted-host mirrors.aliyun.com 

RUN git clone https://gh.llkk.cc/https://github.com/lh3/bwa.git && cd bwa && make

RUN wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download -O ->> hisat2-2.2.1-Linux_x86_64.zip && unzip hisat2-2.2.1-Linux_x86_64.zip

RUN wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2  && bunzip2 samtools-1.17.tar.bz2 && tar -xf samtools-1.17.tar && cd samtools-1.17 && sh configure --prefix=`pwd` && make && make install

RUN wget https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip && unzip gatk-4.4.0.0.zip

RUN apt-get install -y cpanminus

RUN cpanm Cwd
RUN cpanm List::Util
RUN cpanm Getopt::Std
RUN cpanm Parallel::ForkManager
#RUN apt-get install -y r-base r-base-dev

RUN git clone https://gh.llkk.cc/https://github.com/zwycooky/hapBSA.git
RUN chmod +x /hapBSA/scripts/*

# Add the R library path to the PATH variable so that it is available to your applications.
# ENV PATH=$PATH:$R_HOME/lib/R/library

# Add the path of scripts to the PATH variable
ENV PATH=$PATH:/hapBSA/scripts/:/samtools-1.17/:/gatk-4.4.0.0/:/bwa/:/hisat2-2.2.1/

# Start a new shell session when the container runs
CMD ["/bin/bash"]

