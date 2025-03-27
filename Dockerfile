# Use a base image with the desired Linux distribution (e.g., Ubuntu)
FROM ubuntu:latest

# Set non-interactive mode for package installations
ENV DEBIAN_FRONTEND=noninteractive

# Install necessary build tools
RUN apt-get update
RUN apt-get install -y --no-install-recommends make
RUN apt-get install -y --no-install-recommends gcc
RUN apt-get install -y --no-install-recommends wget
RUN apt-get install -y --no-install-recommends curl
RUN apt-get install -y --no-install-recommends ca-certificates
RUN apt-get install -y --no-install-recommends git
RUN apt-get install -y --no-install-recommends libz-dev
RUN apt-get install -y --no-install-recommends libncurses-dev
RUN apt-get install -y --no-install-recommends liblzma-dev
RUN apt-get install -y --no-install-recommends libbz2-dev
RUN apt-get install -y --no-install-recommends libcurl4-openssl-dev
RUN apt-get install -y --no-install-recommends libffi-dev
RUN apt-get install -y --no-install-recommends unzip
RUN apt-get install -y --no-install-recommends parallel
#RUN apt-get install -y --no-install-recommends openjdk-17-jdk
RUN apt update
RUN apt install -y --no-install-recommends openjdk-17-jdk
RUN apt-get install -y --no-install-recommends pkg-config
RUN apt-get install -y --no-install-recommends python3
RUN apt-get install -y --no-install-recommends python3-pip
RUN apt-get install -y --no-install-recommends bwa
RUN apt-get install -y --no-install-recommends samtools
RUN apt-get install -y --no-install-recommends hisat2
RUN apt-get install -y --no-install-recommends cpanminus

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN pip install whatshap --break-system-packages \
    -i http://mirrors.aliyun.com/pypi/simple/ \
    --trusted-host mirrors.aliyun.com

RUN wget https://githubfast.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip && unzip gatk-4.4.0.0.zip

RUN cpanm Cwd
RUN cpanm List::Util
RUN cpanm Math::Random
RUN cpanm Getopt::Std
RUN cpanm Parallel::ForkManager
#RUN apt-get install -y r-base r-base-dev

RUN git clone https://githubfast.com/zwycooky/hapBSA.git
RUN chmod +x /hapBSA/scripts/*

# Add the R library path to the PATH variable so that it is available to your applications.
# ENV PATH=$PATH:$R_HOME/lib/R/library

# Add the path of scripts to the PATH variable
ENV PATH=$PATH:/hapBSA/scripts/

# Start a new shell session when the container runs
CMD ["/bin/bash"]

