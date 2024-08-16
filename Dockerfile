# Filename: PythonVaran
FROM ubuntu
USER root
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends curl nano gcc g++ gnupg unixodbc-dev openssl git \
    software-properties-common ca-certificates build-essential libtiff-dev pandoc \
    libfontconfig1-dev libgdbm-dev libharfbuzz-dev libfribidi-dev zlib1g-dev \
    libncurses5-dev libgdbm-dev libssl-dev libreadline-dev libffi-dev wget libbz2-dev \
    libsqlite3-dev libcurl4-openssl-dev libxml2-dev libfftw3-dev libfftw3-doc && \
    update-ca-certificates && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir /python && cd /python && \
    wget https://www.python.org/ftp/python/3.11.1/Python-3.11.1.tgz && \
    tar -zxvf Python-3.11.1.tgz && \
    cd Python-3.11.1 && \
    ls -lhR && \
    ./configure --enable-optimizations && \
    make install && \
    rm -rf /python

RUN apt-get update && apt-get install r-base -y
COPY requirements.txt /requirements.txt
COPY installation_rpackages.R /installation_rpackages.R
COPY req.txt /req.txt
WORKDIR /
RUN python3 -m pip install -r requirements.txt 
RUN Rscript installation_rpackages.R req.txt

COPY . /

CMD [ "bash"]