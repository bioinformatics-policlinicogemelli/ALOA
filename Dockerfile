# Filename: PythonVaran
FROM ubuntu
USER root
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update &&\
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends curl gcc g++ gnupg unixodbc-dev openssl git &&\
    apt-get install -y software-properties-common ca-certificates &&\
    apt-get install -y build-essential libtiff-dev pandoc libfontconfig1-dev libgdbm-dev libharfbuzz-dev libfribidi-dev zlib1g-dev libncurses5-dev libgdbm-dev libssl-dev libreadline-dev libffi-dev wget libbz2-dev libsqlite3-dev libcurl4-openssl-dev libxml2-dev  && \
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
COPY . /
WORKDIR /
RUN python3 -m pip install -r requirements.txt 
RUN while read -r package; do Rscript -e "install.packages('"$package"')"; done < "req.txt"
RUN R -e "install.packages('remotes')"
RUN Rscript -e "remotes::install_github('r-lib/textshaping')"
RUN Rscript -e "install.packages('pak')"
#RUN Rscript -e "pak::pak('r-lib/ragg')"
RUN Rscript -e "pak::pak('tidyverse')"
RUN Rscript -e "pak::pak('spatstat')"
RUN Rscript -e 'download.file(url = "https://cran.r-project.org/src/contrib/Archive/spatstat.core/spatstat.core_2.4-4.tar.gz", destfile =  "spatstat.core_2.4-4.tar.gz")' 
RUN Rscript -e "install.packages(pkgs='spatstat.core_2.4-4.tar.gz', type='source', repos=NULL)" 

RUN  Rscript -e "remotes::install_github('akoyabio/phenoptr')"
RUN  Rscript -e "pak::pkg_install('akoyabio/phenoptrReports')"


CMD [ "bash"]
