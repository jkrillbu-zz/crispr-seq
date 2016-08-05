FROM python:2

RUN mkdir -p /usr/src/app
WORKDIR /usr/src/app

RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y tcl tcl-dev tk-dev tk 

RUN curl -LO https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2
RUN bunzip2 samtools-0.1.19.tar.bz2 && tar xf samtools-0.1.19.tar && mv samtools-0.1.19 samtools && cd samtools && make
RUN rm samtools-0.1.19.tar
ENV PATH=/usr/src/app/samtools:${PATH}

COPY requirements.txt /usr/src/app/
RUN pip install --no-cache-dir -r requirements.txt

RUN apt-get install -y gfortran
RUN curl -LO https://cran.r-project.org/src/base/R-3/R-3.2.1.tar.gz
RUN tar xf R-3.2.1.tar.gz && cd R-3.2.1 && ./configure --prefix=/usr/local --enable-R-shlib --with-tcltk && make && make install

COPY install.R /usr/src/app/
RUN Rscript install.R

RUN wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.10.tar.bz2
RUN bunzip2 bwa-0.7.10.tar.bz2 && tar xf bwa-0.7.10.tar && mv bwa-0.7.10 bwa && cd bwa && make
RUN rm bwa-0.7.10.tar
ENV PATH=/usr/src/app/bwa:${PATH}

COPY src/  /usr/src/app/src/