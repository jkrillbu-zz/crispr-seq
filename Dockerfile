FROM python:2

RUN mkdir -p /usr/src/app
WORKDIR /usr/src/app

RUN apt-get update
RUN apt-get upgrade -y

RUN curl -LO https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2
RUN bunzip2 samtools-0.1.19.tar.bz2 && tar xf samtools-0.1.19.tar && mv samtools-0.1.19 samtools && cd samtools && make
RUN rm samtools-0.1.19.tar
ENV PATH=/usr/src/app/samtools:${PATH}

COPY requirements.txt /usr/src/app/
RUN pip install --no-cache-dir -r requirements.txt

COPY src/  /usr/src/app/src/