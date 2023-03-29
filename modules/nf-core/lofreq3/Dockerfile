FROM ubuntu:latest

LABEL \
    author="Vojislav Lukic" \
    description="Image for LoFreq3 with bwa, samtools, htslib and nim installed" \
    maintainer="vojislav.lukic@apisassay.com"

RUN apt update && apt install -y \
samtools \
bzip2 \
xz-utils \
wget \
curl \
autoconf \
git \
automake \
make \
gcc \
perl \
zlib1g-dev \
libbz2-dev \
liblzma-dev \
libcurl4-gnutls-dev \
libssl-dev

WORKDIR /opt
RUN git clone https://github.com/lh3/bwa.git
WORKDIR /opt/bwa
RUN make
ENV PATH=/opt/bwa/bin:/opt/bwa:$PATH

WORKDIR /opt
RUN wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2 && tar xf htslib-1.17.tar.bz2
WORKDIR /opt/htslib-1.17
RUN make && make install

ENV PATH=/opt/htslib-1.17/bin:$PATH
ENV LD_LIBRARY_PATH="/opt/htslib-1.17"

WORKDIR /opt
RUN wget https://nim-lang.org/download/nim-1.6.12-linux_x64.tar.xz && tar xf nim-1.6.12-linux_x64.tar.xz
WORKDIR /opt/nim-1.6.12
RUN sh install.sh /usr/local/bin
RUN bin/nim c koch
RUN ./koch boot -d:release
RUN ./koch tools
ENV PATH=/opt/nim-1.6.12/bin:$PATH

WORKDIR /opt
RUN git clone https://github.com/andreas-wilm/lofreq3.git
WORKDIR /opt/lofreq3
RUN nimble build -y

ENV PATH=/opt/lofreq3:$PATH
WORKDIR /opt
