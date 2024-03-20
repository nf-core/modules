FROM python:3.7 as base

ARG SCSPLIT_VERSION="1.0.8"

RUN apt-get update && apt-get install -y python3-pip git

RUN pip install setuptools numpy pandas pysam pyvcf scikit-learn scipy statistics 
RUN wget https://github.com/jon-xu/scSplit/archive/refs/tags/v${SCSPLIT_VERSION}.tar.gz -O /opt/scsplit.gz && \
    tar xvzf /opt/scsplit.gz -C /opt && \
    rm /opt/scsplit.gz

ENV PATH="${PATH}:/opt/scSplit-${SCSPLIT_VERSION}"

CMD [ "bash" ]