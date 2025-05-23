FROM python:3.9 as base

ARG SCSPLIT_VERSION="1.0.8"

RUN apt-get update && apt-get install -y python3-pip git

RUN pip install "setuptools==57.5.0" wheel numpy "pandas<2.0" pysam pyvcf3 "scikit-learn==1.1.3" scipy
RUN wget https://github.com/jon-xu/scSplit/archive/refs/tags/v${SCSPLIT_VERSION}.tar.gz -O /opt/scsplit.gz && \
    tar xvzf /opt/scsplit.gz -C /opt && \
    rm /opt/scsplit.gz

ENV PATH="${PATH}:/opt/scSplit-${SCSPLIT_VERSION}"

CMD [ "bash" ]
