FROM python:3.7 as base

RUN apt-get update && apt-get install -y python3-pip git
RUN pip install setuptools numpy pandas pysam pyvcf scikit-learn scipy statistics
RUN git clone https://github.com/jon-xu/scSplit

CMD [ "bash" ]
