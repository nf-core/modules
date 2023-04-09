FROM tensorflow/tensorflow:1.15.0-py3

RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata
RUN apt-get install -y python3-opencv
RUN apt-get install -y libtiff5-dev git

RUN pip install cython scikit-image==0.14.2 matplotlib tifffile==2020.2.16 scipy==1.1.0 opencv-python==4.3.0.36

RUN pip install git+https://github.com/FZJ-INM1-BDA/pytiff.git@0701f28e5862c26024e8daa34201005b16db4c8f

COPY . /app
