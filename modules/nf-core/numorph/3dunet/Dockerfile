## optimized nvidia base container for older python and tensorflow
FROM nvcr.io/nvidia/tensorflow:20.10-tf1-py3

ARG NUMORPH_3DUNET_VERSION="1.0.9"

LABEL org.opencontainers.image.authors=="Carolin Schwitalla <carolin.schwitalla@uni-tuebingen.de>"
LABEL org.opencontainers.image.description="Docker image for the NuMorph 3D-UNet to perform nuclei quantification on light-sheet microscopy data"
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.version="$NUMORPH_3DUNET_VERSION"



RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    libgl1-mesa-glx \
    libglib2.0-0 \
    libsm6 \
    libxext6 \
    libxrender-dev \
    && rm -rf /var/lib/apt/lists/*


RUN pip install --upgrade pip

# the tool needs this 
RUN pip install git+https://github.com/keras-team/keras-contrib.git

# Install the package
RUN pip install numorph-3dunet

# Create directory for models
RUN mkdir -p /models

CMD ["/bin/bash"]
