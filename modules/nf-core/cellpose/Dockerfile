# GPU image for cellpose segmentation
# Base: PyTorch with CUDA 12.4 (consistent CUDA support)
FROM pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    libgl1 \
    libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip
RUN pip install --no-cache-dir --upgrade pip

# Install cellpose and its dependencies
RUN pip install --no-cache-dir \
    cellpose==4.0.8 \
    numpy \
    scipy \
    matplotlib \
    scikit-image \
    opencv-python-headless

# Set default shell
CMD ["/bin/bash"]
