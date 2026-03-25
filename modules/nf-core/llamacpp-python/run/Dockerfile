FROM nvidia/cuda:12.4.1-devel-ubuntu22.04

RUN apt-get update && apt-get install -y python3 python3-pip
RUN pip3 install llama-cpp-python  \
    --extra-index-url https://abetlen.github.io/llama-cpp-python/whl/cu124
