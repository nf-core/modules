FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

# Install necessary tools
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    && rm -rf /var/lib/apt/lists/*


# Download the zip file
RUN wget "https://github.com/Visceral-Project/EvaluateSegmentation/raw/refs/heads/master/builds/Ubuntu/EvaluateSegmentation-2020.08.28-Ubuntu.zip"

# Unzip the file
RUN unzip EvaluateSegmentation-2020.08.28-Ubuntu.zip -d bin/

RUN chmod +x bin/EvaluateSegmentation

# Optional: delete the zip file to save space
RUN rm EvaluateSegmentation-2020.08.28-Ubuntu.zip
