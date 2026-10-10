# Dockerfile to create container with Cell Ranger
# Push to quay.io/nf-core/universc

FROM continuumio/miniconda3:latest

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Set different version
ENV CELLRANGER_VER=10.1.0
ENV UNIVERSC_VER=1.2.7
ENV MULTIQC_VER=1.35
ENV FASTQ_PAIR_VER=1.0
ENV SAMTOOLS_VER=1.16.1
ENV STAR_VER=2.5.1b

# Install system dependencies
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        gzip \
        pigz \
        rename \
        cpio \
        procps \
        wget \
    && rm -rf /var/lib/apt/lists/*

ENV CONDA_ENV=universc
ENV CELLRANGER_HOME=/opt/cellranger-${CELLRANGER_VER}
ENV UNIVERSC_HOME=/opt/universc

# Install cellranger
WORKDIR /opt
RUN wget -O cellranger-${CELLRANGER_VER}.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-${CELLRANGER_VER}.tar.gz?Expires=1791491800&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=l6e~Yieufc5adjybkCi72fMAS6fswTsGLQ2xFnuV9KCJk-SMiSuPYNOEYP4Cb5YiMDxktMxKH5rtVqpXtVTCPFbpOBbCxViTRDneuaX9g7IZ0UxckG4vl6Aoqc6E0HNVxlr02SB6ZkpPH-C0GwysFQyOeW8UV1NI6wl7vz~oJv6KlesfYtMzcGULwRTk7fYgxmnem5C7qHjbh4L2UVHkGdDe7kC0eElESHY0f9sw4cZCVDVLuLOTGcpzjMOmYDXsDqSyLy67hxiyReO~6lxlvyAkDagUMbmACLFlgUHGIbShvE8UWYxyWFG0DdMaP072SAcrYg38MZG25QXgBUkykg__" \
    && tar -xzf "cellranger-${CELLRANGER_VER}.tar.gz" \
    && test -x "${CELLRANGER_HOME}/cellranger" \
    && rm -f "cellranger-${CELLRANGER_VER}.tar.gz"

ENV PATH="${CELLRANGER_HOME}:${CELLRANGER_HOME}/bin:${UNIVERSC_HOME}:${PATH}"

# Install Conda dependencies
RUN conda create -y \
        -n "${CONDA_ENV}" \
        -c conda-forge \
        -c bioconda \
        python=3.10 \
        "fastq-pair=${FASTQ_PAIR_VER}" \
        "star=${STAR_VER}" \
        "samtools=${SAMTOOLS_VER}" \
        "multiqc=${MULTIQC_VER}" \
    && conda clean -afy

ENV PATH="/opt/conda/envs/${CONDA_ENV}/bin:${PATH}"

# Clone UniverSC, including Git LFS files
RUN wget -q \
        "https://github.com/minoda-lab/universc/archive/refs/tags/${UNIVERSC_VER}.tar.gz" \
        -O "/tmp/universc-${UNIVERSC_VER}.tar.gz" \
    && tar -xzf "/tmp/universc-${UNIVERSC_VER}.tar.gz" -C /tmp \
    && mv "/tmp/universc-${UNIVERSC_VER}" "${UNIVERSC_HOME}" \
    && rm -f "/tmp/universc-${UNIVERSC_VER}.tar.gz"

# The launcher is now supplied by the cloned repository:
# /universc/launch_universc.sh
# Build the small Cell Ranger reference datasets
RUN cd "${UNIVERSC_HOME}/test/cellranger_reference/cellranger-tiny-ref" \
    && rm -rf 3.0.0 1.2.0 \
    && cellranger mkref \
        --genome=3.0.0 \
        --fasta=genome-3.0.0.fa \
        --genes=genes-3.0.0.gtf \
    && cellranger mkref \
        --genome=1.2.0 \
        --fasta=genome-1.2.0.fa \
        --genes=genes-1.2.0.gtf

# Build UniverSC references
RUN cd ${UNIVERSC_HOME} \
    && make reference

# Link test FASTQ directories into the Cell Ranger installation
RUN ln -s \
        "${UNIVERSC_HOME}/test/shared/cellranger-tiny-fastq" \
        "${CELLRANGER_HOME}/cellranger-tiny-fastq" \
    && ln -s \
        "${UNIVERSC_HOME}/test/cellranger_reference/cellranger-tiny-ref" \
        "${CELLRANGER_HOME}/cellranger-tiny-ref"

# Install the UniverSC launcher in the Cell Ranger executable directory
RUN ln -s \
        "${UNIVERSC_HOME}/launch_universc.sh" \
        "${CELLRANGER_HOME}/bin/universc" \
    && ln -s \
        "${UNIVERSC_HOME}/sub" \
        "${CELLRANGER_HOME}/bin/sub" \
    && ln -s \
        "${UNIVERSC_HOME}/whitelists" \
        "${CELLRANGER_HOME}/bin/whitelists"

# Preserve the original UniverSC compatibility workaround
RUN cp \
        "${CELLRANGER_HOME}/lib/python/cellranger/chemistry.py" \
        "${CELLRANGER_HOME}/lib/python/cellranger/check.py"

# UniverSC 1.2.7 compatibility patches:
# - Bash array syntax
# - Cell Ranger 10.1.0 --create-bam syntax
# - Disable obsolete UniverSC Cloupe MRO manipulation for Cell Ranger 10.1.0

# Cell Ranger 10.1.0 changed the internal MRO structure used for
# Cloupe generation. UniverSC 1.2.7 attempts to restore/disable the
# legacy CLOUPE_PREPROCESS stage, which causes MRO preflight/parse
# failures with Cell Ranger 10.1.0.
#
# For Cell Ranger 10.1.0, bypass UniverSC's legacy Cloupe MRO
# manipulation. This does not disable Cell Ranger's core quantification,
# cell calling, matrices, or other primary analysis outputs.
RUN sed -i 's/\[\["$technology" == "vasa-drop" \]\]/\[\[ "$technology" == "vasa-drop" \]\]/' \
        "${UNIVERSC_HOME}/launch_universc.sh" \
    && sed -i 's/$index1(\[0\])/"${index1[0]}"/g' \
        "${UNIVERSC_HOME}/launch_universc.sh" \
    && sed -i 's/${index2}(\[0\])/"${index2[0]}"/g' \
        "${UNIVERSC_HOME}/launch_universc.sh" \
    && sed -i 's/$index2(\[0\])/"${index2[0]}"/g' \
        "${UNIVERSC_HOME}/launch_universc.sh" \
    && sed -i 's/bam="--no-bam"/bam="--create-bam=false"/' \
        "${UNIVERSC_HOME}/launch_universc.sh" \
    && sed -i 's/bam=""/bam="--create-bam=true"/' \
        "${UNIVERSC_HOME}/launch_universc.sh" \
    && sed -i '2640c\        if [[ $technology == "10x" ]] \&\& [[ $barcodefile == "default:10x" ]] \&\& [[ ${cellrangerversion} != "10.1.0" ]]; then' \
        "${UNIVERSC_HOME}/launch_universc.sh" \
    && sed -i '2684c\        elif [[ ${cellrangerversion} != "10.1.0" ]] \&\& { [[ $lastcall_p == "default:10x" ]] || [[ ! -f $lastcallfile ]]; }; then' \
        "${UNIVERSC_HOME}/launch_universc.sh"

# Match the permissions used by the original image
RUN chmod -R a+rwX \
    "${CELLRANGER_HOME}" \
    "${UNIVERSC_HOME}"
