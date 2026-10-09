// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CADDSV_RUN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d2/d2fd85e11209001044c0b5f8cd663e8058cf006b8e6a0d8712a50857677ab664/data'
        : 'community.wave.seqera.io/library/caddsv_snakemake_pulp_bedops_pruned:799f9464989eb8fe'}"

    input:
    tuple val(meta), path(variants)
    path annotations_dir
    path config

    output:
    tuple val(meta), path("caddsv_results/scored/*.tsv"), emit: tsv
    tuple val(meta), path("caddsv_results"), emit: results
    tuple val(meta), path("caddsv_run_*.log"), emit: log, optional: true
    tuple val("${task.process}"), val('caddsv'), eval("caddsv --version"), emit: versions_caddsv, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_gz = variants.extension == 'gz'
    def real_ext = is_gz ? file(variants.baseName).extension : variants.extension
    if (!(real_ext in ['bed', 'tsv'])) {
        error("Unsupported CADD-SV input suffix: ${variants}")
    }
    def run_input = "${prefix}.${real_ext}"

    def config_arg = config ? "--config ${config}" : ""
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export PYTHONNOUSERSITE=1
    export HF_HUB_OFFLINE="\${HF_HUB_OFFLINE:-1}"
    export TRANSFORMERS_OFFLINE="\${TRANSFORMERS_OFFLINE:-1}"

    if ${is_gz}; then
        gzip -cdf "${variants}" > "${run_input}"
    else
        cp -L "${variants}" "${run_input}"
    fi

    # The module environment provides every workflow rule dependency, so the
    # inner Snakemake must not create conda envs or pull rule containers.
    # --drop-metadata: Snakemake 8 otherwise calls `conda info` to hash each
    # rule's conda env when recording job metadata, even with conda disabled.
    caddsv run "${run_input}" \\
        --annotations-dir "${annotations_dir}" \\
        --output-dir caddsv_results \\
        --threads ${task.cpus} \\
        --no-use-conda \\
        --snakemake-args "--drop-metadata" \\
        ${config_arg} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p caddsv_results/scored
    printf 'chr\\tstart\\tend\\ttype\\tCADD-SV_PHRED\\tCADD-SV_score\\n' > caddsv_results/scored/${prefix}_score.tsv
    """
}
