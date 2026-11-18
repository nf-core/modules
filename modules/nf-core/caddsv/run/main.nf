process CADDSV_RUN {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/97/970be9efa82690830164eeabb89b675528d8a2761f844532988246af4a8e6885/data'
:         'community.wave.seqera.io/library/bedops_bedtools_caddsv_htslib_pruned:0e20bb2899a7a4f6' }"
    input:
    tuple val(meta), path(variants)
    path(annotations_dir)
    path(config)

    output:
    tuple val(meta), path("caddsv_results/scored/*.tsv"), emit: tsv
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
        error "Unsupported CADD-SV input suffix: ${variants}"
    }
    def run_input = "${prefix}.${real_ext}"

    // TSV input (REF, ALT, [TYPE], [ID]) is only valid in sequence-only mode
    def seqonly_arg = real_ext == 'tsv' ? '--seqonly' : ''
    def config_arg = config ? "--config ${config}" : ""
    // Snakemake pins OMP_NUM_THREADS to a rule's threads (1 by default), which
    // limits CPU SegmentNT inference to a single thread
    def nt_rules = ['run_nt_script', 'run_nt_script_ref', 'seqonly_run_nt_sample', 'seqonly_run_nt_ref']
    def set_threads = nt_rules.collect { rule -> "${rule}=${task.cpus}" }.join(' ')
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
    caddsv run "${run_input}" \\
        --annotations-dir "${annotations_dir}" \\
        --output-dir caddsv_results \\
        --threads ${task.cpus} \\
        --no-use-conda \\
        --snakemake-args "--set-threads ${set_threads}" \\
        ${seqonly_arg} \\
        ${config_arg} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_gz = variants.extension == 'gz'
    def real_ext = is_gz ? file(variants.baseName).extension : variants.extension
    """
    mkdir -p caddsv_results/scored
    if [ "${real_ext}" = "tsv" ]; then
        printf 'id\\ttype\\tCADD-SV_seqonly_PHRED\\tCADD-SV_seqonly_score\\n' > caddsv_results/scored/${prefix}_seqonly_score.tsv
    else
        printf '#chr\\tstart\\tend\\ttype\\tCADD-SV_PHRED\\tCADD-SV_score\\n' > caddsv_results/scored/${prefix}_score.tsv
    fi
    """
}
