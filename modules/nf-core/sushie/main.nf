process SUSHIE {
    tag "$meta.id"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bb/bb9fae2a9fe86eaafc849e163bc267cbb7a27320fffb938b3d2fc8ba20153f71/data'
        : 'community.wave.seqera.io/library/python_pip_c-compiler_git_pruned:eab3b4c658ec0e54'}"

    input:
    tuple val(meta), path(study_locus_files)
    path(ld_files)
    val(sample_sizes)

    output:
    tuple val(meta), path("*.sushie.corr.tsv.gz")   , emit: corr
    tuple val(meta), path("*.sushie.cs.tsv.gz")     , emit: cs
    tuple val(meta), path("*.sushie.weights.tsv.gz"), emit: weights
    tuple val(meta), path("*.log")                  , emit: log
    tuple val("${task.process}"), val('sushie'), eval('sushie --version'), topic: versions, emit: versions_sushie

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SUSHIE module does not support Conda. Please use Docker instead."
    }
    """
    sushie \\
    finemap \\
    --summary \\
    --gwas ${study_locus_files.join(' ')} \\
    --ld ${ld_files.join(' ')} \\
    --sample-size ${sample_sizes} \\
    --compress \\
    --output ${prefix} \\
    ${args}
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SUSHIE module does not support Conda. Please use Docker instead."
    }
    """
    echo ${args}

    echo "" | gzip > ${prefix}.sushie.corr.tsv.gz
    echo "" | gzip > ${prefix}.sushie.cs.tsv.gz
    echo "" | gzip > ${prefix}.sushie.weights.tsv.gz
    touch ${prefix}.log
    """
}
