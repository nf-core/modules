process SUSHIE {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    container "nf-core/sushie:0.19"

    input:
    tuple val(meta), path(study_locus_files)
    path(ld_files)
    val(sample_sizes)

    output:
    tuple val(meta), path("*.sushie.corr.tsv.gz")   , emit: corr
    tuple val(meta), path("*.sushie.cs.tsv.gz")     , emit: cs
    tuple val(meta), path("*.sushie.weights.tsv.gz"), emit: weights
    tuple val(meta), path("*.log")                  , emit: log
    tuple val("${task.process}"), val('sushie'), val('0.19'), topic: versions, emit: versions_sushie

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // HOME is set to a writable location to avoid pathlib to fail when creating cache files
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SUSHIE module does not support Conda. Please use Docker instead."
    }
    """
    export HOME=\$PWD/nxf_home

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
