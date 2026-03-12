process SUSHIE {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/cameronlloyd/sushie:0.19"

    input:
    tuple val(meta), path(study_locus_files)
    path(ld_files)
    val(sample_sizes)

    output:
    tuple val(meta), path("*.sushie.corr.tsv")   , emit: corr
    tuple val(meta), path("*.sushie.cs.tsv")     , emit: cs
    tuple val(meta), path("*.sushie.weights.tsv"), emit: weights
    tuple val("${task.process}"), val('sushie'), val('0.17'), topic: versions, emit: versions_sushie

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sushie \\
    finemap \\
    --summary \\
    --gwas ${study_locus_files.join(' ')} \\
    --ld ${ld_files.join(' ')} \\
    --sample-size ${sample_sizes} \\
    --output ${prefix} \\
    ${args}
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.sushie.corr.tsv
    touch ${prefix}.sushie.cs.tsv
    touch ${prefix}.sushie.weights.tsv
    """
}
