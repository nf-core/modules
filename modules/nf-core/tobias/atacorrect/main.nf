process TOBIAS_ATACORRECT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tobias:0.17.5--py310h3479294_0':
        'quay.io/biocontainers/tobias:0.17.5--py310h3479294_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(peaks), path(fasta)

    output:
    tuple val(meta), path("*_corrected.bw")  , emit: corrected
    tuple val(meta), path("*_expected.bw")   , emit: expected
    tuple val(meta), path("*_uncorrected.bw"), emit: uncorrected, optional: true
    tuple val(meta), path("*_bias.bw")       , emit: bias, optional: true
    tuple val(meta), path("*_atacorrect.pdf"), emit: report, optional: true
    tuple val("${task.process}"), val('tobias'), eval('TOBIAS --version'), topic: versions, emit: versions_tobias

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p .matplotlib
    export MPLCONFIGDIR="\${PWD}/.matplotlib"

    ln -s $bam ${prefix}.bam
    if [[ "$bai" == *.csi ]]; then
        ln -s $bai ${prefix}.bam.csi
    else
        ln -s $bai ${prefix}.bam.bai
    fi

    TOBIAS ATACorrect \\
        --bam ${prefix}.bam \\
        --genome $fasta \\
        --peaks $peaks \\
        --prefix $prefix \\
        --outdir . \\
        --cores ${task.cpus} \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_corrected.bw
    touch ${prefix}_expected.bw
    touch ${prefix}_uncorrected.bw
    touch ${prefix}_bias.bw
    touch ${prefix}_atacorrect.pdf
    """
}
