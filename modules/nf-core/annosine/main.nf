process ANNOSINE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/annosine2:2.0.7--pyh7cba7a3_0':
        'biocontainers/annosine2:2.0.7--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(fasta)
    val mode

    output:
    tuple val(meta), path("${prefix}.log")  , emit: log
    tuple val(meta), path("${prefix}.fa")   , emit: fa      , optional: true
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}_annosine"
    def VERSION = '2.0.7' // WARN: Manually update when changing Bioconda assets
    if ( "$fasta" == "${prefix}.fa" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    AnnoSINE_v2 \\
        $args \\
        --threads $task.cpus \\
        $mode \\
        $fasta \\
        $prefix \\
        &> >(tee ${prefix}.log 2>&1)

    mv \\
        $prefix/Seed_SINE.fa \\
        ${prefix}.fa \\
        || echo 'AnnoSINE_v2 did not find SINE sequences. See log for details!'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annosine: $VERSION
    END_VERSIONS
    """

    stub:
    prefix      = task.ext.prefix ?: "${meta.id}_annosine"
    def VERSION = '2.0.7' // WARN: Manually update when changing Bioconda assets
    if ( "$fasta" == "${prefix}.fa" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.log
    touch ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annosine: $VERSION
    END_VERSIONS
    """
}
