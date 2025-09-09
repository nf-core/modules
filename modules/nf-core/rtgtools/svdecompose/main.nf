process RTGTOOLS_SVDECOMPOSE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rtg-tools:3.12.1--hdfd78af_0':
        'biocontainers/rtg-tools:3.12.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(input), path(tbi)

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: index
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ""
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def index_vcf = tbi ? "" : "rtg index ${input}"
    def avail_mem = task.memory.toGiga() + "G"
    if ("$input" == "${prefix}.vcf.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    ${index_vcf}

    rtg RTG_MEM=$avail_mem svdecompose \\
        ${args} \\
        --input=${input} \\
        --output=${prefix}.vcf.gz

    rtg index ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rtg-tools: \$(echo \$(rtg version | head -n 1 | awk '{print \$4}'))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$input" == "${prefix}.vcf.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    echo | gzip -n > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rtg-tools: \$(echo \$(rtg version | head -n 1 | awk '{print \$4}'))
    END_VERSIONS
    """
}
