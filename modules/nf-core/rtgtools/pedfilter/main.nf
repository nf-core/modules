process RTGTOOLS_PEDFILTER {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::rtg-tools=3.12.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rtg-tools:3.12.1--hdfd78af_0':
        'biocontainers/rtg-tools:3.12.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.{vcf.gz,ped}") , emit: output
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def extension = args.contains("--vcf") ? "vcf.gz" : "ped"

    if( "${prefix}.${extension}" == "${input}" ) {
        error "The input and output file have the same name, please use another ext.prefix."
    }

    def postprocess = extension == "vcf.gz" ? "| rtg bgzip ${args2} -" : ""

    """
    rtg pedfilter \\
        ${args} \\
        ${input} \\
    ${postprocess} > ${prefix}.${extension}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rtgtools: \$(echo \$(rtg version | head -n 1 | awk '{print \$4}'))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    def extension = args.contains("--vcf") ? "vcf.gz" : "ped"

    if( "${prefix}.${extension}" == "${input}" ) {
        error "The input and output file have the same name, please use another ext.prefix."
    }

    """
    touch ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rtgtools: \$(echo \$(rtg version | head -n 1 | awk '{print \$4}'))
    END_VERSIONS
    """
}
