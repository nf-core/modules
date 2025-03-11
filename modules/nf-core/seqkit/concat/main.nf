process SEQKIT_CONCAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0':
        'biocontainers/seqkit:2.9.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(input, stageAs: 'in/*')

    output:
    tuple val(meta), path("*.{fasta,fastq,fa,fq,fas,fna,faa}"), emit: fastx
    path "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args         ?: ""
    def prefix      = task.ext.prefix       ?: "${meta.id}"
    def file_type   = input instanceof List ? input[0].getExtension() : input.getExtension()
    """
    seqkit \\
        concat \\
        $args \\
        in/* > ${prefix}.${file_type}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    """
    touch ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """
}
