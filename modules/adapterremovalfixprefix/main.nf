def VERSION = '0.05'

process ADAPTERREMOVALFIXPREFIX {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::adapterremovalfixprefix=0.0.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/adapterremovalfixprefix:0.0.5--hdfd78af_2':
        'quay.io/biocontainers/adapterremovalfixprefix:0.0.5--hdfd78af_2' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.fq.gz"), emit: fixed_fastq
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$fastq" == "${prefix}.fq.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    AdapterRemovalFixPrefix \\
        $fastq \\
        $args \\
        | gzip > ${prefix}.fq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        adapterremovalfixprefix: $VERSION
    END_VERSIONS
    """
}
