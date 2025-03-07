process PORTCULLIS_FULL {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::portcullis=1.2.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/portcullis:1.2.4--py38haf070c8_0':
        'biocontainers/portcullis:1.2.4--py38haf070c8_0' }"

    input:
    tuple val(meta) , path(bam)
    tuple val(meta2), path(bed)
    tuple val(meta3), path(fasta)

    output:
    tuple val(meta), path("*.log"), emit: log
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.tab"), emit: tab
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    portcullis \\
        full \\
        ${args} \\
        -t ${task.cpus} \\
        -o $prefix \\
        -r $bed \\
        $fasta \\
        $bam > ${prefix}.portcullis.log

    cp ${prefix}/3-filt/*.bed .
    cp ${prefix}/3-filt/*.tab .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        portcullis: \$(portcullis --version |& sed '1!d ; s/portcullis //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.portcullis.log
    touch ${prefix}.portcullis.bed
    touch ${prefix}.portcullis.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        portcullis: \$(portcullis --version |& sed '1!d ; s/portcullis //')
    END_VERSIONS
    """
}

