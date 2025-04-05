process COPTR_INDEX {
    tag '$indexfasta'
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coptr:1.1.4--pyhdfd78af_3':
        'biocontainers/coptr:1.1.4--pyhdfd78af_3' }"

    input:
    tuple val(meta), path(indexfasta, stageAs: "fastafolder/*")

    output:
    tuple val(meta), path("bowtie2"), emit: index_dir
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir bowtie2
    coptr \
        index \
        $args \
        --bt2-threads $task.cpus \
        fastafolder \
        bowtie2/${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coptr: \$(coptr |& sed -E '11!d ; s/CoPTR.*?\\(v(.*?)\\).*/\\1/')
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir bowtie2
    touch bowtie2/${prefix}.{1..4}.bt2
    touch bowtie2/${prefix}.rev.{1,2}.bt2
    touch bowtie2/${prefix}.genomes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coptr: \$(coptr |& sed -E '11!d ; s/CoPTR.*?\\(v(.*?)\\).*/\\1/')
    END_VERSIONS
    """
}
