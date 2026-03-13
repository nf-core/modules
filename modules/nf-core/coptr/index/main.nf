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
    tuple val("${task.process}"), val('coptr'), eval("coptr |& sed -E '11!d ; s/CoPTR.*?\\(v(.*?)\\).*/\\1/'"), emit: versions_coptr, topic: versions

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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir bowtie2
    touch bowtie2/${prefix}.{1..4}.bt2
    touch bowtie2/${prefix}.rev.{1,2}.bt2
    touch bowtie2/${prefix}.genomes
    """
}
