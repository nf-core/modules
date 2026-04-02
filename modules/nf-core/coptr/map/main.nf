process COPTR_MAP {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coptr:1.1.4--pyhdfd78af_3':
        'biocontainers/coptr:1.1.4--pyhdfd78af_3' }"

    input:
    tuple val(meta), path(fasta, stageAs: "fastafolder/*")
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('coptr'), eval("coptr |& sed -E '11!d ; s/CoPTR.*?\\(v(.*?)\\).*/\\1/'"), emit: versions_coptr, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def prefix2 = task.ext.prefix ?: "${meta2.id}"

    def paired_end = ""
    if ( ! meta.single_end ) {
        paired_end = "--paired"
    }
    """
    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
    [ -z "\$INDEX" ] && INDEX=`find -L ./ -name "*.rev.1.bt2l" | sed "s/\\.rev.1.bt2l\$//"`
    [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

    coptr \
        map \
        $args $paired_end \
        --threads $task.cpus \
        \$INDEX \
        fastafolder \
        .

    mv ${prefix}.bam ${prefix}_${prefix2}.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
