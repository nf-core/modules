process CDHIT_CDHITEST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cd-hit%3A4.8.1--h5b5514e_7':
        'biocontainers/cd-hit:4.8.1--h5b5514e_7' }"

    input:
    tuple val(meta), path(sequences)

    output:
    tuple val(meta), path("*.{fa,fq}")    ,emit: fasta
    tuple val(meta), path("*.clstr")      ,emit: clusters
    tuple val("${task.process}"), val('cdhit'), eval("cd-hit-est -h | sed -n '1s/.*version \\([0-9.]*\\).*/\\1/p'"), topic: versions, emit: versions_cdhitest

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "${sequences}" ==~ /(.*f[astn]*a(.gz)?$)/ ? "fa" : "fq"

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[cd-hit-est] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    cd-hit-est \\
        $args \\
        -i ${sequences} \\
        -o ${prefix}.${suffix} \\
        -M $avail_mem \\
        -T $task.cpus
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "${sequences}" ==~ /(.*f[astn]*a(.gz)?$)/ ? "fa" : "fq"
    """
    echo "${args}"
    touch ${prefix}.${suffix}
    touch ${prefix}.${suffix}.clstr
    """
}
