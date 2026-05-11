process SYLPH_SKETCH {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sylph:0.9.0--ha6fb395_0'
        : 'quay.io/biocontainers/sylph:0.9.0--ha6fb395_0'}"

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path('my_sketches/*.sylsp'), path('database.syldb'), emit: sketch_fastq_genome
    tuple val("${task.process}"), val('sylph'), eval('sylph -V | sed "s/sylph //g"'), topic: versions, emit: versions_sylph

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fastq = meta.single_end ? "-r ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    sylph sketch \\
        ${args} \\
        ${fastq} \\
        -g ${reference} \\
        -S ${prefix} \\
        -d my_sketches \\
        -t ${task.cpus}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p my_sketches
    touch my_sketches/${prefix}.sylsp
    touch database.syldb
    """
}
