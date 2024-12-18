process SYLPH_SKETCH {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sylph:0.7.0--h919a2d8_0' :
        'biocontainers/sylph:0.7.0--h919a2d8_0' }"

    input:
    tuple val(meta), path(reads)
    path(reference)

    output:
    tuple val(meta), path("my_sketches/*.sylsp"), path("database.syldb"), emit: sketch_fastq_genome
    path "versions.yml"                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fastq = meta.single_end ? "-r ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    sylph sketch \\
        $args \\
        $fastq \\
        -g $reference \\
        -S $prefix \\
        -d my_sketches

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: \$(sylph -V|awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p my_sketches
    touch my_sketches/${prefix}.sylsp
    touch database.syldb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: \$(sylph -V|awk '{print \$2}')
    END_VERSIONS
    """
}
