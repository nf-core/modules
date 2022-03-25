process KAIJU_KAIJU {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::kaiju=1.8.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kaiju:1.8.2--h5b5514e_1':
        'quay.io/biocontainers/kaiju:1.8.2--h5b5514e_1' }"

    input:
    tuple val(meta), path(reads)
    tuple path(db), path(dbnodes)

    output:
    tuple val(meta), path('*.tsv'), emit: results
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = meta.single_end ? "-i ${reads}" : "-i ${reads[0]} -j ${reads[1]}"
    """
    kaiju \\
        $args \\
        -z $task.cpus \\
        -t ${dbnodes} \\
        -f ${db} \\
        -o ${prefix}.tsv \\
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kaiju: \$(echo \$( kaiju -h 2>&1 | sed -n 1p | sed 's/^.*Kaiju //' ))
    END_VERSIONS
    """
}
