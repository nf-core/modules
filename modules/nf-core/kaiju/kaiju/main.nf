process KAIJU_KAIJU {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::kaiju=1.8.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kaiju:1.8.2--h5b5514e_1':
        'quay.io/biocontainers/kaiju:1.8.2--h5b5514e_1' }"

    input:
    tuple val(meta), path(reads)
    path(db)

    output:
    tuple val(meta), path('*.tsv'), emit: results
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = meta.single_end ? "-i ${reads}" : "-i ${reads[0]} -j ${reads[1]}"
    """
    dbnodes=`find -L ${db} -name "*nodes.dmp"`
    dbname=`find -L ${db} -name "*.fmi" -not -name "._*"`
    kaiju \\
        $args \\
        -z $task.cpus \\
        -t \$dbnodes \\
        -f \$dbname \\
        -o ${prefix}.tsv \\
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kaiju: \$(echo \$( kaiju -h 2>&1 | sed -n 1p | sed 's/^.*Kaiju //' ))
    END_VERSIONS
    """
}
