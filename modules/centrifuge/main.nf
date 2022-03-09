process CENTRIFUGE {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::centrifuge=1.0.4_beta" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/centrifuge:1.0.4_beta--h9a82719_6' :
        'quay.io/biocontainers/centrifuge:1.0.4_beta--h9a82719_6' }"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path('*report.txt')  , emit: report
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired = meta.single_end ? "-U ${reads}" :  "-1 ${reads[0]} -2 ${reads[1]}"
    """
    tar -xf $db
    centrifuge \\
        -x ${db.toString().replace(".tar.gz","")} \\
        -p $task.cpus \\
        $paired \\
        --report-file ${prefix}.centrifuge.report.txt \\
        $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        centrifuge: \$( centrifuge --version  | sed -n 1p | sed 's/^.*centrifuge-class version //')
    END_VERSIONS
    """
}