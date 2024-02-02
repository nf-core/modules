process CENTRIFUGE_KREPORT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/centrifuge-core:1.0.4.1--h43eeafb_0' :
        'biocontainers/centrifuge-core:1.0.4.1--h43eeafb_0' }"

    input:
    tuple val(meta), path(report)
    path db

    output:
    tuple val(meta), path('*.txt'), emit: kreport
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    db_name=`find -L ${db} -name "*.1.cf" -not -name "._*"  | sed 's/\\.1.cf\$//'`
    centrifuge-kreport -x \$db_name ${report} > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        centrifuge: \$( centrifuge --version  | sed -n 1p | sed 's/^.*centrifuge-class version //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        centrifuge: \$( centrifuge --version  | sed -n 1p | sed 's/^.*centrifuge-class version //')
    END_VERSIONS
    """
}
