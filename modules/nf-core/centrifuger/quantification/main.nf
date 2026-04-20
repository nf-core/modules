process CENTRIFUGER_QUANTIFICATION {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/centrifuger:1.1.0--hf426362_0':
        'biocontainers/centrifuger:1.1.0--hf426362_0' }"

    input:
    tuple val(meta), path(classification_file)
    path db

    output:
    tuple val(meta), path("${meta.id}.tsv"), emit: report_file
    tuple val("${task.process}"), val('centrifuger'), eval("centrifuger -v 2>&1 | sed 's/Centrifuger v//'"), emit: versions_centrifuger,  topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    db_name=`find -L ${db} -name "*.1.cfr" -not -name "._*"  | sed 's/\\.1.cfr\$//'`

    centrifuger-quant \\
        -x \$db_name \\
        -c ${classification_file} \\
        ${args} > ${prefix}.tsv
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    #output
    touch ${prefix}.tsv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        centrifuger -v 2>&1 | sed 's/Centrifuger v//')
    END_VERSIONS
    """

}
