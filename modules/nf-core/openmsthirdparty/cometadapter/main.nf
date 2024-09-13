process OPENMSTHIRDPARTY_COMETADAPTER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms-thirdparty:3.1.0--h9ee0642_4' :
        'biocontainers/openms-thirdparty:3.1.0--h9ee0642_4' }"

    input:
    tuple val(meta), path(mzml), path(fasta)

    output:
    tuple val(meta), path("*.idXML"), emit: idxml
    tuple val(meta), path("*.tsv")  , emit: pin, optional: true
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    CometAdapter \\
        -in $mzml \\
        -database $fasta \\
        -out ${prefix}.idXML \\
        -threads $task.cpus \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CometAdapter: \$(CometAdapter 2>&1 | grep -E '^Version(.*)' | sed 's/Version: //g' | cut -d ' ' -f 1 | cut -d '-' -f 1)
        Comet: \$(comet 2>&1 | grep -E "Comet version.*" | sed 's/Comet version //g' | sed 's/"//g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.idXML
    touch ${prefix}_pin.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CometAdapter: \$(CometAdapter 2>&1 | grep -E '^Version(.*)' | sed 's/Version: //g' | cut -d ' ' -f 1 | cut -d '-' -f 1)
        Comet: \$(comet 2>&1 | grep -E "Comet version.*" | sed 's/Comet version //g' | sed 's/"//g')
    END_VERSIONS
    """
}
