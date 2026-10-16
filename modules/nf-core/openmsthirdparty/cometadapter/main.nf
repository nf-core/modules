process OPENMSTHIRDPARTY_COMETADAPTER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/64/6475e53b3776cd3b71b8ff525f01da9c45dbe9401a00e5630d91ab976f20e8aa/data'
:         'community.wave.seqera.io/library/openms-thirdparty:3.6.0--566f6f290d39203b' }"

    input:
    tuple val(meta), path(mzml), path(fasta)

    output:
    tuple val(meta), path("*.idXML"), emit: idxml
    tuple val(meta), path("*.tsv")  , emit: pin, optional: true
    tuple val("${task.process}"), val('CometAdapter'), eval("CometAdapter --help 2>&1 | sed -nE 's/^Version: ([0-9.]+).*/\\1/p'"), emit: versions_cometadapter, topic: versions
    tuple val("${task.process}"), val('Comet'), eval("comet 2>&1 | sed -n 's/.*Comet version \" *\\(.*\\)\".*/\\1/p'"), emit: versions_comet, topic: versions

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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.idXML
    touch ${prefix}_pin.tsv
    """
}
