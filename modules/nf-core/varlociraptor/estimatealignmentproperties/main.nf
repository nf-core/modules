process VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f1/f1d2b5a7830e21d0407be9958b69810f16fac8008fb899c345f635ff5636de04/data'
:         'community.wave.seqera.io/library/varlociraptor:8.9.5--9445270256470baf' }"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.alignment-properties.json"), emit: alignment_properties_json
    tuple val("${task.process}"), val('varlociraptor'), eval("varlociraptor --version | sed 's/^varlociraptor //'"), topic: versions, emit: versions_varlociraptor

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    varlociraptor estimate alignment-properties \\
        ${fasta} \\
        --bams ${bam} \\
        ${args} \\
        > ${prefix}.alignment-properties.json
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.alignment-properties.json
    """
}
