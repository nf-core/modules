process AGAT_CONVERTBED2GFF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.6.1--pl5321hdfd78af_1' :
        'biocontainers/agat:1.6.1--pl5321hdfd78af_1' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.gff"), emit: gff
    tuple val("${task.process}"), val('agat'), eval('agat_convert_bed2gff.pl --help | grep Version | cut -d" " -f11'), emit: versions_agat, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    agat_convert_bed2gff.pl \\
        --bed ${bed} \\
        --output ${prefix}.gff \\
        ${args}

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gff

    """
}
