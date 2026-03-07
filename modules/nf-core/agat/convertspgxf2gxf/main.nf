process AGAT_CONVERTSPGXF2GXF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/03/033434db0bd6ba28660401e1059286f36641fd8ce55faa11973fe5eaf312adcd/data' :
        'community.wave.seqera.io/library/agat:1.5.1--ae3cd948ce5e9795' }"

    input:
    tuple val(meta), path(gxf)

    output:
    tuple val(meta), path("*.agat.gff"), emit: output_gff
    tuple val(meta), path("*.log")     , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    agat_convert_sp_gxf2gxf.pl \\
        --gxf ${gxf} \\
        --output ${prefix}.agat.gff \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_convert_sp_gxf2gxf.pl --help | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.agat.gff
    touch ${gxf}.agat.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_convert_sp_gxf2gxf.pl --help | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p')
    END_VERSIONS
    """
}
