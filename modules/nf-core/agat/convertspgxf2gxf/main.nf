process AGAT_CONVERTSPGXF2GXF {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::agat=1.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.0.0--pl5321hdfd78af_0' :
        'biocontainers/agat:1.0.0--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("*.agat.gff"), emit: output_gff
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    agat_convert_sp_gxf2gxf.pl \\
        --gff $gff \\
        --output ${prefix}.agat.gff \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_convert_sp_gxf2gxf.pl --help | sed '4!d; s/.*v//')
    END_VERSIONS
    """
}
