process AGAT_SQSTATBASIC {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::agat=1.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.0.0--pl5321hdfd78af_0' :
        'quay.io/biocontainers/agat:1.0.0--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("*.txt"), emit: stats_txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    agat_sq_stat_basic.pl \\
        -i $gff \\
        --output ${prefix}.stats.txt \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_sq_stat_basic.pl --help |head -n4 | tail -n1 | grep -Eo '[0-9.]+')
    END_VERSIONS
    """
}
