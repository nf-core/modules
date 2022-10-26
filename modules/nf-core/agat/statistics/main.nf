process AGAT_STATISTICS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::agat=0.9.2-1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:0.9.2--pl5321hdfd78af_1' :
        'quay.io/biocontainers/agat:0.9.2--pl5321hdfd78af_1 ' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("*.txt"), emit: stats_file
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def program = args.contains("--basic") ? "agat_sq_stat_basic.pl" :
                  "agat_sp_statistics.pl"

    def inputarg = args.contains("--basic") ? "-i" : "--gff"

    args = args.replace("--basic", "")

    """

    ${program} \\
        $args \\
        --output ${prefix}.stats.txt \\
        ${inputarg} $gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_sp_statistics.pl --help |head -n3 | tail -n1 | grep -Eo '[0-9.]+')
    END_VERSIONS
    """
}
