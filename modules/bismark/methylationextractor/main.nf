process BISMARK_METHYLATIONEXTRACTOR {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::bismark=0.23.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bismark:0.23.0--0' :
        'quay.io/biocontainers/bismark:0.23.0--0' }"

    input:
    tuple val(meta), path(bam)
    path index

    output:
    tuple val(meta), path("*.bedGraph.gz")         , emit: bedgraph
    tuple val(meta), path("*.txt.gz")              , emit: methylation_calls
    tuple val(meta), path("*.cov.gz")              , emit: coverage
    tuple val(meta), path("*_splitting_report.txt"), emit: report
    tuple val(meta), path("*.M-bias.txt")          , emit: mbias
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    // Some arguments are only applicable to paired_end data
    def overlap          = params.no_overlap           ? '--no_overlap'                                  : '--include_overlap'
    def ignore_r2        = params.ignore_r2 > 0        ? "--ignore_r2 ${params.ignore_r2}"               : ''
    def ignore_3prime_r2 = params.ignore_3prime_r2 > 0 ? "--ignore_3prime_r2 ${params.ignore_3prime_r2}" : ''

    if (meta.single_end) {
        """
        bismark_methylation_extractor \\
            --bedGraph \\
            --counts \\
            --gzip \\
            --report \\
            -s \\
            $args \\
            $bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
        END_VERSIONS
        """
    } else {
        """
        bismark_methylation_extractor \\
            --bedGraph \\
            --counts \\
            --gzip \\
            --report \\
            -p \\
            $args \\
            $overlap \\
            $ignore_r2 \\
            $ignore_3prime_r2 \\
            $bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
        END_VERSIONS
        """
    }

}
