def VERSION = '1.2.15' // Version information not provided by tool on CLI

process LEEHOM {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::leehom=1.2.15" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/leehom:1.2.15--h29e30f7_1' :
        'quay.io/biocontainers/leehom:1.2.15--h29e30f7_1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.bam")          , optional: true, emit: bam
    tuple val(meta), path("${prefix}.fq.gz")        , optional: true, emit: fq_pass
    tuple val(meta), path("${prefix}.fail.fq.gz")   , optional: true, emit: fq_fail
    tuple val(meta), path("${prefix}_r1.fq.gz")     , optional: true, emit: unmerged_r1_fq_pass
    tuple val(meta), path("${prefix}_r1.fail.fq.gz"), optional: true, emit: unmerged_r1_fq_fail
    tuple val(meta), path("${prefix}_r2.fq.gz")     , optional: true, emit: unmerged_r2_fq_pass
    tuple val(meta), path("${prefix}_r2.fail.fq.gz"), optional: true, emit: unmerged_r2_fq_fail
    tuple val(meta), path("*.log")                                  , emit: log
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    if (reads.toString().endsWith('.bam')) {
        """
        leeHom \\
            $args \\
            -t $task.cpus \\
            -o ${prefix}.bam \\
            --log ${prefix}.log \\
            $reads

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            leehom: $VERSION
        END_VERSIONS
        """
    } else if (meta.single_end) {
        """
        leeHom \\
            $args \\
            -t $task.cpus \\
            -fq1 $reads \\
            -fqo $prefix \\
            --log ${prefix}.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            leehom: $VERSION
        END_VERSIONS
        """
    } else {
        """
        leeHom \\
            $args \\
            -t $task.cpus \\
            -fq1 ${reads[0]} \\
            -fq2 ${reads[1]} \\
            -fqo $prefix \\
            --log ${prefix}.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            leehom: $VERSION
        END_VERSIONS
        """
    }
}
