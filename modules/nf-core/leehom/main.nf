process LEEHOM {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/leehom:1.2.15--h29e30f7_1' :
        'biocontainers/leehom:1.2.15--h29e30f7_1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.bam")          , emit: bam                , optional: true
    tuple val(meta), path("${prefix}.fq.gz")        , emit: fq_pass            , optional: true
    tuple val(meta), path("${prefix}.fail.fq.gz")   , emit: fq_fail            , optional: true
    tuple val(meta), path("${prefix}_r1.fq.gz")     , emit: unmerged_r1_fq_pass, optional: true
    tuple val(meta), path("${prefix}_r1.fail.fq.gz"), emit: unmerged_r1_fq_fail, optional: true
    tuple val(meta), path("${prefix}_r2.fq.gz")     , emit: unmerged_r2_fq_pass, optional: true
    tuple val(meta), path("${prefix}_r2.fail.fq.gz"), emit: unmerged_r2_fq_fail, optional: true
    tuple val(meta), path("*.log")                  , emit: log
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2.15' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    if (reads.toString().endsWith('.bam')) {
        """
        leeHom \\
            ${args} \\
            -t ${task.cpus} \\
            -o ${prefix}.bam \\
            --log ${prefix}.log \\
            ${reads}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            leehom: ${VERSION}
        END_VERSIONS
        """
    } else if (meta.single_end) {
        """
        leeHom \\
            ${args} \\
            -t ${task.cpus} \\
            -fq1 ${reads} \\
            -fqo ${prefix} \\
            --log ${prefix}.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            leehom: ${VERSION}
        END_VERSIONS
        """
    } else {
        """
        leeHom \\
            ${args} \\
            -t ${task.cpus} \\
            -fq1 ${reads[0]} \\
            -fq2 ${reads[1]} \\
            -fqo ${prefix} \\
            --log ${prefix}.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            leehom: ${VERSION}
        END_VERSIONS
        """
    }

    stub:
    prefix        = task.ext.prefix ?: "${meta.id}"
    def VERSION   = '1.2.15' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    is_bam        = reads.toString().endsWith('.bam')
    is_single_end = meta.single_end

    """
    if [[ "${is_bam}" == "true" ]]; then
        touch ${prefix}.bam
    else
        echo "" | gzip > ${prefix}.fq.gz
        echo "" | gzip > ${prefix}.fail.fq.gz
        if [[ "${is_single_end}" == "false" ]]; then
            echo "" | gzip > ${prefix}_r1.fq.gz
            echo "" | gzip > ${prefix}_r1.fail.fq.gz
            echo "" | gzip > ${prefix}_r2.fq.gz
            echo "" | gzip > ${prefix}_r2.fail.fq.gz
        fi
    fi
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        leehom: ${VERSION}
    END_VERSIONS
    """

}
