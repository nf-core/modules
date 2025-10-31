process PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pb-cpg-tools:3.0.0--h9ee0642_0':
        'biocontainers/pb-cpg-tools:3.0.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bed.gz"),      emit: bed
    tuple val(meta), path("*.bed.gz.tbi"),  emit: bed_index
    tuple val(meta), path("*.bw"),          emit: bigwig
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    aligned_bam_to_cpg_scores \\
        --bam ${bam} \\
        --output-prefix ${prefix} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbcpgtools: \$(aligned_bam_to_cpg_scores --version | sed 's/aligned_bam_to_cpg_scores //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

    echo "" | gzip > ${prefix}.combined.bed.gz
    touch ${prefix}.combined.bed.gz.tbi
    touch ${prefix}.combined.bw

    echo "" | gzip > ${prefix}.hap1.bed.gz
    touch ${prefix}.hap1.bed.gz.tbi
    touch ${prefix}.hap1.bw

    echo "" | gzip > ${prefix}.hap2.bed.gz
    touch ${prefix}.hap2.bed.gz.tbi
    touch ${prefix}.hap2.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbcpgtools: \$(aligned_bam_to_cpg_scores --version | sed 's/aligned_bam_to_cpg_scores //')
    END_VERSIONS
    """
}
