process HOMER_MERGEPEAKS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // singularity build url: https://wave.seqera.io/view/builds/bd-9c603739ae7d4fd3_1
    // docker build url: https://wave.seqera.io/view/builds/bd-08c7bb832e96c6bd_1
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-edger_homer_samtools_pruned:9c603739ae7d4fd3'
        : 'community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-edger_homer_samtools_pruned:08c7bb832e96c6bd'}"

    input:
    tuple val(meta), path(peaks)  // peaks can be a list of peak files

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}" + "_merged_peaks"
    def peak_files = peaks instanceof List ? peaks.join(' ') : peaks
    def VERSION = '5.1'

    """
    mergePeaks \\
        $args \\
        -prefix ${prefix} \\
        $peak_files \\
        > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: ${VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_merged_peaks"
    def VERSION = '5.1'
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: ${VERSION}
    END_VERSIONS
    """
}
