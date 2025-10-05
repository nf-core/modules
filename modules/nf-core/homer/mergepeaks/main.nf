process HOMER_MERGEPEAKS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/homer:5.1--586e112615c26ca1':
        'community.wave.seqera.io/library/homer:5.1--851c883b9d1473f9' }"

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
    """
    mergePeaks \\
        $args \\
        -prefix ${prefix} \\
        $peak_files \\
        > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: \$(mergePeaks 2>&1 | grep -i "version" | sed 's/.*version //i' || echo "4.11")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: 4.11
    END_VERSIONS
    """
}
