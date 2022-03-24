process BBMAP_PILEUP {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bbmap=38.92 bioconda::samtools=1.13 pigz=2.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-008daec56b7aaf3f162d7866758142b9f889d690:f5f55fc5623bb7b3f725e8d2f86bedacfd879510-0' :
        'quay.io/biocontainers/mulled-v2-008daec56b7aaf3f162d7866758142b9f889d690:f5f55fc5623bb7b3f725e8d2f86bedacfd879510-0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.stats.txt"), emit: covstats
    tuple val(meta), path("*.hist.txt") , emit: hist
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pileup.sh \\
        -Xmx${task.memory.toGiga()}g \\
        in=${bam} \\
        out=${prefix}.coverage.stats.txt \\
        hist=${prefix}.coverage.hist.txt \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
