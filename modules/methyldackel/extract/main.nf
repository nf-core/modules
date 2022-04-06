process METHYLDACKEL_EXTRACT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::methyldackel=0.6.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/methyldackel:0.6.0--h22771d5_0' :
        'quay.io/biocontainers/methyldackel:0.6.0--h22771d5_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.bedGraph"), emit: bedgraph
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    MethylDackel extract \\
        $args \\
        $fasta \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methyldackel: \$(MethylDackel --version 2>&1 | cut -f1 -d" ")
    END_VERSIONS
    """
}
