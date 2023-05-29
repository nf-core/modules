process METHYLDACKEL_EXTRACT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::methyldackel=0.6.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/methyldackel:0.6.0--h22771d5_0' :
        'biocontainers/methyldackel:0.6.0--h22771d5_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.bedGraph") , optional: true, emit: bedgraph
    tuple val(meta), path("*.methylKit"), optional: true, emit: methylkit
    path  "versions.yml"                                , emit: versions

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
