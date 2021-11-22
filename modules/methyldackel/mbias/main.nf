process METHYLDACKEL_MBIAS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::methyldackel=0.6.0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/methyldackel:0.6.0--h22771d5_0"
    } else {
        container "quay.io/biocontainers/methyldackel:0.6.0--h22771d5_0"
    }

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.mbias.txt"), emit: txt
    path  "versions.yml"                , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    MethylDackel mbias \\
        $args \\
        $fasta \\
        $bam \\
        $prefix \\
        --txt \\
        > ${prefix}.mbias.txt

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(MethylDackel --version 2>&1 | cut -f1 -d" ")
    END_VERSIONS
    """
}
