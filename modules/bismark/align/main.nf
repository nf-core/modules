process BISMARK_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::bismark=0.23.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bismark:0.23.0--0"
    } else {
        container "quay.io/biocontainers/bismark:0.23.0--0"
    }

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*bam")       , emit: bam
    tuple val(meta), path("*report.txt"), emit: report
    tuple val(meta), path("*fq.gz")     , optional:true, emit: unmapped
    path "versions.yml"                 , emit: versions

    script:
    def prefix     = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def fastq      = meta.single_end ? reads : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    bismark \\
        $fastq \\
        $options.args \\
        --genome $index \\
        --bam

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
