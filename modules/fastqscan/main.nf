process FASTQSCAN {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::fastq-scan=0.4.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/fastq-scan:0.4.4--h7d875b9_0"
    } else {
        container "quay.io/biocontainers/fastq-scan:0.4.4--h7d875b9_0"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.json"), emit: json
    path "versions.yml"            , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    zcat $reads | \\
        fastq-scan \\
        $options.args > ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo \$(fastq-scan -v 2>&1) | sed 's/^.*fastq-scan //' )
    END_VERSIONS
    """
}
