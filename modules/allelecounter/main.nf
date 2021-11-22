process ALLELECOUNTER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::cancerit-allelecount=4.3.0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/cancerit-allelecount:4.3.0--h41abebc_0"
    } else {
        container "quay.io/biocontainers/cancerit-allelecount:4.3.0--h41abebc_0"
    }

    input:
    tuple val(meta), path(input), path(input_index)
    path loci
    path fasta

    output:
    tuple val(meta), path("*.alleleCount"), emit: allelecount
    path "versions.yml"                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def reference_options = fasta ? "-r $fasta": ""

    """
    alleleCounter \\
        $args \\
        -l $loci \\
        -b $input \\
        $reference_options \\
        -o ${prefix}.alleleCount

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(alleleCounter --version)
    END_VERSIONS
    """
}
