process COOLER_DIGEST {
    tag "$fasta"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::cooler=0.8.11" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/cooler:0.8.11--pyh3252c3a_0"
    } else {
        container "quay.io/biocontainers/cooler:0.8.11--pyh3252c3a_0"
    }

    input:
    path fasta
    path chromsizes
    val  enzyme

    output:
    path "*.bed"                  , emit: bed
    path "versions.yml"           , emit: versions

    script:
    """
    cooler digest \\
        $options.args \\
        -o "${fasta.baseName}_${enzyme.replaceAll(/[^0-9a-zA-Z]+/, '_')}.bed" \\
        $chromsizes \\
        $fasta \\
        $enzyme

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(cooler --version 2>&1 | sed 's/cooler, version //')
    END_VERSIONS
    """
}
