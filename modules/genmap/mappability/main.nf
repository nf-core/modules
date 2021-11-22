process GENMAP_MAPPABILITY {
    tag '$fasta'
    label 'process_high'

    conda (params.enable_conda ? "bioconda::genmap=1.3.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/genmap:1.3.0--h1b792b2_1"
    } else {
        container "quay.io/biocontainers/genmap:1.3.0--h1b792b2_1"
    }

    input:
    path index

    output:
    path "*.wig"        , optional:true, emit: wig
    path "*.bedgraph"   , optional:true, emit: bedgraph
    path "*.txt"        , optional:true, emit: txt
    path "versions.yml"                , emit: versions

    script:
    """
    genmap \\
        map \\
        $args \\
        -I $index \\
        -O mappability

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(genmap --version 2>&1 | sed 's/GenMap version: //; s/SeqAn.*\$//')
    END_VERSIONS
    """
}
