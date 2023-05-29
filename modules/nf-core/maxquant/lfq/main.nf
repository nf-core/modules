process MAXQUANT_LFQ {
    tag "$meta.id"
    label 'process_long'
    conda "bioconda::maxquant=2.0.3.0=py310hdfd78af_1"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/maxquant:2.0.3.0--py310hdfd78af_1"
    } else {
        container "biocontainers/maxquant:2.0.3.0--py310hdfd78af_1"
    }

    input:
    tuple val(meta), path(fasta), path(paramfile)
    path raw

    output:
    tuple val(meta), path("*.txt"), emit: maxquant_txt
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            maxquant: \$(maxquant --version 2>&1 > /dev/null | cut -f2 -d\" \")
    END_VERSIONS
    sed \"s_<numThreads>.*_<numThreads>$task.cpus</numThreads>_\" ${paramfile} > mqpar_changed.xml
    sed -i \"s|PLACEHOLDER|\$PWD/|g\" mqpar_changed.xml
    mkdir temp
    maxquant mqpar_changed.xml
    mv combined/txt/*.txt .
    """
}
