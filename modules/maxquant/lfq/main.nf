process MAXQUANT_LFQ {
    tag "$meta.id"
    label 'process_long'

    conda (params.enable_conda ? "bioconda::maxquant=2.0.1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/maxquant:2.0.1.0--py39hdfd78af_2"
    } else {
#        container "wombatp/maxquant-pipeline:dev"
	container "quay.io/biocontainers/maxquant:2.0.1.0
    }

    input:
    tuple val(meta), path(fasta), path(paramfile)
    path raw

    output:
    tuple val(meta), path("*.txt"), emit: maxquant_txt
    path "*.version.txt"          , emit: version

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"  
    
    """
	
    maxquant --version | head -n1 - >  maxquant.version.txt
    sed \"s_<numThreads>.*_<numThreads>$task.cpus</numThreads>_\" ${paramfile} > mqpar_changed.xml
    sed -i \"s|PLACEHOLDER|\$PWD/|g\" mqpar_changed.xml
    mkdir temp
    maxquant mqpar_changed.xml
    mv combined/txt/*.txt .
    """
}
